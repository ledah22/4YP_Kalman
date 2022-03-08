run('artifact_plots')

combined_model.phase_model_dimensions = 2;
combined_model.Ts = 0.0005;                 % Sampling Frequency
combined_model.sim_dur = 3;                 % Used by Simulink
combined_model.num_rhythms = 1;             % Denoted as N
combined_model.num_samples = combined_model.sim_dur/combined_model.Ts;
combined_model.window = 5000;               % Length of the sliding window of available input data
combined_model.param_est_length = 2000;     %Number of starting samples used for parameter estimation

%WARNING!!! If you change num_rhythms, make sure to adapt input_data and f

%% Generate data for testing purposes:
t=0:combined_model.Ts:(combined_model.sim_dur-combined_model.Ts);
f_alpha = 8;
f_alpha2 = 12;
ampl = 80*10^(-3);
brain_data = ampl*sin(2*pi*f_alpha*t);
brain_data2 = ampl*4*sin(2*pi*f_alpha2*t);

strong_noise = rand(1, combined_model.num_samples)*sqrt(ampl);
weak_noise = rand(1, combined_model.num_samples)*sqrt(ampl/100);
use_weak_noise = 1;

noise_var = ampl;
noise_data = strong_noise;
if use_weak_noise
    noise_var = ampl/100;
    noise_data = weak_noise;
end

tms_pulses=zeros(combined_model.num_samples, 1);
pulse_start = 500; % set to 250 when input_data uses the air sample
tms_pulses(pulse_start)=0.5;

test_simple = brain_data';
test_noise = (brain_data + noise_data)';
test_multi = (brain_data + brain_data2)';

input_data = test_simple + tms_pulses;

%% Section 4: Noise characterisation

% from artefacts:
phase_sigmaT2 = 1;    %WHY DOES THIS NOT MAKE ANY DIFFERENCE
d = 5;
alpha = 0.5;
d_tot = combined_model.num_samples - pulse_start;

%from phase:
Q2 = noise_var*eye(combined_model.phase_model_dimensions);
R2 = 10;     % >=10 for strong noise, can even be zero for weak noise

%% Artefact removal parameters
%Prekontrolisi!!!!!!!!!!!!!!!!!!!!!>>>>>>
fs = 1/y_arm.interval;
fn = fs/2;
f0 = 50;
f_ratio = f0/fn;

notchWidth = 0.1;
notchZeros = [exp(sqrt(-1)*pi*f_ratio), exp(-sqrt(-1)*pi*f_ratio)];
notchPoles = (1-notchWidth) * notchZeros;

num = poly(notchZeros); %  Get moving average filter coefficients
denum = poly(notchPoles); %  Get autoregressive filter coefficients

y_arm_values = filter(num, denum, y_arm.values);
y_air_values = filter(num, denum, y_air.values);

noise_count = 12500;
meanE = mean(y_arm_values(1:noise_count));
sigmaE = (6.5*0.001/3); %%observed from the data, assume the range where most data can be found is +-3*sigma

meanTMS = mean(y_air_values);
sigmaV = 0.25;

p = 3;
eeg = ar(y_arm_values(1:noise_count), p);

tms_data = iddata(y_arm_values, (double(1:1:y_arm.points)*y_arm.interval)', y_arm.interval);
na = 3; nb = 3; nk = 1;
tms = arx(tms_data, [na nb nk]);

%% Section 3: State space representation - Artefact removal subsystem
Ae = [- eeg.A(2:(p+1)); eye(p-1, p)];
Ge = eye(p, 1);
Ce = eye(1, p);

ro = max(na, nb);
At = zeros(ro);
Bt = (tms.B/tms.A)*eye(ro, 1);
Ct = eye(1, ro);

% Complete system:
A1 = [Ae zeros(p, ro); zeros(ro, p) At];
B1 = [zeros(p, 1); Bt];
C1 = [Ce Ct];
G1 = [Ge zeros(p, ro); zeros(ro, 1) eye(ro)];

%% Phase tracking parameters
combined_model.f = f_alpha;            % Expected frequency of the signal
combined_model.a = 0.99*eye(2*combined_model.num_rhythms);
combined_model.omega = 2*pi*combined_model.Ts*combined_model.f;

% Calculations for the A matrix of this subsystem:
O_cell = cell(combined_model.num_rhythms);
for j = 1:combined_model.num_rhythms
    O_j = [cos(combined_model.omega(j)) -sin(combined_model.omega(j)); sin(combined_model.omega(j)) cos(combined_model.omega(j))];
    O_cell{j} = O_j;
end
O =blkdiag(O_cell{:});
M = repmat([1 0], 1, combined_model.num_rhythms);

%% Section 3: State space representation - Phase tracking subsystem

A2 = combined_model.a*eye(combined_model.phase_model_dimensions*combined_model.num_rhythms)*O;
B2 = [0; 0];
C2 = M;

%% Section 3: State space representation - Concatenated system
combined_model.A = [A1 zeros(ro+p, combined_model.phase_model_dimensions)];
combined_model.B = [B1; B2];
combined_model.C = [zeros(1, ro+p) C2];
combined_model.C = [combined_model.C; combined_model.C];

% Process noise covariance:
combined_model.Q = zeros(combined_model.num_samples*(ro+p+combined_model.phase_model_dimensions), ro+p+combined_model.phase_model_dimensions);
for t=1:combined_model.num_samples
    if t<pulse_start || t>pulse_start+d
        Q1_temp = [sigmaE^2 zeros(1, ro); zeros(ro, 1) zeros(ro)];
    else
        Q1_temp = [sigmaE^2 zeros(1, ro); zeros(ro, 1) phase_sigmaT2*eye(ro)];    
    end
    combined_model.Q((t-1)*(ro+p+combined_model.phase_model_dimensions)+1:t*(ro+p+combined_model.phase_model_dimensions), :) = [G1*Q1_temp*G1' zeros(6, 2); zeros(2, 6) Q2];
end

% IS IT VALID TO SUM R'S?!?!??!
% Measurement noise covariance:
combined_model.R = zeros(1, combined_model.num_samples);
for t=1:combined_model.num_samples
    if t>pulse_start && t<pulse_start + d
        combined_model.R(t) = sigmaV^2 + R2;
    else
        if t>pulse_start + d && t<pulse_start + d_tot
            combined_model.R(t) = sigmaV^2*exp(-alpha*(t-pulse_start-d)) + R2;
        end
    end
end

%% Initialise values:
combined_model.P0 = eye(ro+p+combined_model.phase_model_dimensions);
combined_model.x_0 = zeros((ro+p+combined_model.phase_model_dimensions)*combined_model.num_rhythms, 1);

f_curr = combined_model.f;
a_curr = combined_model.a;
sigmaR_curr = combined_model.R;
f = f_curr;
a = a_curr;
sigmaR = sigmaR_curr;

xhat = zeros(combined_model.window+1, ro+p+combined_model.phase_model_dimensions);
P = zeros((ro+p+combined_model.phase_model_dimensions)*(combined_model.window+1), ro+p+combined_model.phase_model_dimensions);
x_time = zeros(combined_model.window+1, ro+p+combined_model.phase_model_dimensions);
P_time = zeros((ro+p+combined_model.phase_model_dimensions)*(combined_model.window+1), ro+p+combined_model.phase_model_dimensions);