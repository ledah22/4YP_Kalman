run('artifact_plots')

phase_ns_model.phase_model_dimensions = 2;
phase_ns_model.Ts = 0.0005;                 % Sampling Frequency
phase_ns_model.sim_dur = 3;                 % Used by Simulink
phase_ns_model.num_rhythms = 1;             % Denoted as N
phase_ns_model.num_samples = phase_ns_model.sim_dur/phase_ns_model.Ts;
phase_ns_model.window = 5000;               % Length of the sliding window of available input data
phase_ns_model.param_est_length = 2000;     %Number of starting samples used for parameter estimation
phase_ns_model.clk = 1:1:phase_ns_model.num_samples;

%WARNING!!! If you change num_rhythms, make sure to adapt input_data and f

%% Generate data for testing purposes:
t= 0:phase_ns_model.Ts:(phase_ns_model.sim_dur-phase_ns_model.Ts);
f_alpha = 8;
f_alpha2 = 12;
ampl = 80*10^(-3);
brain_data = ampl*sin(2*pi*f_alpha*t);
brain_data2 = ampl*4*sin(2*pi*f_alpha2*t);

strong_noise = rand(1, phase_ns_model.num_samples)*sqrt(ampl);
weak_noise = rand(1, phase_ns_model.num_samples)*sqrt(ampl/100);
use_weak_noise = 1;

noise_var = ampl;
noise_data = strong_noise;
if use_weak_noise
    noise_var = ampl/100;
    noise_data = weak_noise;
end

tms_pulses=zeros(phase_ns_model.num_samples, 1);
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
%d_tot = phase_ns_model.num_samples - pulse_start;
d_tot = 5;

%from phase:
Q2 = noise_var*eye(phase_ns_model.phase_model_dimensions);
R2 = 1000;     % >=10 for strong noise, can even be zero for weak noise

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
phase_ns_model.f = f_alpha;            % Expected frequency of the signal
phase_ns_model.a = 0.99*eye(2*phase_ns_model.num_rhythms);
phase_ns_model.omega = 2*pi*phase_ns_model.Ts*phase_ns_model.f;

% Calculations for the A matrix of this subsystem:
O_cell = cell(phase_ns_model.num_rhythms);
for j = 1:phase_ns_model.num_rhythms
    O_j = [cos(phase_ns_model.omega(j)) -sin(phase_ns_model.omega(j)); sin(phase_ns_model.omega(j)) cos(phase_ns_model.omega(j))];
    O_cell{j} = O_j;
end
O =blkdiag(O_cell{:});
M = repmat([1 0], 1, phase_ns_model.num_rhythms);

%% Section 3: State space representation - Phase tracking subsystem

A2 = phase_ns_model.a*O;
B2 = [0; 0];
C2 = M;

%% Section 3: State space representation - Concatenated system
phase_ns_model.A = A2;
phase_ns_model.B = [1; 0];
phase_ns_model.C = C2;
phase_ns_model.C = [phase_ns_model.C; phase_ns_model.C];

model1.A = A1;
model1.B = B1;
model1.C = C1;
model1.G = G1;

model2.A = A2;
model2.B = B2;
model2.C = C2;

% Process noise covariance:
model1.Q = zeros(phase_ns_model.num_samples*(ro+1), ro+1);
model2.Q = Q2;
phase_ns_model.Q = Q2;
for t=1:phase_ns_model.num_samples
    if t<pulse_start || t>pulse_start+d
        Q1_temp = [sigmaE^2 zeros(1, ro); zeros(ro, 1) zeros(ro)];
    else
        Q1_temp = [sigmaE^2 zeros(1, ro); zeros(ro, 1) phase_sigmaT2*eye(ro)];    
    end
    model1.Q((t-1)*(ro+1)+1:t*(ro+1), :) = Q1_temp;
    %phase_ns_model.Q((t-1)*(ro+p+phase_ns_model.phase_model_dimensions)+1:t*(ro+p+phase_ns_model.phase_model_dimensions), :) = [G1*Q1_temp*G1' zeros(6, 2); zeros(2, 6) Q2];
end

% IS IT VALID TO SUM R'S?!?!??!
% Measurement noise covariance:
phase_ns_model.R = zeros(1, phase_ns_model.num_samples);
model1.R = zeros(1, phase_ns_model.num_samples);
model2.R = R2;
for t=1:phase_ns_model.num_samples
    if t>pulse_start && t<pulse_start + d
        model1.R(t) = sigmaV^2;
        phase_ns_model.R(t) = model1.R(t) + R2;
    else
        if t>pulse_start + d && t<pulse_start + d_tot
            model1.R(t) = sigmaV^2*exp(-alpha*(t-pulse_start-d));
            phase_ns_model.R(t) = model1.R(t) + R2;
        end
    end
end

%% Initialise values:
phase_ns_model.P0 = eye(phase_ns_model.phase_model_dimensions);
phase_ns_model.x_0 = zeros(phase_ns_model.phase_model_dimensions*phase_ns_model.num_rhythms, 1);

model1.P0 = eye(ro+p);
model1.x_0 = zeros((ro+p)*phase_ns_model.num_rhythms, 1);
model2.P0 = eye(phase_ns_model.phase_model_dimensions);
model2.x_0 = zeros(phase_ns_model.phase_model_dimensions*phase_ns_model.num_rhythms,1);

% f_curr = phase_ns_model.f;
% a_curr = phase_ns_model.a;
% sigmaR_curr = phase_ns_model.R;
% f = f_curr;
% a = a_curr;
% sigmaR = sigmaR_curr;

xhat = zeros(phase_ns_model.window+1, phase_ns_model.phase_model_dimensions);
P = zeros(phase_ns_model.phase_model_dimensions*(phase_ns_model.window+1), phase_ns_model.phase_model_dimensions);
x_time = zeros(phase_ns_model.window+1, phase_ns_model.phase_model_dimensions);
P_time = zeros(phase_ns_model.phase_model_dimensions*(phase_ns_model.window+1), phase_ns_model.phase_model_dimensions);