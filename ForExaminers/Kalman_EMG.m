run('artifact_plots')

%% Section 1: Preprocessing of EMG data to filter out 50Hz
%Source: https://dsp.stackexchange.com/a/1090/59989
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


%% Section 2: Define the system properties and testing data:

using_real_data = 1; %Set to 0 if testing on synthesised data, set to 1 if using the EMG measurements
end_samples_values = [3000, 30000];
end_samples = end_samples_values(using_real_data+1);

Ts_values = [0.001, y_arm.interval];
Ts = Ts_values(using_real_data+1);

noise_count = 1.25/Ts;
sigmaE = 0.025;
sigmaV = 0.25;

% Input to the system is the TMS pulse:
pulses_values = [500, 14080];
pulses = pulses_values(using_real_data+1);
intensities = 0.5;
t=0:Ts:((end_samples*Ts)-Ts);
tms_pulses=zeros(size(t));
tms_pulses(pulses)=intensities;

% Generate data for testing purposes:
f_alpha = 8;
ampl = 80*10^(-3);
noise_data = normrnd(0, ampl/10, [1, end_samples]);
brain_data = ampl*sin(2*pi*f_alpha*(0:Ts:(3-Ts)));
test_data = brain_data + noise_data + tms_pulses;

%Define the input data:
if using_real_data
    input_data = y_arm_values';
else
    input_data = test_data;
end

p = 3;
eeg = ar(input_data(1:noise_count), p);
tms_data = iddata(input_data', (double(1:1:end_samples)*Ts)', Ts);

na = 3; nb = 3; nk = 1;
tms = arx(tms_data, [na nb nk])

%% Section 3: State-space representation
Ae = [- eeg.A(2:(p+1)); eye(p-1, p)];
Ge = eye(p, 1);
Ce = eye(1, p);

ro = max(na, nb);
At = zeros(ro);
Bt = (tms.B/tms.A)*eye(ro, 1);
Ct = eye(1, ro);
% Non-stationary behaviour:
Gm = [Ge zeros(p, ro); zeros(ro, 1) eye(ro)];

% Complete system:
A = [Ae zeros(p, ro); zeros(ro, p) At];
B = [zeros(p, 1); Bt];
G = Gm;
C = [Ce Ct];

%Initial values:
x0 = [0.3 0 0 0 0 0];
P0 = eye(ro+p);

%% Section 4: Noise and sampling time

% Tunning parameters:
sigmaT2 = 0.1;
pulse_starts = [pulses, 14080];
pulse_start = pulse_starts(using_real_data+1)-1;

d = 5;
alphas = [0.5, 0.1];
alpha = alphas(using_real_data+1);

d_tot_values = [15, 250];
d_tot = d_tot_values(using_real_data+1);


% Process noise covariance:
Q = zeros(end_samples*(ro+1), ro+1);
for t=1:end_samples
    if t<pulse_start || t>pulse_start+d
        Q((t-1)*(ro+1)+1:t*(ro+1), :) = [sigmaE^2 zeros(1, ro); zeros(ro, 1) zeros(ro)];
    else
        Q((t-1)*(ro+1)+1:t*(ro+1), :) = [sigmaE^2 zeros(1, ro); zeros(ro, 1) sigmaT2*eye(ro)];
    end
end

% Measurement noise covariance:
R = zeros(1, end_samples)+0.005;

for t=1:30000
    if t>pulse_start && t<pulse_start + d
        R(t) = sigmaV^2;
    else
        if t>pulse_start + d && t<pulse_start + d_tot
            R(t) = sigmaV^2*exp(-alpha*(t-pulse_start-d));
        end
    end
end
