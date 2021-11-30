%% Section 1: Experiment settings and general values
Ts = 0.0001;
sim_dur = 3;
num_rhythms = 1;
num_samples = sim_dur/Ts;

% Generate clean data for testing purposes:
f_alpha = 8;
brain_data = 80*10^(-3)*sin(2*pi*f_alpha*(0:Ts:(sim_dur - Ts)));
%brain_data = 2*pi*mod((0:Ts:(num_samples-1)*Ts)/(1/f_alpha), 1/f_alpha);
noise_data = rand(1, num_samples)*sqrt(0.001);

input_data = [(brain_data + noise_data)'; zeros(length(brain_data), 1)];
tms_pulses = zeros(length(input_data), 2);

%% Section 2: preprocessing to filter out 50Hz
%Source: https://dsp.stackexchange.com/a/1090/59989
fs = 1/Ts;
fn = fs/2;
f0 = 50;
f_ratio = f0/fn;

notchWidth = 0.1;
notchZeros = [exp(sqrt(-1)*pi*f_ratio), exp(-sqrt(-1)*pi*f_ratio)];
notchPoles = (1-notchWidth) * notchZeros;

num = poly(notchZeros); %  Get moving average filter coefficients
denum = poly(notchPoles); %  Get autoregressive filter coefficients

%input_data = filter(num, denum, input_data);

%% Section 3: Parameter estimation
f = f_alpha;
a = 1;
omega = 2*pi*f*Ts;
Q_const = 0.001;
sigma2_R = 0;

O =repmat([cos(omega) -sin(omega); sin(omega) cos(omega)], num_rhythms);
M = repmat([1 0], 1, num_rhythms);

%% Section 4: Kalman filter state space
% Initialisation
x_0 = zeros(2*num_rhythms, 1);
Q = Q_const*eye(2*num_rhythms);

%State space
A = a*eye(2*num_rhythms)*O;
B=zeros(2*num_rhythms);    %% this is where we could add tms?
Gm = Q;
C = M;

R = sigma2_R*eye(1, num_samples);
%Q = zeros(num_samples*2*num_rhythms, 2*num_rhythms);