%% Section 1: Experiment settings and general values
Ts = 0.0001;                % Sampling Frequency
sim_dur = 1.5;                % Used by Simulink
num_rhythms = 1;            % Denoted as N
num_samples = sim_dur/Ts;

%WARNING!!! If you change num_rhythms, make sure to adapt input_data and f
%to this!!!

% Generate clean data for testing purposes:
f_alpha = 6;                % Rhythm that we want to track
f_alpha2 = 12;
ampl = 80*10^(-3);
brain_data = ampl*sin(2*pi*f_alpha*(0:Ts:(sim_dur - Ts)));
brain_data2 = ampl*4*sin(2*pi*f_alpha2*(0:Ts:(sim_dur - Ts)));
noise_var = 0.001*80;
noise_data = rand(1, num_samples)*sqrt(noise_var);

test_simple = brain_data';
test_noise = (brain_data + noise_data)';
test_multi = (brain_data + brain_data2)';

input_data = test_simple;
tms_pulses = zeros(length(input_data), 2*num_rhythms);  % NO INPUT FOR NOW

%% POTENTIALLY NOT NEEDED ANYMORE BECAUSE EM METHOD PICKS THE FREQ OF INTEREST?
% CHECK ON REAL DATA AFTERWARDS TO CONFIRM
% Section 2: preprocessing to filter out 50Hz
%Source: https://dsp.stackexchange.com/a/1090/59989
% fs = 1/Ts;
% fn = fs/2;
% f0 = 50;
% f_ratio = f0/fn;
% 
% notchWidth = 0.1;
% notchZeros = [exp(sqrt(-1)*pi*f_ratio), exp(-sqrt(-1)*pi*f_ratio)];
% notchPoles = (1-notchWidth) * notchZeros;
% 
% num = poly(notchZeros); %  Get moving average filter coefficients
% denum = poly(notchPoles); %  Get autoregressive filter coefficients

%input_data = filter(num, denum, input_data);

%% Section 3: Parameter estimation
f_single = [f_alpha];
f_multi = [f_alpha, f_alpha2];

f = f_single;
a = 0.99;
omega = 2*pi*Ts*f;
Q_const = noise_var;
sigma2_R = 1;

O_cell = cell(num_rhythms);

for j = 1:num_rhythms
    O_j = [cos(omega(j)) -sin(omega(j)); sin(omega(j)) cos(omega(j))];
    O_cell{j} = O_j;
end
O =blkdiag(O_cell{:});
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

%FOR NON-STATIONARY Q = zeros(num_samples*2*num_rhythms, 2*num_rhythms);
%R = sigma2_R*eye(1, num_samples);
R = sigma2_R;
%Q = zeros(2*num_rhythms, 2*num_rhythms);