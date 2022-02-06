%% Section 1: Experiment settings and general values
model.state_dimensions = 2;
model.Ts = 0.0005;                % Sampling Frequency
model.sim_dur = 1.5;                % Used by Simulink
model.num_rhythms = 1;            % Denoted as N
model.num_samples = model.sim_dur/model.Ts;

model.param_est_length = 2000; %Number of starting samples used for parameter estimation

%WARNING!!! If you change num_rhythms, make sure to adapt input_data and f
%to this!!!

% Generate clean data for testing purposes:
f_alpha = 6;                % Rhythm that we want to track
f_alpha2 = 12;
ampl = 80*10^(-3);
brain_data = ampl*sin(2*pi*f_alpha*(0:model.Ts:(model.sim_dur - model.Ts)));
brain_data2 = ampl*4*sin(2*pi*f_alpha2*(0:model.Ts:(model.sim_dur - model.Ts)));
noise_var = 0.001*80;
noise_data = rand(1, model.num_samples)*sqrt(noise_var);

test_simple = brain_data';
test_noise = (brain_data + noise_data)';
test_multi = (brain_data + brain_data2)';

input_data = test_simple;
tms_pulses = zeros(length(input_data), model.state_dimensions*model.num_rhythms);  % NO INPUT FOR NOW

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

model.f = f_single;
model.a = 0.99*eye(2*model.num_rhythms);
model.omega = 2*pi*model.Ts*model.f;
model.Q_const = noise_var;
model.sigma2_R = 1*eye(2*model.num_rhythms);

O_cell = cell(model.num_rhythms);

for j = 1:model.num_rhythms
    O_j = [cos(model.omega(j)) -sin(model.omega(j)); sin(model.omega(j)) cos(model.omega(j))];
    O_cell{j} = O_j;
end
model.O =blkdiag(O_cell{:});
model.M = repmat([1 0], 1, model.num_rhythms);

%% Section 4: Kalman filter state space
% Initialisation
model.x_0 = zeros(model.state_dimensions*model.num_rhythms, 1);
model.Q = Q_const*eye(model.state_dimensions*model.num_rhythms);

%State space
model.A = model.a*eye(model.state_dimensions*model.num_rhythms)*model.O;
model.B=zeros(model.state_dimensions*model.num_rhythms);    %% this is where we could add tms?
model.Gm = model.Q;
model.C = model.M;

%FOR NON-STATIONARY Q = zeros(num_samples*model.state_dimensions*num_rhythms, model.state_dimensions*num_rhythms);
%R = sigma2_R*eye(1, num_samples);
model.R = sigma2_R;
%Q = zeros(model.state_dimensions*num_rhythms, 2*num_rhythms);