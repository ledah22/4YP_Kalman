%% Section 1: Experiment settings and general values
model.state_dimensions = 2;
model.Ts = 0.001;                % Sampling Frequency
model.sim_dur = 3;                % Used by Simulink
model.num_rhythms = 1;            % Denoted as N
model.num_samples = model.sim_dur/model.Ts;
model.window = 5000;        % Length of the sliding window of available input data
model.param_est_length = 2000;    %Number of starting samples used for parameter estimation

%WARNING!!! If you change num_rhythms, make sure to adapt input_data and f
%to this!!!

% Generate clean data for testing purposes:
f_alpha = 8;                % Rhythm that we want to track
f_alpha2 = 12;
ampl = 80*10^(-3);
brain_data = ampl*sin(2*pi*f_alpha*(0:model.Ts:(model.sim_dur - model.Ts)));
brain_data2 = ampl*4*sin(2*pi*f_alpha2*(0:model.Ts:(model.sim_dur - model.Ts)));
brain_data3 = brain_data;
t= 0:model.Ts:(model.sim_dur-model.Ts);
reset_start = 450;
burst_start = 70;
brain_data3(reset_start+1:end) = ampl*sin(2*pi*f_alpha*t(reset_start+1:end)+pi/2);
%brain_data3(1:burst_start) = zeros(1, burst_start);

strong_noise = normrnd(0, 5*ampl, [1, model.num_samples]);
weak_noise = normrnd(0, ampl/10, [1, model.num_samples]);
noise_var = (ampl/100); %ampl when using strong noise, ampl/100 for weak noise
noise_var = 10e-4;

test_simple = brain_data';
test_noise_strong = (brain_data + strong_noise)';
test_noise_weak = (brain_data + weak_noise)';
test_multi = (brain_data + brain_data2)';
test_reset = brain_data3;
test_burst = test_noise_strong;
%test_burst = test_simple;
test_burst(1:burst_start) = zeros(1, burst_start);

input_data = test_noise_strong;
%test_simple(1:burst_start) = zeros(1, burst_start);
tms_pulses = zeros(1, length(input_data));

%% Section 3: Parameter estimation
f_single = [f_alpha];
f_multi = [f_alpha, f_alpha2];

model.f = f_single;
model.a = 0.99*eye(2*model.num_rhythms);
model.omega = 2*pi*model.Ts*model.f;
model.Q_const = 0.00001;
model.sigma2_R = noise_var; %>=10 for strong noise, can even be zero for weak noise

f_curr = model.f;
a_curr = model.a;
sigmaR_curr = model.sigma2_R;
f = f_curr;
a = a_curr;
sigmaR = sigmaR_curr;

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
model.P_0 = zeros(model.state_dimensions*model.num_rhythms);
model.Q = model.Q_const*eye(model.state_dimensions*model.num_rhythms);

%State space
model.A = model.a*eye(model.state_dimensions*model.num_rhythms)*model.O;
model.B = [0; 0];
model.Gm = model.Q;
model.C = model.M;
model.C = [[1 0]; [1 0]];
%model.C = [1 0];
model.D = 0;

%FOR NON-STATIONARY Q = zeros(num_samples*model.state_dimensions*num_rhythms, model.state_dimensions*num_rhythms);
%R = sigma2_R*eye(1, num_samples);
model.R = model.sigma2_R;
%Q = zeros(model.state_dimensions*num_rhythms, 2*num_rhythms);

xhat = zeros(model.window+1, 2);
P = zeros(2*(model.window+1), 2);
x_time = zeros(model.window+1, 2);
P_time = zeros(2*(model.window+1), 2);