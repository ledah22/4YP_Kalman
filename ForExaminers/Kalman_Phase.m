%% Section 1: Experimental param settings and testing data
model.state_dimensions = 2;
model.Ts = 0.001;                 % Sampling Frequency
model.sim_dur = 3;                % Used by Simulink
model.num_rhythms = 1;            % Denoted as N
model.num_samples = model.sim_dur/model.Ts;
model.window = 5000;
model.param_est_length = 2000;    %Number of starting samples used for parameter estimation

% Generate clean data for testing purposes:
f_alpha = 8;                % Frequency to be tracked
f_alpha2 = 12;
ampl = 80*10^(-3);
brain_data = ampl*sin(2*pi*f_alpha*(0:model.Ts:(model.sim_dur - model.Ts)));
brain_data2 = ampl*4*sin(2*pi*f_alpha2*(0:model.Ts:(model.sim_dur - model.Ts)));
brain_data3 = brain_data;

t= 0:model.Ts:(model.sim_dur-model.Ts);
reset_start = 450;
burst_start = 70;
brain_data3(reset_start+1:end) = ampl*sin(2*pi*f_alpha*t(reset_start+1:end)+pi/2);

noise_ampl = [ampl/100, ampl/10, ampl, 5*ampl];
strong_noise = normrnd(0, noise_ampl(3), [1, model.num_samples]);
weak_noise = normrnd(0, noise_ampl(2), [1, model.num_samples]);
noise_var = 10e-4;

test_simple = brain_data';
test_noise_strong = (brain_data + strong_noise)';
test_noise_weak = (brain_data + weak_noise)';
test_multi = (brain_data + brain_data2)';
test_reset = brain_data3;
test_burst = test_noise_strong;
test_burst(1:burst_start) = zeros(1, burst_start);

input_data = test_noise_strong;
tms_pulses = zeros(1, length(input_data)); %only phase is tracked in this module, artefacts are not considered

%% Section 3: Model and noise parameters
model.f = f_alpha;
model.a = 0.99*eye(2*model.num_rhythms);
model.omega = 2*pi*model.Ts*model.f;
model.Q_const = 0.00001;
model.sigma2_R = noise_var;

f_curr = model.f;
a_curr = model.a;
sigmaR_curr = model.sigma2_R;
f = f_curr;
a = a_curr;

model.Q = model.Q_const*eye(model.state_dimensions*model.num_rhythms);
model.R = model.sigma2_R;

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

%State space
model.A = model.a*eye(model.state_dimensions*model.num_rhythms)*model.O;
model.B = [0; 0];
model.Gm = model.Q;
model.D = 0;
model.C = [model.M; model.M]; % Mathematicaly model.C = model.M, but the
% form written here is to satisfy Simulink dimensional constraints as it
% won't pass 1x2 matrix otherwise.

xhat = zeros(model.window+1, 2);
P = zeros(2*(model.window+1), 2);
x_time = zeros(model.window+1, 2);
P_time = zeros(2*(model.window+1), 2);