%% Section 1: Experimental param settings
phase_model.phase_model_dimensions = 2;
phase_model.Ts = 0.001;                  % Sampling Frequency
phase_model.sim_dur = 3;                 % Used by Simulink
phase_model.num_rhythms = 1;             % Denoted as N
phase_model.num_samples = phase_model.sim_dur/phase_model.Ts;
phase_model.param_est_window = 2000;

%% Section 2: Generate data for testing purposes:
t= 0:phase_model.Ts:(phase_model.sim_dur-phase_model.Ts);
f_alpha = 6;
f_alpha2 = 12;
ampl = 80*10^(-3);
burst_start = 250;
brain_data = ampl*sin(2*pi*f_alpha*t);
brain_data(1:burst_start) = zeros(1, burst_start);
brain_data2 = brain_data;
brain_data2(3000:end) = ampl*sin(2*pi*f_alpha*t(3000:end)-pi/2);

noise_ratios = [0, 0.0001, 0.01, 0.1, 1, 10];
noise_ratio = noise_ratios(5);
noise_var = ampl*noise_ratio;
noise_data = normrnd(0, noise_var, [1, phase_model.num_samples]);
noise_data(1:burst_start) = zeros(1, burst_start);

tms_pulses=zeros(phase_model.num_samples, 1);
one_pulse = 1; %Set to 0 to test multiple pulses. Set to 1 to test with a
% single TMS pulse. To exclude TMS pulses, set to 1, but set the pulse
% intensity to 0.

if one_pulse
    pulse_start = 450;
    pulse_height = 0.5;
    %pulse_height = 0;       %for testing on accelerometer tremor data
else
    pulse_start = [500, 900, 1240];
    pulse_height = [0.5 -0.25 0.2];
end

brain_data2(pulse_start(1):end) = ampl*sin(2*pi*f_alpha*t(pulse_start(1):end)+pi/2);
brain_data2(1:burst_start) = zeros(1, burst_start);
tms_pulses(pulse_start) = pulse_height;

test_simple = brain_data';
test_noise = (brain_data + noise_data)';
test_reset = (brain_data2 + noise_data)';
test_multi = (brain_data + brain_data2)';
test_data = test_noise;

using_real_data = 0; % 0 for synthesised data, 1 for Parkinsonian tremor
% accelerometer data, 2 for essential tremor accelerometer data

if using_real_data == 0
    input_data = test_data + tms_pulses;
else
    all_data = combine_data();
    if using_real_data == 1
        input_data = all_data.data1(2).x;
        clean_data = all_data.data1(2).xf;
    else
        input_data = reshape(all.data2(2).x, [], 1);
        clean_data = reshape(all.data2(2).xf, [], 1);
    end
end

%% Section 3: Noise characterisation

% from artefact removal:
d = 5;
alpha = 0.5;
d_tot = 15;
sigmaV = 0.25;                  % used for R1
sigmaE = 0.025;                 % used for Q1

%from phase tracking:
Q_phase = 0.0001*eye(phase_model.phase_model_dimensions);
R_phase = noise_var+0.025;

%% Section 4: Phase tracking parameters
phase_model.f = f_alpha;            % Expected frequency of the signal
phase_model.a = 0.99999*eye(2*phase_model.num_rhythms);
phase_model.omega = 2*pi*phase_model.Ts*phase_model.f;

% Calculations for the A matrix of this subsystem:
O_cell = cell(phase_model.num_rhythms);
for j = 1:phase_model.num_rhythms
    O_j = [cos(phase_model.omega(j)) -sin(phase_model.omega(j)); sin(phase_model.omega(j)) cos(phase_model.omega(j))];
    O_cell{j} = O_j;
end
O =blkdiag(O_cell{:});
M = repmat([1 0], 1, phase_model.num_rhythms);

%% Section 5: Artefact Removal parameters

input_data_analytical = hilbert(test_simple);
input_data_phase = angle(input_data_analytical);
p = 2;
syst = ar(input_data_phase, p);

%% Section 6: Concatenated state-space, noise matrices and initialisation
phase_model.A = phase_model.a*O;
phase_model.B = [1; 0];
phase_model.C = M;
artefact_model.A = [- syst.A(2:(p+1)); eye(p-1, p)];
artefact_model.B = [0; 0];
artefact_model.C = [1, 0];
artefact_model.D = 1;

% Process noise covariance:
phase_model.Q = Q_phase*eye(phase_model.phase_model_dimensions);
artefact_model.Q = eye(phase_model.num_samples*phase_model.phase_model_dimensions, phase_model.phase_model_dimensions);
for t=1:phase_model.num_samples
    for count = 1:length(pulse_start)
        if t<pulse_start(count) || t>pulse_start(count)+d
            Q_temp = sigmaE^2*eye(phase_model.phase_model_dimensions);
        else
            Q_temp = (sigmaE^2+pulse_height(count))*eye(phase_model.phase_model_dimensions);
        end
        artefact_model.Q((t-1)*(phase_model.phase_model_dimensions)+1:t*(phase_model.phase_model_dimensions), :) = Q_temp;
    end
end

% Measurement noise covariance:
artefact_model.R = zeros(1, phase_model.num_samples);
phase_model.R = R_phase;

% Initialise values:
phase_model.P0 = 0.001*eye(phase_model.phase_model_dimensions);
phase_model.x_0 = zeros(phase_model.phase_model_dimensions*phase_model.num_rhythms, 1);

%% Section 7: Tune Parameters:
[phase_model.A_test, phase_model.a, Q_phase, phase_model.R] = modular_parameter_estimation(test_simple, phase_model);
phase_model.Q = Q_phase*eye(phase_model.phase_model_dimensions);
phase_model.A = phase_model.a*O;