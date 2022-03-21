phase_ns_model.phase_model_dimensions = 2;
phase_ns_model.Ts = 0.0005;                  % Sampling Frequency
phase_ns_model.sim_dur = 3;                 % Used by Simulink
phase_ns_model.num_rhythms = 1;             % Denoted as N
phase_ns_model.num_samples = phase_ns_model.sim_dur/phase_ns_model.Ts;
phase_ns_model.param_est_window = 2000;

%WARNING!!! If you change num_rhythms, make sure to adapt input_data and f

%% Generate data for testing purposes:
t= 0:phase_ns_model.Ts:(phase_ns_model.sim_dur-phase_ns_model.Ts);
f_alpha = 8;
f_alpha2 = 12;
ampl = 80*10^(-3);
brain_data = ampl*sin(2*pi*f_alpha*t);
brain_data2 = brain_data;
brain_data2(3000:end) = ampl*sin(2*pi*f_alpha*t(3000:end)-pi/2);

noise_ratio = [0, 0.0001, 0.01, 0.1, 1, 10];
noise_var = ampl*noise_ratio(1);
noise_data = rand(1, phase_ns_model.num_samples)*sqrt(noise_var);

tms_pulses=zeros(phase_ns_model.num_samples, 1);
% pulse_start = [500, 900, 1740];
% pulse_height = [0.5 0.05 0.2];
pulse_start = 600;
pulse_height = 0.5;
brain_data2(pulse_start:end) = ampl*sin(2*pi*f_alpha*t(pulse_start:end)+pi/2);
tms_pulses(pulse_start) = pulse_height;

test_simple = brain_data';
test_noise = (brain_data + noise_data)';
test_reset = brain_data2';
test_multi = (brain_data + brain_data2)';

input_data = test_reset + tms_pulses;

%% Section 4: Noise characterisation

% from artefacts:
d = 5;
alpha = 0.5;
d_tot = 15;
sigmaV = 0.25;                  % used for R1
sigmaE = 0.025;                 % used for Q1

%from phase:
Q2 = 0*eye(phase_ns_model.phase_model_dimensions);
R2 = noise_var+0.025;

%% Phase tracking parameters
phase_ns_model.f = f_alpha;            % Expected frequency of the signal
phase_ns_model.a = 0.99*eye(2*phase_ns_model.num_rhythms);
phase_ns_model.omega = 2*pi*phase_ns_model.Ts*phase_ns_model.f;

%DODAJ NAMESTANJE PARAM

% Calculations for the A matrix of this subsystem:
O_cell = cell(phase_ns_model.num_rhythms);
for j = 1:phase_ns_model.num_rhythms
    O_j = [cos(phase_ns_model.omega(j)) -sin(phase_ns_model.omega(j)); sin(phase_ns_model.omega(j)) cos(phase_ns_model.omega(j))];
    O_cell{j} = O_j;
end
O =blkdiag(O_cell{:});
M = repmat([1 0], 1, phase_ns_model.num_rhythms);

%% Section 3: State space representation - Concatenated system
phase_ns_model.A = phase_ns_model.a*O;
phase_ns_model.B = [1; 0];
phase_ns_model.C = M;
%phase_ns_model.C = [phase_ns_model.C; phase_ns_model.C];

input_data_analytical = hilbert(test_simple);
input_data_phase = angle(input_data_analytical);
p = 2;
syst = ar(input_data_phase, p);
phase_double_model.A = [- syst.A(2:(p+1)); eye(p-1, p)];
phase_double_model.B = [0; 0];
phase_double_model.C = [1, 0];
phase_double_model.D = 1;
%phase_double_model.C = [phase_double_model.C; phase_double_model.C];

% Process noise covariance:
% model1.Q = zeros(phase_ns_model.num_samples*(ro+1), ro+1);
% model2.Q = Q2;
phase_ns_model.Q2 = Q2;
phase_ns_model.Q1 = eye(phase_ns_model.num_samples*phase_ns_model.phase_model_dimensions, phase_ns_model.phase_model_dimensions);
for t=1:phase_ns_model.num_samples
    for count = 1:length(pulse_start)
        if t<pulse_start(count) || t>pulse_start(count)+d
            Q1_temp = sigmaE^2*eye(phase_ns_model.phase_model_dimensions);
            %Q1_temp = 0*eye(phase_ns_model.phase_model_dimensions);
            %         Q1_temp = [sigmaE^2 zeros(1, ro); zeros(ro, 1) zeros(ro)];
        else
            Q1_temp = (sigmaE^2+pulse_height(count))*eye(phase_ns_model.phase_model_dimensions);
            %Q1_temp = 0*eye(phase_ns_model.phase_model_dimensions);
            %         Q1_temp = [sigmaE^2 zeros(1, ro); zeros(ro, 1) phase_sigmaT2*eye(ro)];
        end
        %     model1.Q((t-1)*(ro+1)+1:t*(ro+1), :) = Q1_temp;
        phase_ns_model.Q1((t-1)*(phase_ns_model.phase_model_dimensions)+1:t*(phase_ns_model.phase_model_dimensions), :) = Q1_temp;
    end
end

% IS IT VALID TO SUM R'S?!?!??!
% Measurement noise covariance:
phase_ns_model.R1 = zeros(1, phase_ns_model.num_samples);
phase_ns_model.R2 = R2;
for t=1:phase_ns_model.num_samples
    for pulse = pulse_start
        if t>pulse && t<pulse + d
            phase_ns_model.R1(t) = (sigmaV+pulse_height(count));
            phase_ns_model.R1(t) = 0;
        else
            if t>=pulse + d && t<pulse + d_tot
                phase_ns_model.R1(t) = (sigmaV+pulse_height(count))*exp(-alpha*(t-pulse-d));
                phase_ns_model.R1(t) = 0;
            end
        end
    end
end

%% Initialise values:
phase_ns_model.P0 = 0.001*eye(phase_ns_model.phase_model_dimensions);
phase_ns_model.x_0 = zeros(phase_ns_model.phase_model_dimensions*phase_ns_model.num_rhythms, 1);

%% Tune Parameters:
[phase_ns_model.A_test, phase_ns_model.a, phase_ns_model.Q2, phase_ns_model.R2] = modular_parameter_estimation(test_noise, phase_ns_model);
phase_ns_model.A = phase_ns_model.a*O;