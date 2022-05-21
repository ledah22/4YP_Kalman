function [phase_model, artefact_model, input_data, tms_pulses] = tremor_kalman_test(acc_data, Ts)
phase_model.phase_model_dimensions = 2;
phase_model.Ts = Ts;                  % Sampling Frequency
phase_model.num_rhythms = 1;             % Denoted as N
%phase_model.num_samples = phase_model.sim_dur/phase_model.Ts;
data_length = length(acc_data);
phase_model.num_samples = floor(0.8*data_length);
phase_model.sim_dur = phase_model.num_samples*phase_model.Ts;                 % Used by Simulink
phase_model.param_est_window = 2000;

%WARNING!!! If you change num_rhythms, make sure to adapt input_data and f

%% Generate data for testing purposes:
tms_pulses=zeros(phase_model.num_samples, 1);
pulse_start = 450;
pulse_height = 0;       %for testing acc tremor data
tms_pulses(pulse_start) = pulse_height;

input_data = acc_data(1:phase_model.num_samples) - mean(acc_data(1:phase_model.num_samples));

%% Section 4: Noise characterisation

% from artefacts:
d = 5;
alpha = 0.5;
d_tot = 15;
sigmaV = 0.25;                  % used for R1
sigmaE = 0.025;                 % used for Q1

%% Phase tracking parameters
phase_model.f = 2;            % Expected frequency of the signal
phase_model.a = 0.99*eye(2*phase_model.num_rhythms);
phase_model.omega = 2*pi*phase_model.Ts*phase_model.f;

% Calculations for the A matrix of this subsystem:
O_cell = cell(phase_model.num_rhythms);
for j = 1:phase_model.num_rhythms
    O_j = [cos(phase_model.omega(j)) -sin(phase_model.omega(j)); sin(phase_model.omega(j)) cos(phase_model.omega(j))];
    O_cell{j} = O_j;
end
O =blkdiag(O_cell{:});
M = repmat([1 0], 1, phase_model.num_rhythms);

%% Section 3: State space representation - Concatenated system
phase_model.A = phase_model.a*O;
phase_model.B = [1; 0];
phase_model.C = M;

artefact_model.A = [1 0; 0 1];
artefact_model.B = [0; 0];
artefact_model.C = [1, 0];
artefact_model.D = 1;

% Process noise covariance:
phase_model.Q = 0*eye(phase_model.phase_model_dimensions);
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
phase_model.R = 0.8;
artefact_model.R = zeros(1, phase_model.num_samples)+std(input_data);
for t=1:phase_model.num_samples
    for pulse = pulse_start
        if t>pulse && t<pulse + d
            %phase_model.R1(t) = sigmaV+pulse_height(count);
            %artefact_model.R(t) = 0;
        else
            if t>=pulse + d && t<pulse + d_tot
                %phase_model.R1(t) = (sigmaV+pulse_height(count))*exp(-alpha*(t-pulse-d));
                %artefact_model.R(t) = 0;
            end
        end
    end
end

%% Initialise values:
phase_model.P0 = 0.001*eye(phase_model.phase_model_dimensions);
phase_model.x_0 = zeros(phase_model.phase_model_dimensions*phase_model.num_rhythms, 1);
phase_model.x_0 = [0; 1];
%% Tune Parameters:
 [phase_model.A_test, phase_model.a, Q_phase, phase_model.R] = modular_parameter_estimation(input_data', phase_model);
 phase_model.Q = Q_phase*eye(phase_model.phase_model_dimensions);
phase_model.A = phase_model.a*O;
%phase_model.R = 1.075;
end