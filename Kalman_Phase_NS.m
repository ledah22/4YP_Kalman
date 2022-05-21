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

noise_ratio = [0, 0.0001, 0.01, 0.1, 1, 10];
noise_var = ampl*noise_ratio(3);
noise_data = rand(1, phase_ns_model.num_samples)*sqrt(noise_var);

tms_pulses=zeros(phase_ns_model.num_samples, 1);
pulse_start = 500; % set to 250 when input_data uses the air sample
pulse_height = 0.5;
tms_pulses(pulse_start) = pulse_height;

test_simple = brain_data';
test_noise = (brain_data + noise_data)';
test_multi = (brain_data + brain_data2)';

input_data = test_noise + tms_pulses;

%% Section 4: Noise characterisation

% from artefacts:
d = 5;
alpha = 0.5;
d_tot = 50;
%sigmaV = 0.25; % observed from the data, assume the range where most data can be found is +-3*sigma
sigmaV = 0.25;               % used for R1
sigmaE = (6.5*0.001/3);     % used for Q1

%from phase:
Q2 = noise_var*eye(phase_ns_model.phase_model_dimensions);
R2 = 0;     % >=10 for strong noise, can even be zero for weak noise

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

%% Section 3: State space representation - Concatenated system
phase_ns_model.A = phase_ns_model.a*O;
phase_ns_model.B = [1; 0];
phase_ns_model.C = M;
phase_ns_model.C = [phase_ns_model.C; phase_ns_model.C];


input_data_analytical = hilbert(test_simple);
input_data_phase = angle(input_data_analytical);
p = 2;
syst = ar(input_data_phase, p);
phase_double_model.A = [- syst.A(2:(p+1)); eye(p-1, p)];
phase_double_model.B = [0; 0];
phase_double_model.C = [1, 0];
phase_double_model.C = [phase_double_model.C; phase_double_model.C];
phase_double_model.Q = 0;
phase_double_model.R = 0.05;

% Process noise covariance:
% model1.Q = zeros(phase_ns_model.num_samples*(ro+1), ro+1);
% model2.Q = Q2;
%phase_ns_model.Q = Q2;
phase_ns_model.Q = eye(phase_ns_model.num_samples*phase_ns_model.phase_model_dimensions, phase_ns_model.phase_model_dimensions);
for t=1:phase_ns_model.num_samples
    if t<pulse_start || t>pulse_start+d
        Q1_temp = sigmaE^2*eye(phase_ns_model.phase_model_dimensions);
        %         Q1_temp = [sigmaE^2 zeros(1, ro); zeros(ro, 1) zeros(ro)];
    else
        Q1_temp = (sigmaE^2+pulse_height)*eye(phase_ns_model.phase_model_dimensions);
        %         Q1_temp = [sigmaE^2 zeros(1, ro); zeros(ro, 1) phase_sigmaT2*eye(ro)];
    end
    %     model1.Q((t-1)*(ro+1)+1:t*(ro+1), :) = Q1_temp;
    phase_ns_model.Q((t-1)*(phase_ns_model.phase_model_dimensions)+1:t*(phase_ns_model.phase_model_dimensions), :) = Q1_temp+Q2;
end

% IS IT VALID TO SUM R'S?!?!??!
% Measurement noise covariance:
phase_ns_model.R = R2* eye(1, phase_ns_model.num_samples);
for t=1:phase_ns_model.num_samples
    if t>pulse_start && t<pulse_start + d
        %phase_ns_model.R(t) = sigmaV + R2;
        phase_ns_model.R(t) = sigmaV;
    else
        if t>=pulse_start + d && t<pulse_start + d_tot
            %phase_ns_model.R(t) = sigmaV*exp(-alpha*(t-pulse_start-d)) + R2;
            phase_ns_model.R(t) = sigmaV*exp(-alpha*(t-pulse_start-d));
        end
    end
end

%% Initialise values:
phase_ns_model.P0 = eye(phase_ns_model.phase_model_dimensions);
phase_ns_model.x_0 = zeros(phase_ns_model.phase_model_dimensions*phase_ns_model.num_rhythms, 1);