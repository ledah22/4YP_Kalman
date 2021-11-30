run('artifact_plots')

%% Section 1: preprocessing to filter out 50Hz
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


%% Section 2: define the system properties:

end_samples = 30000;
Ts = y_arm.interval;

noise_count = 1.25/Ts;
meanE = mean(y_arm_values(1:noise_count));
%sigmaE = sqrt(sum((y_arm_values(1:noise_count)-meanE).^2)/(noise_count-1));
sigmaE = (6.5*0.001/3); %%observed from the data, assume the range where most data can be found is +-3*sigma

meanTMS = mean(y_air_values);       %%CHECK
%sigmaV = double(sum((y_air_values - meanTMS).^2)/y_air.points);
sigmaV = 0.25;

p = 3;
eeg = ar(y_arm_values(1:noise_count), p)

tms_data = iddata(y_arm_values, (double(1:1:end_samples)*Ts)', y_arm.interval);
na = 3; nb = 3; nk = 1;
tms = arx(tms_data, [na nb nk])

%% Section 3: State space representation
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

%% Section 4: Noise and sampling time

% Tunning parameters:
sigmaT2 = 1;    %WHY DOES THIS NOT MAKE ANY DIFFERENCE
pulse_start = 14080; % set to 250 when input_data uses the air sample
d = 5;
alpha = 0.05;
%d_tot = end_samples - pulse_start;
d_tot = 250;

% Process noise covariance:
Q = zeros(30000*(ro+1), ro+1);
for t=1:30000
    if t<pulse_start || t>pulse_start+d
        Q((t-1)*(ro+1)+1:t*(ro+1), :) = [sigmaE^2 zeros(1, ro); zeros(ro, 1) zeros(ro)];
    else
        Q((t-1)*(ro+1)+1:t*(ro+1), :) = [sigmaE^2 zeros(1, ro); zeros(ro, 1) sigmaT2*eye(ro)];
    end
end
%Q = [sigmaE zeros(1, ro); zeros(ro, 1) sigmaT2*eye(ro)];

% Measurement noise covariance:
R = zeros(1, 30000);
for t=1:30000
    if t>pulse_start && t<pulse_start + d
        R(t) = sigmaV^2;
    else
        if t>pulse_start + d && t<pulse_start + d_tot
            R(t) = sigmaV^2*exp(-alpha*(t-pulse_start-d));
        end
    end
end

% Input to the system is the TMS pulse:
t=0:Ts:(3-Ts);
tms_pulses=zeros(size(t));
tms_pulses(pulse_start)=0.5;

% Generate data for testing purposes:
f_alpha = 10;
noise_data = rand(1, end_samples)*sqrt(sigmaE);
brain_data = 80*10^(-3)*sin(2*pi*f_alpha*(0:Ts:(3-Ts)));
test_data = brain_data + noise_data + tms_pulses;

%Define the input data:
input_data = y_arm_values;
input_data = test_data;
