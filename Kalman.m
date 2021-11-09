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
% plot(x_arm, y_arm.values)
% hold on
% plot(x_arm, y_arm_values)

y_air_values = filter(num, denum, y_air.values);
% plot(x_arm, y_air.values)
% hold on
% plot(x_arm, y_air_values)

%% Section 2: define the system properties:

% Original code:
% Copyright 2018 The MathWorks, Inc.
% Pendulum model
% Gravity
g = 9.81; % [m/s^2]
% Pendulum mass
m = 1; % [kg]
% Pendulum length
l = 0.5; % [m]

% My code:
noise_count = 1.25/y_arm.interval;
meanE = mean(y_arm_values(1:noise_count));
sigmaE = sqrt(sum((y_arm_values(1:noise_count)-meanE).^2)/(noise_count-1));
sigmaE = (6.5*0.001/3); %%observed from the data, assume the range where most data can be found is +-3*sigma

meanTMS = mean(y_air_values);
sigmaV = double(sum((y_air_values - meanTMS).^2)/y_air.points);

p = 3;
D = 0; q = 0;
eeq = arima(p, D, q)
tms_data = iddata(y_arm.values, (double(1:1:30000)/10000)', y_arm.interval);
na = 3; nb = 3; nk = 1;
tms = arx(tms_data, [na nb nk])

%% Section 3: State space representation
Ae = [-cell2mat(eeq.AR); eye(p-1, p)];
Ge = eye(p, 1);
Ce = eye(1, p);

ro = max(na, nb);
At = zeros(ro);
Bt = tms.B/tms.A;
Ct = eye(1, ro);
% Non-stationary behaviour:
Gm = [Ge zeros(p, ro); zeros(ro, 1) eye(ro)];

A = [Ae zeros(p, ro); zeros(ro, p) At];
B = [zeros(p, 1); Bt];
G = Gm;
C = [Ce Ct];

%% Section 4: Noise and sampling time
sigmaT2 = zeros()
% Process noise covariance
Q = [sigmaE zeros(1, ro); zeros(ro, 1) sigmaT2*eye(ro)];
% Measurement noise covariance
R = sigmaV;
% Sampling time
Ts = y_arm.interval; % [s] 