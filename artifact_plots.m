load('artefact_air_isidora.mat')
load('artefact_arm_isidora.mat')

y_arm = MEP_RC7_20211028_125143_000_wave_data;
x_arm = y_arm.start:1:(y_arm.points-1);
x_arm = double(x_arm);
x_arm = x_arm*0.0001;

y_air = MEP_RC7_20211028_124023_000_wave_data;
x_air = y_air.start:1:(y_air.points-1);
x_air = double(x_air);
x_air = (x_air*1.0*y_air.interval)';


subplot(2, 1, 1)
plot(x_air, y_air.values)
title('Electrical equipment response')
xlabel('Time         [s]')
ylabel('Voltage          [mV]')

subplot(2, 1, 2)
plot(x_arm, y_arm.values)
title('Combined response (muscles and artifacts)')
xlabel('Time         [s]')
ylabel('Voltage          [mV]')