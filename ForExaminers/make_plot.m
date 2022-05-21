function [] = make_plot(input, output, fs, num_samples)
close all
time = (0:1:(num_samples-1))/fs;
plot(time, input)
hold on
plot(time, output(2:end))
title('Artefact removal on EMG data')
xlabel('Time (s)')
ylabel('Voltage (mV)')
legend('Input Data', 'Output Data')
set(gca,'FontSize',14)
set(groot,'defaultLineLineWidth',0.01)
%set(groot,'defaultLineLineWidth',2)
end