function [] = phaseplot(thetas, Ts, input_data)
[l, num_rhythms] = size(thetas);
input_data_analytical = hilbert(input_data);
input_data_phase = angle(input_data_analytical);

figure
for j = 1:num_rhythms
    subplot(num_rhythms, 1, j);
    plot(1:1:l, thetas(:, j));
end
title("Phase tracked at sampling time of "+ 1/Ts + "Hz")
xlabel('Time (num. of samples)')
ylabel('Amplitude (Radians and Volts)')
hold on
plot(input_data)
hold on
plot(input_data_phase)
legend('Phase', 'Input signal', 'Hilbert Transform')