function [] = phaseplot(thetas, Ts, input_data)
[l, num_rhythms] = size(thetas);
figure
for j = 1:num_rhythms
    subplot(num_rhythms, 1, j);
    %plot(1:1:1000, thetas(1:1000, j));
    plot(1:1:l, thetas(:, j));
end
title("Phase tracked at sampling time of "+ 1/Ts + "Hz")
xlabel('Time (num. of samples)')
ylabel('Amplitude (Radians and Volts)')
hold on
plot(input_data)
legend('Phase', 'Input signal')