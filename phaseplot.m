function [] = phaseplot(output_data, Ts, input_data, plotHilbert)
[l, num_rhythms] = size(output_data);

figure
plot(input_data)
xlabel('Time (num. of samples)')
ylabel('Amplitude (Radians and Volts)')
hold on
for j = 1:num_rhythms
    subplot(num_rhythms, 1, j);
    plot(1:1:l, output_data(:, j));
end

if plotHilbert
    input_data_analytical = hilbert(input_data);
    input_data_phase = angle(input_data_analytical);
    
    hold on
    plot(input_data_phase)
    legend('Input signal', 'Phase', 'Hilbert Transform')
    title("Phase tracked at sampling time of "+ 1/Ts + "Hz")
else
    legend('Input signal', 'Processed Signal')
    title("Amplitude removal at sampling time of "+ 1/Ts + "Hz")
end
% 
% figure
% for j = 1:num_rhythms
%     subplot(num_rhythms, 1, j);
%     plot(1:1:l, output_data(:, j)-input_data);
%     title("Estimation error")
% end