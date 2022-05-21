function [] = phaseplot(output_data, f, Ts, input_data, plotHilbert, clean_data, burst_start, reset_start)
if reset_start == -1
    reset_start = length(clean_data);
end

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
% set(gca,'FontSize',14)
% set(groot,'defaultLineLineWidth',0.4)

if plotHilbert
    input_data_phase = angle(1i.*hilbert(clean_data((burst_start+1):end)));
    hold on
    plot([zeros(burst_start, 1); input_data_phase])
    legend('Input signal', 'Phase', 'Hilbert Transform')
    title("Tracking " + f + "Hz sampled at "+ 1/Ts + "Hz.")
else
    legend('Input signal', 'Processed Signal')
    title("Amplitude removal at sampling time of "+ 1/Ts + "Hz.")
end

end