function [] = tremorplot(output_data, f, Ts, input_data, plotHilbert, clean_data, title_str, input_on)
[l, num_rhythms] = size(output_data);
legend_entries = strings(1, 1+plotHilbert+input_on);
legend_entries(1, 1) = "Kalman phase";

figure
for j = 1:num_rhythms
    subplot(num_rhythms, 1, j);
    plot(1:1:l, output_data(:, j));
end
set(gca,'FontSize',14)
set(groot,'defaultLineLineWidth',0.00001)
xlabel('Time (num. of samples)')
ylabel('Amplitude (Radians and Volts)')
title(title_str)

if input_on
    hold on
    plot(input_data)
    legend_entries(1, 2) = "Input data";
end

if plotHilbert
    input_data_analytical = hilbert(clean_data(1:length(output_data)));
    input_data_phase = angle(i*input_data_analytical);
    hold on
    plot(input_data_phase)
    legend_entries(1, 2+input_on) = "Hilbert phase";
end

legend(legend_entries')

end
% 
% figure
% for j = 1:num_rhythms
%     subplot(num_rhythms, 1, j);
%     plot(1:1:l, output_data(:, j)-input_data);
%     title("Estimation error")
% end