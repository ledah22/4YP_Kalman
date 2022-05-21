function [SNR_i, SNR_o, NF, CircErr, CircSD] = combinedplot(output_data, output_states, f, Ts, input_data, plotHilbert, clean_data, burst_start, no_tms_clean_data)
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
SNR_i = snr(no_tms_clean_data(:, (burst_start+1):end)', clean_data(burst_start+1:end)-no_tms_clean_data(:, (burst_start+1):end)');
%SNR_o = snr(no_tms_clean_data(:, (burst_start+1):end)', clean_states(burst_start+1:(end-1), 1)-no_tms_clean_data(:, (burst_start+1):end)');
SNR_o = 1;
NF = SNR_i/SNR_o;
CircSD = 0;
CircErr = 0;

if plotHilbert
    input_data_analytical = hilbert(no_tms_clean_data((burst_start+1):end)');
    input_data_phase = angle(1i*input_data_analytical);

    CircErr = mean(abs(input_data_phase-output_data((burst_start+2):end, 1)))/(2*pi);
    CircSD = sqrt(-2*log(abs(mean(exp(1i*(input_data_phase-output_data((burst_start+2):end, 1)))))));
    %SNR_o = snr(input_data_phase, output_data((burst_start+2):end) - input_data_phase);
    NF = SNR_i/SNR_o;
    hold on
    plot([pi*eye(burst_start, 1); input_data_phase])
    legend('Input signal', 'Phase', 'Hilbert Transform')
    title("Tracking " + f + "Hz sampled at "+ 1/Ts + "Hz.")
else
    legend('Input signal', 'Processed Signal')
    title("Amplitude removal at sampling time of "+ 1/Ts + "Hz.")
end
% 
% figure
% for j = 1:num_rhythms
%     subplot(num_rhythms, 1, j);
%     plot(1:1:l, output_data(:, j)-input_data);
%     title("Estimation error")
% end