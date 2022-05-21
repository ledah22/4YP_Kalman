function [SNR_i, SNR_o, NF, CircErr, CircSD] = phaseplot(output_data, output_states, f, Ts, input_data, plotHilbert, clean_data, burst_start, reset_start)
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
set(gca,'FontSize',14)
set(groot,'defaultLineLineWidth',0.4)
%SNR_i = snr(clean_data, input_data-clean_data);
SNR_i = 1;
%SNR_o = snr(clean_data, output_states(1:(end-1), 1)-clean_data);
SNR_o = 0;
NF = SNR_i/SNR_o;
CircSD = 0;
CircErr = 0;

if plotHilbert
    input_data_analytical1 = hilbert(clean_data((burst_start+1):reset_start));
    try
        input_data_analytical2 = hilbert(clean_data((reset_start+1):end));
    catch
        input_data_analytical2 = [];
    end
    %input_data_analytical = [input_data_analytical1; input_data_analytical2];
    input_data_phase1 = angle(1i.*input_data_analytical1);
    input_data_phase2 = angle(1i.*input_data_analytical2);
    input_data_phase = [input_data_phase1, input_data_phase2]';
    input_data_phase = angle(1i.*hilbert(clean_data((burst_start+1):end)));

    %CircErr = mean(abs(input_data_phase-output_data(2:end, 1)))/(2*pi);
    %CircSD = sqrt(-2*log(abs(sum(exp(1i*(input_data_phase-output_data(2:end, 1))))/length(input_data_phase))));
    %SNR_o = snr(input_data_phase, output_data(2:end) - input_data_phase);
    %NF = SNR_i/SNR_o;
    hold on
    plot([zeros(burst_start, 1); input_data_phase])
    %plot(input_data_phase)
    legend('Input signal', 'Phase', 'Hilbert Transform')
    title("Tracking " + f + "Hz sampled at "+ 1/Ts + "Hz.")
else
    legend('Input signal', 'Processed Signal')
    title("Amplitude removal at sampling time of "+ 1/Ts + "Hz.")
end

end
% 
% figure
% for j = 1:num_rhythms
%     subplot(num_rhythms, 1, j);
%     plot(1:1:l, output_data(:, j)-input_data);
%     title("Estimation error")
% end