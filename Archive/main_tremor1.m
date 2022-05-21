function [data] = main_tremor1()
% Process data without stimulation recorded with study participant at rest
cohort = [1 4 5 6 7 8 10 11 12 14];
dom_axis = [3, 3, 3, 3, 1, 1, 3, 1, 3, 1]; % 1 for x, 2 for y, 3 for z
% sample_length = 373229;
% data = zeros(length(cohort)*2*3, sample_length);
% % Data saves all processing outputs, tremor data in all axes and filtered
% % data for reference in the following order: raw x, filtered x, raw y,
% % filtered y, raw z, filtered z.

for iii = 1:length(cohort)
    
    load(strcat('data1/P0',num2str(cohort(iii)),'_baseline.mat'));
    
    data = SmrData.WvData;
    
    % Old sampling rate: 10417 Hz
    samplerateold = SmrData.SR;
    
    % Downsample data to 1 kHz
    samplerate = 1000;
    ts = timeseries(data, 0 : (1 / samplerateold):((size(data, 2)-1) / samplerateold));
    ts1 = resample(ts, 0 : (1 / samplerate):((size(data, 2)-1)/samplerateold), 'linear');
    all_data(1:size(ts1.data, 1), 1:size(ts1.data, 3)) = ts1.data;
    
    % Tremor data: channels 3, 5 and 6
    % X: sideways; Y: forward; Z: up and down
    tremor = all_data([3 5 6], :);
    tremorz = tremor(1, 1:end);
    tremory = tremor(2, 1:end);
    tremorx = tremor(3, 1:end);
    
    % Find the dominant tremor frequency for bandpass filtering
    for ax = 1 : 3
        [Pxx, F] = pwelch(tremor(ax, :), samplerate, [], samplerate, samplerate)
        frange = F(3 : 10);
        Pxxrange = Pxx(3 : 10);
        Freqpeak(ax, :) = frange(find(Pxxrange == max(Pxxrange)));
        Ppeak(ax, :) = max(Pxxrange);
        ps_curves(ax, :) = Pxx;
    end
    peak_ax = [(Freqpeak(find(Ppeak == max(Ppeak)))) (find(Ppeak == max(Ppeak)))];
    Fpeak_ns = peak_ax(1);
    
    % Set up bandpass filter +/- 2Hz around dominant frequency
    if (Fpeak_ns - 2) >= 1
        [b, a] = butter(2, [(Fpeak_ns - 2)/(0.5 * samplerate) (Fpeak_ns + 2)/(0.5 * samplerate)], 'bandpass'); %15
    else
        [b, a] = butter(2, [(1)/(0.5 * samplerate) (Fpeak_ns + 2)/(0.5 * samplerate)], 'bandpass'); %15
    end
    
    % Filter tremor data and apply accelerometer conversion factor
    tremorxf = filtfilt(b, a, tremorx) * 10 * 9.81 / 0.5;
    tremoryf = filtfilt(b, a, tremory) * 10 * 9.81 / 0.5;
    tremorzf = filtfilt(b, a, tremorz) * 10 * 9.81 / 0.5;
    
%     data((iii-1)*6+1, :) = tremorx;
%     data((iii-1)*6+2, :) = tremorxf;
%     data((iii-1)*6+3, :) = tremory;
%     data((iii-1)*6+4, :) = tremoryf;
%     data((iii-1)*6+5, :) = tremorz;
%     data((iii-1)*6+6, :) = tremorzf;
    
    [phase_model, artefact_model, input_data, tms_pulses] = tremor_kalman_test(tremorx);
    sim('Kalman_Modular_Simple');
    title_str = 'Sample '+ num2str(cohort(iii)) + ' x axis';
    tremorplot(out.thetas(1:phase_model.sim_dur), phase_model.f, phase_model.Ts, input_data(1:phase_model.sim_dur), 1, tremorxf(1:phase_model.sim_dur), title_str);
    
    [phase_model, artefact_model, input_data, tms_pulses] = tremor_kalman_test(tremorz);
    sim('Kalman_Modular_Simple');
    title_str = 'Sample '+ num2str(cohort(iii)) + ' z axis';
    tremorplot(out.thetas(1:phase_model.sim_dur), phase_model.f, phase_model.Ts, input_data(1:phase_model.sim_dur), 1, tremorzf(1:phase_model.sim_dur), title_str);
    
end


end