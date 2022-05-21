function combined_data = combine_data()

addpath '/Users/isidoraradenkovic/Desktop/Education/4YP/4YP_Kalman/data1'
% This will be different on other devices. However, this is just an example
% of how to extract the tremor data for processing. To test on
% accelerometer data, plase obtain data from the respective authors first.

%% Adaptation of the script created by Beatriz Silveira de Arruda

% Process data without stimulation recorded with study participant at rest
cohort = [1 4 5 6 7 8 10 11 12 14];
% cohort = [8, 14]; %To test on smaller batches, less memory exhausting
% dom_axis = [3, 3, 3, 3, 1, 1, 3, 1, 3, 1]; % Dominant tremor axis, denote
% 1 for x, 2 for y, 3 for z
% data1 = []; use to only focus on the other dataset

for iii = 1:length(cohort)
    
    load(strcat('P0',num2str(cohort(iii)),'_baseline.mat'))
    
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
    
    [d, c] = butter(2, 1/(0.5 * samplerate), 'high');
    tremorx = filtfilt(d, c, tremorx) * 10 * 9.81 / 0.5;
    tremory = filtfilt(d, c, tremory) * 10 * 9.81 / 0.5;
    tremorz = filtfilt(d, c, tremorz) * 10 * 9.81 / 0.5;
    
    data_temp.name = "Baseline " + num2str(cohort(iii));
    data_temp.x = tremorx;
    data_temp.xf = tremorxf;
    data_temp.y = tremory;
    data_temp.yf = tremoryf;
    data_temp.z = tremorz;
    data_temp.zf = tremorzf;
    
    data1(iii) = data_temp;
    
end

clearvars data_temp;
%% Adaptation of the script written by Carolina Reis
%%% data files are names SmrData where:
% SR: sampling rate
% WvTits: legends for data rows in WvData
% Signal: acc data Z axis; acc2: acc data Y axis; acc3: acc data X axis;
% FltSig: "Signal" channel band-passed filtered 2-8Hz
% Ignore Trg and codes (stimualtion related signals)

%%% indices of begining of posture (hu) and end of posture (hd)

hu={[160001 358233 518806 660400 817804];...
    [13119 110602 229200 349501 469200];...
    [129703 233503 427738 591202 702603];...
    [558301 775910 933460 1062297 1178351];...
    [21159 147999 288926 438422 568188];...
    [5926 146868 278294 425572 552266];...
    [104219 288259];...
    [7242 141494 279017 403458 555862];...
    [18225 163365 390458 585048 735348];...
    [71526 263629 444666]};

hd={[242500 447100 603892 750900 902099];...
    [49229 175200 289095 416900 533019];...
    [165400 303387 510893 645594 778100];...
    [638400 844595 996600 1131877 1239913];...
    [80208 209078 353536 495625 628761];...
    [72343 215572 339246 487600 618855];...
    [208160 394933];...
    [81747 211850 342088 476237 621414];...
    [75565 250362 502353 648016 812846];...
    [126495 333852 504802]};

cohort = [2 3 4 5 8 10 11 13 16 17];
% cohort = [5];

for iii = 1:length(cohort)
    load(strcat('P0',num2str(cohort(iii)),'_NS.mat'))
    sample = SmrData.WvData;
    tremorf = sample(1, :);
    
    [d, c] = butter(2, 2.5/(0.5 * SmrData.SR), 'high');
    tremor = filtfilt(d, c, sample(3, :)) * 10 * 9.81 / 0.5;
    tremorf = filtfilt(d, c, sample(1, :)) * 10 * 9.81 / 0.5;
    
    [f, e] = butter(2, 10/(0.5 * SmrData.SR), 'low');
    tremorf = filtfilt(f, e, tremorf) * 10 * 9.81 / 0.5;
    
    row_u = cell2mat(hu(iii));
    row_d = cell2mat(hd(iii));
    for ii = 1:length(row_u)
        data_temp.name = "NS " + num2str(cohort(iii))+" from " + num2str(row_u(ii)) + " to " + num2str(row_d(ii));
        data_temp.SR = SmrData.SR;
        data_temp.x = tremor(row_u(ii):row_d(ii));
        data_temp.xf = tremorf(row_u(ii):row_d(ii));
        
        % Downsampling
        samplerateold = data_temp.SR;
        samplerate = 1000;
        ts = timeseries(data_temp.x, 0 : (1 / samplerateold):((size(data_temp.x, 2)-1) / samplerateold));
        ts1 = resample(ts, 0 : (1 / samplerate):((size(data_temp.x, 2)-1)/samplerateold), 'linear');
        l = length(ts1.data);
        data_temp.x = ts1.data(1, 1, 1:l);
        
        ts = timeseries(data_temp.xf, 0 : (1 / samplerateold):((size(data_temp.xf, 2)-1) / samplerateold));
        ts2 = resample(ts, 0 : (1 / samplerate):((size(data_temp.xf, 2)-1)/samplerateold), 'linear');
        l = length(ts2.data);
        data_temp.xf = ts2.data(1, 1, :);
        data_temp.SR = 1000;
        
        data2(iii, ii) = data_temp;
    end
end

%% Assemble data
combined_data.data1 = data1;
combined_data.data2 = data2;
end