%% Fiber Photometry Pre-processing
% Input: CSC file with your fiber photometry data
% Output: TSD structure with your processed fiber photometry data 

% Pre-processing steps: 
% 1. Filtering: to reduce noise and electrical artifacts 
% 2. Detrending: to correct for photobleaching 
% 3. Normalization: conversion to dF/F or Z-scoring or both
% 4. Structuring: save as a structure compatable with the van der Meer lab code base
% optional sanity check: save as a structure compatable with other code bases to verify
% your preprocessing results in the same output (here I used the Tritsch
% lab code base ✓ ).

% written by Mimi Janssen, 4/08/2023
% Citations: 
% This code draws heavily from Thomas Akams photometry preprocessing code.
% Thomas Akam: https://github.com/ThomasAkam/photometry_preprocessing/blob/master/Photometry%20data%20preprocessing.ipynb
% Tritsch lab code base: https://github.com/TritschLab/TLab_Toolbox
% van der Meer lab code base: https://github.com/mimijanssen/vandermeerlab-replay-da

%% Add vanderMeer Code Base to Your Path 
% restoredefaultpath; clear classes; % start with a clean slate
% cd('C:\Users\mimia\Documents\Toolboxes\vandermeerlab-replay-da\code-matlab\shared'); % or, wherever your code is located -- NOTE \shared subfolder!
% p = genpath(pwd); % create list of all folders from here
% addpath(p);
 
%% Load Data
clear; clc;

cd 'D:\M556\M556_2025_03_06_recording6' % path with your csc fiber data file
file_name = 'M556_2025_03_06'; % file name that your processed data will be saved as

cfg.fc = {'CSC30.ncs'};
csc_photo = LoadCSC(cfg);

% extracts FP signal, time, and sampling rate
FP_data = [];
FP_data.acq.Fs = csc_photo.cfg.hdr{1}.SamplingFrequency; % set FP_data.acq.Fs to sampling frequency rate (5000 points per second) 
FP_data.acq.time = csc_photo.tvec - csc_photo.tvec(1); % set FP_data.acq.Fs time to the time vector subtracted by the first time point
% this initializes the time vector to start at 0 then 2.0 x 10 ^-4 or 0.0002 and 4.0 x
%x 10 ^ -4 0.0004 ... 0.0002 second time interval means 5,000 points per
%second
FP_data.acq.FP{1} = csc_photo.data';

% rename variables and store in a tsd struct
FP = tsd;
FP.data = (FP_data.acq.FP{1})*(1239/39); %voltage multiplier
FP.tvec = FP_data.acq.time; 
FP.cfg.hdr{1}.SamplingFrequency = FP_data.acq.Fs;

%% Raw Signal
% To plot raw data use FP and multiply by 31.77 (multiplier b/c of voltage
% divider)
sessionTitle = 'CW_';
last_time = length(FP.tvec); %the value of the last time point is how many seconds the recording was
timerange1 = 10*FP.cfg.hdr{1}.SamplingFrequency; 
timerange2= 100*FP.cfg.hdr{1}.SamplingFrequency; 
timerange3= length(FP.tvec); 
time_ranges = [timerange1, timerange2, timerange3];

figure(1)
for t_i = 1:length(time_ranges)
    t_range = 1:time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(FP.tvec(t_range), FP.data(t_range), 'Color', [0 0.5 0])
    %title([num2str(time_ranges(t_i)),' samples'], 'Interpreter','none')
    ylabel('Fiber Signal (V)'); xlabel('Time (s)');
end

%% Downsampled to 1000 Hz 
% FS is currently 5000 and I'm changing it to 1000 Hz;
% The GFP mice I sampled at 1,600 Hz and I'm going to keep it that way -
% not decimate... 

%dsf = FP.cfg.hdr{1}.SamplingFrequency/1000; 
%FP.data = decimate(FP.data,dsf);
%FP.tvec = downsample(FP.tvec,dsf);
%FP.cfg.hdr{1}.SamplingFrequency = FP.cfg.hdr{1}.SamplingFrequency/dsf;

% note: decimate uses filtfilt to filter "in both directions to make sure 
% the filtered data has zero phase. Make a data vector properly prepended
% and appended to filter forwards and back so end effects can be
% obliterated."

%% Extract Events
% Extracts photobeam events.
LoadExpKeys();
cfg_evt = [];
cfg_evt.eventList = ExpKeys.eventList;
cfg_evt.eventLabel = ExpKeys.eventLabel;
evt = LoadEvents(cfg_evt);

%% ✧･ﾟ: *✧･ﾟ:*     STEP 1: FILTERING　　 *:･ﾟ✧*:･ﾟ✧

%% Denoised
% Why: to filter out electrical noise greater than 20 Hz 
%       "recording has large electrical noise artifacts, likely due to the
%       high gain amplifiers in the photodetectors picking up signals from
%       nearby mobile phone. The artifacts are very short pulses and can be
%       greatly reduced by running a median filter before the standard low
%       pass filter. 
% Method: Median filter and butterworth Low pass filter
% Note: Temporal dynamics of the biosensor are on the level of subseconds
% (X), so filtering out > 20 Hz signals should be ok. 

% It is important to use filtfilt as this is a way to try and preserve any
% phase relationships. Forward filtering methods may shift the phase 
% of your signal and obscure relationships between the LFP and behavioral 
% or neural events.(theta phase precession, or cross-frequency coupling, 
% will be affected). For neural data, where even small phase 
% shifts can be problematic, we therefore take an alternative approach: we 
% filter the signal forwards and backwards, such that the net phase 
% response is zero: no phase shift! This is accomplished using the 
% ''filtfilt()'' function. 

% Medfilt1 has no phase shift because it computes the median in a window
% centered on the current sample.

% Median filter
FP_denoised = medfilt1(FP.data);

figure(2)
plot(FP.tvec,FP.data);
hold on
plot(FP.tvec,FP_denoised);
title('Median Filtered FP Signal')
ylabel('Signal (V)')
xlabel('Time (s)')

% Butterworth filter
fc = 20; % cut off frequency
fs = FP.cfg.hdr{1}.SamplingFrequency; % sampling rate 
[b,a] = butter(8,fc/(fs/2)); % 8th order filter
FP.FP_denoised = filtfilt(b, a, FP_denoised);

figure(3)
plot(FP.tvec,FP.data);
hold on
plot(FP.tvec,FP.FP_denoised)
hold off
title('20 Hz Butterworth Filtered FP Signal over Raw Data')
ylabel('Signal (V)')
xlabel('Time (s)')
legend('raw','butterworth','Location','northeast')

% pwelch to check if butterworth worked: 
% pwelch of raw data
wsize = FP.cfg.hdr{1}.SamplingFrequency * 4; % sampling frequency * 4s is the amount of samples (the window size needed)
[Pxx,F] = pwelch(FP.data,hanning(wsize),wsize/2,[],FP.cfg.hdr{1}.SamplingFrequency);
figure(4)
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 100]);
title('PSD for Raw Fiber')

[Pxx_butter,F_butter] =  pwelch(FP.FP_denoised,hanning(wsize),wsize/2,[],FP.cfg.hdr{1}.SamplingFrequency);
figure(5)
plot(F_butter,10*log10(Pxx_butter),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 100]);
title('PSD for Denoised Fiber')
% pwelch of butterworth (20Hz noise should be cut) 

%% ✧･ﾟ: *✧･ﾟ:*     STEP 2: DETREND　　 *:･ﾟ✧*:･ﾟ✧

%% Global Detrend: (double exponential) for 3 windows
% Why: to account for photobleaching 
% Method: fit an exponential decay to the data and subtract this
% exponential fit from the signal . 
% Adam used double exponential fit because of the multiple sources of
% fluorescence that contributed to the bleaching. (autofluorescence from
% fiber, brain tissue, and flurophore which may bleach at different rates)

% I used detrended in this way for three windows since my data has three
% distinct parts where the animal was picked up and moved, causing the
% fiber patch to also shift. 
 
% time that pre task period ends; time that post task period starts
pre = ExpKeys.prerecord(2) - csc_photo.tvec(1); % time of pre sleep period ends, initialized  
post = ExpKeys.postrecord(1) - csc_photo.tvec(1); % time of post sleep period, initialized  

pre_idx = nearest_idx3(pre,FP.tvec);
post_idx = nearest_idx3(post,FP.tvec);

windows = [1,pre_idx,post_idx]; % determines how many samples there are in each window (s)
% windows chosen: 
% start to sleep session1 ends 
% sleep session1 ends to sleep session 2 starts 
% sleep session 2 ends to end.

FP_detrend_win = [];
F_expfit_win =[];

t = FP.tvec;
y = FP.FP_denoised;

for iter_win = 1:1:3 % iterate through 3 time periods
    F = @(x,t)x(5)+x(1)*exp(-x(2)*t) + x(3)*exp(-x(4)*t);

    if iter_win == 1
        max_sig = max(FP.FP_denoised(1:windows(2)));
        x0 = [max_sig 0.005 max_sig/2 0.005 0];
    elseif iter_win ==2 
        max_sig = max(FP.FP_denoised(windows(2)+1:windows(3)));
        x0 = [max_sig 0.005 max_sig/2 0.005 0];
    else 
        max_sig = max(FP.FP_denoised(windows(3)+1:end));
        x0 = [max_sig 0.005 max_sig/2 0.005 0];
    end

    if iter_win == 1
        xunc = lsqcurvefit(F,x0,t(1:windows(2)),y(1:windows(2)));
        F_expfit = F(xunc,t(1:windows(2)));
        % subtract the fit from signal
        FP_detrended_win = FP.FP_denoised(1:windows(2)) - F_expfit;
        % used for further processing
        F_expfit_win = [F_expfit_win F_expfit'];
        FP_detrend_win = [FP_detrend_win FP_detrended_win'];
    elseif iter_win == 2 
        xunc = lsqcurvefit(F,x0,t(windows(2)+1:windows(3)),y(windows(2)+1:windows(3)));
        F_expfit = F(xunc,t(windows(2)+1:windows(3)));
        FP_detrended_win = FP.FP_denoised(windows(2)+1:windows(3)) - F_expfit;
        F_expfit_win = [F_expfit_win F_expfit'];
        FP_detrend_win = [FP_detrend_win FP_detrended_win'];
    else
        xunc = lsqcurvefit(F,x0,t(windows(3)+1:end),y(windows(3)+1:end));
        F_expfit = F(xunc,t(windows(3)+1:end));
        FP_detrended_win = FP.FP_denoised(windows(3)+1:end) - F_expfit;
        F_expfit_win = [F_expfit_win F_expfit'];
        FP_detrend_win = [FP_detrend_win FP_detrended_win'];
    end
end

FP.FP_detrend_win = FP_detrend_win; 

figure(6)
plot(FP.tvec,FP.FP_detrend_win);
hold on
plot(FP.tvec,F_expfit_win)
hold off
title('3 Window Double Exponential Detrend')
ylabel('Signal (V)')
xlabel('Time (s)')
legend('detrended signal','expfit','Location','northeast')

%% ✧･ﾟ: *✧･ﾟ:*     STEP 3: NORMALIZATION　　 *:･ﾟ✧*:･ﾟ✧

%% Normalization for Windowed Detrend Based on Task (double exponential)
% dF = (FP-baseline)./baseline;
dF_win = 100.*FP_detrend_win./F_expfit_win; % delta F/F
% This is the same as doing: dF = (FP-baseline)./baseline;
% dF_win2 = 100.*(FP_denoised' - F_expfit_win )./F_expfit_win;

% Z-score
% Alternatively, we can normalize by z-scoring each session 
% subtracting the mean and dividing by standard deviation.

%F_zscored_win = (FP_detrend_win - mean(FP_detrend_win))./std(FP_detrend_win); %just detrended and z-scored
%zdF_win = (dF_win - mean(dF_win))./std(dF_win); % delta F , z-scored 
F_zscored_win = zscore(FP_detrend_win);
zdF_win = zscore(dF_win);

FP.dF_win = dF_win;  % dF/F 
FP.F_zscored_win = F_zscored_win; % detrended and z-scored
FP.zdF_win = zdF_win; % dF/F, detrended, and z-scored 

figure(7)
plot(FP.tvec,FP.zdF_win)
hold on
plot(FP.tvec,FP.dF_win) % transients seem to be larger here 
title('dF/F (z-scored)')
ylabel('dF/F (z-scored')
xlabel('Time (s)')

%% Local Detrend: Windowed Detrend using Locdetrend
% another way to detrend the signal is a simple linear detrend. 
addpath('C:\Users\mimia\Documents\GitHub\vandermeerlab-replay-da\chronux-master\spectral_analysis\helper')
addpath('C:\Users\mimia\Documents\GitHub\vandermeerlab-replay-da\chronux-master\spectral_analysis\continuous')

% 1 min detrend, 1 sample stepsize  
FP_detrend_60s = locdetrend(FP.FP_denoised,FP.cfg.hdr{1}.SamplingFrequency,[60 1]);
% 10 s detrend, 1 sample stepsize 
FP_detrend_10s = locdetrend(FP.FP_denoised,FP.cfg.hdr{1}.SamplingFrequency,[10 1]); 

% 1 min detrend, almost no overlap
% if your step size is the size of your window, then you will have no
% overlap
FP_detrend_60s_no = locdetrend(FP.FP_denoised,FP.cfg.hdr{1}.SamplingFrequency,[61 60]);
% There is slight overlap on the edges since I had a edge artifact when I
% did 60, 60. 

FP.detrend_60s = FP_detrend_60s; 
FP.detrend_10s = FP_detrend_10s;
FP.detrend_60s_no = FP_detrend_60s; 

% detrended and filtered signal (using 10s window for detrend)
figure(8)
plot(FP.tvec, FP.FP_denoised)
hold on
plot(FP.tvec,FP_detrend_60s_no)
title('Detrended and Filtered FP Signal (60s)')
ylabel('Signal (V)')
xlabel('Time (s)')

%% Normalization for Windowed Detrend (locdetrend)
% z-scored locdetrended signal
%zF_win_60s_no = (FP_detrend_60s_no - mean(FP_detrend_60s_no))./std(FP_detrend_60s_no);
%zF_win_60s = (FP_detrend_60s - mean(FP_detrend_60s))./std(FP_detrend_60s); % you don't need to apply the . becase it should find the standard deviation for everything... try just doing zscore... and see if you get a different output
%zF_win_10s = (FP_detrend_10s - mean(FP_detrend_10s))./std(FP_detrend_10s);

zF_win_60s_no = zscore(FP_detrend_60s_no);
zF_win_60s = zscore(FP_detrend_60s);
zF_win_10s = zscore(FP_detrend_10s);

FP.zF_win_60s = zF_win_60s; 
FP.zF_win_10s = zF_win_10s;
FP.zF_win_60s_no = zF_win_60s_no; 

% filtered, detrended, and normalized signal (using 60s window for detrend)
figure(9)
subplot(2,1,1)
plot(t,zF_win_60s)
title('Filtered, Detrended, and Normalized FP Signal (60s)')
ylabel('Signal z-scored (V)')
xlabel('Time (s)')
subplot(2,1,2)
plot(t,zF_win_60s_no)
title('Filtered, Detrended, Normalized FP Signal (60s, almost no overlap)')
ylabel('Signal z-scored (V)')
xlabel('Time (s)')
 
%% Preprocess using Tritsch code base 

% Parameters
% params.dsRate = 100; % Downsampling rate if you want to downsample the signal
% %This dsRate will also be applied to all signals during the analysis
% %pipeline
% 
% % Filter Parameters
% params.FP.lpCut = 20; % Cut-off frequency for filter
% params.FP.filtOrder = 8; % Order of the filter
% 
% % Baseline Parameters
% params.FP.basePrc = 5; % Percentile value from 1 - 100 to use when finding baseline points
% %Note: Lower percentiles are used because the mean of signal is not true
% %baseline
% params.FP.winSize = 10; % Window size for baselining in seconds
% params.FP.winOv = 0; %Window overlap size in seconds
% params.FP.interpType = 'linear'; % 'linear' 'spline' 
% params.FP.fitType = 'interp'; % Fit method 'interp' , 'exp' , 'line'
% 
% % Demodulation Parameters
% %When demodulating signals, the filter creates edge artifacts. We record
% %for a few seconds longer, so we can remove x seconds from the beginning
% %and end
% %Adjust the variable to "0" if it's a normal photometry recording
% params.FP.sigEdge = 15; %15; %Time in seconds of data to be removed from beginning and end of signal
%params.FP.modFreq = [319 217];

% Preprocess CW mode
%addpath('C:\Users\mimia\Documents\Replay-DA\Analysis\const')

%FP_data_Tlab = processFP(params, FP_data);
%FP_Tlab = tsd(FP_data_Tlab.final.time', FP_data_Tlab.final.FP{1}', 'FP');
% for M534 rec 2 or M547 rec 1... this did not work... need to reshape
% something

% Save TSD Structure- TLAB
% filename = append(file_name, "Tlabprocessed.mat");
% save(filename, '-struct','FP_Tlab')

%% ✧･ﾟ: *✧･ﾟ:*     STEP 4: SAVING　　 *:･ﾟ✧*:･ﾟ✧

%% Save TSD Structure - VLAB
filename = append(file_name, "processed.mat");
save(filename, '-struct','FP')
