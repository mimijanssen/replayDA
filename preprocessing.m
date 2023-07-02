%% Fiber Photometry Pre-processing
% Author: Miriam Janssen
% Date: 04/08/2023
% Citations: The pipeline is inspired by Thomas Akam's code from https://github.com/ThomasAkam/photometry_preprocessing/blob/master/Photometry%20data%20preprocessing.ipynb

% Pre-processing steps: 
% 1. Filtering: to reduce noise and electrical artifacts 
% 2. Detrending: to correct for photobleaching 
% 3. (Omitted) Movement correction: to remove movement artifacts
% 4. Normalization: conversion to dF/F or Z-scoring or both

%% Load Data
clear; clc;
cd 'C:\Data\M437\M437-2023-06-16_track';  LoadExpKeys;
file_name = 'M437-2023-06-16';

addpath(genpath('Users\mimia\Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\mimia\Documents\GitHub\replay_DA\analysis\photometry'));

cfg.fc = {'CSC30.ncs'};
csc_photo = LoadCSC(cfg);

% extracts FP signal, time, and sampling rate
FP_data = [];
FP_data.acq.Fs = csc_photo.cfg.hdr{1}.SamplingFrequency;
FP_data.acq.time = csc_photo.tvec - csc_photo.tvec(1); 
FP_data.acq.FP{1} = csc_photo.data';

% rename variables
FP = (FP_data.acq.FP{1})*(1239/39); %1239/39 is the voltage multiplier
time = FP_data.acq.time; 
FS = FP_data.acq.Fs;

%% Raw Signal
% To plot raw data use FP and multiply by 31.77 (multiplier b/c of voltage
% divider)
sessionTitle = 'CW_';
last_time = length(time);
timerange1 = 10/0.0002; %datapoint range for 10 s
timerange2= 100/0.0002; %datapoint range for 100 s 
timerange3= length(time); % datapoint range for all data
time_ranges = [timerange1, timerange2, timerange3]; %in seconds 

figure(1)
for t_i = 1:length(time_ranges)
    t_range = 1:time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(time(t_range), FP(t_range), 'Color', [0 0.5 0])
    title([sessionTitle, num2str(time_ranges(t_i)),'samples'], 'Interpreter','none')
    ylabel('Fiber Signal (V)'); xlabel('Time (s)');
end

%% Downsampled to 1000 Hz 
% FS is currently 5000 and I'm changing it to 1000 Hz; consider 250 Hz
dsf = FS/1000;
FP = decimate(FP,dsf);
time = downsample(time,dsf);
FS = FS/dsf;

%% Extract Events
LoadExpKeys();
cfg_evt = [];
cfg_evt.eventList = ExpKeys.eventList;
cfg_evt.eventLabel = ExpKeys.eventLabel;
evt = LoadEvents(cfg_evt);

% Other Events if you have them
%sleep_sess1_starts = evt.t{1,5} - csc_photo.tvec(1); % initialze time of events
%sleep_sess1_ends = evt.t{1,6} - csc_photo.tvec(1); 
%sleep_sess2_starts = evt.t{1,7} - csc_photo.tvec(1);
%sleep_sess2_ends = evt.t{1,8} - csc_photo.tvec(1); 
%sleep1_starts = find(abs(time-sleep_sess1_starts) < 0.0005);
%sleep1_ends = find(abs(time-sleep_sess1_ends) < 0.0005);
%sleep2_starts = find(abs(time-sleep_sess2_starts) < 0.0005);
%sleep2_ends = find(abs(time-sleep_sess2_ends) < 0.0005);
%event_time = [time(sleep1_starts), time(sleep1_ends), time(sleep2_starts), time(sleep2_ends)];
%event_label = {'sleep session1 starts', 'sleep session1 ends','sleep session2 starts','sleep session2 ends'};

%% Denoised
% Why: to filter out electrical noise greater than 10 Hz 
%       "recording has large electrical noise artifacts, likely due to the
%       high gain amplifiers in the photodetectors picking up signals from
%       nearby mobile phone. The artifacts are very short pulses and can be
%       greatly reduced by running a median filter before the standard low
%       pass filter. 
% Method: Median filter & Low pass filter
% Note: Temporal dynamics of the biosensor are on the level of subseconds
% (X), so filtering out > 10 Hz signals should be ok. 

% Median filter: remove electrical artifacts 
FP_denoised = medfilt1(FP);

% check once
figure(2)
plot(time,FP);
title('Median Filtered FP Signal')
ylabel('Signal (V)')
xlabel('Time (s)')

% Butterworth Low pass filter 
% Note: Trisch lab used 20 Hz
fc = 20; % frequency
[b,a] = butter(2,fc/(FS/2)); % 2nd order
%freqz(b,a,[],FP_data.acq.Fs)
FP_denoised= filter(b,a,FP_denoised);
xlim([1 inf])

figure(3)
plot(time,FP_denoised)
title('20 Hz Butterworth Filtered FP Signal')
ylabel('Signal (V)')
xlabel('Time (s)')
xlim([1 inf])

figure(4)
plot(time,FP);
hold on
plot(time,FP_denoised)
hold off
title('20 Hz Butterworth Filtered FP Signal over Median Filtered FP Signal')
ylabel('Signal (V)')
xlabel('Time (s)')
legend('median','butterworth','Location','northeast')
xlim([1 inf]) 

%% Start Buffer
% Photobleaching is exponential and often greatest in the first few
% seconds of recording. This paper recommends removing 2-5 seconds from the
% beginning and end of the recording file 
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7853640/

% also removing the first 5 seconds : This is when I move the mouse from
% the pedestal to the track and might cause some major movement artifacts. 

% Start buffer
% remove first 3 second; number of samples to remove = 0.001;
FP_denoised(1:3/0.001) = []; 
time(1:3/0.001) = [];

%% Detrend 
% Why: to account for photobleaching 
% Method: fit an exponential decay to the data and subtract this
% exponential fit from the signal . 
% Adam used double exponential fit because of the multiple sources of
% fluorescence that contributed to the bleaching. (autofluorescence from
% fiber, brain tissue, and flurophore which may bleach at different rates) 

% minimize least square error with lsqcurvefit()
% input time, FP_denoised 
t = time;
y = FP_denoised;

% create a model
F = @(x,time)x(5)+x(1)*exp(-x(2)*time) + x(3)*exp(-x(4)*time);

% initial parameter guess
max_sig = max(FP_denoised);
x0 = [max_sig 0.005 max_sig/2 0.005 0];

% solve least squares
xunc = lsqcurvefit(F,x0,t,y);
F_expfit = F(xunc,time);

% subtract the fit from signal 
FP_detrended = FP_denoised - F_expfit;

figure(5)
plot(t,FP_denoised)
hold on
plot(t,F_expfit)
hold off
title('Exponential Fit to Filtered FP Signal')
ylabel('Signal (V)')
xlabel('Time (s)')
xlim([0 inf])

figure(6)
plot(t,FP_detrended)
title('Detrended and Filtered FP Signal')
ylabel('Signal (V)')
xlabel('Time (s)')
xlim([0 inf])

%% Motion Correction 
% see: https://github.com/ThomasAkam/photometry_preprocessing/blob/master/Photometry%20data%20preprocessing.ipynb

%% Normalization 
% Why: to combine data across sessions and/or subjects 
% "Different sessions may have different levels of fluorophores expression, 
%  excitation light, and autofluorescence" (Thomas Akam).
% Method: dF/F or Z-score

% F(t) - F0 / F0; 
dF = 100.*FP_detrended./F_expfit;

figure(7)
plot(t,dF)
title('dF Signal')
ylabel('Signal dF/F (%)')
xlabel('Time (s)')
xlim([5 inf])

% Z-score
% Alternatively, we can normalize by z-scoring each session 
% subtracting the mean and deviding by standard deviation. 
% x = current data
% avg_x = mean of the population ; or alternatively the median 
% stddev_x = standard deviation of the population 

% over whole session: 

F_zscored = (FP_detrended - mean(FP_detrended))./std(FP_detrended);
zdF = (dF - mean(dF))./std(dF);

figure(8)
plot(t,F_zscored)
title('Signal z-scored')
ylabel('Signal (z-scored)')
xlabel('Time (s)')
xlim([0 inf])

figure(9)
plot(t,zdF)
title('Signal dFz-scored')
ylabel('Signal (dFz-scored)')
xlabel('Time (s)')
xlim([0 inf])

%% Moving Window
% with moving windows: (no overlap)
window = 0.5 ; 
window_points = 0.5/0.001; % for 1,000 Hz.

F_zscored_moving = [];

for iter = 1:window_points:length(FP_detrended)
    if iter < length(FP_detrended)-window_points % less than FP_detrended window points, use end.... 
        F_zscored_moving = [F_zscored_moving ; FP_detrended(iter:iter+window_points-1,1) - mean(FP_detrended(iter:iter+window_points-1,1))./std(FP_detrended(iter:iter+window_points-1,1))];
    else 
        F_zscored_moving = [F_zscored_moving ; FP_detrended(iter:end,1) - mean(FP_detrended(iter:end,1))./std(FP_detrended(iter:end,1))];

    end
end

figure(9)
plot(t,F_zscored_moving)
title('Signal dFz-scored moving average')
ylabel('Signal (dFz-scored)')
xlabel('Time (s)')
xlim([0 inf])

%% Save Variables 
filename = append(file_name, "processed.mat");
data.t = t;
data.F_detrend = FP_detrended;
data.FP_z = F_zscored;
data.dF = dF; 
data.zdF = zdF;
data.zDF_mov = F_zscored_moving;
data.window = window; %seconds 
data.evtt = event_time;
data.evtlabel = event_label;
save(filename, '-struct','data')
disp('saved')