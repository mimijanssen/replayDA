%% Fiber Photometry Pre-processing for Drug Sessions
% Input: CSC file with your fiber photometry data
% Output: TSD structure with your processed fiber photometry data 

% Pre-processing steps: 
% 1. Filtering: to reduce noise and electrical artifacts 
% 2. Detrending: to correct for photobleaching - using a line fit through
% the baseline period 
% 3. Normalization: conversion to dF/F or Z-scoring 
% 4. Structuring: save as a structure compatable with the van der Meer lab code base
% optional sanity check: save as a structure compatable with other code bases to verify
% your preprocessing results in the same output (here I used the Tritsch
% lab code base ✓ ).

% written by Mimi Janssen, 4/08/2023; edit 11/15/2024
% Citations: 
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
cd 'F:\M533\M533_2024_08_29_amp2';  LoadExpKeys; %LoadMetadata;
file_name = 'M533_2024_08_29';

addpath(genpath('Users\mimia\Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\mimia\Documents\GitHub\replay_DA\analysis\photometry'));
addpath(genpath('C:\Users\mimia\Documents\GitHub\std_shade'));

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
FP.data = (FP_data.acq.FP{1})*(1239/39); % voltage multiplier
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
dsf = FP.cfg.hdr{1}.SamplingFrequency/1000; % This is hard coded so make sure you change this if you collect at a different sampling rate. Sorry will fix later.
FP.data = decimate(FP.data,dsf);
FP.tvec = downsample(FP.tvec,dsf);
FP.cfg.hdr{1}.SamplingFrequency = FP.cfg.hdr{1}.SamplingFrequency/dsf;

% note: decimate uses filtfilt to filter "in both directions to make sure 
% the filtered data has zero phase. Make a data vector properly prepended
% and appended to filter forwards and back so end effects can be
% obliterated."

%% Extract Events
LoadExpKeys();
cfg_evt = [];
cfg_evt.eventList = ExpKeys.eventList;
cfg_evt.eventLabel = ExpKeys.eventLabel;
evt = LoadEvents(cfg_evt);

% Events
pedestal_starts = evt.t{1,5} - csc_photo.tvec(1); 
drug_admin = evt.t{1,7} - csc_photo.tvec(1); % initialzied time of drug administration
admin_done = evt.t{1,8} - csc_photo.tvec(1);
pedestal2_ends = evt.t{1,10} - csc_photo.tvec(1); 

% tolerance -- change to nearest
ind_ped_starts = nearest_idx(pedestal_starts,FP.tvec); %find(abs(time-pedestal_starts) < 0.0005);
ind_drug = nearest_idx(drug_admin,FP.tvec); %find(abs(time-drug_admin) < 0.0005);
ind_admin = nearest_idx(admin_done,FP.tvec); %find(abs(time-admin_done) < 0.0005);
ind_ped_ends = nearest_idx(pedestal2_ends,FP.tvec); %find(abs(time-pedestal2_ends) < 0.0005);
event_time = [FP.tvec(ind_ped_starts), FP.tvec(ind_drug), FP.tvec(ind_admin), FP.tvec(ind_ped_ends)]; % idk why i did it like this
event_label = {'pedestal session starts', 'drug admin','admin done','pedestal session2 ends'};

%% ✧･ﾟ: *✧･ﾟ:*     STEP 1: FILTERING　　 *:･ﾟ✧*:･ﾟ✧

% Median filter: remove electrical artifacts 

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
% pwelch of filtered data
[Pxx_butter,F_butter] =  pwelch(FP.FP_denoised,hanning(wsize),wsize/2,[],FP.cfg.hdr{1}.SamplingFrequency);
figure(5)
plot(F_butter,10*log10(Pxx_butter),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 100]);
title('PSD for Denoised Fiber')

%% ✧･ﾟ: *✧･ﾟ:*     STEP 2: DETREND　　 *:･ﾟ✧*:･ﾟ✧

% Detrend - line fit on baseline period

% time from pedestal start to drug admin. 

t = FP.tvec(ind_ped_starts:ind_drug);
y = FP_denoised(ind_ped_starts:ind_drug);

% FOR LINEAR
% create a model
%F = @(x,t)x(2)+x(1)*t;

% initial parameter guess
%max_sig = max(y);
%x0 = [max_sig/2 0];

% FOR DOUBLE EXP 
% create a model
F = @(x,t)x(5)+x(1)*exp(-x(2)*t) + x(3)*exp(-x(4)*t);

% initial parameter guess
max_sig = max(y);
x0 = [max_sig 0.005 max_sig/2 0.005 0];

% solve least squares
xunc = lsqcurvefit(F,x0,t,y); % mode fit to baseline data
F_expfit = F(xunc,FP.tvec); % now fit to whole time

% subtract the fit from signal 
FP_detrended = FP_denoised - F_expfit;
% used for further processing 

FP.FP_detrend_base = FP_detrended; % same in the tsd structure.

figure(5)
plot(FP.tvec,FP_denoised)
hold on
plot(FP.tvec,F_expfit)
hold off
title('Exponential Fit to Filtered FP Signal')
ylabel('Signal (V)')
xlabel('Time (s)')
xlim([5 inf])

figure(6)
plot(FP.tvec,FP_detrended)
title('Detrended and Filtered FP Signal')
ylabel('Signal (V)')
xlabel('Time (s)')
xlim([5 inf])

%% ✧･ﾟ: *✧･ﾟ:*     STEP 3: NORMALIZATION　　 *:･ﾟ✧*:･ﾟ✧
% Method: dF/F or Z-score

% dF = (FP-baseline)./baseline;
dF_base = 100.*FP_detrended./F_expfit; % delta F
% This is the same as doing: dF = (FP-baseline)./baseline;;

% Z-score
F_zscored = (FP_detrended - mean(FP_detrended))./std(FP_detrended); %just detrended and z-scored
zdF_base = (dF_base - mean(dF_base))./std(dF_base); % delta F , z-scored 

FP.dF_base = dF_base;  % dF/F 
FP.F_zscored_base = F_zscored; % detrended and z-scored
FP.zdF_base = zdF_base; % dF/F and z-scored 

figure(7)
plot(FP.tvec,F_zscored)
title('dF/F (z-scored)')
ylabel('dF/F (z-scored')
xlabel('Time (s)')

% Note to self: use other axis to plot saline data~ 

%% ✧･ﾟ: *✧･ﾟ:*     STEP 4: SAVING　　 *:･ﾟ✧*:･ﾟ✧
% Save TSD Structure - VLAB
filename = append(file_name, "processed.mat");
save(filename, '-struct','FP')
