%% Plot Raw Fiber Data
% example: cd = C:\Data\M406\2022-10-30_M406_L_sleep_track_train
%folder = 'C:\Data\M406\2022-10-30_M406_L_sleep_track_train' ;
%cd folder 

%% Path
addpath(genpath('Users\mimia\Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\mimia\Documents\GitHub\replay_DA\analysis\photometry'));
addpath(genpath('C:\Users\mimia\Documents\replay-DA\analysis'));

%% Load Data
cfg.fc = {'CSC30.ncs'};
csc_photo = LoadCSC(cfg);

%% Raw Data
FP_data = []; % initializing fiber photometry data
FP_data.acq.Fs = csc_photo.cfg.hdr{1}.SamplingFrequency; %set FP_data.acq.Fs to sampling frequency rate (5000 points per second) 
% is this what the sampling rate does? 
FP_data.acq.time = csc_photo.tvec - csc_photo.tvec(1); %set FP_data.acq.Fs time to the time vector subtracted by the first time point
%this initializes the time vector to start at 0 then 2.0 x 10 ^-4 or 0.0002 and 4.0 x
%x 10 ^ -4 0.0004 ... 0.0002 second time interval means 5,000 points per
%second
FP_data.acq.FP{1} = csc_photo.data'; %fiber data 

%rename for plot
FP = FP_data.acq.FP{1};
time = FP_data.acq.time; 

%% Plot Data
sessionTitle = 'CW_';
last_time = length(time); %the value of the last time point is how many seconds the recording was
timerange1 = 10/0.0002; %datapoint range for 10 s
timerange2= 100/0.0002; %datapoint range for 100 s 
timerange3= length(time); % datapoint range for all data
time_ranges = [timerange1, timerange2, timerange3]; %in seconds 
%Fs = FP_data.final.Fs;

for t_i = 1:length(time_ranges)
    t_range = 1:time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(time(t_range), FP(t_range), 'Color', [0 0.5 0])
    title([sessionTitle, num2str(time_ranges(t_i))], 'Interpreter','none')
    ylabel('Fiber Signal (V)'); xlabel('Time (s)');
end

%% Detrend data 
%parameters
params.FP.winSize = 10; % 10 second window size for baselining in seconds
params.FP.winOv = 0; % 0 window overlap size in seconds
params.FP.fitType = 'line'; % fit a line 
params.FP.interpType = 'linear'; %linear interpolation
params.FP.dsRate = 1; % downsample rate 
params.FP.basePrc = 5; % Percentile value from 1 - 100 to use when finding baseline points
%Note: Lower percentiles are used because the mean of signal is not true
%baseline

rawFS = FP_data.acq.Fs;
Fs = rawFS/params.FP.dsRate;

[dF,baseline] = detrendFP(FP,params.FP.interpType,params.FP.fitType,params.FP.basePrc,params.FP.winSize,params.FP.winOv,Fs);


%% Plot dF Data
%rename for plot
dF_FP = dF; %scaled, normalized, and detrended
dF_time = time; 
FP_norm = (FP-baseline(1))./baseline(1); %normalized
FP_scaled = FP_norm*100;

sessionTitle = 'Samples';
last_time = length(time); %the value of the last time point is how many seconds the recording was
timerange1 = 10/0.0002; %datapoint range for 10 s (5,000 points per second)
timerange2= 100/0.0002; %datapoint range for 100 s 
timerange3= length(time); % datapoint range for all data
time_ranges = [timerange1, timerange2, timerange3]; %in seconds 
%Fs = FP_data.final.Fs;

for t_i = 1:length(time_ranges)
    t_range = 1:time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(time(t_range), FP_scaled(t_range), 'Color', [0 0.7 0])
    hold on
    plot(time(t_range), dF_FP(t_range), 'Color', [0 0.5 0])
    hold on
    if t_i == 3
    yline(baseline(1));
    else
    end
    hold on
    legend('F', 'dF','init baseline','Box','off','Orientation','horizontal')
    title([sessionTitle, num2str(time_ranges(t_i))], 'Interpreter','none')
    ylabel('Fiber Signal (V)'); xlabel('Time (s)');
end




