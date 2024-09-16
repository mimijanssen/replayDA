%% Plot Raw Fiber Data
% example: cd = C:\Data\M406\2022-10-30_M406_L_sleep_track_train
%folder = 'C:\Data\M406\2022-10-30_M406_L_sleep_track_train' ;
%cd folder 

%% Path
addpath(genpath('Users\mimia\Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\mimia\Documents\GitHub\replay_DA\analysis\photometry'));
%% Load Data
cfg.fc = {'CSC30.ncs'};
csc_photo = LoadCSC(cfg);

%% Plot 
FP_data = []; % initializing fiber photometry data
FP_data.acq.Fs = csc_photo.cfg.hdr{1}.SamplingFrequency; %set FP_data.acq.Fs to sampling frequency rate (5000 points per second) 
FP_data.acq.time = csc_photo.tvec - csc_photo.tvec(1); %set FP_data.acq.Fs time to the time vector subtracted by the first time point
%this initializes the time vector to start at 0 then 2.0 x 10 ^-4 or 0.0002 and 4.0 x
%x 10 ^ -4 0.0004 ... 0.0002 second time interval means 5,000 points per
%second
FP_data.acq.FP{1} = csc_photo.data'; %fiber data 

%rename for plot
FP = (FP_data.acq.FP{1})*(1239/39); %can use 32I need to verify what the voltage divider is. but I'm using this for now. 
time = FP_data.acq.time; 

%% Detrend raw data with locdetrend
FP_detrended = locdetrend(FP);

%% Load events
LoadExpKeys();
cfg_evt = [];
cfg_evt.eventList = ExpKeys.eventList;
cfg_evt.eventLabel = ExpKeys.eventLabel;
evt = LoadEvents(cfg_evt);


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

%% Plot Detrended data
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
    plot(time(t_range), FP_detrended(t_range), 'Color', [0 0.5 0])
    title([sessionTitle, num2str(time_ranges(t_i)),'samples'], 'Interpreter','none')
    ylabel('Fiber Signal (V)'); xlabel('Time (s)');
end

% plot a figure around those x times
selected_x_end = nearest_idx3(selected_x, time) + (5/0.0002);
selected_x_start = nearest_idx3(selected_x, time) - (5/0.0002);

figure(2)
for t_i2 = 1:4
    subplot(4, 1, t_i2);
    plot((time(selected_x_start(t_i2):selected_x_end(t_i2)))-time(selected_x_start(t_i2)), FP_detrended(selected_x_start(t_i2):selected_x_end(t_i2)), 'Color', [0 0.5 0])
title(['10 second window around selected points'])
ylabel('Fiber Signal (V)'); xlabel('Time (s)');
end


%% Plot Trial Data 
first_time = evt.t{1,1}(1) - csc_photo.tvec(1);
last_time = evt.t{1,2}(20) - csc_photo.tvec(1); 

%indinit = find(time==first_time);
% tolerance 
indinit = find(abs(time-first_time) < 0.0001);
%indfin = find(time ==last_time);
indfin = find(abs(time-last_time) < 0.0001);

plot(time(indinit:indfin), FP(indinit:indfin), 'Color', [0 0.5 0])

title('Linear Track Signal')
ylabel('Fiber Signal (V)'); xlabel('Time (s)');


%% Plot Wavesurfer over NYX
load('M406_20221030_0001.mat')

% for this you should subtract from all values the initial time 
plot(time(indinit:indfin)-time(indinit), FP(indinit:indfin), 'Color', [0 0.5 0])
hold on
plot(data.acq.time,data.acq.FP{1,1})

title('Linear Track Signal')
ylabel('Fiber Signal (V)'); xlabel('Time (s)');
legend('Neuralynx','Wavesurfer')
hold off

