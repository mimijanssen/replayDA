%% Parameters
clear all; clc
%IF NOT RAW
params.dsRate = 100; % Downsampling rate if you want to downsample the signal
%This dsRate will also be applied to all signals during the analysis
%pipeline

% Filter Parameters
params.FP.lpCut = 10; % Cut-off frequency for filter
params.FP.filtOrder = 8; % Order of the filter

% Baseline Parameters
params.FP.basePrc = 5; % Percentile value from 1 - 100 to use when finding baseline points
%Note: Lower percentiles are used because the mean of signal is not true
%baseline
params.FP.winSize = 10; % Window size for baselining in seconds
params.FP.winOv = 0; %Window overlap size in seconds
params.FP.interpType = 'linear'; % 'linear' 'spline' 
params.FP.fitType = 'interp'; % Fit method 'interp' , 'exp' , 'line'

% Demodulation Parameters
%When demodulating signals, the filter creates edge artifacts. We record
%for a few seconds longer, so we can remove x seconds from the beginning
%and end
%Adjust the variable to "0" if it's a normal photometry recording
params.FP.sigEdge = 15; %15; %Time in seconds of data to be removed from beginning and end of signal
params.FP.modFreq = [319 217];


%% Load Neuralynx CSC photometry data
cd C:\Data\M406\2022-10-30_M406_L_sleep_track_train;  LoadExpKeys; LoadMetadata;
cfg.fc = {'CSC30.ncs'};
csc_photo = LoadCSC(cfg);

FP_data = [];
FP_data.acq.Fs = csc_photo.cfg.hdr{1}.SamplingFrequency;
FP_data.acq.time = csc_photo.tvec - csc_photo.tvec(1); % initialized time
FP_data.acq.FP{1} = csc_photo.data';

FP = FP_data.acq.FP{1};
fsig = tsd(FP_data.acq.time, FP);

%% Process data recorded in CW mode
% IF NOT RAW
% The process detrends, demodulates, removes artifacts, and finds delta F
%FP_data = processFP(params, FP_data);
%FP = tsd(FP_data.final.time', FP_data.final.FP{1}', 'FP');

%% Path
addpath(genpath('Users\mimia\Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\mimia\Documents\GitHub\matlab-utils\chronux\spectral_analysis\continuous'));
addpath(genpath('C:\Users\mimia\Documents\Replay-DA\analysis\photometry'));

%% Load events
LoadExpKeys();
cfg_evt = [];
cfg_evt.eventList = ExpKeys.eventList;
cfg_evt.eventLabel = ExpKeys.eventLabel;
evt = LoadEvents(cfg_evt);

%% Detrend data
FP_detrend = locdetrend(FP); %global detrend

%make a new structure with detrend 
csc_photo2 = csc_photo;
csc_photo2.data = FP_detrend';

%% Seperate Data
% Early trials 
early_csc_photo = csc_photo2.data(1:(length(csc_photo2.data)/3));
% Middle trials 

% Late trials 

%% Align data based on maximum value within that time window. 
% within a certain window (5s before and after an event) specified by
% evt.t, find the maximum value of signal and save that time stamp into a
% new evt.t 
cfg_def.dt = 0.01; % time step, specify this for 'interp' mode
cfg_def.mode = 'interp'; % 'raw' or 'interp'; if 'interp', need to specify cfg_def.dt
cfg_def.interp_mode = 'linear';
cfg_def.window = [-5 5]; % start and end times of window (in s) % is this 5 seconds? if it is adding 5 to the time? 

% time ordered from least to greatest 
evt_ordered = sort([evt.t{1},evt.t{2}]); 

% gets an interval from raw times 
cfg = ProcessConfig2(cfg_def,csc_photo2);
t_in = evt_ordered; 
tsd_in = cfg;
if ~isfield(t_in,'type') % assume raw times
    this_iv = iv(t_in + cfg.window(1), t_in + cfg.window(2));
end
        % make big idx matrix
        start_idx = nearest_idx3(this_iv.tstart,tsd_in.tvec);
        end_idx = nearest_idx3(this_iv.tend,tsd_in.tvec);

for triali = 1:1:length(this_iv.tstart)
    plot([cfg.data(start_idx(triali):end_idx(triali))])
    max_signal = max([cfg.data(start_idx(triali):end_idx(triali))]); %max signal value
    signal_index = find([cfg.data(start_idx(triali):end_idx(triali))] == max_signal); % index of max signal value 
    max_signal_time = cfg.tvec(signal_index); %value of time
    evt2_ordered(triali) = max_signal_time; 
    hold on
end
    max_signal_time2 = nearest_idx3(evt2_ordered,tsd_in.tvec);
evt2_ordered = max_signal_time2';


%% PSTH
% Color matrix 
count = 0.9 ;
jet = zeros(40,3);
for zi = 1:1:39
    count =  count - 0.020; %gets darker over trials 
    jet(zi,2) = count;
end

% for df
cfg_def.dt = 0.01; % time step, specify this for 'interp' mode
cfg_def.mode = 'interp'; % 'raw' or 'interp'; if 'interp', need to specify cfg_def.dt
cfg_def.interp_mode = 'linear';
cfg_def.window = [-5 5]; % start and end times of window (in s)

% Order time events from least to greatest 
%evt_ordered = sort([evt.t{1},evt.t{2}]); % this is all left and then all right but we want left right left right
% This changes the input to TSDpeth so you need to change that. 

% IF MAX - MIN in that time is > 0.2 V then delete that trial 
if evt.t{1}(1) < evt.t{2}(1) % if the left is smaller than the right time: 
for triali = 1:1:40 %(min(length(evt.t{1}),length(evt.t{2})))*2 %length(evt.t{1}) 40 for the color but 20 for the trial 
if mod (triali,2) == 0 
 right_end = evt_ordered(triali); 
    right_end_peth = TSDpeth(cfg_def, csc_photo2, right_end); %raw = cscphoto
    if (max(right_end_peth.data*32) - min(right_end_peth.data*32)) < 0.13
    max_signal =  max(right_end_peth.data);
    max_signal_index = find(right_end_peth.data == max_signal);
    time_plot = right_end_peth.tvec(max_signal_index);
    time_plotstart = right_end.tvec(max_signal_index)-2;
    time_plotend = right_end.tvec(max_signal_index)+2;
    plot(time_plotstart:time_plotend,right_end_peth.data*31.77,'Color', [jet(triali,:)])    
    else 
        continue
    end
else 
    left_end = evt_ordered(triali); 
    left_end_peth = TSDpeth(cfg_def, csc_photo2, left_end); %RAW = csc_photo; dF = FP
    if (max(left_end_peth.data*32) - min(left_end_peth.data*32)) < 0.13
    max_signal =  max(left_end_peth.data);
    max_signal_index = find(left_end_peth.data == max_signal);
    time_plot = left_end_peth.tvec(max_signal_index);
    time_plotstart = left_end.tvec(max_signal_index)-2;
    time_plotend = left_end.tvec(max_signal_index)+2;
    plot(time_plotstart:time_plotend,left_end_peth.data*31.77,'Color', [jet(triali,:)])    
    else
        continue
    end
    hold on
end
   % legend('Left end', 'Right end')
end

colormap(jet);
cbh = colorbar;
cbh.Ticks = linspace (0,1,9); % 9 - need to change this based on length
cbh.TickLabels = num2cell(0:5:40); % 40 - need to change this based on length
cbh.Label.String = 'Trial';
xlabel('time (s)')
ylabel('F')
title('F PET with 10s window')

else 
for triali = 1:1:(min(length(evt.t{1}),length(evt.t{2})))*2;
if mod (triali,2) == 0 
 left_end = evt_ordered(triali); 
    left_end_peth = TSDpeth(cfg_def, csc_photo2, left_end); %RAW = csc_photo; dF = FP
    if (max(left_end_peth.data*32) - min(left_end_peth.data*32)) < 0.13
    max_signal =  max(left_end_peth.data);
    max_signal_index = find(left_end_peth.data == max_signal);
    time_plot = left_end_peth.tvec(max_signal_index);
    time_plotstart = left_end.tvec(max_signal_index)-2;
    time_plotend = left_end.tvec(max_signal_index)+2;
    plot(time_plotstart:time_plotend,left_end_peth.data*31.77,'Color', [jet(triali,:)])
    %plot(left_end_peth.tvec, (left_end_peth.data)*31.77, 'Color', [jet(triali,:)]);
    else 
        continue
    end
else 
    right_end = evt_ordered(triali); 
    right_end_peth = TSDpeth(cfg_def, csc_photo2, right_end);
    if (max(right_end_peth.data*32) - min(right_end_peth.data*32)) < 0.13
    max_signal =  max(right_end_peth.data);
    max_signal_index = find(right_end_peth.data == max_signal);
    time_plot = right_end_peth.tvec(max_signal_index);
    time_plotstart = right_end.tvec(max_signal_index)-2;
    time_plotend = right_end.tvec(max_signal_index)+2;
    plot(time_plotstart:time_plotend,right_end_peth.data*31.77,'Color', [jet(triali,:)])
    else
        continue
    end
    hold on
end
   % legend('Left end', 'Right end')
end

colormap(jet);
cbh = colorbar;
cbh.Ticks = linspace (0,1,9); % 9 - need to change this based on length
cbh.TickLabels = num2cell(0:5:40); % 40 - need to change this based on length
cbh.Label.String = 'Trial';
xlabel('time (s)')
ylabel('F')
title('F PET with 10s window')

end

%m409 10/28 won't plot with FP - why is that? 
% 10/27 won't plot with FP
% 10/30 won't plot with FP 


%% PSTH
% Color matrix 
count = 0.9 ;
jet = zeros(40,3);
for zi = 1:1:39
    count =  count - 0.020; %gets darker over trials 
    jet(zi,2) = count;
end

% for df
cfg_def.dt = 0.01; % time step, specify this for 'interp' mode
cfg_def.mode = 'interp'; % 'raw' or 'interp'; if 'interp', need to specify cfg_def.dt
cfg_def.interp_mode = 'linear';
cfg_def.window = [-5 5]; % start and end times of window (in s)

% Order time events from least to greatest 
%evt_ordered = sort([evt.t{1},evt.t{2}]); % this is all left and then all right but we want left right left right
% This changes the input to TSDpeth so you need to change that. 

% IF MAX - MIN in that time is > 0.2 V then delete that trial 
if evt.t{1}(1) < evt.t{2}(1) % if the left is smaller than the right time: 
for triali = 1:1:40 %(min(length(evt.t{1}),length(evt.t{2})))*2 %length(evt.t{1}) 40 for the color but 20 for the trial 
if mod (triali,2) == 0 
 right_end = evt2_ordered(triali); 
    right_end_peth = TSDpeth(cfg_def, csc_photo2, right_end); %raw = cscphoto
    if (max(right_end_peth.data*32) - min(right_end_peth.data*32)) < 0.13
    plot(right_end_peth.tvec, (right_end_peth.data)*31.77, 'Color',[jet(triali,:)]); % some sessions will not plot because data is NaN
    else 
        continue
    end
else 
    left_end = evt2_ordered(triali); 
    left_end_peth = TSDpeth(cfg_def, csc_photo2, left_end); %RAW = csc_photo; dF = FP
    if (max(left_end_peth.data*32) - min(left_end_peth.data*32)) < 0.13
    plot(left_end_peth.tvec, (left_end_peth.data)*31.77, 'Color', [jet(triali,:)]);
    else
        continue
    end
    hold on
end
   % legend('Left end', 'Right end')
end

colormap(jet);
cbh = colorbar;
cbh.Ticks = linspace (0,1,9); % 9 - need to change this based on length
cbh.TickLabels = num2cell(0:5:40); % 40 - need to change this based on length
cbh.Label.String = 'Trial';
xlabel('time (s)')
ylabel('F')
title('F PET with 10s window')

else 
for triali = 1:1:(min(length(evt.t{1}),length(evt.t{2})))*2;
if mod (triali,2) == 0 
 left_end = evt2_ordered(triali); 
    left_end_peth = TSDpeth(cfg_def, csc_photo2, left_end); %RAW = csc_photo; dF = FP
    if (max(left_end_peth.data*32) - min(left_end_peth.data*32)) < 0.13
    plot(left_end_peth.tvec, (left_end_peth.data)*31.77, 'Color', [jet(triali,:)]);
    else 
        continue
    end
else 
    right_end = evt2_ordered(triali); 
    right_end_peth = TSDpeth(cfg_def, csc_photo2, right_end);
    if (max(right_end_peth.data*32) - min(right_end_peth.data*32)) < 0.13
    plot(right_end_peth.tvec, (right_end_peth.data)*31.77, 'Color',[jet(triali,:)]);
    else
        continue
    end
    hold on
end
   % legend('Left end', 'Right end')
end

colormap(jet);
cbh = colorbar;
cbh.Ticks = linspace (0,1,9); % 9 - need to change this based on length
cbh.TickLabels = num2cell(0:5:40); % 40 - need to change this based on length
cbh.Label.String = 'Trial';
xlabel('time (s)')
ylabel('F')
title('F PET with 10s window')

end

%m409 10/28 won't plot with FP - why is that? 
% 10/27 won't plot with FP
% 10/30 won't plot with FP 


%% Look at average time between triggers
% microseconds to seconds
left_event_s = evt.t{1}(1:20) ;
right_event_s = evt.t{2}(1:20);

between_beams = left_event_s - right_event_s;
max(between_beams) %86, 107
min(between_beams) %21 , 12 % so you wouldn't
% see the other beams signal unless I make the PSTH over maybe 15 seconds

% Max
%406
% 10/27: 138
% 10/28: 86
% 10/30: 28
% 10/31: 25
% 11/10: 16

%409
% 10/27: 128
% 10/28: 107
% 10:30: 28
% 10/31: 19


% Min
% 406
% 10/27: 61
% 10/28: 21
% 10/30: 11
% 10/31: 11
% 11/10: 8

% M409
% 10/27: 8.9
% 10/28: 12
% 10/30: 12
% 10/31: 12







