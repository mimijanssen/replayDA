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
cd C:\Data\M411\2023-01-23_M411_rac_track;  LoadExpKeys; LoadMetadata;
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
FP_data = processFP(params, FP_data);
FP = tsd(FP_data.final.time', FP_data.final.FP{1}', 'FP');

%% Path
addpath(genpath('Users\mimia\Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\mimia\Documents\Replay-DA\analysis\photometry'));

%% Load events
LoadExpKeys();
cfg_evt = [];
cfg_evt.eventList = ExpKeys.eventList;
cfg_evt.eventLabel = ExpKeys.eventLabel;
evt = LoadEvents(cfg_evt);

%% PSTH
% color 
% increments of 0.05
% increments of 0.025 
% 40 row matrix x 3 colomn
% color matrix 
count = 0 ;

jet = zeros(40,3);
for zi = 1:1:40

    count =  count + 0.025; 
    jet(zi,2) = count;
end

jet1 = [0 0.05 0; 0 0.10 0; 0 0.15 0;0 0.20 0; 0 0.25 0; 0 0.30 0; 0 0.35 0; 0 0.4 0; 0 0.45 0; 0 0.50 0; 0 0.55 0; 0 0.60 0; 0 0.65 0; 0 0.70 0; 0 0.75 0; 0 0.80 0; 0 0.85 0; 0 0.90 0 ; 0 0.95 0; 0 1 0];
jet2 = [0 0 0.05; 0 0 0.10; 0 0 0.15;0 0 0.20 ; 0 0 0.25; 0 0 0.30 ; 0 0 0.35; 0 0 0.4 ; 0 0 0.45 ; 0 0 0.50 ; 0 0 0.55; 0 0 0.60; 0 0 0.65; 0 0 0.70; 0 0 0.75; 0 0 0.80; 0 0 0.85; 0 0 0.90 ; 0 0 0.95; 0 0 1];

for triali = 1:1:20;%length(evt.t{1})

cfg_def.dt = 0.01; % time step, specify this for 'interp' mode
cfg_def.mode = 'interp'; % 'raw' or 'interp'; if 'interp', need to specify cfg_def.dt
cfg_def.interp_mode = 'linear';
cfg_def.window = [-2 2]; % start and end times of window (in s)


left_end = evt.t{1}(triali); %(20); % 3
right_end = evt.t{2}(triali); %(20); % 4 


left_end_peth = TSDpeth(cfg_def, csc_photo, left_end); %RAW = csc_photo; dF = FP
right_end_peth = TSDpeth(cfg_def, csc_photo, right_end);
    plot(left_end_peth.tvec, (left_end_peth.data)*32, 'Color', [jet1(triali,:)]); 
    hold on;
    plot(right_end_peth.tvec, (right_end_peth.data)*32, 'Color',[jet1(triali,:)]); % some sessions will not plot because data is NaN 
    hold on
    %legend('Left end', 'Right end')
end

colormap(jet1);
cbh = colorbar;
cbh.Ticks = linspace (0,1,20);
cbh.TickLabels = num2cell(1:1:21);
cbh.Label.String = 'Trial';
xlabel('time (s)')
ylabel('F')
title('F PET with 10s window')

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







