%% 
clear, clc;
cd 'F:\M453\M453-2024-01-14_recording2'
%FP = load('M453_2024_01_15processed.mat'); % CHANGE THIS PER SESSION 

%%
% extract events 
LoadExpKeys
cfg_evt = [];
evt2 = LoadEvents(cfg_evt);

please = [];
please.fc = ExpKeys.goodSpikes;
S = LoadNST(please);

% extract LFP 
csc_name = [];
csc_name.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(csc_name); % csc with good ripples

% initialize LFP 
lfp_time = csc.tvec- csc.tvec(1); % lfp time 
lfp = csc.data; 

cfg_fiber.fc = {'CSC30.ncs'};
csc_photo_fc = LoadCSC(cfg_fiber);

%% rescale fiber 

%photo.data = rescale(csc_photo_fc.data,-0.1,1);
%%

% plot that
cfg_plot = [];
cfg_plot.lfp(1) = csc;
cfg_plot.lfp(2) = csc_photo_fc;
cfg_plot.lfpWidth = 1.5;
%cfg_plot.lfpHeight = 50;
h = MultiRaster(cfg_plot,S);