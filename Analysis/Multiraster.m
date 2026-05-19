%% 
cd C:\Data\M433\2023-09-21_M433_recording3
%% 
load('M433-2023-09-21processed.mat'); % CHANGE THIS PER SESSION 
%%
LoadExpKeys();

please = [];
please.fc = ExpKeys.goodSpikes;
S = LoadNST(please);

%% 
S = restrict(S,ExpKeys.postrecord(1),ExpKeys.postrecord(2));

%%
csc_name = [];
csc_name.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(csc_name); % csc with good ripples

csc = restrict(csc,ExpKeys.postrecord(1),ExpKeys.postrecord(2));
%%
%
cfg_plot = [];
cfg_plot.lfp(1) = csc;
cfg_plot.lfp(2) = theta_data.power_post;
%cfg_plot.evt = metadata.taskvars.trial_iv_L; % "left" trials on the T-maze
cfg.lfpHeight = 30;
h = MultiRaster(cfg_plot,S);