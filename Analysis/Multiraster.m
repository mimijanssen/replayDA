%% 
cd C:\Data\M433\2023-09-21_M433_recording3
%% 
load('M433-2023-09-21processed.mat'); % CHANGE THIS PER SESSION 
%%
LoadExpKeys();

please = [];
please.fc = ExpKeys.goodSpikes;
S = LoadNST(please);

please = []; please.fc = {'CSC15.ncs'};
csc_s = LoadCSC(please);

please = []; please.fc = {'CSC30.ncs'};
csc = LoadCSC(please);
csc.data = dF_win; 
csc.tvec = t;
csc.cfg.hdr{1,1}.SamplingFrequency = 1000; 

%
cfg_plot = [];
cfg_plot.lfp(1) = csc_s;
cfg_plot.lfp(2) = csc;
%cfg_plot.evt = metadata.taskvars.trial_iv_L; % "left" trials on the T-maze
cfg.lfpHeight = 30;
h = MultiRaster(cfg_plot,S);