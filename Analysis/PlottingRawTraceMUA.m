%% Pick the recording you want to look at. 
clear, clc;
cd 'D:\M453\M453-2024-01-15_recording3'
FP = load('M453_2024_01_15processed'); % fiber data processed with my pipeline
iv_in = load('2024-01-15_M453_recording3detectedSWRs');

%% extract events 
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

%% Restore variables for plotting on same time scale with preprocessed files. 
S2 = S; 
S2.t{1,1} = S.t{1,1} - csc.tvec(1); 

csc2 = tsd;
csc2.tvec = csc.tvec - csc.tvec(1);
csc2.data = csc.data;
csc2.cfg = FP.cfg; 

FP2 = tsd;
FP2.tvec = FP.tvec;
FP2.data = FP.zF_win_60s'; 
FP2.cfg = FP.cfg; 

%% Ducktrap time
cfg_plot = [];
%cfg_plot.evt = iv_in;% actually you can just check if these were detected
%SWRs by checking the interval.... 
cfg_plot.lfp(1) = csc2;
cfg_plot.lfp(2) = FP2; %csc_photo_fc;
h = MultiRaster(cfg_plot,S2);

%% Check
iv_in.evt.tstart2 = iv_in.evt.tstart - csc.tvec(1); 
iv_in.evt.tend2 = iv_in.evt.tend - csc.tvec(1); 
iv_in.evt.center = (iv_in.evt.tstart2 + iv_in.evt.tend2)/2 ; 
%idx = nearest_idx3(4920,iv_in.evt.tstart); 
%time = (iv_in.evt.tend(idx));
% OK was not a detected one...