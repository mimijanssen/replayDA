% SWR peaks vs. Frequency Bands and Movement Speed 

% 1) find SWRs 
% 2) find SWR-DA peaks within one second
% 3) fit a spline
% 4) run a spectrogram
% 5) save movement speed
% 6) run a correlation between signals 

%%  Load Data
clear; clc;
rng(pi)
cd 'F:\M533\M533_2024_08_20_recording2'; 
file_name = 'M533_2024_08_20'; 

FP_file=dir('*processed*');
FP = load(FP_file.name); % preprocessed data. 

SWR_file = dir('*detectedSWRs*');
load(SWR_file.name) % SWR intervals.

Track_file = dir('*_track.mat*'); 
load(Track_file.name); % for pseudo_outcomes

DLC_file = dir('*convertedDLC*'); 
P = readtable(DLC_file.name,'PreserveVariableNames',true);

% Display a message if everything loaded correctly 
% meaning the right number of files were found : 
if length(FP_file) == 1 && length(SWR_file)==1 && length(Track_file) ==1 && length(DLC_file)==1
    disp('loading looks good')
else
    disp('check files')
end

%% 1) Find SWRs 
% extract events (times of sleep sessions)
LoadExpKeys
cfg_evt = [];
evt2 = LoadEvents(cfg_evt);

% load raw fiber data 
%cfg_fiber.fc = {'CSC30.ncs'};
%raw_fiber = LoadCSC(cfg_fiber);
%raw_fiber_time = raw_fiber.tvec - raw_fiber.tvec(1); 

% extract LFP 
csc_name = [];
csc_name.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(csc_name); % csc with good ripples
% initialize LFP 
lfp_time = csc.tvec- csc.tvec(1); % lfp time 
lfp = csc.data; 

%time = FP.tvec; % fiber time - processed

% time of post recording 
post = ExpKeys.postrecord(1) - csc.tvec(1); % time of post sleep period, initialized  

% initialize SWR interval
SWR_start = evt.tstart- csc.tvec(1);
SWR_end = evt.tend- csc.tvec(1);
SWR_iv = [SWR_start SWR_end];
SWR_ind_start = nearest_idx3(SWR_iv(:,1),lfp_time); %find(abs(lfp-SWR_iv(1,1)) < 0.0005); % only for one but can we extend this to everything??
SWR_ind_end = nearest_idx3(SWR_iv(:,2),lfp_time); %find(abs(lfp-SWR_iv(1,2)) < 0.0005); % only for one but can we extend this to everything??

SWR_ind_mid = (SWR_ind_start + SWR_ind_end)/2;  %middle timepoint for each SWR

% find that corresponding SWR time index closest to that 
SWR_ind_start_post = nearest_idx3(post,SWR_iv(:,1)); %find(abs(lfp-SWR_iv(1,1)) < 0.0005); % only for one but can we extend this to everything??
SWR_ind_end_post = nearest_idx3(post,SWR_iv(:,2)); %find(abs(lfp-SWR_iv(1,2)) < 0.0005); % only for one but can we extend this to everything??
SWR_ind_mid_post = (SWR_ind_start_post + SWR_ind_end_post)/2;  %middle index 

% keep all SWR after that time -- should be 654
post_SWR_ind = SWR_ind_mid(SWR_ind_mid > SWR_ind_mid(round(SWR_ind_mid_post))); %SWR index that is greater than first post time point

% Pre track rest ^.^ 
pre_SWR_ind = SWR_ind_mid(SWR_ind_mid < SWR_ind_mid(round(SWR_ind_mid_post))); %all middle SWR timepoints that are less then the 

% for fiber after swr
SWR_time_mid = zeros(length(SWR_ind_mid),1); 
SWR_fiber_ind = zeros(length(SWR_ind_mid),1); 
for i = 1:1:size(SWR_ind_mid,1)
    SWR_time_mid(i) = lfp_time(round(SWR_ind_mid(i))); % lfp time in (s) for a swr
    SWR_fiber_ind(i) = nearest_idx3(SWR_time_mid(i),FP.tvec); % find the corresponidng time in fibFPer and saves the index.
end

pre_count = length(pre_SWR_ind);
post_count = length(post_SWR_ind);

prerecord_init = ExpKeys.prerecord(1)-csc.tvec(1);
prerecord_end = ExpKeys.prerecord(2)-csc.tvec(1);

postrecord_init = ExpKeys.postrecord(1)-csc.tvec(1);
postrecord_end = ExpKeys.postrecord(2)-csc.tvec(1);


%% FIND SWR-DA peaks within One Second. 
% dF from 2 seconds 
seconds = 2; 
samples = (seconds*FP.cfg.hdr{1,1}.SamplingFrequency)/2; % samples I will take before and after swrs

x1 = 1:1:1000;%1001:1:3000; %1:1:2000;% % two seconds before for preproc data
x2 = 1001:1:2001; %3001:1:5000; %2001:1:4000;%  % two seconds after for preproc data

for i = 1:height(matrix_sess) % iterate through each swr. 
    swr_time = lfp_time(round(SWR_ind_mid(i))); % swr lfp time  
    fiber_index = nearest_idx3(swr_time, FP.tvec); % fiber index closest to middle swr_time 
    signal = FP.zF_win_60s(fiber_index-samples:fiber_index+samples); 
    signal_raw = FP.data(fiber_index-samples:fiber_index+samples); 
    matrix_sess.('OnesBeforePeak')(i) = max(signal(x1)); %-min(signal(x1)); % maybe the average signal might be better than the lowest signal?? 
    [matrix_sess.('OnesAfterPeak')(i),I] = max(signal(x2)); %-min(signal(x2)); 
    matrix_sess.('TimeAfterPeak')(i) = FP.tvec(I); % time of the peak post swr
    matrix_sess.('OnesBeforePeakRAW')(i) = max(signal_raw(x1)); %-min(signal(x1)); % maybe the average signal might be better than the lowest signal?? 
    matrix_sess.('OnesAfterPeakRAW')(i) = max(signal_raw(x2)); %-min(signal(x2)); 
end

% I changed how I did this so it is just max and not max - min values... 
% ask matt if this is ok ~

% populate AUC information 
for i = 1:height(matrix_sess) % iterate through each swr. 
    swr_time = lfp_time(round(SWR_ind_mid(i))); % swr lfp time  
    fiber_index = nearest_idx3(swr_time, FP.tvec); % fiber index closest to middle swr_time 
    signal = FP.zF_win_60s(fiber_index-samples:fiber_index+samples); 
    signal_raw = FP.data(fiber_index-samples:fiber_index+samples); 
    matrix_sess.('OnesBeforeAUC')(i) = trapz(x1,signal(x1)); 
    matrix_sess.('OnesAfterAUC')(i) = trapz(x2,signal(x2)); 
    matrix_sess.('OnesBeforeAUCRAW')(i) = trapz(x1,signal_raw(x1)); 
    matrix_sess.('OnesAfterAUCRAW')(i) = trapz(x2,signal_raw(x2)); 
end
