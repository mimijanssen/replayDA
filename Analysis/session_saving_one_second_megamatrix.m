%% session saving -- HUGE MATRIX 
% method: 1) for each session, go to that directory and save everything you
% need in a table 2) load all session tables and join them. 

% This is step one. 
addpath('C:\Users\mimia\Documents\GitHub\replayDA\Analysis')

% I might want to try this with F_zscored_win (z-scored windowed dtrend)
% and with zdF_win (z-scored and dF/F over windows) and dF_win (dF/F)  
%% Set Directory
% input information 
clear; clc;
rng(pi)
cd 'D:\M600\M600_2025_01_22_recording8'; 
file_name = 'M600_2025_01_22'; 
mouseID = ['M600'];
session = 8; 
mouse = convertMouse(mouseID); % converted mouse number 

%% Load Files
FP_file=dir('*processed*');
FP = load(FP_file.name); % preprocessed data. 

SWR_file = dir('*detectedSWRs*');
load(SWR_file.name) % SWR intervals.

Track_file = dir('*_track.mat*'); 
load(Track_file.name); % for pseudo_outcomes

%DLC_file = dir('*convertedDLC*'); 
%P = readtable(DLC_file.name,'PreserveVariableNames',true);

% Display a message if everything loaded correctly 
% meaning the right number of files were found : 
if length(FP_file) == 1 && length(SWR_file)==1 && length(Track_file) ==1 %&& length(DLC_file)==1
    disp('loading looks good')
else
    disp('check files')
end

%% Determine SWR Counts for Pre and Post Sleep Sessions
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

swr_des.swr_count = [pre_count post_count];
swr_des.swr_label = ["Pre" "Post"]; 

%% save frequency and save avg duration in the same way as pre_count...
% ~~ Frequency ~~
% pre-task sleep time (min)
prerecord_init = ExpKeys.prerecord(1)-csc.tvec(1);
prerecord_end = ExpKeys.prerecord(2)-csc.tvec(1);

time_pre = (prerecord_end - prerecord_init)/60; % in minutes
freq_pre = pre_count/time_pre; 

% post-task sleep time (min) 
postrecord_init = ExpKeys.postrecord(1)-csc.tvec(1);
postrecord_end = ExpKeys.postrecord(2)-csc.tvec(1);

time_post = (postrecord_end- postrecord_init)/60;
freq_post = post_count/time_post; 

swr_des.freq = [freq_pre freq_post];

% ~~ Duration ~~
% in miliseconds 

% find SWR indicies that are before the posts 
SWR_start_pre = lfp_time(round(SWR_ind_start(SWR_ind_start < SWR_ind_start(round(SWR_ind_mid_post))))); % in sec
SWR_end_pre = lfp_time(round(SWR_ind_end(SWR_ind_end < SWR_ind_end(round(SWR_ind_mid_post))))); % find all timepoints that are before the middle of the first post swr
pre_avg_dur = (mean(SWR_end_pre - SWR_start_pre))*1000; 
% 0.0677 seconds -- about 70 ms. checks out!

SWR_start_post = lfp_time(round(SWR_ind_start(SWR_ind_start > SWR_ind_start(round(SWR_ind_mid_post))))); % find all timepoints that are before the middle of the first post swr
SWR_end_post = lfp_time(round(SWR_ind_end(SWR_ind_end > SWR_ind_end(round(SWR_ind_mid_post))))); % find all timepoints that are before the middle of the first post swr
post_avg_dur = (mean(SWR_end_post - SWR_start_post))*1000; 
% 0.0723 seconds  

swr_des.dur = [pre_avg_dur post_avg_dur];

% ~~ Save Variables ~~ 
%filename = append(file_name, "swr_des.mat");
%save(filename, '-struct','swr_des')

%% start matrix
% each swr has it's own row 
matrix_sess = array2table(zeros(length(SWR_ind_mid),17),'VariableNames',{'mouseID','sess','swrID','PrePost','OnesPreRaw','OnesPostRaw','OnesPreProc','OnesPostProc','OnesBeforePeak','OnesAfterPeak','OnesBeforeAUC','OnesAfterAUC','TimeAfterPeak','OnesBeforePeakRAW','OnesAfterPeakRAW','OnesBeforeAUCRAW','OnesAfterAUCRAW'});

%% input mouse/session identity information 
matrix_sess.('mouseID')(:,1) = mouse; 
matrix_sess.('sess')(:,1) = session;
matrix_sess.('swrID')(:) = linspace(0,1,length(matrix_sess.('swrID')))';
% 1 for pre
% everything before SWR_ind_mid_post index gets a 1. 
matrix_sess.('PrePost')(1:round(SWR_ind_mid_post)-1,1) = 1;
% 2 for post. 
matrix_sess.('PrePost')(round(SWR_ind_mid_post):end,1) = 2;

%% populate matrix with structure information. 
seconds = 1; 
samples = (seconds*FP.cfg.hdr{1,1}.SamplingFrequency)/2; % samples I will take before and after swrs

% raw data 
matrix_sess.('OnesPreRaw') = cell(height(matrix_sess), 1);
matrix_sess.('OnesPostRaw') = cell(height(matrix_sess), 1);
% preprocessed data 
matrix_sess.('OnesPreProc') = cell(height(matrix_sess), 1);
matrix_sess.('OnesPostProc') = cell(height(matrix_sess), 1);

for i = 1:height(matrix_sess) % iterate through each swr. 
    raw_data_struct = struct(); % structure of the data for an individual swr
    data_struct = struct(); % structure of the data for an individual swr
    swr_time = lfp_time(round(SWR_ind_mid(i))); % swr lfp time (initialized) 
    fiber_index = nearest_idx3(swr_time, FP.tvec); % fiber index closest to middle swr_time -- ok make sure this time is initialized 
    if matrix_sess.('PrePost')(i) == 1 % if pre-track rest
        raw_data_struct.signal = FP.data(fiber_index-samples:fiber_index+samples); % saving signal
        data_struct.signal = FP.zF_win_60s(fiber_index-samples:fiber_index+samples); % saving signal
        raw_data_struct.tvec = FP.tvec(fiber_index-samples:fiber_index+samples); % saving time as well- even though it should be the same for each
        data_struct.tvec = raw_data_struct.tvec;     % saving time as well- even though it should be the same for each
        matrix_sess.('OnesPreRaw'){i} = raw_data_struct;
        matrix_sess.('OnesPreProc'){i} = data_struct;
    else % else - post-track rest
        raw_data_struct.signal = FP.data(fiber_index-samples:fiber_index+samples); % saving signal
        data_struct.signal = FP.zF_win_60s(fiber_index-samples:fiber_index+samples); % saving signal
        raw_data_struct.tvec = FP.tvec(fiber_index-samples:fiber_index+samples);     % saving time as well- even though it should be the same for each
        data_struct.tvec = raw_data_struct.tvec;     % saving time as well- even though it should be the same for each
        matrix_sess.('OnesPostRaw'){i} = raw_data_struct;
        matrix_sess.('OnesPostProc'){i} = data_struct;
    end
end


%% populate dF information on from preproc data 
% dF from 2 seconds 
x1 = 1:1:500;%1001:1:3000; %1:1:2000;% % two seconds before for preproc data
x2 = 501:1:1000; %3001:1:5000; %2001:1:4000;%  % two seconds after for preproc data

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

%% Save everything
cd 'D:\SWR_DA_MegaMatrix_1s'
filename = append(file_name, "mega1.mat");
save(filename,'matrix_sess')
   