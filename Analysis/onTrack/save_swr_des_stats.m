%% session saving -- HUGE MATRIX 
% method: 1) for each session, go to that directory and save everything you
% need in a table 2) load all session tables and join them. 

% 6/30/26 added sleep to the matrix

addpath('C:\Users\mimia\Documents\GitHub\replayDA\Analysis') % my analysis code exists here
addpath('C:\Users\mimia\Documents\GitHub\vandermeerlab-replay-da\code-matlab\tasks\Alyssa_Tmaze\beta'); % SWR amplitude detection exists here

%% Set Directory
% input information 
clear; clc;
rng(pi)
cd 'F:\M654\M654_2026_02_08_recording8' 
file_name = 'M654_2026_02_08'; 
mouseID = ['M654'];
session = 8; 
mouse = convertMouse(mouseID); % converted mouse number 


%% Load Files
FP_file=dir('*processed*');
FP = load(FP_file.name); % preprocessed data. 

SWR_file = dir('*detectedSWRs*'); %_reg %basic
load(SWR_file.name) % SWR intervals.

Track_file = dir('*_track.mat*'); 
load(Track_file.name); % for pseudo_outcomes

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

% extract LFP 
csc_name = [];
csc_name.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(csc_name); % csc with good ripples

% initialize LFP 
lfp_time = csc.tvec- csc.tvec(1); % lfp time 
lfp = csc.data; 

post = ExpKeys.postrecord(1) - csc.tvec(1); % time of post sleep period, initialized  
pre_end = ExpKeys.prerecord(2) - csc.tvec(1); % time of post sleep period, initialized  

% initialize event times
% sleep events 
post = ExpKeys.postrecord(1) - csc.tvec(1);
% swr events
SWR_start = evt.tstart- csc.tvec(1); %evt.tstart; just tstart for basic
SWR_end = evt.tend- csc.tvec(1);
SWR_iv = [SWR_start SWR_end];
SWR_ind_start = nearest_idx3(SWR_iv(:,1),lfp_time); % could have just done evt.tstart and csc.time .... 
SWR_ind_end = nearest_idx3(SWR_iv(:,2),lfp_time);

SWR_ind_mid = (SWR_ind_start + SWR_ind_end)/2;  %middle timepoint for each SWR

% track sleep
SWR_ind_start_post = nearest_idx3(post,SWR_iv(:,1)); %find(abs(lfp-SWR_iv(1,1)) < 0.0005); % only for one but can we extend this to everything??
SWR_ind_end_post = nearest_idx3(post,SWR_iv(:,2)); %find(abs(lfp-SWR_iv(1,2)) < 0.0005); % only for one but can we extend this to everything??
SWR_ind_mid_post = (SWR_ind_start_post + SWR_ind_end_post)/2;  %middle index 
% this in the middle index of the first post SWR

% find the corresponding SWR time index closest to the end of pre track
% sleep
SWR_ind_start_pre_end = nearest_idx3(pre_end,SWR_iv(:,1)); 
SWR_ind_end_pre_end = nearest_idx3(pre_end,SWR_iv(:,2)); 
SWR_ind_mid_pre_end = (SWR_ind_start_pre_end + SWR_ind_end_pre_end)/2;  %middle index 

% all post-track rest swrs
post_SWR_ind = SWR_ind_mid(SWR_ind_mid > SWR_ind_mid(round(SWR_ind_mid_post))); %SWR index that is greater than first post time point

% all pre-track rest swrs 
pre_SWR_ind = SWR_ind_mid(SWR_ind_mid < SWR_ind_mid(round(SWR_ind_mid_pre_end))); %all middle SWR timepoints that are less then the 

% track 
track_SWR_ind = SWR_ind_mid(SWR_ind_mid_pre_end:SWR_ind_mid_post); %all middle SWR timepoints

% for fiber after swr
SWR_time_mid = zeros(length(SWR_ind_mid),1); 
SWR_fiber_ind = zeros(length(SWR_ind_mid),1); 
for i = 1:1:size(SWR_ind_mid,1)
    SWR_time_mid(i) = lfp_time(round(SWR_ind_mid(i))); % lfp time in (s) for a swr
    SWR_fiber_ind(i) = nearest_idx3(SWR_time_mid(i),FP.tvec); % find the corresponidng time in fiber and saves the index. FP.tvec is in seconds 
end

pre_count = length(pre_SWR_ind);
post_count = length(post_SWR_ind);
track_count = length(track_SWR_ind);

swr_des.swr_count = [pre_count track_count post_count];
swr_des.swr_label = ["Pre" "Track" "Post"]; 

%% save frequency and save avg duration in the same way as pre_count...
% ~~ Frequency ~~
% pre-task sleep time (min)
prerecord_init = ExpKeys.prerecord(1)-csc.tvec(1);
prerecord_end = ExpKeys.prerecord(2)-csc.tvec(1);

time_pre = (prerecord_end - prerecord_init)/60; % in minutes
freq_pre = pre_count/time_pre; 

% track time (min) 
track_init = ExpKeys.task(1)-csc.tvec(1);
track_end = ExpKeys.task(2)-csc.tvec(1);

time_track = (track_end - track_init)/60; % in minutes
freq_track = track_count/time_track; 

% post-task sleep time (min) 
postrecord_init = ExpKeys.postrecord(1)-csc.tvec(1);
postrecord_end = ExpKeys.postrecord(2)-csc.tvec(1);

time_post = (postrecord_end- postrecord_init)/60;
freq_post = post_count/time_post; 

swr_des.freq = [freq_pre freq_track freq_post];

% ~~ Duration ~~
% in miliseconds 

% find SWR indicies that are before the posts 
SWR_start_pre = lfp_time(round(SWR_ind_start(SWR_ind_start < SWR_ind_start(round(SWR_ind_mid_pre_end))))); % in sec
SWR_end_pre = lfp_time(round(SWR_ind_end(SWR_ind_end < SWR_ind_end(round(SWR_ind_mid_pre_end))))); % find all timepoints that are before the middle of the first post swr
pre_avg_dur = (mean(SWR_end_pre - SWR_start_pre))*1000; 
% 0.0677 seconds -- about 70 ms. checks out!

SWR_start_track = lfp_time(round(SWR_ind_start(SWR_ind_mid_pre_end:SWR_ind_mid_post))); % in sec
SWR_end_track = lfp_time(round(SWR_ind_end(SWR_ind_mid_pre_end:SWR_ind_mid_post))); % find all timepoints that are before the middle of the first post swr
track_avg_dur = (mean(SWR_end_track - SWR_start_track))*1000;

SWR_start_post = lfp_time(round(SWR_ind_start(SWR_ind_start > SWR_ind_start(round(SWR_ind_mid_post))))); % find all timepoints that are before the middle of the first post swr
SWR_end_post = lfp_time(round(SWR_ind_end(SWR_ind_end > SWR_ind_end(round(SWR_ind_mid_post))))); % find all timepoints that are before the middle of the first post swr
post_avg_dur = (mean(SWR_end_post - SWR_start_post))*1000; 
% 0.0723 seconds  

swr_des.dur = [pre_avg_dur track_avg_dur post_avg_dur];

% ~~ Save Variables ~~ 
cd ('D:\trackcounts')
filename = append(file_name, "swr_des_track.mat");
save(filename, '-struct','swr_des')   