%% session saving -- HUGE MATRIX 
% method: 1) for each session, go to that directory and save everything you
% need in a table 2) load all session tables and join them. 

% This is step one. 

%% Set Directory
clear; clc;
rng(pi)
cd 'D:\M433\M433_2023_09_19_recording1'; 
file_name = 'M453_2024_01_20'; 

%% Load Files
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
time = FP.tvec; % fiber time 
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
for i = 1:1:size(SWR_ind_mid)
    SWR_time_mid(i) = lfp_time(round(SWR_ind_mid(i))); % lfp time in (s) for a swr
    SWR_fiber_ind(i) = nearest_idx3(SWR_time_mid(i),time); % find the corresponidng time in fibFPer and saves the index.
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
filename = append(file_name, "swr_des.mat");
save(filename, '-struct','swr_des')

%% start matrix
% each swr is a row so it is length(swr_indx, x) 
matrix_sess = array2table(zeros(length(SWR_ind_mid),13),'VariableNames',{'test1','test2','test3','test4','test5','test6','test7','test8','test9','test10','test11','test12','test13'});
% working on this here...
test_struct = {[1,2,3,4,5,6,7,8,98,9]}; 

matrix_sess.test2 = test_struct;

%% Extract fiber after swrs 
prepros_signal = [];
prepros_signal = FP.zF_win_60s; 

SWR_ind_mid_post = (SWR_ind_start_post + SWR_ind_end_post)/2;  %middle index 

post_SWR_ind = SWR_ind_mid(SWR_ind_mid > SWR_ind_mid(round(SWR_ind_mid_post))); %all middle SWR timepoints

% Pre track rest ^.^ 

pre_SWR_ind = SWR_ind_mid(SWR_ind_mid < SWR_ind_mid(round(SWR_ind_mid_post))); %all middle SWR timepoints

% for fiber after swr
SWR_time_mid = zeros(length(SWR_ind_mid),1); 
SWR_fiber_ind = zeros(length(SWR_ind_mid),1); 
for i = 1:1:size(SWR_ind_mid)
    SWR_time_mid(i) = lfp_time(round(SWR_ind_mid(i)));
    SWR_fiber_ind(i) = nearest_idx3(SWR_time_mid(i),time);
end

% output: avg_circ_zdF_extract_pre, avg_fiber_pre  & post
seconds = 8; % edit this
samples = (seconds*FP.cfg.hdr{1,1}.SamplingFrequency)/2; %(seconds*1000)/2;  % divide by two because you want 4s + and - directions 

%% post 
% keep all SWR after that time -- should be 654
post_SWR_ind = SWR_ind_mid(SWR_ind_mid > SWR_ind_mid(round(SWR_ind_mid_post))); %all middle SWR timepoints

% find lfp time index and fiber signal index 
SWR_time_mid_post = zeros(length(post_SWR_ind),1); 
SWR_fiber_ind_post = zeros(length(post_SWR_ind),1); 
for i = 1:1:size(post_SWR_ind)
    SWR_time_mid_post(i) = lfp_time(round(post_SWR_ind(i)));
    SWR_fiber_ind_post(i) = nearest_idx3(SWR_time_mid_post(i),time);
end

% PETH -----------------------------------------------------------------
% z_10s - detrended, denoised/filtered, normalized (z-scored)
zdF_extract_post = zeros(length(post_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
time_extract_post =  zeros(length(post_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for ievt = 1:1:length(post_SWR_ind)-1 
    timeset = time((SWR_fiber_ind_post(ievt)-samples):(SWR_fiber_ind_post(ievt)+samples)); % pick fiber events that are 4 seconds each way
    time_extract_post(ievt,:) = time((SWR_fiber_ind_post(ievt)-samples):(SWR_fiber_ind_post(ievt)+samples))-timeset(1); 
    zdF_extract_post(ievt,:) = (prepros_signal((SWR_fiber_ind_post(ievt)-samples):(SWR_fiber_ind_post(ievt)+samples)));
end

avg_fiber_post = nanmean(zdF_extract_post);
std_fiber_post = 2*std(zdF_extract_post);

X = prepros_signal; 
N = 1000; % number of circshifts 
K=randi([1 length(prepros_signal)],1, N); % pick a random number between 1 and number of samples ... 100 times 
events_num = length(post_SWR_ind);
% initialize 
circ_zdF_extract = zeros(events_num, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
avg_circ_zdF_extract_post = zeros(N,seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for iter_circ = 1:1:N % 1 through number N (for 1000 circshifts) 
    Y = circshift(X,K(iter_circ)); % circshift the entire fiber signal based on the random number 
    for ievt = 1:1:length(post_SWR_ind)-1 % for each SWR event pick out 1-978, pick out that subset of the fiber signal 
    % find time for x axis
       timeset = time((SWR_fiber_ind_post(ievt)-samples):(SWR_fiber_ind_post(ievt)+samples)); % pick fiber events that are 4 seconds each way
       circ_zdF_extract(ievt,:) = (Y((SWR_fiber_ind_post(ievt)-samples):(SWR_fiber_ind_post(ievt)+samples))); %resets every time
    end
    avg_circ_zdF_extract_post(iter_circ,:) = nanmean(circ_zdF_extract); % this is a mean over all the trials .... 
end

circ_avg_fiber_post = nanmean(avg_circ_zdF_extract_post);
circ_std_fiber_post = 2*std(avg_circ_zdF_extract_post);

%% Pre 
% Pre-sleep session only plot 
pre_SWR_ind = SWR_ind_mid(SWR_ind_mid < SWR_ind_mid(round(SWR_ind_mid_post))); %all middle SWR timepoints

% find lfp time index and fiber signal index 
SWR_time_mid_pre = zeros(length(pre_SWR_ind),1); 
SWR_fiber_ind_pre = zeros(length(pre_SWR_ind),1); 
for i = 1:1:size(pre_SWR_ind)
    SWR_time_mid_pre(i) = lfp_time(round(pre_SWR_ind(i)));
    SWR_fiber_ind_pre(i) = nearest_idx3(SWR_time_mid_pre(i),time);
end

% PETH ----------------------------------------------------------------_
zdF_extract_pre = zeros(length(pre_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
time_extract_pre =  zeros(length(pre_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for ievt = 1:1:length(pre_SWR_ind)-1
    % find time for x axis
    timeset = time((SWR_fiber_ind_pre(ievt)-samples):(SWR_fiber_ind_pre(ievt)+samples)); % pick fiber events that are 4 seconds each way
    time_extract_pre(ievt,:) = time((SWR_fiber_ind_pre(ievt)-samples):(SWR_fiber_ind_pre(ievt)+samples))-timeset(1); 
    zdF_extract_pre(ievt,:) = (prepros_signal((SWR_fiber_ind_pre(ievt)-samples):(SWR_fiber_ind_pre(ievt)+samples)));
end
% last row is all zeros... 

avg_fiber_pre = mean(zdF_extract_pre);
%SEM_fiber = std([avg_fiber],0,1)/sqrt(length([]));
std_fiber_pre = 2*std(zdF_extract_pre);

std_top_pre = avg_fiber_pre + 2*std(zdF_extract_pre);
std_bot_pre = avg_fiber_pre - 2*std(zdF_extract_pre);

% CIRCSHIFT ------------------------------------------------------------
% elements in the array X by K positions. shifts by [m,n] n dimension. 
X = prepros_signal; 
K=randi([1 length(prepros_signal)],1, N); % pick a random number between 1 and number of samples ... 100 times 
events_num = length(pre_SWR_ind);
% initialize 
circ_zdF_extract = zeros(events_num, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
avg_circ_zdF_extract_pre = zeros(N,seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for iter_circ = 1:1:N % 1 through number N
    Y = circshift(X,K(iter_circ)); % circshift the entire fiber signal based on the random number 
    for ievt = 1:1:length(pre_SWR_ind)-1 % for each SWR event pick out 1-978, pick out that subset of the fiber signal 
    % find time for x axis
       timeset = time((SWR_fiber_ind_pre(ievt)-samples):(SWR_fiber_ind_pre(ievt)+samples)); % pick fiber events that are 4 seconds each way
       circ_zdF_extract(ievt,:) = (Y((SWR_fiber_ind_pre(ievt)-samples):(SWR_fiber_ind_pre(ievt)+samples))); %resets every time
    end
    avg_circ_zdF_extract_pre(iter_circ,:) = nanmean(circ_zdF_extract);
end

circ_avg_fiber_pre = mean(avg_circ_zdF_extract_pre);
%SEM_fiber = std([avg_fiber],0,1)/sqrt(length([]));
circ_std_fiber_pre = 2*std(avg_circ_zdF_extract_pre);
std_top_pre_circ = circ_avg_fiber_pre + 2*std(avg_circ_zdF_extract_pre);
std_bot_pre_circ = circ_avg_fiber_pre - 2*std(avg_circ_zdF_extract_pre);

%% Extract fiber after random points 

% CHANGED THIS to exp instead of prerecord
lfp_pre = restrict(csc,ExpKeys.prerecord(1), ExpKeys.prerecord(2)); % restricted lfps
lfp_post = restrict(csc,ExpKeys.postrecord(1),ExpKeys.postrecord(2)); 

% extract the same number of random time points as there were SWRs
samples_lfp = FS*4; 
index_pre = randi([samples_lfp,length(lfp_pre.data)-samples_lfp],1,length(pre_SWR_ind)); 
index_post = randi([samples_lfp,length(lfp_post.data)-samples_lfp],1,length(post_SWR_ind)); 
% pick out an index between the 4s from the start and end of the lfp signal


% find lfp time index and fiber signal index 
SWR_time_mid_pre_con = zeros(length(pre_SWR_ind),1); 
SWR_fiber_ind_pre_con = zeros(length(pre_SWR_ind),1); 
%lfp_pre.tvec_init = lfp_pre.tvec - lfp_pre.tvec(1);
for i = 1:1:size(pre_SWR_ind)
    SWR_time_mid_pre_con(i) = lfp_pre.tvec(index_pre(i))-lfp_pre.tvec(1) ; % subtracted time ... % find time corresponding to index 
    % lfp_time will be different for post since it is initialized 
    SWR_fiber_ind_pre_con(i) = nearest_idx3(SWR_time_mid_pre_con(i),time); % closest fiber time to swr time...
end

% PETH ----------------------------------------------------------------_
zdF_extract_pre_con = zeros(length(pre_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
time_extract_pre_con =  zeros(length(pre_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for ievt = 1:1:length(pre_SWR_ind)-1
    % find time for x axis
    timeset_con = time((SWR_fiber_ind_pre_con(ievt)-samples):(SWR_fiber_ind_pre_con(ievt)+samples)); % pick fiber events that are 4 seconds each way
    time_extract_pre_con(ievt,:) = time((SWR_fiber_ind_pre_con(ievt)-samples):(SWR_fiber_ind_pre_con(ievt)+samples))-timeset(1); 
    zdF_extract_pre_con(ievt,:) = (prepros_signal((SWR_fiber_ind_pre_con(ievt)-samples):(SWR_fiber_ind_pre_con(ievt)+samples)));
end


% PETH post  ---------------------------------------------------------------
% find lfp time index and fiber signal index 
SWR_time_mid_post_con = zeros(length(post_SWR_ind),1); 
SWR_fiber_ind_post_con = zeros(length(post_SWR_ind),1); 
for i = 1:1:size(post_SWR_ind)
    SWR_time_mid_post_con(i) = lfp_post.tvec(round(index_post(i)))-lfp_post.tvec(1); % find time 
    SWR_fiber_ind_post_con(i) = nearest_idx3(SWR_time_mid_post_con(i),time); % closest fiber time to swr time...
end

zdF_extract_post_con = zeros(length(post_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
time_extract_post_con =  zeros(length(post_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for ievt = 1:1:length(post_SWR_ind)
    % find time for x axis
    if SWR_fiber_ind_post_con(ievt)+samples <= length(time) % if ripple has sufficient samples (not towards the boarder of the recording
        timeset_con = time((SWR_fiber_ind_post_con(ievt)-samples):(SWR_fiber_ind_post_con(ievt)+samples)); % pick fiber events that are 4 seconds each way
        time_extract_post_con(ievt,:) = time((SWR_fiber_ind_post_con(ievt)-samples):(SWR_fiber_ind_post_con(ievt)+samples))-timeset(1);
        zdF_extract_post_con(ievt,:) = (prepros_signal((SWR_fiber_ind_post_con(ievt)-samples):(SWR_fiber_ind_post_con(ievt)+samples)));
    end
end
