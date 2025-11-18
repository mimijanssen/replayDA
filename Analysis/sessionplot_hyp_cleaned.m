%% SESSION PLOT : Hypothesis-Driven Plots 
% input: preprocessed fiber data, detected SWRs, good LFP csc
% output: session plot with: 
% ------ fiber signal after swrs (ask Wolford about this)
% ------ heat map of all trials fiber signal - 1SD shuffle

%% Load data 
% load fiber data 
clear; clc; 
%rng(pi);

cd 'D:\M595\M595_2025_10_05_recording14';
FP = load('M595_2025_10_05processed.mat');

% extract SWR intervals
load('M595_2025_10_05_recording14-manualIV.mat')

file_name = 'M595_2025_10_05'; 

addpath('C:\Users\mimia\Documents\Toolboxes\shadedErrorBar')

%%

% extract events 
LoadExpKeys
cfg_evt = [];
evt2 = LoadEvents(cfg_evt);

% extract LFP 
csc_name = [];
csc_name.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(csc_name);
time = FP.tvec; % fiber time 


% initialize LFP 
lfp_time = csc.tvec- csc.tvec(1);
lfp = csc.data;
%time = t; 

FS = csc.cfg.hdr{1}.SamplingFrequency; % set FP_data.acq.Fs to sampling frequency rate (5000 points per second) 

% time of post recording 
post = ExpKeys.postrecord(1) - csc.tvec(1); % time of post sleep period, initialized  
pre_end = ExpKeys.prerecord(2) - csc.tvec(1); % time of post sleep period, initialized  

% initialize SWR interval
SWR_start = evt.tstart- csc.tvec(1);
SWR_end = evt.tend- csc.tvec(1);
SWR_iv = [SWR_start SWR_end];
SWR_ind_start = nearest_idx3(SWR_iv(:,1),lfp_time);
SWR_ind_end = nearest_idx3(SWR_iv(:,2),lfp_time); 

SWR_ind_mid = (SWR_ind_start + SWR_ind_end)/2;  %middle LFP timepoint for each SWR

% find the SWR start and end time closest to the time of post sleep
% session. 
SWR_ind_start_post = nearest_idx3(post,SWR_iv(:,1)); 
SWR_ind_end_post = nearest_idx3(post,SWR_iv(:,2)); 

SWR_ind_start_pre_end = nearest_idx3(pre_end,SWR_iv(:,1)); 
SWR_ind_end_pre_end = nearest_idx3(pre_end,SWR_iv(:,2)); 

%% Extract fiber after swrs 
prepros_signal = [];
prepros_signal = FP.zF_win_60s; 

% Middle of the post SWR
SWR_ind_mid_post = (SWR_ind_start_post + SWR_ind_end_post)/2;  %middle index 
% keep all SWR after that time

% Pre track rest (anything before end of pretrack rest) 
SWR_ind_mid_pre_end = (SWR_ind_start_pre_end + SWR_ind_end_pre_end)/2;  %middle index 

% for fiber after swr
SWR_time_mid = zeros(length(SWR_ind_mid),1); 
SWR_fiber_ind = zeros(length(SWR_ind_mid),1); 
for i = 1:1:size(SWR_ind_mid) 
    SWR_time_mid(i) = lfp_time(round(SWR_ind_mid(i)));
    SWR_fiber_ind(i) = nearest_idx3(SWR_time_mid(i),time);
end

seconds = 8; % edit this
samples = (seconds*FP.cfg.hdr{1,1}.SamplingFrequency)/2;

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
zdF_extract_post = zeros(length(post_SWR_ind)-1, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
time_extract_post =  zeros(length(post_SWR_ind)-1, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for ievt = 1:1:length(post_SWR_ind)-1 %-3 for M600 recording 1; -1 for everyone else. 
    timeset = time((SWR_fiber_ind_post(ievt)-samples):(SWR_fiber_ind_post(ievt)+samples)); % pick fiber events that are 4 seconds each way
    time_extract_post(ievt,:) = time((SWR_fiber_ind_post(ievt)-samples):(SWR_fiber_ind_post(ievt)+samples))-timeset(1); 
    zdF_extract_post(ievt,:) = (prepros_signal((SWR_fiber_ind_post(ievt)-samples):(SWR_fiber_ind_post(ievt)+samples)));
end

avg_fiber_post = nanmean(zdF_extract_post);
std_fiber_post = 2*std(zdF_extract_post);

X = prepros_signal; 
N = 1000; % number of circshifts 
K=randi([1 length(prepros_signal)],1, N); % pick a random number between 1 and number of samples ... 100 times 
events_num = length(post_SWR_ind)-1;
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
pre_SWR_ind = SWR_ind_mid(SWR_ind_mid < SWR_ind_mid(round(SWR_ind_mid_pre_end))); %all middle SWR timepoints

% find lfp time index and fiber signal index 
SWR_time_mid_pre = zeros(length(pre_SWR_ind),1); 
SWR_fiber_ind_pre = zeros(length(pre_SWR_ind),1); 
for i = 1:1:size(pre_SWR_ind)
    SWR_time_mid_pre(i) = lfp_time(round(pre_SWR_ind(i)));
    SWR_fiber_ind_pre(i) = nearest_idx3(SWR_time_mid_pre(i),time);
end

% PETH ----------------------------------------------------------------_
zdF_extract_pre = zeros(length(pre_SWR_ind)-1, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
time_extract_pre =  zeros(length(pre_SWR_ind)-1, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for ievt = 1:1:length(pre_SWR_ind)-1 
    timeset = time((SWR_fiber_ind_pre(ievt)-samples):(SWR_fiber_ind_pre(ievt)+samples)); % pick fiber events that are 4 seconds each way
    time_extract_pre(ievt,:) = time((SWR_fiber_ind_pre(ievt)-samples):(SWR_fiber_ind_pre(ievt)+samples))-timeset(1); 
    zdF_extract_pre(ievt,:) = (prepros_signal((SWR_fiber_ind_pre(ievt)-samples):(SWR_fiber_ind_pre(ievt)+samples)));
end

avg_fiber_pre = mean(zdF_extract_pre);
std_fiber_pre = 2*std(zdF_extract_pre);

std_top_pre = avg_fiber_pre + 2*std(zdF_extract_pre);
std_bot_pre = avg_fiber_pre - 2*std(zdF_extract_pre);

% CIRCSHIFT ------------------------------------------------------------
% elements in the array X by K positions. shifts by [m,n] n dimension. 
X = prepros_signal; 
K=randi([1 length(prepros_signal)],1, N); % pick a random number between 1 and number of samples ... 100 times 
events_num = length(pre_SWR_ind)-1;
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

%% Track 
track_SWR_ind = SWR_ind_mid(SWR_ind_mid_pre_end:SWR_ind_mid_post); %all middle SWR timepoints

% find lfp time index and fiber signal index 
SWR_time_mid_track = zeros(length(track_SWR_ind),1); 
SWR_fiber_ind_track = zeros(length(track_SWR_ind),1); 
for i = 1:1:size(track_SWR_ind)
    SWR_time_mid_track(i) = lfp_time(round(track_SWR_ind(i)));
    SWR_fiber_ind_track(i) = nearest_idx3(SWR_time_mid_track(i),time);
end

% PETH ----------------------------------------------------------------_
zdF_extract_track = zeros(length(track_SWR_ind)-1, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
time_extract_track =  zeros(length(track_SWR_ind)-1, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for ievt = 1:1:length(track_SWR_ind)-1 
    timeset = time((SWR_fiber_ind_track(ievt)-samples):(SWR_fiber_ind_track(ievt)+samples)); % pick fiber events that are 4 seconds each way
    time_extract_track(ievt,:) = time((SWR_fiber_ind_track(ievt)-samples):(SWR_fiber_ind_track(ievt)+samples))-timeset(1); 
    zdF_extract_track(ievt,:) = (prepros_signal((SWR_fiber_ind_track(ievt)-samples):(SWR_fiber_ind_track(ievt)+samples)));
end

avg_fiber_track = mean(zdF_extract_track);
std_fiber_track = 2*std(zdF_extract_track);

std_top_track = avg_fiber_track + 2*std(zdF_extract_track);
std_bot_track = avg_fiber_track - 2*std(zdF_extract_track);

% CIRCSHIFT ------------------------------------------------------------
% elements in the array X by K positions. shifts by [m,n] n dimension. 
X = prepros_signal; 
K=randi([1 length(prepros_signal)],1, N); % pick a random number between 1 and number of samples ... 100 times 
events_num = length(track_SWR_ind)-1;
% initialize 
circ_zdF_extract = zeros(events_num, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
avg_circ_zdF_extract_track = zeros(N,seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for iter_circ = 1:1:N % 1 through number N
    Y = circshift(X,K(iter_circ)); % circshift the entire fiber signal based on the random number 
    for ievt = 1:1:length(track_SWR_ind)-1 % for each SWR event pick out 1-978, pick out that subset of the fiber signal 
    % find time for x axis
       timeset = time((SWR_fiber_ind_track(ievt)-samples):(SWR_fiber_ind_track(ievt)+samples)); % pick fiber events that are 4 seconds each way
       circ_zdF_extract(ievt,:) = (Y((SWR_fiber_ind_track(ievt)-samples):(SWR_fiber_ind_track(ievt)+samples))); %resets every time
    end
    avg_circ_zdF_extract_track(iter_circ,:) = nanmean(circ_zdF_extract);
end

circ_avg_fiber_track = mean(avg_circ_zdF_extract_track);
%SEM_fiber = std([avg_fiber],0,1)/sqrt(length([]));
circ_std_fiber_track = 2*std(avg_circ_zdF_extract_track);
std_top_track_circ = circ_avg_fiber_track + 2*std(avg_circ_zdF_extract_track);
std_bot_track_circ = circ_avg_fiber_track - 2*std(avg_circ_zdF_extract_track);

%% HEAT MAP: Extract fiber after random points or SWR points

% restrict lfp to time intervals 
prerecord_init = ExpKeys.prerecord(1)-csc.tvec(1);
prerecord_end = ExpKeys.prerecord(2)-csc.tvec(1);

postrecord_init = ExpKeys.postrecord(1)-csc.tvec(1);
postrecord_end = ExpKeys.postrecord(2)-csc.tvec(1);

track_init = prerecord_end;
track_end = postrecord_init;

lfp_pre = restrict(csc,ExpKeys.prerecord(1), ExpKeys.prerecord(2)); 
lfp_post = restrict(csc,ExpKeys.postrecord(1),ExpKeys.postrecord(2));
lfp_track = restrict(csc,ExpKeys.prerecord(2),ExpKeys.postrecord(1));

% extract the same number of random time points as there were SWRs
samples_lfp = FS*4; 
index_pre = randi([samples_lfp,length(lfp_pre.data)-samples_lfp],1,length(pre_SWR_ind)); 
index_post = randi([samples_lfp,length(lfp_post.data)-samples_lfp],1,length(post_SWR_ind)); 
index_track = randi([samples_lfp,length(lfp_track.data)-samples_lfp],1,length(track_SWR_ind)); 

% PRE:
SWR_time_mid_pre_con = zeros(length(pre_SWR_ind),1); 
SWR_fiber_ind_pre_con = zeros(length(pre_SWR_ind),1); 
for i = 1:1:size(pre_SWR_ind)
    SWR_time_mid_pre_con(i) = lfp_pre.tvec(index_pre(i))-lfp_pre.tvec(1) ; 
    SWR_fiber_ind_pre_con(i) = nearest_idx3(SWR_time_mid_pre_con(i),time);
end

zdF_extract_pre_con = zeros(length(pre_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
time_extract_pre_con =  zeros(length(pre_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for ievt = 1:1:length(pre_SWR_ind)-1
    timeset_con = time((SWR_fiber_ind_pre_con(ievt)-samples):(SWR_fiber_ind_pre_con(ievt)+samples)); % pick fiber events that are 4 seconds each way
    time_extract_pre_con(ievt,:) = time((SWR_fiber_ind_pre_con(ievt)-samples):(SWR_fiber_ind_pre_con(ievt)+samples))-timeset(1); 
    zdF_extract_pre_con(ievt,:) = (prepros_signal((SWR_fiber_ind_pre_con(ievt)-samples):(SWR_fiber_ind_pre_con(ievt)+samples)));
end


% POST: 
SWR_time_mid_post_con = zeros(length(post_SWR_ind),1); 
SWR_fiber_ind_post_con = zeros(length(post_SWR_ind),1); 
for i = 1:1:size(post_SWR_ind)
    SWR_time_mid_post_con(i) = lfp_post.tvec(round(index_post(i)))-lfp_post.tvec(1); % find time 
    SWR_fiber_ind_post_con(i) = nearest_idx3(SWR_time_mid_post_con(i),time); % closest fiber time to swr time...
end

zdF_extract_post_con = zeros(length(post_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
time_extract_post_con =  zeros(length(post_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for ievt = 1:1:length(post_SWR_ind)
    if SWR_fiber_ind_post_con(ievt)+samples <= length(time) % if ripple has sufficient samples (not towards the boarder of the recording
        timeset_con = time((SWR_fiber_ind_post_con(ievt)-samples):(SWR_fiber_ind_post_con(ievt)+samples)); % pick fiber events that are 4 seconds each way
        time_extract_post_con(ievt,:) = time((SWR_fiber_ind_post_con(ievt)-samples):(SWR_fiber_ind_post_con(ievt)+samples))-timeset(1);
        zdF_extract_post_con(ievt,:) = (prepros_signal((SWR_fiber_ind_post_con(ievt)-samples):(SWR_fiber_ind_post_con(ievt)+samples)));
    end
end

% TRACK: 
SWR_time_mid_track_con = zeros(length(track_SWR_ind),1); 
SWR_fiber_ind_track_con = zeros(length(track_SWR_ind),1); 
for i = 1:1:size(track_SWR_ind)
    SWR_time_mid_track_con(i) = lfp_track.tvec(round(index_track(i)))-lfp_track.tvec(1); % find time 
    SWR_fiber_ind_track_con(i) = nearest_idx3(SWR_time_mid_track_con(i),time); % closest fiber time to swr time...
end

zdF_extract_track_con = zeros(length(track_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
time_extract_track_con =  zeros(length(track_SWR_ind), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for ievt = 1:1:length(track_SWR_ind)
    if SWR_fiber_ind_track_con(ievt)+samples <= length(time) % if ripple has sufficient samples (not towards the boarder of the recording
        timeset_con = time((SWR_fiber_ind_track_con(ievt)-samples):(SWR_fiber_ind_track_con(ievt)+samples)); % pick fiber events that are 4 seconds each way
        time_extract_track_con(ievt,:) = time((SWR_fiber_ind_track_con(ievt)-samples):(SWR_fiber_ind_track_con(ievt)+samples))-timeset(1);
        zdF_extract_track_con(ievt,:) = (prepros_signal((SWR_fiber_ind_track_con(ievt)-samples):(SWR_fiber_ind_track_con(ievt)+samples)));
    end
end

%% 
pre_heatmap = zdF_extract_pre - circ_std_fiber_pre; 
post_heatmap = zdF_extract_post - circ_std_fiber_post; 
track_heatmap = zdF_extract_track - circ_std_fiber_track; 


updated_map_pre = zeros(size(pre_heatmap));
updated_map_post = zeros(size(post_heatmap));
updated_map_track = zeros(size(track_heatmap));


for iter1 = 1:1: length(circ_std_fiber_pre)
    for iter_rows = 1:1:length(pre_SWR_ind)-1
        if zdF_extract_pre_con(iter_rows,iter1) > 0
            updated_map_pre(iter_rows,iter1) = zdF_extract_pre(iter_rows,iter1) - circ_std_fiber_pre(iter1);
        else
            updated_map_pre(iter_rows,iter1) = zdF_extract_pre(iter_rows,iter1) + circ_std_fiber_pre(iter1);
        end
    end
    for iter_rows2 = 1:1:length(post_SWR_ind)-1
        if zdF_extract_post(iter_rows2,iter1) > 0
            updated_map_post(iter_rows2,iter1) = zdF_extract_post(iter_rows2,iter1) - circ_std_fiber_post(iter1);
        else
            updated_map_post(iter_rows2,iter1) = zdF_extract_post(iter_rows2,iter1) + circ_std_fiber_post(iter1);
        end
    end
    for iter_rows3 = 1:1:length(track_SWR_ind)-1
        if zdF_extract_track(iter_rows3,iter1) > 0
            updated_map_track(iter_rows3,iter1) = zdF_extract_track(iter_rows3,iter1) - circ_std_fiber_track(iter1);
        else
            updated_map_track(iter_rows3,iter1) = zdF_extract_track(iter_rows3,iter1) + circ_std_fiber_track(iter1);
        end
    end
end


pre_heatmap_sorted = sortrows(pre_heatmap, 7201, 'descend');
post_heatmap_sorted = sortrows(post_heatmap, 7201, 'descend');
track_heatmap_sorted = sortrows(track_heatmap, 7201, 'descend');

updated_pre_heatmap_sorted = sortrows(updated_map_pre, 7201, 'descend');
updated_post_heatmap_sorted = sortrows(updated_map_post, 7201, 'descend');
updated_track_heatmap_sorted = sortrows(updated_map_track, 7201, 'descend');

%% shuffled heatmaps 
% getting a constant value here... 
pre_heatmap_con = zdF_extract_pre_con  - circ_std_fiber_pre; % need to edit this 
updated_map_pre_con = zeros(size(pre_heatmap_con));

post_heatmap_con = zdF_extract_post_con - circ_std_fiber_post; 
updated_map_post_con = zeros(size(post_heatmap_con));

track_heatmap_con = zdF_extract_track_con - circ_std_fiber_track; 
updated_map_track_con = zeros(size(track_heatmap_con));

for iter1 = 1:1: length(circ_std_fiber_pre)
    for iter_rows = 1:1:length(pre_SWR_ind)-1
        if zdF_extract_pre_con(iter_rows,iter1) > 0
            updated_map_pre_con(iter_rows,iter1) = zdF_extract_pre_con(iter_rows,iter1) - circ_std_fiber_pre(iter1);
        else
            updated_map_pre_con(iter_rows,iter1) = zdF_extract_pre_con(iter_rows,iter1) + circ_std_fiber_pre(iter1);
        end
    end
    for iter_rows2 = 1:1:length(post_SWR_ind)-1
        if zdF_extract_post_con(iter_rows2,iter1) > 0
            updated_map_post_con(iter_rows2,iter1) = zdF_extract_post_con(iter_rows2,iter1) - circ_std_fiber_post(iter1);
        else
            updated_map_post_con(iter_rows2,iter1) = zdF_extract_post_con(iter_rows2,iter1) + circ_std_fiber_post(iter1);
        end
    end
    for iter_rows3 = 1:1:length(track_SWR_ind)-1
        if zdF_extract_track_con(iter_rows3,iter1) > 0
            updated_map_track_con(iter_rows3,iter1) = zdF_extract_track_con(iter_rows3,iter1) - circ_std_fiber_track(iter1);
        else
            updated_map_track_con(iter_rows3,iter1) = zdF_extract_track_con(iter_rows3,iter1) + circ_std_fiber_track(iter1);
        end
    end
end

pre_heatmap_sorted_con = sortrows(pre_heatmap, 7201, 'descend');
post_heatmap_sorted_con = sortrows(post_heatmap, 7201, 'descend');
track_heatmap_sorted_con = sortrows(track_heatmap, 7201, 'descend');

updated_pre_heatmap_sorted_con = sortrows(updated_map_pre_con, 7201, 'descend');
updated_post_heatmap_sorted_con = sortrows(updated_map_post_con, 7201, 'descend');
updated_track_heatmap_sorted_con = sortrows(updated_map_track_con, 7201, 'descend');


%% Figures together 
fig2 = figure(2);

med_c = [104,187,225]./255; % rgb(167, 199, 231) rgb(255, 165, 0)  blue: rgb(104,187,227)
low_c = [78,178,101]./255; % color

subplot(3,4,1:2)
shadedErrorBar(time_extract_pre(1,:),circ_avg_fiber_pre,circ_std_fiber_pre,'lineProps','-k','transparent',1) % subtract the circ mean here 
hold on
plot(time_extract_pre(1,:),circ_avg_fiber_pre,'LineWidth',2,'Color','k') % subtract the circ mean here 
% plot average on top with larger line
hold on
plot(time_extract_pre(1,:),avg_fiber_pre,'LineWidth',2,'Color',low_c)
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
ylim([-1.5 1.5])
xlim([0 8])
xticks([0 4 8])
xticklabels({'-4','0','4'})
title('Pre Track Rest PETH','FontSize', 20)
ylabel('Mean [DA] (z-score)','FontSize', 16)
xlabel('Time from SWR (s)','FontSize', 16)
legend('','shuffle','signal','Location','northwest')
legend boxoff

subplot(3,4,9:10)
shadedErrorBar(time_extract_post(1,:),circ_avg_fiber_post,circ_std_fiber_post,'lineProps','-k','transparent',1) % subtracted the circ mean here
hold on
plot(time_extract_post(1,:),circ_avg_fiber_post,'LineWidth',2,'Color','k') %subtracted the circ mean here
hold on
plot(time_extract_post(1,:),avg_fiber_post,'LineWidth',2,'Color',low_c)
hold on
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
ylim([-1.5 1.5])
xticks([0 4 8])
xticklabels({'-4','0','4'})
title('Post Track Rest PETH','FontSize', 20)
ylabel('Mean [DA] (z-score)','FontSize', 16)
xlabel('Time from SWR (s)','FontSize', 16)
legend('','shuffle','signal','Location','northwest')
legend boxoff 

subplot(3,4,5:6)
shadedErrorBar(time_extract_track(1,:),circ_avg_fiber_track,circ_std_fiber_track,'lineProps','-k','transparent',1) % subtract the circ mean here 
hold on
plot(time_extract_track(1,:),circ_avg_fiber_track,'LineWidth',2,'Color','k') % subtract the circ mean here 
% plot average on top with larger line
hold on
plot(time_extract_track(1,:),avg_fiber_track,'LineWidth',2,'Color',low_c)
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
ylim([-1.5 1.5])
xlim([0 8])
xticks([0 4 8])
xticklabels({'-4','0','4'})
title('On Track PETH','FontSize', 20)
ylabel('Mean [DA] (z-score)','FontSize', 16)
xlabel('Time from SWR (s)','FontSize', 16)
legend('','shuffle','signal','Location','northwest')
legend boxoff

subplot(3,4,3)
imagesc(updated_pre_heatmap_sorted) 
colorbar
xticks([0 6400 12800])
xticklabels({'-4','0','4'})
xlabel('Time from SWR (s)')
ylabel('zdF - 1 sd of shuffle')
title('Pre-task fiber after SWR')

subplot(3,4,11)
imagesc(updated_post_heatmap_sorted) 
colorbar
xticks([0 6400 12800])
xticklabels({'-4','0','4'})
xlabel('Time from SWR (s)')
ylabel('zdF - 1 sd of shuffle')
title('Post-task fiber after SWR')

subplot(3,4,7)
imagesc(updated_track_heatmap_sorted) 
colorbar
xticks([0 6400 12800])
xticklabels({'-4','0','4'})
xlabel('Time from SWR (s)')
ylabel('zdF - 1 sd of shuffle')
title('On-track fiber after SWR')

subplot(3,4,4)
imagesc(updated_pre_heatmap_sorted_con) 
colorbar
xticks([0 6400 12800])
xticklabels({'-4','0','4'})
xlabel('Time from random event (s)')
ylabel('zdF - 1 sd of shuffle')
title('Pre-task fiber after random time-point')

subplot(3,4,12)
imagesc(updated_post_heatmap_sorted_con) 
colorbar
xticks([0 6400 12800])
xticklabels({'-4','0','4'})
xlabel('Time from random event (s)')
ylabel('zdF - 1 sd of shuffle')
title('Post-task fiber after random time-point')

subplot(3,4,8)
imagesc(updated_track_heatmap_sorted_con) 
colorbar
xticks([0 6400 12800])
xticklabels({'-4','0','4'})
xlabel('Time from random event (s)')
ylabel('zdF - 1 sd of shuffle')
title('On-task fiber after random time-point')

set(gcf,'Color',[1,1,1])
shg

hold off

txt = {'Session Plot: Hypothesis Suplots'};
sgtitle(txt)

fig2.WindowState = 'maximized';
 %% Save figure

cd 'C:\Users\mimia\Documents\ReplayDA Figures\M595'
saveas(fig2,'M595_ontrack14_hypothesis.png') % CHANGE THIS 

avg_SWR_DA.circ_avg_pre = circ_avg_fiber_pre;
avg_SWR_DA.circ_avg_post = circ_avg_fiber_post;
avg_SWR_DA.circ_avg_track = circ_avg_fiber_track;

avg_SWR_DA.circ_std_pre = circ_std_fiber_pre;
avg_SWR_DA.circ_std_post = circ_std_fiber_post;
avg_SWR_DA.circ_std_track = circ_std_fiber_track;

avg_SWR_DA.avg_fiber_pre = avg_fiber_pre;
avg_SWR_DA.avg_fiber_post = avg_fiber_post;
avg_SWR_DA.avg_fiber_track = avg_fiber_track;

avg_SWR_DA.time = time_extract_pre;

%% Save data
cd 'D:\M595\'
filename = append(file_name, "avg.mat");
save(filename, '-struct','avg_SWR_DA')

