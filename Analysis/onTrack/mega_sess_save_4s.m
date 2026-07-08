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
cd 'D:\M654\M654_2026_02_08_recording8' 
file_name = 'M654_2026_02_08'; 
mouseID = ['M654'];
session = 7; 
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
%filename = append(file_name, "swr_des_basic.mat");
%save(filename, '-struct','swr_des')


%% start matrix
% each swr has it's own row 
matrix_sess = array2table(zeros(length(SWR_ind_mid),27),'VariableNames',{'mouseID','sess','swrID','PrePost','TwosPreRaw','TwosTrackRaw','TwosPostRaw','TwosPreProc','TwosTrackProc','TwosPostProc','OnesBeforePeak','OnesAfterPeak','OnesBeforeAUC','OnesAfterAUC','TimeAfterPeak','OnesBeforePeakRAW','OnesAfterPeakRAW','OnesBeforeAUCRAW','OnesAfterAUCRAW','SWRdur','SWRamp','SWRpower','SWRtimestart','SWRtimeend','SWRindmid','SWRindstart','SWRindend'});

%% input mouse/session identity information 
matrix_sess.('mouseID')(:,1) = mouse; 
matrix_sess.('sess')(:,1) = session;
%matrix_sess.('history')(:,1) = training;
matrix_sess.('swrID')(:) = linspace(0,1,length(matrix_sess.('swrID')))';
matrix_sess.('SWRindmid')(:) = SWR_ind_mid;
matrix_sess.('SWRtimestart')(:) = SWR_start;
matrix_sess.('SWRtimeend')(:) = SWR_end;
matrix_sess.('SWRdur')(:) = SWR_end(:) - SWR_start(:); 
matrix_sess.('SWRindstart')(:) = SWR_ind_start;
matrix_sess.('SWRindend')(:) = SWR_ind_end;

% IF RUNNING TRACK SESSIONS USE THIS: 
matrix_sess.('PrePost')(1:round(SWR_ind_mid_pre_end)-1,1) = 1;
% 2 for track
matrix_sess.('PrePost')(round(SWR_ind_mid_pre_end):round(SWR_ind_mid_post)-1,1) = 2;
% 3 for post. 
matrix_sess.('PrePost')(round(SWR_ind_mid_post):end,1) = 3;

%% Populate matrix with structure information: Save 4 s before and after SWR

seconds = 4;
samples = seconds * FP.cfg.hdr{1,1}.SamplingFrequency;

badSWR = false(height(matrix_sess),1);

for i = 1:height(matrix_sess)

    swr_time = lfp_time(round(matrix_sess.SWRindmid(i)));
    fiber_index = nearest_idx3(swr_time, FP.tvec);

    % Require a complete window on both sides of the SWR
    if (fiber_index - samples < 1) || ...
       (fiber_index + samples > length(FP.data))
        badSWR(i) = true;
    end
end

% Remove incomplete SWRs from the table
matrix_sess(badSWR,:) = [];

matrix_sess.TwosPreRaw   = cell(height(matrix_sess),1);
matrix_sess.TwosPostRaw  = cell(height(matrix_sess),1);
matrix_sess.TwosTrackRaw   = cell(height(matrix_sess),1);
matrix_sess.TwosPreProc  = cell(height(matrix_sess),1);
matrix_sess.TwosTrackProc  = cell(height(matrix_sess),1);
matrix_sess.TwosPostProc = cell(height(matrix_sess),1);

for i = 1:height(matrix_sess)
    swr_time = lfp_time(round(matrix_sess.SWRindmid(i)));
    fiber_index = nearest_idx3(swr_time, FP.tvec);
    idx = (fiber_index-samples):(fiber_index+samples);
    raw_data_struct = struct();
    data_struct = struct();
    raw_data_struct.signal = FP.data(idx);
    raw_data_struct.tvec   = FP.tvec(idx);
    data_struct.signal = FP.zF_win_60s(idx);
    data_struct.tvec   = raw_data_struct.tvec;
    if matrix_sess.PrePost(i) == 1
        matrix_sess.TwosPreRaw{i}  = raw_data_struct;
        matrix_sess.TwosPreProc{i} = data_struct;
    elseif matrix_sess.PrePost(i) == 2
        matrix_sess.TwosTrackRaw{i}  = raw_data_struct;
        matrix_sess.TwosTrackProc{i} = data_struct;
    elseif matrix_sess.PrePost(i) == 3
        matrix_sess.TwosPostRaw{i}  = raw_data_struct;
        matrix_sess.TwosPostProc{i} = data_struct;
    end
end

%% populate SWR details
% 18-23 :
% 'SWRdur','SWRamp','SWRpower','SWRf100ms','SWRtimestart','SWRtimeend';

samples_swr = (0.1*csc.cfg.hdr{1,1}.SamplingFrequency); % samples I will take before and after swrs

% filter the LFP band 
cfg = [];
cfg.f = [120 250]; % 140-220
cfg.display_filter = 0; 

% zscored LFP 
SWRz = zscore_tsd(csc); 

SWRf = FilterLFP(cfg,csc);

% obtain power and z-score it 
SWRp= LFPpower([],SWRf);
SWRp_z = zscore_tsd(SWRp);

% obtain amplitude and z-score it ? 
LoadMetadata % for freqs
SWRa = amSWR([],metadata.SWRfreqs,SWRf);
SWRa_z = zscore_tsd(SWRa); % should be proportional to power!

for i = 1:height(matrix_sess) % iterate through each swr. 
    % swrduration = SWR_end - SWR_start 
    % swramp = mean(SWRa.data(SWR_ind_start(i):SWR_ind_end(i)))
    matrix_sess.('SWRamp')(i) = mean(SWRa.data(matrix_sess.SWRindstart(i):matrix_sess.SWRindend(i)));
    matrix_sess.('SWRpower')(i) = mean(SWRp.data(matrix_sess.SWRindstart(i):matrix_sess.SWRindend(i)));
    alt_index = round(matrix_sess.SWRindmid(i));
    matrix_sess.('SWR100ms'){i} = SWRz.data(alt_index-samples_swr:alt_index+samples_swr);
end

%% populate dF information on from preproc data 
% dF from 2 seconds 
x1 = 1:1:FP.cfg.hdr{1,1}.SamplingFrequency; %1:1:500;%1001:1:3000; %1:1:2000;% % two seconds before for preproc data
x2 = FP.cfg.hdr{1,1}.SamplingFrequency + 1:1:2*FP.cfg.hdr{1,1}.SamplingFrequency + 1; % 1601:1:3200; %501:1:1000; %3001:1:5000; %2001:1:4000;%  % two seconds after for preproc data
% Is this right for the GFP Mice?

for i = 1:height(matrix_sess) % iterate through each swr. 
    swr_time = lfp_time(round(matrix_sess.SWRindmid(i))); % swr lfp time  
    fiber_index = nearest_idx3(swr_time, FP.tvec); % fiber index closest to middle swr_time 
    signal = FP.zF_win_60s(fiber_index-samples:fiber_index+samples); 
    signal_raw = FP.data(fiber_index-samples:fiber_index+samples); 
    matrix_sess.('OnesBeforePeak')(i) = max(signal(x1)); %-min(signal(x1)); % maybe the average signal might be better than the lowest signal?? 
    [matrix_sess.('OnesAfterPeak')(i),I] = max(signal(x2)); %-min(signal(x2)); 
    matrix_sess.('TimeAfterPeak')(i) = FP.tvec(I); % time of the peak post swr
    matrix_sess.('OnesBeforePeakRAW')(i) = max(signal_raw(x1)); %-min(signal(x1)); % maybe the average signal might be better than the lowest signal?? 
    matrix_sess.('OnesAfterPeakRAW')(i) = max(signal_raw(x2)); %-min(signal(x2)); 
end

% populate AUC information 
for i = 1:height(matrix_sess) % iterate through each swr. 
    swr_time = lfp_time(round(matrix_sess.SWRindmid(i))); % swr lfp time  
    fiber_index = nearest_idx3(swr_time, FP.tvec); % fiber index closest to middle swr_time 
    signal = FP.zF_win_60s(fiber_index-samples:fiber_index+samples); 
    signal_raw = FP.data(fiber_index-samples:fiber_index+samples); 
    matrix_sess.('OnesBeforeAUC')(i) = trapz(x1,signal(x1)); 
    matrix_sess.('OnesAfterAUC')(i) = trapz(x2,signal(x2)); 
    matrix_sess.('OnesBeforeAUCRAW')(i) = trapz(x1,signal_raw(x1)); 
    matrix_sess.('OnesAfterAUCRAW')(i) = trapz(x2,signal_raw(x2)); 
end

%% Save everything
cd 'D:\SWR_DA_MegaMatrix_4s_track'
filename = append(file_name, "mega1_4.mat");
save(filename,'matrix_sess')
   