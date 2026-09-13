%% session saving -- HUGE MATRIX for NREM and Wake 
% method: 1) for each session, go to that directory and save everything you
% need in a table 2) load all session tables and join them. 
addpath('C:\Users\mimia\Documents\GitHub\replayDA\Analysis') % my analysis code exists here
addpath('C:\Users\mimia\Documents\GitHub\vandermeerlab-replay-da\code-matlab\tasks\Alyssa_Tmaze\beta'); % SWR amplitude detection exists here
addpath('cd C:\Users\mimia\Documents\GitHub\sleep')

% https://www.sciencedirect.com/science/article/pii/S0306452202006693

%% Set Directory
% input information 
clear; clc;
rng(pi)
cd 'F:\M433\M433_2023_09_26_recording8'; 
file_name = 'M433_2023_09_26'; 
load('labels.mat'); % 2 is wake, 3 is nrem, 4 is undefined. 
mouseID = ['M433'];
session = 8; 
mouse = convertMouse(mouseID); % converted mouse number 

%% Load Files
FP_file=dir('*processed*');
FP = load(FP_file.name); % preprocessed data. 

SWR_file = dir('*detectedSWRs*');
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

seconds = 1; % one second before and after 
samples = (seconds*FP.cfg.hdr{1,1}.SamplingFrequency); % samples I will take before and after swrs
time = FP.tvec; % fiber time 

%% Determine SWR Counts for Sleep Sessions
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

% initialize event times
% sleep events 
post = ExpKeys.postrecord(1) - csc.tvec(1);
% swr events
SWR_start = evt.tstart- csc.tvec(1);
SWR_end = evt.tend- csc.tvec(1);
SWR_iv = [SWR_start SWR_end];
SWR_ind_start = nearest_idx3(SWR_iv(:,1),lfp_time); % could have just done evt.tstart and csc.time .... 
SWR_ind_end = nearest_idx3(SWR_iv(:,2),lfp_time); % index 

SWR_ind_mid = (SWR_ind_start + SWR_ind_end)/2;  %middle index for SWR

% for fiber after swr
SWR_time_mid = zeros(length(SWR_ind_mid),1); 
SWR_time_start = zeros(length(SWR_ind_mid),1); 
SWR_time_end = zeros(length(SWR_ind_mid),1); 

SWR_fiber_ind = zeros(length(SWR_ind_mid),1); 
for i = 1:1:size(SWR_ind_mid,1)
    SWR_time_mid(i) = lfp_time(round(SWR_ind_mid(i))); % lfp time in (s) for a swr
    SWR_time_start(i) = lfp_time(round(SWR_ind_start(i))); % lfp time in (s) for a swr
    SWR_time_end(i) = lfp_time(round(SWR_ind_end(i))); % lfp time in (s) for a swr
    SWR_fiber_ind(i) = nearest_idx3(SWR_time_mid(i),FP.tvec); % find the corresponidng time in fiber and saves the index. FP.tvec is in seconds 
end


%% find SWR times that are in nrem sleep stage and wake sleep stage
% ~ for every swr time, check if that is an nrem stage in labels (==3) 
% ~ for every swr time, check if that is an nrem stage in labels (==3) 

labels_time_center = (0:length(labels)-1)' * 2 + 1;

SWR_time_nrem = [];
SWR_nrem_start = []; % list of start times for these guys 
SWR_nrem_end = []; % list of end times for duration

SWR_time_wake = [];
SWR_wake_start = []; % list of start times for these guys 
SWR_wake_end = []; % list of end times for duration

sleep_labels = zeros(length(SWR_ind_mid),1);

for swri = 1:1:length(SWR_time_mid)
    label_ind = nearest_idx3(SWR_time_mid(swri),labels_time_center); % index of labels_time 
    if labels(label_ind) == 3
        SWR_time_nrem = [SWR_time_nrem; SWR_time_mid(swri)];
        SWR_nrem_start = [SWR_nrem_start; SWR_time_start(swri)];
        SWR_nrem_end = [SWR_nrem_end; SWR_time_end(swri)];
        sleep_labels(swri) = 2; % nrem is 2
    end
    if labels(label_ind) == 2
        SWR_time_wake = [SWR_time_wake; SWR_time_mid(swri)];
        SWR_wake_start = [SWR_wake_start; SWR_time_start(swri)];
        SWR_wake_end = [SWR_wake_end; SWR_time_end(swri)];
        sleep_labels(swri) = 1; % wake is 1 
    end
end

% ~ convert to fiber time index
SWR_fiber_ind_nrem = zeros(length(SWR_time_nrem),1); 
for i = 1:1:length(SWR_time_nrem)
    SWR_fiber_ind_nrem(i) = nearest_idx3(SWR_time_nrem(i),time);
end

% ~ convert to fiber time index
SWR_fiber_ind_wake = zeros(length(SWR_time_wake),1); 
for i = 1:1:length(SWR_time_wake)
    SWR_fiber_ind_wake(i) = nearest_idx3(SWR_time_wake(i),time);
end


%% save nrem and wake swr count 
wake_count = length(SWR_fiber_ind_wake);
nrem_count = length(SWR_fiber_ind_nrem);

swr_sleep.swr_count = [wake_count nrem_count];
swr_sleep.swr_label = ["wake" "nrem"]; 

%% save frequency and save avg duration in the same way as pre_count...
% ~~ Frequency ~~
% find the time for nrem and wake. Each label is 2 seconds
num_wake = sum(labels == 2)*2; % seconds
num_nrem = sum(labels == 3)*2;

freq_nrem = nrem_count/num_nrem; 
freq_wake = wake_count/num_wake; 

swr_sleep.freq = [freq_wake freq_nrem];

% ~~ Duration ~~
% in miliseconds 
wake_avg_dur = (mean(SWR_wake_end - SWR_wake_start))*1000; 
% 66 seconds (just like pre)

nrem_avg_dur = (mean(SWR_nrem_end - SWR_nrem_start))*1000; 
% 72 seconds (just like post)

swr_sleep.dur = [wake_avg_dur nrem_avg_dur];

%% start matrix (1 second matrix) -- SAVING ALL SWRS so also save the frequency of the theta phase in the lfp. also the power of the delta phase.
% each swr has it's own row 
matrix_sess = array2table(zeros(length(SWR_ind_mid),26),'VariableNames',{'mouseID','sess','swrID','sleep','TwosWakeRaw','TwosNremRaw','TwosWakeProc','TwosNremProc','OnesBeforePeak','OnesAfterPeak','OnesBeforeAUC','OnesAfterAUC','TimeAfterPeak','OnesBeforePeakRAW','OnesAfterPeakRAW','OnesBeforeAUCRAW','OnesAfterAUCRAW','SWRdur','SWRamp','SWRpower','SWR100ms','SWRtimestart','SWRtimeend','thetap','deltap','swrp'});


% maybe i want to add rpe strength to this too.... high - med rpe.
% and track speed. 
%% input mouse/session identity information 
matrix_sess.('mouseID')(:,1) = mouse; 
matrix_sess.('sess')(:,1) = session;
matrix_sess.('swrID')(:) = linspace(0,1,length(matrix_sess.('swrID')))';
matrix_sess.('sleep')(:)= sleep_labels;
% 1 for pre
% everything before SWR_ind_mid_post index gets a 1. 
%matrix_sess.('PrePost')(1:round(SWR_ind_mid_post)-1,1) = 1;
% 2 for post. 
%matrix_sess.('PrePost')(round(SWR_ind_mid_post):end,1) = 2;

%% populate matrix with structure information: Save one second before and after. 
seconds = 1; % one second before and after 
samples = (seconds*FP.cfg.hdr{1,1}.SamplingFrequency); % samples I will take before and after swrs

% raw data 
matrix_sess.('TwosWakeRaw') = cell(height(matrix_sess), 1);
matrix_sess.('TwosNremRaw') = cell(height(matrix_sess), 1);
% preprocessed data 
matrix_sess.('TwosWakeProc') = cell(height(matrix_sess), 1);
matrix_sess.('TwosNremProc') = cell(height(matrix_sess), 1);

for i = 1:height(matrix_sess) % iterate through each swr. 
    raw_data_struct = struct(); % structure of the data for an individual swr
    data_struct = struct(); % structure of the data for an individual swr
    swr_time = lfp_time(round(SWR_ind_mid(i))); % swr lfp time (initialized) 
    fiber_index = nearest_idx3(swr_time, FP.tvec); % fiber index closest to middle swr_time -- ok make sure this time is initialized 
    if matrix_sess.('sleep')(i) == 1 % if wake
        raw_data_struct.signal = FP.data(fiber_index-samples:fiber_index+samples); % saving signal
        data_struct.signal = FP.zF_win_60s(fiber_index-samples:fiber_index+samples); % saving signal
        raw_data_struct.tvec = FP.tvec(fiber_index-samples:fiber_index+samples); % saving time as well- even though it should be the same for each
        data_struct.tvec = raw_data_struct.tvec;     % saving time as well- even though it should be the same for each
        matrix_sess.('TwosWakeRaw'){i} = raw_data_struct;
        matrix_sess.('TwosWakeProc'){i} = data_struct;
    else % nrem
        raw_data_struct.signal = FP.data(fiber_index-samples:fiber_index+samples); % saving signal
        data_struct.signal = FP.zF_win_60s(fiber_index-samples:fiber_index+samples); % saving signal
        raw_data_struct.tvec = FP.tvec(fiber_index-samples:fiber_index+samples);     % saving time as well- even though it should be the same for each
        data_struct.tvec = raw_data_struct.tvec;     % saving time as well- even though it should be the same for each
        matrix_sess.('TwosNremRaw'){i} = raw_data_struct;
        matrix_sess.('TwosNremProc'){i} = data_struct;
    end
end

%% populate SWR details
% 18-23 :
% 'SWRdur','SWRamp','SWRpower','SWRf100ms','SWRtimestart','SWRtimeend';
%cd 'F:\M433\M433_2023_09_19_recording1'; 
LoadMetadata

samples_swr = (0.1*csc.cfg.hdr{1,1}.SamplingFrequency); % samples I will take before and after swrs

% filter the LFP band 
cfg = [];
cfg.f = [140 220];
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

matrix_sess.('SWR100ms') = cell(height(matrix_sess), 1);

for i = 1:height(matrix_sess) % iterate through each swr. 
    % swrduration = SWR_end - SWR_start 
    data_struct = struct(); % structure of the data for an individual swr
    matrix_sess.('SWRdur')(i) = SWR_end(i) - SWR_start(i); 
    % swramp = mean(SWRa.data(SWR_ind_start(i):SWR_ind_end(i)))
    matrix_sess.('SWRamp')(i) = mean(SWRa.data(SWR_ind_start(i):SWR_ind_end(i)));
    matrix_sess.('SWRpower')(i) = mean(SWRp.data(SWR_ind_start(i):SWR_ind_end(i)));
    matrix_sess.('SWRtimestart')(i) = SWR_start(i);
    matrix_sess.('SWRtimeend')(i) = SWR_end(i);
    alt_index = round(SWR_ind_mid(i));
    data_struct.data = SWRz.data(alt_index-samples_swr:alt_index+samples_swr);
    matrix_sess.('SWR100ms'){i} = data_struct;
end

%% populate dF information on from preproc data 
% dF from 2 seconds 
x1 = 1:1:FP.cfg.hdr{1,1}.SamplingFrequency; %1:1:500;%1001:1:3000; %1:1:2000;% % two seconds before for preproc data
x2 = FP.cfg.hdr{1,1}.SamplingFrequency + 1:1:2*FP.cfg.hdr{1,1}.SamplingFrequency + 1; % 1601:1:3200; %501:1:1000; %3001:1:5000; %2001:1:4000;%  % two seconds after for preproc data
% Is this right for the GFP Mice?

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

%% theta and delta power 
% restrict lfp to time intervals 
lfp_pre = restrict(csc,ExpKeys.prerecord(1), ExpKeys.prerecord(2)); 
lfp_post = restrict(csc,ExpKeys.postrecord(1),ExpKeys.postrecord(2));

% ~ Theta 
% Filter LFP 
params = get_theta_params_forMimi();

cfg = [];
cfg.f = params.theta.band ;
cfg.type = 'fdesign';
theta_lfp_pre = FilterLFP(cfg, lfp_pre);
theta_lfp_post = FilterLFP(cfg,lfp_post);
theta_data.lfp_tsd_pre = theta_lfp_pre;
theta_data.lfp_tsd_post = theta_lfp_post; 

% hilbert transform to get power 
cfg_zsc = [];
cfg_zsc.output = 'power';
theta_data.power_pre = LFPpower(cfg_zsc, theta_data.lfp_tsd_pre);
theta_data.power_post = LFPpower(cfg_zsc, theta_data.lfp_tsd_post);

% zscore the power
theta_data.zpower_pre = zscore_tsd(theta_data.power_pre);
theta_data.zpower_post = zscore_tsd(theta_data.power_post);

% ~ SWR (for 2 s instead of centered on SWR... instantaneous ... can be for sanity check.)
cfg = [];
cfg.f = [140 250];
cfg.display_filter = 0; 
swr_lfp = FilterLFP(cfg, lfp_pre);
swr_data.lfp_tsd_pre = swr_lfp;
swr_lfp = FilterLFP(cfg, lfp_post);
swr_data.lfp_tsd_post = swr_lfp;

% hilbert transform to get power  
cfg_zsc = [];
cfg_zsc.output = 'power';
swr_data.power_pre = LFPpower(cfg_zsc, swr_data.lfp_tsd_pre);
swr_data.power_post = LFPpower(cfg_zsc, swr_data.lfp_tsd_post);

% zscore the power
swr_data.zpower_pre = zscore_tsd(swr_data.power_pre);
swr_data.zpower_post = zscore_tsd(swr_data.power_post);

% ~ Delta 
% For Delta you need to decimate 
cfg = [];
cfg.decimate = 10;   % 1500 → 150 Hz
csc_ds_pre = decimate_tsd(cfg, lfp_pre);
csc_ds_post = decimate_tsd(cfg, lfp_post);

cfg = [];
cfg.f = [0.5 4];
cfg.display_filter = 0;
delta_lfp = FilterLFP(cfg, csc_ds_pre);
delta_data.lfp_tsd_pre = delta_lfp;
delta_lfp = FilterLFP(cfg, csc_ds_post);
delta_data.lfp_tsd_post = delta_lfp;

% hilbert transform to get power 
cfg_zsc = [];
cfg_zsc.output = 'power';
delta_data.power_pre = LFPpower(cfg_zsc, delta_data.lfp_tsd_pre);
delta_data.power_post = LFPpower(cfg_zsc, delta_data.lfp_tsd_post);

% zscore the power
delta_data.zpower_pre = zscore_tsd(delta_data.power_pre);
delta_data.zpower_post = zscore_tsd(delta_data.power_post);

% initialized lfp 
%delta_data.zpower.lfp = delta_data.zpower.tvec - delta_data.zpower.tvec(1); 

%% Add the power to the matrix
% originally I took the mean over the period of SWR but now I'm taking an
% instantenous point at the center of the SWR. 

SWR_time_center = (evt.tstart + evt.tend)./2 ; 
SWR_ind_mid_pre = nearest_idx3(SWR_time_center,theta_data.lfp_tsd_pre.tvec); 
SWR_ind_mid_post = nearest_idx3(SWR_time_center,theta_data.lfp_tsd_post.tvec); 
SWR_ind_mid_pre_d = nearest_idx3(SWR_time_center,delta_data.lfp_tsd_pre.tvec); 
SWR_ind_mid_post_d = nearest_idx3(SWR_time_center,delta_data.lfp_tsd_post.tvec); 

i_pre = 0 ; 
i_post = 0; 
for i = 1:height(matrix_sess) % iterate through each swr. 
    if SWR_time_mid(i) < post % if pre
        i_pre = i_pre + 1; 
        matrix_sess.('thetap')(i) = theta_data.zpower_pre.data(SWR_ind_mid_pre(i_pre));  %mean(theta_data.zpower_pre.data(SWR_ind_start_pre(i):SWR_ind_end_pre(i)));
        matrix_sess.('swrp')(i) = swr_data.zpower_pre.data(SWR_ind_mid_pre(i_pre)); %mean(swr_data.zpower_pre.data(SWR_ind_start_pre(i):SWR_ind_end_pre(i)));
        matrix_sess.('deltap')(i) = delta_data.zpower_pre.data(SWR_ind_mid_pre_d(i_pre)); %mean(delta_data.zpower_pre.data(SWR_ind_start_pre_d(i):SWR_ind_end_pre_d(i)));
    elseif SWR_time_mid(i) > post
        i_post = i_post + 1; 
        matrix_sess.('thetap')(i) = theta_data.zpower_post.data(SWR_ind_mid_post(i_post)); %mean(theta_data.zpower_post.data(SWR_ind_start_post(i):SWR_ind_end_post(i))); % something needs to be initialized here. 
        matrix_sess.('swrp')(i) = swr_data.zpower_post.data(SWR_ind_mid_post(i_post)); %mean(swr_data.zpower_post.data(SWR_ind_start_post(i):SWR_ind_end_post(i)));
        matrix_sess.('deltap')(i) = delta_data.zpower_post.data(SWR_ind_mid_post_d(i_post)); %mean(delta_data.zpower_post.data(SWR_ind_start_post_d(i):SWR_ind_end_post_d(i)));
    end
end


%% Save average theta, delta and ratio for NREM and WAKE. 
% And also save the average theta and delta and theta/delta ratio for the
% sleep stage. 
fs = 1500;                 % sampling rate (Hz)
epoch_length = 2;         % seconds
samples_per_epoch = fs * epoch_length;

fs2 = 150; 
samples2 = fs2*epoch_length;

undefined = find(labels==4); 

theta_data.zpower_pre.tvec_init= theta_data.zpower_pre.tvec - theta_data.zpower_pre.tvec(1);
delta_data.zpower_pre.tvec_init= delta_data.zpower_pre.tvec - delta_data.zpower_pre.tvec(1);

theta_data.zpower_post.tvec_init= theta_data.zpower_post.tvec - theta_data.zpower_pre.tvec(1);
delta_data.zpower_post.tvec_init= delta_data.zpower_post.tvec - delta_data.zpower_pre.tvec(1);


theta_avg_wake = [];
swr_avg_wake = [];
delta_avg_wake = [];
theta_avg_nrem = [];
swr_avg_nrem = [];
delta_avg_nrem = [];

for i = 1:1:length(labels) % for each 2 second period 
    % based on i,find the corresponding time interval. then find the
    % corresponding index and the corresponding data information. 

    seconds_t = i*2 -  2; % first sample 
    
    if seconds_t < 1800
        ind_t = nearest_idx3(seconds_t, theta_data.zpower_pre.tvec_init); % can be used for swr too. 
        ind_d = nearest_idx3(seconds_t, delta_data.zpower_pre.tvec_init);
        theta_temp = nanmean(theta_data.zpower_pre.data(ind_t:ind_t+samples_per_epoch));
        swr_temp = nanmean(swr_data.zpower_pre.data(ind_t:ind_t+samples_per_epoch));
        delta_temp = nanmean(delta_data.zpower_pre.data(ind_d:ind_d+samples2));
    elseif seconds_t > 1800 & labels(i) ~= 4 & seconds_t+2 < delta_data.zpower_post.tvec_init(end)
        ind_t = nearest_idx3(seconds_t, theta_data.zpower_post.tvec_init); % can be used for swr too. 
        ind_d = nearest_idx3(seconds_t, delta_data.zpower_post.tvec_init);
        theta_temp = nanmean(theta_data.zpower_post.data(ind_t:ind_t+samples_per_epoch));
        swr_temp = nanmean(swr_data.zpower_post.data(ind_t:ind_t+samples_per_epoch));
        delta_temp = nanmean(delta_data.zpower_post.data(ind_d:ind_d+samples2));
    end

    if labels(i) == 2 
        theta_avg_wake = [theta_avg_wake; theta_temp];
        swr_avg_wake = [swr_avg_wake; swr_temp];
        delta_avg_wake = [delta_avg_wake; delta_temp];
    elseif labels(i) == 3 
        theta_avg_nrem = [theta_avg_nrem; theta_temp];
        swr_avg_nrem = [swr_avg_nrem; swr_temp];
        delta_avg_nrem = [delta_avg_nrem; delta_temp];  
    end
end


theta_wake = nanmean(theta_avg_wake);
delta_wake = nanmean(delta_avg_wake);
swr_wake = nanmean(swr_avg_wake);

theta_nrem = nanmean(theta_avg_nrem);
delta_nrem = nanmean(delta_avg_nrem);
swr_nrem = nanmean(swr_avg_nrem);

swr_sleep.theta = [theta_wake theta_nrem];
swr_sleep.swr = [swr_wake swr_nrem];
swr_sleep.delta = [delta_wake delta_nrem];


%% Plot these powers with the spectrogram??

%zlfpr_pre = restrict(SWRz,ExpKeys.prerecord(1), ExpKeys.prerecord(2)); % if you don't have this, implement it (it's one line of code!)
%theta_pre = restrict(theta_data.zpower, ExpKeys.prerecord(1), ExpKeys.prerecord(2));
%%
%[S,F,T,P] = spectrogram(zlfpr_pre.data,hanning(500),50,0:40,csc.cfg.hdr{1,1}.SamplingFrequency);

%%
% figure(1); clf
% 
% % ---- Spectrogram ----
% yyaxis left
% imagesc(T, F, 10*log10(P));
% axis xy
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% colormap jet
% colorbar
% %ylim([0 20]);
% hold on
% 
% % ---- Theta power (right axis) ----
% yyaxis right
% plot(theta_pre.tvec - theta_pre.tvec(1), ...
%      theta_pre.data, ...
%      'b', 'LineWidth', 2)
% 
% ylabel('Theta power (z-score)')
% ylim([-4 4])   % adjust as needed
% 
% title('Theta Spectrogram with Theta Power Overlay')
%% RPE variable and track to huge matrix? 

%%
% ~~ Save Variables ~~ 
cd F:\sleepcounts
filename = append(file_name, "swr_sleep.mat");
save(filename, '-struct','swr_sleep')

%% Save everything
cd 'F:\SWR_DA_Mega_1s_theta'
filename = append(file_name, "mega1_new.mat");
save(filename,'matrix_sess')
   