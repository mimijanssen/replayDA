%% SESSION PLOT : Descriptive Plots 
% input: preprocessed and raw fiber data, detected SWRs, good LFP csc, and DLC points 
% output: 
% ~ session plot with: 
% ------ raw fiber signal at 2 different timescales
% ------ fiber RPEs 
% ------ swr counts (pre and post task) 
% ------ LFP spectrogram with speed : to see if coupled with theta
% ------ position
% ------ average speed (pre and post task)
% ~ processed data: 
% ------ average rpe data 
% ------ swr counts

% written by Mimi Janssen, 8/16/2024

%% establish paths 
% restoredefaultpath; clear classes; % start with a clean slate
% 
% cd('C:\Users\mimia\Documents\Toolboxes\vandermeerlab-replay-da\code-matlab\shared'); % or, wherever your code is located -- NOTE \shared subfolder!
% p = genpath(pwd); % create list of all folders from here
% addpath(p);
 
%% Load Data 
clear; clc;
rng(pi)
cd 'D:\M548\M548_2024_08_25_recording1'; 
FP = load('M548_2024_08_25processed.mat'); % fiber data processed with my pipeline
load('M548_2024_08_25detectedSWRs.mat') % SWR intervals
load('M548_2024-08-25_track.mat') % for pseudo_outcomes % CHANGE THIS 

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
% for M452, I ran into a problem where swr start was closer to a later
% index (83) and swr end was closer to a earlier index (82). I'll include
% the later because maybe more of it was included? 

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
swr_count = [pre_count post_count];
swr_label = ["Pre" "Post"]; 

%% load DLC points & process the points to get speed 
% output: keepbodyx, keepbodyy, speed_pre, speed_post

P = readtable('M548_2024_08_25-convertedDLC_resnet50_Linear TrackApr5shuffle1_100000.csv','PreserveVariableNames',true); % CHANGE THIS 
A = table2array(P);
% midbody = [A(:,1), A(:,17) A(:,18), A(:,19)];
% body_ind = find(midbody(:,4) >= 0.99);
% frames = A(:,1); 
% frame count 
frames = A(:,1); 
% BODY PARTS: 
midbody = [A(:,1), A(:,17) A(:,18), A(:,19)];
[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, [] );
bodyx = midbody(:,2);
bodyy = midbody(:,3);
% Remove points that are not well predicted 
keep_bodyx = bodyx(midbody(:,4)>=0.98);
keep_bodyy = bodyy(midbody(:,4)>=0.98); % I want to pad with nan 

% pad with nan that are not well predicted 
padbodyx = bodyx;
padbodyy = bodyy; 
padbodyx(midbody(:,4)<=0.98)=NaN;
padbodyy(midbody(:,4)<=0.98)=NaN;

%bad_body_pred = (length(keep_bodyx)/length(bodyx))*100 ; % 99.3% accurate
fpos = {};
fpos.type = {'tsd'};
fpos.tvec = (frames*(1/30))'; 
fpos.label= {'x','y'}; 
fpos.data = [(medfilt1(padbodyx,'omitnan'))'; (medfilt1(padbodyy,'omitnan'))'];
fpos.cfg = {};
linspd = getLinSpd([],fpos);

% edited to initialize expcode times

prerecord_init = ExpKeys.prerecord(1)-csc.tvec(1);
prerecord_end = ExpKeys.prerecord(2)-csc.tvec(1);
track_init = ExpKeys.task(1)-csc.tvec(1);
track_end = ExpKeys.task(2)-csc.tvec(1);
postrecord_init = ExpKeys.postrecord(1)-csc.tvec(1);
postrecord_end = ExpKeys.postrecord(2)-csc.tvec(1);
% restrict linspd to pre and post...

spd_pre = restrict(linspd,prerecord_init, prerecord_end); % if you don't have this, implement it (it's one line of code!)
spd_post = restrict(linspd,postrecord_init, postrecord_end); % if you don't have this, implement it (it's one line of code!)

% initialize spd_pre and spd_post
spd_post.tvec = spd_post.tvec - spd_post.tvec(1);
spd_pre.tvec = spd_pre.tvec - spd_pre.tvec(1);

sample_video = median(diff(Timestamps)); % 1 s / average time between samples 
fs_video = 1/sample_video; % frames per ms... %29.97 frames/s ... 

fposx = {};
fposx.type = {'tsd'};
fposx.tvec = fpos.tvec; 
fposx.label= {'x'}; 
fposx.data = [fpos.data(1,:)];
fposx.cfg = {};
fposx.cfg.history.mfun = {};
fposx.cfg.history.cfg = {};

% restrict pre rest
pos_pre = restrict(fposx,prerecord_init, prerecord_end); % if you don't have this, implement it (it's one line of code!)

% restrict post rest
pos_post = restrict(fposx,postrecord_init, postrecord_end); 

%
posx_pre_start = nearest_idx3(prerecord_init,fpos.tvec); % start time for pre this is fine because didn't initialzie either yet
posx_pre_end = nearest_idx3(prerecord_end,fpos.tvec); % end time for pre
posx_post_start = nearest_idx3(postrecord_init,fpos.tvec);
posx_post_end = nearest_idx3(postrecord_end,fpos.tvec);

spdx = diff(padbodyx); % has the same points as linspd which is nice 
spdy = diff(padbodyy);

% speed is one point off maybe restrict pos and then do linspeed
speed_pre = median(linspd.data(posx_pre_start:posx_pre_end));
speed_post = median(linspd.data(posx_post_start:posx_post_end));

%% extract fiber after swrs 
% output: avg_circ_zdF_extract_pre, avg_fiber_pre  & post

% Establish what type of preprocessed signal you want to use: 
prepros_signal = FP.zF_win_60s; 

seconds = 8; % seconds you want to show in plot
samples = (seconds*FP.cfg.hdr{1,1}.SamplingFrequency)/2;  % divide by two because you want 4s + and - directions 

zdF_win1_extract = zeros(length(SWR_ind_mid), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
time_win1_extract =  zeros(length(SWR_ind_mid), seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);
% the time is not aligning
for  ievt = 1:1:length(SWR_ind_mid)-1
    %if SWR_fiber_ind(ievt)+samples > 
    timeset = time((SWR_fiber_ind(ievt)-samples):(SWR_fiber_ind(ievt)+samples)); % pick fiber events that are 4 seconds each way
    time_win1_extract(ievt,:) = time((SWR_fiber_ind(ievt)-samples):(SWR_fiber_ind(ievt)+samples))-timeset(1); 
    zdF_win1_extract(ievt,:) = (prepros_signal((SWR_fiber_ind(ievt)-samples):(SWR_fiber_ind(ievt)+samples)));
end

avg_fiber_win1 = nanmean(zdF_win1_extract);
std_fiber_win1 = 2*std(zdF_win1_extract);

% post 
% find time where post-sleep starts 
post = ExpKeys.postrecord(1) - csc.tvec(1); % time of post sleep period, initialized  

% find that corresponding SWR time index closest to that 
SWR_ind_start_post = nearest_idx3(post,SWR_iv(:,1)); %find(abs(lfp-SWR_iv(1,1)) < 0.0005); % only for one but can we extend this to everything??
SWR_ind_end_post = nearest_idx3(post,SWR_iv(:,2)); %find(abs(lfp-SWR_iv(1,2)) < 0.0005); % only for one but can we extend this to everything??

% find midpoint 
SWR_ind_mid_post = (SWR_ind_start_post + SWR_ind_end_post)/2;  %middle index 

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
addpath('C:\Users\mimia\Documents\Toolboxes\raacampbell-shadedErrorBar')

avg_fiber_post = nanmean(zdF_extract_post);
std_fiber_post = 2*std(zdF_extract_post);

X = prepros_signal; 
N = 1000; % number of circshifts 
K=randi([1 length(prepros_signal)],1, N); % pick a random number between 1 and number of samples ... 100 times 
events_num = length(post_SWR_ind);
% initialize 
circ_zdF_extract = zeros(events_num, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
avg_circ_zdF_extract_post = zeros(N,seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for iter_circ = 1:1:N % 1 through number N
    Y = circshift(X,K(iter_circ)); % circshift the entire fiber signal based on the random number 
    for ievt = 1:1:length(post_SWR_ind)-1 % for each SWR event pick out 1-978, pick out that subset of the fiber signal 
    % find time for x axis
       timeset = time((SWR_fiber_ind_post(ievt)-samples):(SWR_fiber_ind_post(ievt)+samples)); % pick fiber events that are 4 seconds each way
       circ_zdF_extract(ievt,:) = (Y((SWR_fiber_ind_post(ievt)-samples):(SWR_fiber_ind_post(ievt)+samples))); %resets every time
    end
    avg_circ_zdF_extract_post(iter_circ,:) = nanmean(circ_zdF_extract);
end
% last row is all zeros... 
circ_avg_fiber_post = nanmean(avg_circ_zdF_extract_post);
circ_std_fiber_post = 2*std(avg_circ_zdF_extract_post);


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

avg_fiber_pre = nanmean(zdF_extract_pre);
%SEM_fiber = std([avg_fiber],0,1)/sqrt(length([]));
std_fiber_pre = 2*std(zdF_extract_pre);

%std_top_pre = avg_fiber_pre + 2*std(zdF_extract_pre);
%std_bot_pre = avg_fiber_pre - 2*std(zdF_extract_pre);

% CIRCSHIFT ------------------------------------------------------------
% elements in the array X by K positions. shifts by [m,n] n dimension. 
X = prepros_signal; 
N = 1000; % number of circshifts 
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

circ_avg_fiber_pre = nanmean(avg_circ_zdF_extract_pre);
%SEM_fiber = std([avg_fiber],0,1)/sqrt(length([]));
circ_std_fiber_pre = 2*std(avg_circ_zdF_extract_pre);

%% preprocess csc for spectrogram
zlfp = zscore_tsd(csc);

% restrict pre rest
zlfpr_pre = restrict(zlfp,ExpKeys.prerecord(1), ExpKeys.prerecord(2)); % if you don't have this, implement it (it's one line of code!)

% restrict post rest
zlfpr_post = restrict(zlfp,ExpKeys.postrecord(1), ExpKeys.postrecord(2)); % if you don't have this, implement it (it's one line of code!)

% spectrogram function pre rest
[S_pre,F_pre,T_pre,P_pre] = spectrogram(zlfpr_pre.data,hanning(7500),3750,1:20,csc.cfg.hdr{1,1}.SamplingFrequency);

% spectrogram function post rest
[S_post,F_post,T_post,P_post] = spectrogram(zlfpr_post.data,hanning(7500),3750,1:20,csc.cfg.hdr{1,1}.SamplingFrequency);

% figure(2)
subplot(2,1,1)
yyaxis left
imagesc(T_pre,F_pre,10*log10(P_pre)); % converting to dB as usual
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
hold on
subplot(2,1,2)
yyaxis left
imagesc(T_post,F_post,10*log10(P_post)); % converting to dB as usual
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
hold on

%% Extract track codes 
% Load the pseudo prob data 
LoadExpKeys();
cfg_evt = [];
cfg_evt.eventList = ExpKeys.eventList;
cfg_evt.eventLabel = ExpKeys.eventLabel;
evt3 = LoadEvents(cfg_evt);

cfg_fiber.fc = {'CSC30.ncs'};
csc_photo_fc = LoadCSC(cfg_fiber);

evt_ordered = sort([evt3.t{1},evt3.t{2}]); % don't need to sort anymore but ok % NEED code to determine which is first... ??
photobeam_times = evt_ordered - csc_photo_fc.tvec(1); % initialize time 

high_t = photobeam_times(pseudo_outcomes(1:length(pseudo_outcomes)) == 3); % s
med_t = photobeam_times(pseudo_outcomes(1:length(pseudo_outcomes)) == 2); % should be 36 not 52
low_t = photobeam_times(pseudo_outcomes(1:length(pseudo_outcomes)) == 1); % 6 

% Color
med_c = [104,187,225]./255; % rgb(167, 199, 231) rgb(255, 165, 0)  blue: rgb(104,187,227)
low_c = [0,104,87]./255; 
high_c = [255, 165,0]./255; % rgb(204, 204, 255) periwinkle

all_events = [{high_t}, {med_t}, {low_t}];
all_colors = [{high_c},{med_c},{low_c}];

window = 8*FP.cfg.hdr{1,1}.SamplingFrequency; % automate based on sampling frequency 8seconds/ 
t = [];
t = FP.tvec; %FP_Tlab.tvec- FP_Tlab.tvec(1);
% initialize values for average matrix 
high = zeros(length(high_t),window*2); %rows = trial, columns = signal 
med = zeros(length(med_t),window*2);
low = zeros(length(low_t),window*2);

for iter = 1:1:3
    hold on
    for ptime = 1:1:length(all_events{1,iter})
        indxpb = find(abs(t-all_events{1,iter}(ptime)) < 0.0005); % find the time of the event. 0.0005
        init_trial = (indxpb - (window)); 
        end_trial = (indxpb + (window));
            if iter ==1 
                high(ptime,:) = prepros_signal(init_trial:end_trial-1); 
                %h1 = plot((t(init_trial:end_trial))-t(init_trial), zdF(init_trial:end_trial), 'Color',[all_colors{1,iter}]);
            elseif iter == 2
                med(ptime,:) = prepros_signal(init_trial:end_trial-1); 
                %h2 = plot(t(init_trial:end_trial)-t(init_trial), zdF(init_trial:end_trial), 'Color',[all_colors{1,iter}]); % set t back to 0
            else
                low(ptime,:) = prepros_signal(init_trial:end_trial-1); 
                %h3 = plot((t(init_trial:end_trial))-t(init_trial), zdF(init_trial:end_trial), 'Color',[all_colors{1,iter}]);
            end
                t_shared = t(init_trial:end_trial-1)-t(init_trial); % double the window
    end
end

% average by columns
avg_high = mean(high); 
avg_med = mean(med);
avg_low = mean(low);

% std by columns
std_high = std(high); 
std_med = std(med);
std_low = std(low);

%% Pwelch 
% welch on raw fiber and raw lfp 
wsize = csc.cfg.hdr{1,1}.SamplingFrequency * 4; % sampling frequency * 4s is the amount of samples (the window size needed)
wsize2 = csc_photo_fc.cfg.hdr{1,1}.SamplingFrequency*4;
[Pxx,F] = pwelch(csc.data,hanning(wsize),wsize/2,[],csc.cfg.hdr{1,1}.SamplingFrequency);
[Pxx_f,F_f] = pwelch(csc_photo_fc.data,hanning(wsize2),wsize2/2,[],csc_photo_fc.cfg.hdr{1,1}.SamplingFrequency);

%% restrict signal for pwelch
fiber_pre = restrict(csc_photo_fc,ExpKeys.prerecord(1), ExpKeys.prerecord(2)); % if you don't have this, implement it (it's one line of code!)
fiber_track = restrict(csc_photo_fc,ExpKeys.task(1), ExpKeys.task(2)); % if you don't have this, implement it (it's one line of code!)
fiber_post = restrict(csc_photo_fc,ExpKeys.postrecord(1), ExpKeys.postrecord(2)); % if you don't have this, implement it (it's one line of code!)

[Pxx_fpre,F_fpre] = pwelch(fiber_pre.data,hanning(wsize2),wsize2/2,[],csc_photo_fc.cfg.hdr{1,1}.SamplingFrequency);
[Pxx_ftrack,F_ftrack] = pwelch(fiber_track.data,hanning(wsize2),wsize2/2,[],csc_photo_fc.cfg.hdr{1,1}.SamplingFrequency);
[Pxx_fpost,F_fpost] = pwelch(fiber_post.data,hanning(wsize2),wsize2/2,[],csc_photo_fc.cfg.hdr{1,1}.SamplingFrequency);

%% stitch the subplots together 
fig = figure(1);
% fiber 
% 1-3 
%subplot(4,8,1:3)
FP = csc_photo_fc;
subnum= [1,5];
sessionTitle = {'10s','Whole Session'};
last_time = length(FP.tvec); %the value of the last time point is how many seconds the recording was
timerange1 = 10/0.0002; %datapoint range for 10 s
timerange3= length(FP.tvec); % datapoint range for all data
time_ranges = [timerange1, timerange3]; %in seconds 
for t_i = 1:length(time_ranges)
    t_range = 1:time_ranges((t_i));
    subplot(3,4,subnum(t_i));
    plot(FP.tvec(t_range), FP.data(t_range), 'Color', [0 0.5 0])
    title([num2str(sessionTitle{1,t_i})], 'Interpreter','none')
    ylabel('Fiber Signal (V)'); xlabel('Time (s)');
end
hold on

subplot(3,4,9)
shadedErrorBar(t_shared,avg_high,std_high,'lineprops',{'-','color',high_c,'MarkerFaceColor',high_c});
hold on
plot(t_shared,avg_high,'LineWidth',2,'Color',high_c)
hold on
shadedErrorBar(t_shared,avg_med,std_med,'lineprops',{'-','color',med_c,'MarkerFaceColor',med_c});
hold on
plot(t_shared,avg_med,'LineWidth',2,'Color',med_c)

shadedErrorBar(t_shared,avg_low,std_low,'lineprops',{'-','color',low_c,'MarkerFaceColor',low_c});
plot(t_shared,avg_low,'LineWidth',2,'Color',low_c)

xlabel('Time from photobeam break (s)')
xlim([0 16])
xticks([0 8 16])
xticklabels({'-8','0','8'})
legend('high','','medium','','low');
legend boxoff
ylabel('Mean Signal (detrended & z-scored)')
title('Fiber Signal RPE')

subplot(3,4,3)
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 100]);
title('PSD for LFP')

subplot(3,4,2)
plot(F_fpre,10*log10(Pxx_fpre)); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
hold on
plot(F_fpost,10*log10(Pxx_fpost)); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
leg = legend('pre','post');
set(leg,'Box','off')
xlim([0 10]);
title('PSD for Fiber (pre & post track)')

subplot(3,4,6)
plot(F_ftrack,10*log10(Pxx_ftrack),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 20]);
hold on
title('PSD for Fiber (track)')

subplot(3,4,10)
plot(F_ftrack,10*log10(Pxx_ftrack),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 100]);
title('PSD for Fiber (track)- noise check')

subplot(3,4,[7,8])
yyaxis left
imagesc(T_pre,F_pre,10*log10(P_pre)); % converting to dB as usual
%set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
hold on
yyaxis right 
plot(spd_pre,'.k') %linspd.data(posx_pre_start:posx_pre_end)
title('Pre-Track Rest')
xlabel('Time (s)')
ylabel('Speed (pixels/s)')

% 28
subplot(3,4,[11,12])
yyaxis left
imagesc(T_post,F_post,10*log10(P_post)); % converting to dB as usual
%set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
hold on
yyaxis right 
plot(spd_post,'.k') % linspd.data(posx_post_start:posx_post_end)
title('Post-Track Rest')
xlabel('Time (s)')
ylabel('Speed (pixels/s)')

% position subplot
%32
subplot(3,4,4)
plot(keep_bodyx,keep_bodyy,'k.')
title('Midbody Position')
xlabel('X position')
ylabel('Y position')

set(gcf,'Color',[1,1,1])
shg

hold off

%txt = {'Session Plot: Descriptive Suplots',['SWR count: pre-', num2str(swr_count(1)), ', post-', num2str(swr_count(2))]};
txt = {'Session Plot: Descriptive Suplots',['SWR count: pre-', num2str(swr_count(1)), ', post-', num2str(swr_count(2))],[ 'Speed (pixels/s): pre-', num2str(speed_pre), ', post-', num2str(speed_post)]};
sgtitle(txt)

fig.WindowState = 'maximized';
% 
% --------------------------------------------------------------------------

  cd 'C:\Users\mimia\Documents\ReplayDA Figures\M548\recording 1'
  % save descriptive plot
  saveas(fig,'M548_recording1_descriptive.png') % CHANGE THIS 

avg_RPE.t_shared = t_shared;
avg_RPE.avg_high = avg_high;
avg_RPE.std_high = std_high;
avg_RPE.avg_low = avg_low;
avg_RPE.std_low = std_low;
avg_RPE.avg_med = avg_med;
avg_RPE.std_med = std_med;
avg_RPE.swr_count = swr_count;
avg_RPE.swr_label = ['pre','post'];


% other variables to save: 
% fpos
% linspd 
% spd_post 
% spd_pre 

%%
% --------------------------------------------------------------------------
   cd 'D:\M548\avg_data'
   file_name = 'M548_2024_08_25'; 
   filename = append(file_name, "avgRPE.mat");
   save(filename, '-struct','avg_RPE')
   
  % filename = append(file_name, "pos.mat");
  % save(filename, 'fpos','linspd','spd_post','spd_pre')

