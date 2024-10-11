%% SESSION PLOT : Hypothesis-Driven Plots 
% input: preprocessed fiber data, detected SWRs, good LFP csc
% output: session plot with: 
% ------ fiber signal after swrs (ask Wolford about this)
% ------ heat map of all trials fiber signal - 1SD shuffle

%% Load data 
% load fiber data 
clear; clc; 
%rng(pi);

cd 'D:\M453\M453-2024-01-20_recording8';
FP = load('M453_2024_01_20processed.mat');

% extract SWR intervals
load('2024-01-20_M453_recording8detectedSWRs.mat')

file_name = 'M453_2024_01_20'; 

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

%% Extract fiber after swrs 
prepros_signal = [];
prepros_signal = FP.zF_win_60s; 

SWR_ind_mid_post = (SWR_ind_start_post + SWR_ind_end_post)/2;  %middle index 

% keep all SWR after that time -- should be 654
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

% I can't just subtract the SD because then the negative ones will be even
% more negative... 
% it is also can be flipping the sign... if less than SD make 0... 
% or just accept that this analysis will only show you things that are
% larger than SD...

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

% restrict lfp to time intervals 
prerecord_init = ExpKeys.prerecord(1)-csc.tvec(1);
prerecord_end = ExpKeys.prerecord(2)-csc.tvec(1);

postrecord_init = ExpKeys.postrecord(1)-csc.tvec(1);
postrecord_end = ExpKeys.postrecord(2)-csc.tvec(1);

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

%% 
pre_heatmap = zdF_extract_pre - circ_std_fiber_pre; 
post_heatmap = zdF_extract_post - circ_std_fiber_post; 

updated_map_pre = zeros(size(pre_heatmap));
updated_map_post = zeros(size(post_heatmap));


for iter1 = 1:1: length(circ_std_fiber_pre)
    for iter_rows = 1:1:length(pre_SWR_ind)
        if zdF_extract_pre_con(iter_rows,iter1) > 0
            updated_map_pre(iter_rows,iter1) = zdF_extract_pre(iter_rows,iter1) - circ_std_fiber_pre(iter1);
        else
            updated_map_pre(iter_rows,iter1) = zdF_extract_pre(iter_rows,iter1) + circ_std_fiber_pre(iter1);
        end
    end
    for iter_rows2 = 1:1:length(post_SWR_ind)
        if zdF_extract_post(iter_rows2,iter1) > 0
            updated_map_post(iter_rows2,iter1) = zdF_extract_post(iter_rows2,iter1) - circ_std_fiber_post(iter1);
        else
            updated_map_post(iter_rows2,iter1) = zdF_extract_post(iter_rows2,iter1) + circ_std_fiber_post(iter1);
        end
    end
end


pre_heatmap_sorted = sortrows(pre_heatmap, 5001, 'descend');
post_heatmap_sorted = sortrows(post_heatmap, 5001, 'descend');

updated_pre_heatmap_sorted = sortrows(updated_map_pre, 4501, 'descend');
updated_post_heatmap_sorted = sortrows(updated_map_post, 4501, 'descend');

% order heatmap by the largest value at 1s? Should start increasing by then
% 0.5? 

%% shuffled heatmaps 
% getting a constant value here... 
pre_heatmap_con = zdF_extract_pre_con  - circ_std_fiber_pre; % need to edit this 
updated_map_pre_con = zeros(size(pre_heatmap_con));

post_heatmap_con = zdF_extract_post_con - circ_std_fiber_post; 
updated_map_post_con = zeros(size(post_heatmap_con));

for iter1 = 1:1: length(circ_std_fiber_pre)
    for iter_rows = 1:1:length(pre_SWR_ind)
        if zdF_extract_pre_con(iter_rows,iter1) > 0
            updated_map_pre_con(iter_rows,iter1) = zdF_extract_pre_con(iter_rows,iter1) - circ_std_fiber_pre(iter1);
        else
            updated_map_pre_con(iter_rows,iter1) = zdF_extract_pre_con(iter_rows,iter1) + circ_std_fiber_pre(iter1);
        end
    end
    for iter_rows2 = 1:1:length(post_SWR_ind)
        if zdF_extract_post_con(iter_rows2,iter1) > 0
            updated_map_post_con(iter_rows2,iter1) = zdF_extract_post_con(iter_rows2,iter1) - circ_std_fiber_post(iter1);
        else
            updated_map_post_con(iter_rows2,iter1) = zdF_extract_post_con(iter_rows2,iter1) + circ_std_fiber_post(iter1);
        end
    end
end
% only row 8001 has it... why is that? 

% need to also iterate through all rows...

pre_heatmap_sorted_con = sortrows(pre_heatmap, 5001, 'descend');
post_heatmap_sorted_con = sortrows(post_heatmap, 5001, 'descend');

updated_pre_heatmap_sorted_con = sortrows(updated_map_pre_con, 5001, 'descend');
updated_post_heatmap_sorted_con = sortrows(updated_map_post_con, 5001, 'descend');

% order heatmap by the largest value at 1s? Should start increasing by then
% 0.5? 

%% Figures together 
fig2 = figure(2);

med_c = [104,187,225]./255; % rgb(167, 199, 231) rgb(255, 165, 0)  blue: rgb(104,187,227)
low_c = [78,178,101]./255; % color

subplot(2,4,1:2)
shadedErrorBar(time_extract_pre(1,:),circ_avg_fiber_pre,circ_std_fiber_pre,'lineProps','-k','transparent',1)
hold on
plot(time_extract_pre(1,:),circ_avg_fiber_pre,'LineWidth',2,'Color','k')
% plot average on top with larger line
hold on
plot(time_extract_pre(1,:),avg_fiber_pre,'LineWidth',2,'Color',low_c)
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
ylim([-0.8 0.8])
xlim([0 8])
xticks([0 4 8])
xticklabels({'-4','0','4'})
title('Pre Track Rest PETH','FontSize', 20)
ylabel('Mean [DA] (z-score)','FontSize', 16)
xlabel('Time from SWR (s)','FontSize', 16)
legend('','shuffle','signal','Location','northwest')
legend boxoff

subplot(2,4,5:6)
shadedErrorBar(time_extract_post(1,:),circ_avg_fiber_post,circ_std_fiber_post,'lineProps','-k','transparent',1)
hold on
plot(time_extract_post(1,:),circ_avg_fiber_post,'LineWidth',2,'Color','k')
hold on
plot(time_extract_post(1,:),avg_fiber_post,'LineWidth',2,'Color',low_c)
hold on
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
ylim([-0.8 0.8])
xticks([0 4 8])
xticklabels({'-4','0','4'})
title('Post Track Rest PETH','FontSize', 20)
ylabel('Mean [DA] (z-score)','FontSize', 16)
xlabel('Time from SWR (s)','FontSize', 16)
legend('','shuffle','signal','Location','northwest')
legend boxoff 

subplot(2,4,3)
imagesc(updated_pre_heatmap_sorted) 
colorbar
xticks([0 4000 8000])
xticklabels({'-4','0','4'})
xlabel('Time from SWR (s)')
ylabel('zdF - 1 sd of shuffle')
title('Pre-task fiber after SWR')

subplot(2,4,7)
imagesc(updated_post_heatmap_sorted) 
colorbar
xticks([0 4000 8000])
xticklabels({'-4','0','4'})
xlabel('Time from SWR (s)')
ylabel('zdF - 1 sd of shuffle')
title('Post-task fiber after SWR')

subplot(2,4,4)
imagesc(updated_pre_heatmap_sorted_con) 
colorbar
xticks([0 4000 8000])
xticklabels({'-4','0','4'})
xlabel('Time from random event (s)')
ylabel('zdF - 1 sd of shuffle')
title('Pre-task fiber after random time-point')

subplot(2,4,8)
imagesc(updated_post_heatmap_sorted_con) 
colorbar
xticks([0 4000 8000])
xticklabels({'-4','0','4'})
xlabel('Time from random event (s)')
ylabel('zdF - 1 sd of shuffle')
title('Post-task fiber after random time-point')

set(gcf,'Color',[1,1,1])
shg

hold off

txt = {'Session Plot: Hypothesis Suplots'};
sgtitle(txt)

fig2.WindowState = 'maximized';
 %%

  cd 'C:\Users\mimia\Documents\Replay-DA\Figures\M453\recording 8'
  % save descriptive plot
  saveas(fig2,'M453_recording8_hypothesis.png') % CHANGE THIS 

avg_SWR_DA.circ_avg_pre = circ_avg_fiber_pre;
avg_SWR_DA.circ_avg_post = circ_avg_fiber_post;

avg_SWR_DA.circ_std_pre = circ_std_fiber_pre;
avg_SWR_DA.circ_std_post = circ_std_fiber_post;

avg_SWR_DA.avg_fiber_pre = avg_fiber_pre;
avg_SWR_DA.avg_fiber_post = avg_fiber_post;

avg_SWR_DA.time = time_extract_pre;

% circ_avg_fiber_pre + post
% circ_std_fiber_pre + post
% avg_fiber_pre
% time_extract_pre for one session (should all be the same) 
%%
 cd 'D:\M453\avg_data'
 filename = append(file_name, "avg.mat");
 save(filename, '-struct','avg_SWR_DA')

% STDS are correct!
