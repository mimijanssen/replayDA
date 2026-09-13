%% PETH separated by sleep stage? 
%% SESSION PLOT : Hypothesis-Driven Plots 
% input: preprocessed fiber data, detected SWRs, good LFP csc
% output: session plot with: 
% ------ fiber signal after swrs (ask Wolford about this)
% ------ heat map of all trials fiber signal - 1SD shuffle

%% Load data 
% load fiber data 
clear; clc; 
rng(pi);
cd 'F:\M548\M548_2024_08_31_recording7';
FP = load('M548_2024_08_31processed.mat');
load('M548_2024_08_31detectedSWRs.mat')
load('labels.mat'); % 2 is wake, 3 is nrem, 4 is undefined. 
file_name = 'M548_2024_08_31'; 
addpath('C:\Users\mimia\Docments\Toolboxes\shadedErrorBar')

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

%% find SWR times that are in nrem sleep stage 
% ~ for every swr time, check if that is an nrem stage in labels (==3) 
labels_time_center = (0:length(labels)-1)' * 2 + 1;

SWR_time_nrem = [];
for swri = 1:1:length(SWR_time_mid)
    label_ind = nearest_idx3(SWR_time_mid(swri),labels_time_center); % index of labels_time 
    if labels(label_ind) == 3
        SWR_time_nrem = [SWR_time_nrem; SWR_time_mid(swri)];
    end
end

% ~ convert to fiber time index
SWR_fiber_ind_nrem = zeros(length(SWR_time_nrem),1); 
for i = 1:1:length(SWR_time_nrem)
    SWR_fiber_ind_nrem(i) = nearest_idx3(SWR_time_nrem(i),time);
end

% ~ PETH 
zdF_extract_nrem = zeros(length(SWR_time_nrem)-1, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
time_extract_nrem =  zeros(length(SWR_time_nrem)-1, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for ievt = 1:1:length(SWR_time_nrem)-1 %-3 for M600 recording 1; -1 for everyone else because the last SWR sometimes is too close to the end recording
    timeset = time((SWR_fiber_ind_nrem(ievt)-samples):(SWR_fiber_ind_nrem(ievt)+samples)); % pick fiber events that are 4 seconds each way
    time_extract_nrem(ievt,:) = time((SWR_fiber_ind_nrem(ievt)-samples):(SWR_fiber_ind_nrem(ievt)+samples))-timeset(1); 
    zdF_extract_nrem(ievt,:) = (prepros_signal((SWR_fiber_ind_nrem(ievt)-samples):(SWR_fiber_ind_nrem(ievt)+samples)));
end

% ~ circshift 
avg_fiber_nrem = nanmean(zdF_extract_nrem);
std_fiber_nrem = 2*std(zdF_extract_nrem);

X = prepros_signal; 
N = 1000; % number of circshifts 
K=randi([1 length(prepros_signal)],1, N); % pick a random number between 1 and number of samples ... 100 times 
events_num = length(SWR_time_nrem)-1;

circ_zdF_extract_nrem = zeros(events_num, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
avg_circ_zdF_extract_nrem = zeros(N,seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for iter_circ = 1:1:N % 1 through number N (for 1000 circshifts) 
    Y = circshift(X,K(iter_circ)); % circshift the entire fiber signal based on the random number 
    for ievt = 1:1:length(SWR_time_nrem)-1 % for each SWR event pick out 1-978, pick out that subset of the fiber signal 
       timeset = time((SWR_fiber_ind_nrem(ievt)-samples):(SWR_fiber_ind_nrem(ievt)+samples)); % pick fiber events that are 4 seconds each way
       circ_zdF_extract_nrem(ievt,:) = (Y((SWR_fiber_ind_nrem(ievt)-samples):(SWR_fiber_ind_nrem(ievt)+samples))); %resets every time
    end
    avg_circ_zdF_extract_nrem(iter_circ,:) = nanmean(circ_zdF_extract_nrem); % this is a mean over all the trials .... 
end

circ_avg_fiber_nrem = nanmean(avg_circ_zdF_extract_nrem);
circ_std_fiber_nrem = 2*std(avg_circ_zdF_extract_nrem);

%% find SWR times that are in quiet wakefullness
% ~ for every swr time, check if that is an nrem stage in labels (==3) 
labels_time_center = (0:length(labels)-1)' * 2 + 1;

SWR_time_wake = [];
for swri = 1:1:length(SWR_time_mid)
    label_ind = nearest_idx3(SWR_time_mid(swri),labels_time_center); % index of labels_time 
    if labels(label_ind) == 2
        SWR_time_wake = [SWR_time_wake; SWR_time_mid(swri)];
    end
end

% ~ convert to fiber time index
SWR_fiber_ind_wake = zeros(length(SWR_time_wake),1); 
for i = 1:1:length(SWR_time_wake)
    SWR_fiber_ind_wake(i) = nearest_idx3(SWR_time_wake(i),time);
end

% ~ PETH 
zdF_extract_wake = zeros(length(SWR_time_wake)-1, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
time_extract_nrem =  zeros(length(SWR_time_wake)-1, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for ievt = 1:1:length(SWR_time_wake)-1 %-3 for M600 recording 1; -1 for everyone else because the last SWR sometimes is too close to the end recording
    timeset = time((SWR_fiber_ind_wake(ievt)-samples):(SWR_fiber_ind_wake(ievt)+samples)); % pick fiber events that are 4 seconds each way
    time_extract_wake(ievt,:) = time((SWR_fiber_ind_wake(ievt)-samples):(SWR_fiber_ind_wake(ievt)+samples))-timeset(1); 
    zdF_extract_wake(ievt,:) = (prepros_signal((SWR_fiber_ind_wake(ievt)-samples):(SWR_fiber_ind_wake(ievt)+samples)));
end

% ~ circshift 
avg_fiber_wake = nanmean(zdF_extract_wake);
std_fiber_wake = 2*std(zdF_extract_wake);

X = prepros_signal; 
N = 1000; % number of circshifts 
K=randi([1 length(prepros_signal)],1, N); % pick a random number between 1 and number of samples ... 100 times 
events_num = length(SWR_time_wake)-1;

circ_zdF_extract_wake = zeros(events_num, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
avg_circ_zdF_extract_wake = zeros(N,seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for iter_circ = 1:1:N % 1 through number N (for 1000 circshifts) 
    Y = circshift(X,K(iter_circ)); % circshift the entire fiber signal based on the random number 
    for ievt = 1:1:length(SWR_time_wake)-1 % for each SWR event pick out 1-978, pick out that subset of the fiber signal 
       timeset = time((SWR_fiber_ind_wake(ievt)-samples):(SWR_fiber_ind_wake(ievt)+samples)); % pick fiber events that are 4 seconds each way
       circ_zdF_extract_wake(ievt,:) = (Y((SWR_fiber_ind_wake(ievt)-samples):(SWR_fiber_ind_wake(ievt)+samples))); %resets every time
    end
    avg_circ_zdF_extract_wake(iter_circ,:) = nanmean(circ_zdF_extract_wake); % this is a mean over all the trials .... 
end

circ_avg_fiber_wake = nanmean(avg_circ_zdF_extract_wake);
circ_std_fiber_wake = 2*std(avg_circ_zdF_extract_wake);


%% Figures together 
fig2 = figure(2);

med_c = [104,187,225]./255; % rgb(167, 199, 231) rgb(255, 165, 0)  blue: rgb(104,187,227)
low_c = [78,178,101]./255; % color

subplot(1,2,1)
shadedErrorBar(time_extract_wake(1,:),circ_avg_fiber_wake,circ_std_fiber_wake,'lineProps','-k','transparent',1) % subtract the circ mean here 
hold on
plot(time_extract_wake(1,:),circ_avg_fiber_wake,'LineWidth',2,'Color','k') % subtract the circ mean here 
% plot average on top with larger line
hold on
plot(time_extract_wake(1,:),avg_fiber_wake,'LineWidth',2,'Color',low_c)
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
ylim([-0.5 0.5])
xlim([0 8])
xticks([0 4 8])
xticklabels({'-4','0','4'})
title('Quiet Wakefullness PETH','FontSize', 20)
ylabel('Mean [DA] (z-score)','FontSize', 16)
xlabel('Time from SWR (s)','FontSize', 16)
legend('','shuffle','signal','Location','northwest')
legend boxoff

subplot(1,2,2)
shadedErrorBar(time_extract_wake(1,:),circ_avg_fiber_nrem,circ_std_fiber_nrem,'lineProps','-k','transparent',1) % subtract the circ mean here 
hold on
plot(time_extract_wake(1,:),circ_avg_fiber_nrem,'LineWidth',2,'Color','k') % subtract the circ mean here 
% plot average on top with larger line
hold on
plot(time_extract_wake(1,:),avg_fiber_nrem,'LineWidth',2,'Color',low_c)
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
ylim([-0.5 0.5])
xlim([0 8])
xticks([0 4 8])
xticklabels({'-4','0','4'})
title('NREM Rest PETH','FontSize', 20)
ylabel('Mean [DA] (z-score)','FontSize', 16)
xlabel('Time from SWR (s)','FontSize', 16)
legend('','shuffle','signal','Location','northwest')
legend boxoff

set(gcf,'Color',[1,1,1])
shg

hold off

txt = {'Session Plot: Hypothesis Suplots for Sleep Wake Stage'};
sgtitle(txt)

fig2.WindowState = 'maximized';
 %% Save figure
cd 'C:\Users\mimia\Documents\ReplayDA Figures\M548\recording 7'
saveas(fig2,'M548_sleep_wake_hypothesis7.png') % CHANGE THIS 

avg_SWR_DA.circ_avg_nrem = circ_avg_fiber_nrem;
avg_SWR_DA.circ_avg_wake = circ_avg_fiber_wake;

avg_SWR_DA.circ_std_nrem = circ_std_fiber_nrem;
avg_SWR_DA.circ_std_wake = circ_std_fiber_wake;

avg_SWR_DA.avg_fiber_nrem = avg_fiber_nrem;
avg_SWR_DA.avg_fiber_wake = avg_fiber_wake;

%% Save data
cd 'D:\'
filename = append(file_name, "sleepwake.mat");
save(filename, '-struct','avg_SWR_DA')

