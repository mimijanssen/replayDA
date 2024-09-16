%% Plot Raw Fiber Data
% example: cd = C:\Data\M406\2022-10-30_M406_L_sleep_track_train
%folder = 'C:\Data\M406\2022-10-30_M406_L_sleep_track_train' ;
%cd folder 
clear all; clc
cd C:\Data\2023-02-01_M411_saline_track\2023-02-01_M383_rack_track;  LoadExpKeys; %LoadMetadata;

%% Path
addpath(genpath('Users\mimia\Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\mimia\Documents\GitHub\replay_DA\analysis\photometry'));
addpath(genpath('C:\Users\mimia\Documents\GitHub\std_shade'));
%% Load Data
cfg.fc = {'CSC30.ncs'};
csc_photo = LoadCSC(cfg);

%% Plot 
FP_data = [];
FP_data.acq.Fs = csc_photo.cfg.hdr{1}.SamplingFrequency; % set FP_data.acq.Fs to sampling frequency rate (5000 points per second) 
FP_data.acq.time = csc_photo.tvec - csc_photo.tvec(1); %set FP_data.acq.Fs time to the time vector subtracted by the first time point
% this initializes the time vector to start at 0 then 2.0 x 10 ^-4 or 0.0002 and 4.0 x
%x 10 ^ -4 0.0004 ... 0.0002 second time interval means 5,000 points per
%second
FP_data.acq.FP{1} = csc_photo.data'; %fiber data 

% rename variables
FP = (FP_data.acq.FP{1})*(1239/39); %1239/39 is the voltage multiplier
time = FP_data.acq.time; 

%% Detrend Raw Data with Locdetrend
FP_detrended = locdetrend(FP);

%% Load Events
LoadExpKeys();
cfg_evt = [];
cfg_evt.eventList = ExpKeys.eventList;
cfg_evt.eventLabel = ExpKeys.eventLabel;
evt = LoadEvents(cfg_evt);

%% Plot Detrended Data IF YOU ADMINISTERED A DRUG 
% To plot raw data use FP and multiply by 31.77 (multiplier b/c of voltage
% divider)
sessionTitle = 'CW_';
last_time = length(time); %the value of the last time point is how many seconds the recording was
timerange1 = 10/0.0002; %datapoint range for 10 s
timerange2= 100/0.0002; %datapoint range for 100 s 
timerange3= length(time); % datapoint range for all data
time_ranges = [timerange1, timerange2, timerange3]; %in seconds 
%Fs = FP_data.final.Fs;

% INDEX FOR DRUG ADMIN
drug_admin = evt.t{1,5} - csc_photo.tvec(1); % initialzied time of drug administration
admin_done = evt.t{1,6} - csc_photo.tvec(1);
%indinit = find(time==first_time);
% tolerance 
ind_drug = find(abs(time-drug_admin) < 0.0001);
%indfin = find(time ==last_time);
ind_admin = find(abs(time-admin_done) < 0.0001);
c = [0.8 0.7 0.8];
b = [0.3010 0.7450 0.9330];

% INDEX FOR TRACK 
% and last time. 
%first_time = evt.t{1,1}(1) - csc_photo.tvec(1);
%last_time = evt.t{1,2}(20) - csc_photo.tvec(1); 
% tolerance 
%indinit = find(abs(time-first_time) < 0.0001);
%indfin = find(time ==last_time);
%indfin = find(abs(time-last_time) < 0.0001);

% time(indinit:indfin))-time(indinit)
figure(1)
for t_i = 1:length(time_ranges)
    t_range = 1:time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(time(t_range), FP_detrended(t_range), 'Color', [0 0.5 0])
    hold on
    %set(gcf,'color','none');
    if t_i ==3 
        y = ylim;
        %plot([time(ind_drug) time(ind_drug)],[y(1) y(2)])
        %plot([time(ind_admin) time(ind_admin)],[y(1) y(2)])
        %L2 = patch([linspace(time(indinit),time(indfin),40) fliplr(linspace(time(indinit),time(indfin),40))],[(linspace(y(1),y(1),40)) fliplr((linspace(y(2),y(2),40)))],b,'EdgeColor','none');
        %uistack(L2,'bottom')
        %hold on
        L = patch([linspace(time(ind_drug),time(ind_admin),40) fliplr(linspace(time(ind_drug),time(ind_admin),40))],[(linspace(y(1),y(1),40)) fliplr((linspace(y(2),y(2),40)))],c,'EdgeColor','none');
        uistack(L,'bottom')
        legend('Injection','Location','southeast')
        legend('boxoff')
        hold off
    end
    title([sessionTitle, num2str(time_ranges(t_i)),'samples'], 'Interpreter','none')
    ylabel('Fiber Signal (V)'); xlabel('Time (s)');
end

%[selected_x,selected_y] = ginput(4); % array of selected x times. 

% plot a figure around those x times
%selected_x_end = nearest_idx3(selected_x, time) + (5/0.0002);
%selected_x_start = nearest_idx3(selected_x, time) - (5/0.0002);

%figure(2)
%for t_i2 = 1:4
%    subplot(4, 1, t_i2);
%    plot((time(selected_x_start(t_i2):selected_x_end(t_i2)))-time(selected_x_start(t_i2)), FP_detrended(selected_x_start(t_i2):selected_x_end(t_i2)), 'Color', [0 0.5 0])
%title(['10 second window around selected points (x=5s is the selected time)'])
%ylabel('Fiber Signal (V)'); xlabel('Time (s)');
%end

%% Plot Detrended Data PRE-EXP 
% To plot raw data use FP and multiply by 31.77 (multiplier b/c of voltage
% divider)
sessionTitle = 'CW_';
last_time = length(time); %the value of the last time point is how many seconds the recording was
timerange1 = 10/0.0002; %datapoint range for 10 s
timerange2= 100/0.0002; %datapoint range for 100 s 
timerange3= length(time); % datapoint range for all data
time_ranges = [timerange1, timerange2, timerange3]; %in seconds 
%Fs = FP_data.final.Fs;


% time(indinit:indfin))-time(indinit)
figure(1)
for t_i = 1:length(time_ranges)
    t_range = 1:time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(time(t_range), FP_detrended(t_range), 'Color', [0 0.5 0])
    hold on
    %set(gcf,'color','none');
    title([sessionTitle, num2str(time_ranges(t_i)),'samples'], 'Interpreter','none')
    ylabel('Fiber Signal (V)'); xlabel('Time (s)');
end

%[selected_x,selected_y] = ginput(4); % array of selected x times. 



%% Plot RAW Data
% To plot raw data use FP and multiply by 31.77 (multiplier b/c of voltage
% divider)
sessionTitle = 'CW_';
last_time = length(time); %the value of the last time point is how many seconds the recording was
timerange1 = 10/0.0002; %datapoint range for 10 s
timerange2= 100/0.0002; %datapoint range for 100 s 
timerange3= length(time); % datapoint range for all data
time_ranges = [timerange1, timerange2, timerange3]; %in seconds 
%Fs = FP_data.final.Fs;

figure(1)
for t_i = 1:length(time_ranges)
    t_range = 1:time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(time(t_range), FP(t_range), 'Color', [0 0.5 0])
    title([sessionTitle, num2str(time_ranges(t_i)),'samples'], 'Interpreter','none')
    ylabel('Fiber Signal (V)'); xlabel('Time (s)');
end

[selected_x,selected_y] = ginput(4); % array of selected x times. 

% plot a figure around those x times
selected_x_end = nearest_idx3(selected_x, time) + (5/0.0002);
selected_x_start = nearest_idx3(selected_x, time) - (5/0.0002);

figure(2)
for t_i2 = 1:4
    subplot(4, 1, t_i2);
    plot((time(selected_x_start(t_i2):selected_x_end(t_i2)))-time(selected_x_start(t_i2)), FP(selected_x_start(t_i2):selected_x_end(t_i2)), 'Color', [0 0.5 0])
title(['10 second window around selected points (x=5s is the selected time)'])
ylabel('Fiber Signal (V)'); xlabel('Time (s)');
end
%% Plot Detrended Linear Track Data 
% Need an if then statement here to make sure left is first and if not to
% switch the time . alternatively you can sort this and then take the first
% and last time. 
first_time = evt.t{1,1}(1) - csc_photo.tvec(1);
last_time = evt.t{1,2}(20) - csc_photo.tvec(1); 

%indinit = find(time==first_time);
% tolerance 
indinit = find(abs(time-first_time) < 0.0001);
%indfin = find(time ==last_time);
indfin = find(abs(time-last_time) < 0.0001);

plot((time(indinit:indfin))-time(indinit), FP_detrended(indinit:indfin), 'Color', [0 0.5 0])

title('Fiber Photometry Signal on the Linear Track')
ylabel('Fiber Signal (V)'); xlabel('Time (s)');


%% Plot Individual Trial Data on Top of Each Other (DETRENDED) 20 trials
clear indxpb init_trial end_trial ptime

% Color matrix 
count = 0.9 ;
jet = zeros(20,3); % 30
for zi = 1:1:20 %39
    count =  count - 0.04; %0.040; %gets darker over trials %0.020 for 40 trials 
    jet(zi,2) = count;
end

% Plot
evt_ordered = sort([evt.t{1},evt.t{2}]); 
photobeam_times = evt_ordered - csc_photo.tvec(1); % initialize time 

figure(1)
for ptime = 1:1: length(photobeam_times)
    % tolerance 
    indxpb = find(abs(time-photobeam_times(ptime)) < 0.0001); % why does this keep chaning???
    init_trial = (indxpb - (8/0.0002)); % 8 seconds; remember 0.0002 is the time in seconds for each voltage point 
    end_trial = (indxpb + (8/0.0002));
    %if (max(left_end_peth.data*32) - min(left_end_peth.data*32)) < 0.15
    if abs(max(FP_detrended(init_trial:end_trial))) + abs(min(FP_detrended(init_trial:end_trial))) < 0.4

    plot((time(init_trial:end_trial))-time(init_trial), FP_detrended(init_trial:end_trial), 'Color', [jet(ptime,:)])
    else 
    disp('removed artifact')
    end
    hold on
colormap(jet);
cbh = colorbar;
cbh.Ticks = linspace (0,1,5); % 9 - need to change this based on length
cbh.TickLabels = num2cell(0:5:20); % 30; 40 - need to change this based on length
cbh.Label.String = 'Trial';
xlabel('time (s)')
xlim([0 16])
xticks([0 8 16])
xticklabels({'-8','0','8'})
ylabel('F')
title('F PET with 16s window')
end

% Difference between running times: 
running_t = diff(evt_ordered);
avg_run_t = mean(running_t);
med_run_t = median(running_t);
min(running_t) 
max(running_t)
%% Plot Individual Trial Data on Top of Each Other (DETRENDED)
clear indxpb init_trial end_trial ptime

% Color matrix 
count = 0.9 ;
jet = zeros(50,3); 
for zi = 1:1:50
    count =  count - 0.015;%0.040; %gets darker over trials %0.020 for 40 trials 
    jet(zi,2) = count;
end

% Plot
evt_ordered = sort([evt.t{1},evt.t{2}]); 
photobeam_times = evt_ordered - csc_photo.tvec(1); % initialize time 

figure(1)
for ptime = 1:1: length(photobeam_times)
    % tolerance 
    indxpb = find(abs(time-photobeam_times(ptime)) < 0.00005);
    init_trial = (indxpb - (8/0.0002));
    end_trial = (indxpb + (8/0.0002)); % WAS I CHANING THIS????.... 
    %if (max(left_end_peth.data*32) - min(left_end_peth.data*32)) < 0.15
    if abs(max(FP_detrended(init_trial:end_trial))) + abs(min(FP_detrended(init_trial:end_trial))) < 0.25

    plot((time(init_trial:end_trial))-time(init_trial), FP_detrended(init_trial:end_trial), 'Color', [jet(ptime,:)])
    else 
    disp('removed artifact')
    end
    hold on
colormap(jet);
cbh = colorbar;
cbh.Ticks = linspace (0,1,11); % 9 - need to change this based on length
cbh.TickLabels = num2cell(0:5:50); % 30; 40 - need to change this based on length
cbh.Label.String = 'Trial';
xlabel('time (s)')
xlim([0 16])
xticks([0 8 16])
xticklabels({'-8','0','8'})
ylabel('F')
title('F PET with 16s window')
end

% Difference between running times: 
running_t = diff(evt_ordered);
avg_run_t = mean(running_t);
med_run_t = median(running_t);
min(running_t) 
max(running_t)


%% Plot Individual Trial Data on Top of Each Other (RAW)
clear indxpb init_trial end_trial ptime

% Color matrix 
count = 0.9 ;
jet = zeros(40,3);
for zi = 1:1:39
    count =  count - 0.020; %gets darker over trials 
    jet(zi,2) = count;
end

% Plot
evt_ordered = sort([evt.t{1},evt.t{2}]); 
photobeam_times = evt_ordered - csc_photo.tvec(1); % initialize time 

figure(1)
for ptime = 1:1: length(photobeam_times)
    % tolerance 
    indxpb = find(abs(time-photobeam_times(ptime)) < 0.0001);
    init_trial = (indxpb- (2/0.0002));
    end_trial = (indxpb + (2/0.0002));
    %if (max(left_end_peth.data*32) - min(left_end_peth.data*32)) < 0.15
   % if abs(max(FP(init_trial:end_trial))) + abs(min(FP(init_trial:end_trial))) < 0.3

    plot((time(init_trial:end_trial))-time(init_trial), FP(init_trial:end_trial), 'Color', [jet(ptime,:)])
   %else 
   %continue
   % end
    hold on
colormap(jet);
cbh = colorbar;
cbh.Ticks = linspace (0,1,9); % 9 - need to change this based on length
cbh.TickLabels = num2cell(0:5:40); % 40 - need to change this based on length
cbh.Label.String = 'Trial';
xlabel('time (s)')
ylabel('F')
xlim([0 4])
xticks([0 2 4])
xticklabels({'-2','0','2'})
title('F PET with 4s window')
end
%% Plot Individual Trial Data Aligned to Largest Peak
% Color matrix 
clear indxpb init_trial end_trial ptime

figure(2)
for ptime = 1:1:length(photobeam_times)
    % tolerance 
    hold on
    indxpb = find(abs(time-photobeam_times(ptime)) < 0.0001);
    init_trial = (indxpb- (5/0.0002));
    end_trial = (indxpb + (5/0.0002));
    signal_max = max(FP_detrended(init_trial:end_trial)); 
    signal_int = FP_detrended(init_trial:end_trial);
    signal_max_index = find(abs(signal_int - signal_max) < 0.0001); %index within the interval
    % add this index to init_trial 
    %establish a new index around max index 
    init_signal = init_trial + signal_max_index - (5/0.0002);
    end_signal = init_trial + signal_max_index + (5/0.0002);
    if abs(max(FP_detrended(init_signal:end_signal))) + abs(min(FP_detrended(init_signal:end_signal))) < 0.4
    plot((time(init_signal:end_signal)-time(init_trial)), FP_detrended(init_signal:end_signal), 'Color', [jet(ptime,:)])
    else
        %nothing
    end
    hold on
    colormap(jet);
cbh = colorbar;
cbh.Ticks = linspace (0,1,9); % 9 - need to change this based on length
cbh.TickLabels = num2cell(0:5:40); % 40 - need to change this based on length
cbh.Label.String = 'Trial';
xlabel('time (s)')
ylabel('F')
xlim([0 10])
title('F PET with 10s window')
end

%% Split data and plot average of one half and average of the other half
% Plot
time_split1 = [];
time_split2 = [];
FP_detrended_split1 =[];
FP_detrended_split2 = [];

evt_ordered = sort([evt.t{1},evt.t{2}]); 
num_trials = length(evt_ordered);

evt_split1 = evt_ordered(1:num_trials/2);
evt_split2 = evt_ordered(num_trials/2 + 1 : num_trials);

photobeam_times_split1 = evt_split1 - csc_photo.tvec(1); % initialize time 
photobeam_times_split2 = evt_split2 - csc_photo.tvec(1); % initialize time 

for ptime = 1:1: length(evt_ordered)-1
    % tolerance 
    indxpb = find(abs(time-photobeam_times(ptime)) < 0.0001);
    init_trial = (indxpb- (5/0.0002));
    end_trial = (indxpb + (5/0.0002));
    %if (max(left_end_peth.data*32) - min(left_end_peth.data*32)) < 0.15
    if ptime < num_trials/2
    time_split1(ptime,:) = (time(init_trial:end_trial))-time(init_trial);
    FP_detrended_split1(ptime,:) = FP_detrended(init_trial:end_trial);
    else 
    time_split2(ptime-20,:) = (time(init_trial:end_trial))-time(init_trial);
    FP_detrended_split2(ptime-20,:) = FP_detrended(init_trial:end_trial); 
   %else 
   %continue
    end
end

figure(3)
%plot(mean(FP_detrended_split1))
stdshade(FP_detrended_split1,0.5,'g');
hold on 
stdshade(FP_detrended_split2,0.5,'r');
%plot(mean(FP_detrended_split2))
xlabel('time (s)')
ylabel('averaged F')
%xlim([0 10])
legend('first half','second half')
title('F PET with 10s window')
hold off

SEM_split1 = std((FP_detrended_split1))/sqrt(20);
SEM_split2 = std((FP_detrended_split2))/sqrt(20);
y = linspace(-5,5,11);

figure (4)
plot(mean(FP_detrended_split1))
hold on
plot(mean(FP_detrended_split2))
hold on
shadedErrorBar([],mean(FP_detrended_split1),SEM_split1','lineProps','g'); 
hold on
shadedErrorBar([],mean(FP_detrended_split2),SEM_split2'); 

xlabel('time (s)')
ylabel('averaged F')
legend('','','first half','second half')
xticks([0 25000 50001])
xticklabels({'-5','0','5'})
title('F PET with 10s window')
hold off

%SEM = std(x)/sqrt(length(x));

%% Plot Wavesurfer over NYX
load('M406_20221030_0001.mat')

% for this you should subtract from all values the initial time 
plot(time(indinit:indfin)-time(indinit), FP(indinit:indfin), 'Color', [0 0.5 0])
hold on
plot(data.acq.time,data.acq.FP{1,1})

title('Linear Track Signal')
ylabel('Fiber Signal (V)'); xlabel('Time (s)');
legend('Neuralynx','Wavesurfer')
hold off

