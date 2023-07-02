%% Load preprocessed data
clear; clc;
cd 'C:\Data\M437\M437-2023-06-24_track';
load('M437-2023-06-24processed');


%% Extract linear track events
LoadExpKeys();
cfg_evt = [];
cfg_evt.eventList = ExpKeys.eventList;
cfg_evt.eventLabel = ExpKeys.eventLabel;
evt = LoadEvents(cfg_evt);

cfg.fc = {'CSC30.ncs'};
csc_photo = LoadCSC(cfg);

% Load the pseudo prob data 
load('M437_2023-06-24_track.mat') % for pseudo_outcomes

%% Plot Detrended Linear Track Data 
% Need an if then statement here to make sure left is first and if not to
% switch the time . alternatively you can sort this and then take the first
% and last time. 

% order left and right reward data to one row 
if evt.t{1,1}(1) > evt.t{1,2}(1)
    first_time = evt.t{1,1}(1) - csc_photo.tvec(1);
    last_time = evt.t{1,2}(30) - csc_photo.tvec(1);
    left = true ; 
    right = false;
else 
    first_time = evt.t{1,2}(1) - csc_photo.tvec(1);
    last_time = evt.t{1,1}(30) - csc_photo.tvec(1); 
    right = true; 
    left = true; 
end

% extract data of which trial is high, low, medium reward  

%indinit = find(time==first_time);
% tolerance 
indinit = find(abs(t-first_time) < 0.0005);
%indfin = find(time ==last_time);
indfin = find(abs(t-last_time) < 0.0005);

% Gray : [128 133 133]./255
% Dark Green = [0,104,87]./225;
% Light Green = color_s2 = [0,104,87]./225;

plot((t(indinit:indfin))-t(indinit), zdF(indinit:indfin), 'Color', [0,104,87]./225)

title('Fiber Photometry Signal on the Linear Track')
ylabel('Signal (dF z-scored)'); xlabel('Time (s)');


%% Plot Individual Trial Data on Top of Each Other (DETRENDED) 60 trials
clear indxpb init_trial end_trial ptime

% Color matrix 
count = 0.9 ;
jet = zeros(20,3); % 30
for zi = 1:1:60 %39
    count =  count - 0.01; %0.040; %gets darker over trials %0.020 for 40 trials 
    jet(zi,2) = count;
end

% Plot
evt_ordered = sort([evt.t{1},evt.t{2}]); % don't need to sort anymore but ok % NEED code to determine which is first... ??
photobeam_times = evt_ordered - csc_photo.tvec(1); % initialize time 

figure(1)
for ptime = 1:1:length(photobeam_times)
    % tolerance 
    indxpb = find(abs(t-photobeam_times(ptime)) < 0.0005); % why does this keep chaning???
    init_trial = (indxpb - (8/0.001)); % 8 seconds; remember 0.0002 is the time in seconds for each voltage point 
    end_trial = (indxpb + (8/0.001));
    %if (max(left_end_peth.data*32) - min(left_end_peth.data*32)) < 0.15
    %if abs(max(zdF(init_trial:end_trial))) + abs(min(zdF(init_trial:end_trial))) < 30
        plot((t(init_trial:end_trial))-t(init_trial), zdF(init_trial:end_trial), 'Color', [jet(ptime,:)])
   % else 
   %     disp('removed artifact')
   % end
    hold on
colormap(jet);
cbh = colorbar;
cbh.Ticks = linspace (0,1,5); % 9 - need to change this based on length
cbh.TickLabels = num2cell(0:5:60); % 30; 40 - need to change this based on length
cbh.Label.String = 'Trial';
xlabel('time (s)')
xlim([0 16])
xticks([0 8 16])
xticklabels({'-8','0','8'})
ylabel('F')
title('F PET with 16s window')
end

%% Break out linear track data by probability
high_t = photobeam_times(pseudo_outcomes == 3); % should be more than 2... I messed this up i think I messed up high reward I guess...
med_t = photobeam_times(pseudo_outcomes == 2); % should be 36 not 52
low_t = photobeam_times(pseudo_outcomes == 1); % 6 

%% Plot values 
% Color
med_c = [204,204,225]./255; % rgb(167, 199, 231) rgb(255, 165, 0)
low_c = [0,104,87]./255; 
high_c = [255, 165,0]./255; % rgb(204, 204, 255) periwinkle

all_events = [{high_t}, {med_t}, {low_t}];
all_colors = [{high_c},{med_c},{low_c}];

figure(2) 
for iter = 1:1:3
    hold on
    for ptime = 1:1:length(all_events{1,iter})
        indxpb = find(abs(t-all_events{1,iter}(ptime)) < 0.0005); % find the time of the event. 
        init_trial = (indxpb - (8/0.001)); % 8 seconds; remember 0.0002 is the time in seconds for each voltage point
        end_trial = (indxpb + (8/0.001));
        %if abs(max(zdF(init_trial:end_trial))) + abs(min(zdF(init_trial:end_trial))) < 30
            if iter ==1 
                h1 = plot((t(init_trial:end_trial))-t(init_trial), zdF(init_trial:end_trial), 'Color',[all_colors{1,iter}]);
            elseif iter == 2
                h2 = plot(t(init_trial:end_trial)-t(init_trial), zdF(init_trial:end_trial), 'Color',[all_colors{1,iter}]); % set t back to 0
            else
                h3 = plot((t(init_trial:end_trial))-t(init_trial), zdF(init_trial:end_trial), 'Color',[all_colors{1,iter}]);
            end
        %else
            disp('removed artifact')
        %end
    end
end
xlabel('time (s)')
xlim([0 16])
xticks([0 8 16])
xticklabels({'-8','0','8'})
legend([h1, h2, h3],'high','medium','low');
legend boxoff
ylabel('Signal (dF z-scored)')
title('16s time window of fiber signal centered around photobeam break')

%% Plot Averaged Values 
addpath(genpath('C:\Users\mimia\Documents\GitHub\shadedErrorBar'));

% Color
med_c = [104,187,225]./255; % rgb(167, 199, 231) rgb(255, 165, 0)  blue: rgb(104,187,227)
low_c = [0,104,87]./255; 
high_c = [255, 165,0]./255; % rgb(204, 204, 255) periwinkle

all_events = [{high_t}, {med_t}, {low_t}];
all_colors = [{high_c},{med_c},{low_c}];

window = 8/0.001; 

% initialize values for average matrix 
high = zeros(length(high_t),window*2); %rows = trial, columns = signal 
med = zeros(length(med_t),window*2);
low = zeros(length(low_t),window*2);

figure(2)
for iter = 1:1:3
    hold on
    for ptime = 1:1:length(all_events{1,iter})
        indxpb = find(abs(t-all_events{1,iter}(ptime)) < 0.0005); % find the time of the event. 
        init_trial = (indxpb - (window)); % 8 seconds; remember 0.0002 is the time in seconds for each voltage point
        end_trial = (indxpb + (window));
        if abs(max(zdF(init_trial:end_trial))) + abs(min(zdF(init_trial:end_trial))) < 30
            % save each value here to average 
            if iter ==1 
                high(ptime,:) = zdF(init_trial:end_trial-1); 
                %h1 = plot((t(init_trial:end_trial))-t(init_trial), zdF(init_trial:end_trial), 'Color',[all_colors{1,iter}]);
            elseif iter == 2
                med(ptime,:) = zdF(init_trial:end_trial-1); 
                %h2 = plot(t(init_trial:end_trial)-t(init_trial), zdF(init_trial:end_trial), 'Color',[all_colors{1,iter}]); % set t back to 0
            else
                low(ptime,:) = zdF(init_trial:end_trial-1); 
                %h3 = plot((t(init_trial:end_trial))-t(init_trial), zdF(init_trial:end_trial), 'Color',[all_colors{1,iter}]);
            end
        else
            disp('removed artifact')
        end
    end
    t_shared = t(init_trial:end_trial-1)-t(init_trial); % double the window
end

% average by columns
avg_high = mean(high); 
avg_med = mean(med);
avg_low = mean(low);

% std by columns
std_high = std(high); 
std_med = std(med);
std_low = std(low);

shadedErrorBar(t_shared,avg_high,std_high,'lineprops',{'-','color',high_c,'MarkerFaceColor',high_c});
hold on
shadedErrorBar(t_shared,avg_med,std_med,'lineprops',{'-','color',med_c,'MarkerFaceColor',med_c});
shadedErrorBar(t_shared,avg_low,std_low,'lineprops',{'-','color',low_c,'MarkerFaceColor',low_c});

xlabel('Time (s)')
xlim([0 16])
xticks([0 8 16])
xticklabels({'-8','0','8'})
legend('high','medium','low');
legend boxoff
ylabel('Mean Signal (dF z-scored)')
title('16s Time Window of Fiber Signal Centered Around Photobeam Break')

set(gcf,'Color',[1,1,1])
shg

orient(gcf,'landscape')

%% SAVE DATA
file_name = 'M437-2023-06-24';
filename = append(file_name, "4avgsess.mat");
data.t = t;
data.zdF = zdF;
data.evtt = all_events;
data.F_detrend = F_detrend;
data.FP_z = FP_z;
data.dF = dF; 
data.zDF_mov = zDF_mov;
data.window = window; %seconds 
save(filename, '-struct','data')