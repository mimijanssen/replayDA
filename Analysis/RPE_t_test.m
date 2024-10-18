%% Save RPE t-tests 
clear; clc;
cd 'F:\M433\M433_2023_09_19_recording1'; 
FP = load('M433_2023_09_19processed'); % fiber data processed with my pipeline
load('M433_2023-09-19_track') % for pseudo_outcomes % CHANGE THIS 
file_name = 'M433_2023_09_19'; 

%% Extract track codes 

prepros_signal = FP.zF_win_60s; 

LoadExpKeys();
cfg_evt = [];
cfg_evt.eventList = ExpKeys.eventList;
cfg_evt.eventLabel = ExpKeys.eventLabel;
evt3 = LoadEvents(cfg_evt);

cfg_fiber.fc = {'CSC30.ncs'};
csc_photo_fc = LoadCSC(cfg_fiber);

evt_ordered = (sort([evt3.t{1},evt3.t{2}])); % only want 60 of them
evt_ordered = evt_ordered(1:60);
photobeam_times = evt_ordered - csc_photo_fc.tvec(1); % initialize time 

high_log = pseudo_outcomes(1:length(pseudo_outcomes)) == 3; 
med_log = pseudo_outcomes(1:length(pseudo_outcomes)) == 2; 
low_log = pseudo_outcomes(1:length(pseudo_outcomes)) == 1; 

high_t = photobeam_times(high_log(1:length(photobeam_times))); % s
med_t = photobeam_times(med_log(1:length(photobeam_times))); % should be 36 not 52
low_t = photobeam_times(low_log(1:length(photobeam_times))); % 6

% Color
med_c = [104,187,225]./255; % rgb(167, 199, 231) rgb(255, 165, 0)  blue: rgb(104,187,227)
low_c = [0,104,87]./255; % I actually want to change this to purple...
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
   
        %indxpb = find(abs(t-all_events{1,iter}(ptime)) < 0.0005); 
        % for some reason this failed for some points in M460. So I tried a better
        % method: nearest_idx3 
        indxpb = nearest_idx(all_events{1,iter}(ptime),t);

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


%% Find Max values and AUC values. 
x_values = 8001:1:16000;
n = size(low,1);

% area under the curve for low trials
AUC_low = zeros(n,1);
for i_low = 1:1:n
    AUC_low(i_low,:) = trapz(x_values, low(i_low,x_values));
end

% area under the curve for high trials
AUC_high = zeros(n,1);
for i_high = 1:1:n
    AUC_low(i_high,:) = trapz(x_values, high(i_high,x_values));
end

% max values for high trials 
dF_low = zeros(n,1);
for i_low = 1:1:n
    dF_low(i_low,:) = max(min(i_low,(x_values))); 
end

% min values for low trials 
dF_high = zeros(n,1);
for i_high = 1:1:n
    dF_low(i_high,:) = max(high(i_high,(x_values))); 
end

% paired t-test
AUC_low_ttest = 

%%
cd 'F:\M433\avg_data\RPE_t-test';
filename = append(file_name, "RPE_Ttest.mat");
save(filename, 'fpos','linspd','spd_post','spd_pre','spd_track','speed_track')

