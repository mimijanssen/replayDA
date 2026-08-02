%% Save RPE permutation scores (1,000) and means 

clear; clc;
cd 'F:\M548\M548_2024_08_31_recording7'; 
FP = load('M548_2024_08_31processed'); % fiber data processed with my pipeline
load('M548_2024-08-31_track') % for pseudo_outcomes % CHANGE THIS 
file_name = 'M548_2024_08_31'; 

%% Extract track codes 

prepros_signal = FP.zF_win_60s; % 60 s detrend

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
        indxpb = nearest_idx(all_events{1,iter}(ptime),t);
        init_trial = (indxpb - (window)); 
        end_trial = (indxpb + (window));
            if iter == 1 
                high(ptime,:) = prepros_signal(init_trial:end_trial-1); 
            elseif iter == 3
                low(ptime,:) = prepros_signal(init_trial:end_trial-1); 
            end
    end
end

% average by columns
avg_high = mean(high); 
avg_low = mean(low);

%% Mean values
x_values = 8001:1:12000; % 4 seconds after photobeam break
n = size(low,1);

% mean values for low trials
mean_low = zeros(n,1);
for i_low = 1:1:n
    mean_low(i_low,:) = mean(low(i_low,(x_values))); 
end

sess_mean_low = mean(mean_low);

% mean values for high trials 
mean_high = zeros(n,1);
for i_high = 1:1:n
    mean_high(i_high,:) = mean(high(i_high,(x_values))); 
end

sess_mean_high = mean(mean_high);

% %% Permuations 
% n_perms = 1000;
% rng(42);
% 
% perm_stats_high = zeros(n_perms,1);
% perm_stats_low = zeros(n_perms,1);
% 
% 
% for perm = 1:n_perms
%     k = 1; % to iterate through all permutations? 
%     trial_stats_high = zeros(n,1);
% 
%     % ================= high =================
%     for j = 1:size(high,1) % go through each trial 
% 
%         x = high(j,:);
%         n = length(x); % length of data
% 
%         % Random circular shift
%         shift = randi(n);
%         x = circshift(x,shift);
% 
%         % Random timepoint to average
%         max_idx = n - size(x_values,2);
% 
%         inj = randi([1 max_idx]);
% 
%         trial_stats_high(k) = ...
%             mean(x(inj+1:inj+size(x_values,2))) ...
%              - x(inj);
%             k = k+1;
%         end
% 
%         trial_stats_low = zeros(n,1);
%         k = 1;
%         % ================= Low Trials  =================
%         for j = 1:size(low,1)
% 
%             x = low(j,:);
% 
%             n = length(x);
% 
%             shift = randi(n);
%             x = circshift(x,shift);
% 
%             min_idx = 1;
%             max_idx = n - size(x_values,2);
% 
%             if max_idx <= min_idx
%                 continue
%             end
% 
%             inj = randi([min_idx max_idx]);
% 
%             trial_stats_low(k) = ...
%                 mean(x(inj+1:inj+size(x_values,2))) ...
%                 - x(inj);
%             k = k+1;
%         end
% 
%         % Remove unused entries (only happens if a session was skipped)
%         perm_stats_high(perm) = mean(trial_stats_high);
%         perm_stats_low(perm) = mean(trial_stats_low); 
%     end

%% -------------------- Permutation test --------------------

n_perms = 1000;
rng(42);

window = 8 * FP.cfg.hdr{1,1}.SamplingFrequency;
post_idx = 8001:12000;          % same window you currently average

real_high = zeros(length(high_t),1);
real_low  = zeros(length(low_t),1);

% -------------------- Real statistic --------------------

for j = 1:length(high_t)
    idx = nearest_idx(high_t(j),FP.tvec);
    x = FP.zF_win_60s(idx-window:idx+window-1);
    real_high(j) = mean(x(post_idx));
end

for j = 1:length(low_t)
    idx = nearest_idx(low_t(j),FP.tvec);
    x = FP.zF_win_60s(idx-window:idx+window-1);
    real_low(j) = mean(x(post_idx));
end

real_stat = mean(real_high) - mean(real_low);

% -------------------- Null distribution --------------------

perm_stats = zeros(n_perms,1);

signal = FP.zF_win_60s;
N = length(signal);

for perm = 1:n_perms

    % Circularly shift ENTIRE recording
    shift = randi(N);
    shifted = circshift(signal,shift);

    % High trials
    high_vals = zeros(length(high_t),1);

    for j = 1:length(high_t)
        idx = nearest_idx(high_t(j),FP.tvec);
        if idx-window < 1 || idx+window-1 > N
            continue
        end
        x = shifted(idx-window:idx+window-1);
        high_vals(j) = mean(x(post_idx));
    end

    % Low trials
    low_vals = zeros(length(low_t),1);

    for j = 1:length(low_t)
        idx = nearest_idx(low_t(j),FP.tvec);
        if idx-window < 1 || idx+window-1 > N
            continue
        end
        x = shifted(idx-window:idx+window-1);
        low_vals(j) = mean(x(post_idx));
    end

    %% Statistic

    perm_stats(perm) = mean(high_vals) - mean(low_vals);

end

% -------------------- Empirical p-value --------------------

p_one = (sum(perm_stats >= real_stat)+1)/(n_perms+1);

fprintf('\nObserved statistic = %.3f\n',real_stat)
fprintf('Permutation mean = %.3f\n',mean(perm_stats))
fprintf('One-tailed p = %.4f\n',p_one)

% -------------------- Plot --------------------

figure

histogram(perm_stats,50,...
    'FaceColor',[0.7 0.7 0.7],...
    'EdgeColor','none')

hold on

xline(real_stat,'r','LineWidth',2)

xline(mean(perm_stats),'k--','LineWidth',2)

xlabel('High - Low response')
ylabel('Count')

title(sprintf('Permutation test   p = %.4f',p_one))

box off
set(gca,'FontSize',12)

%%
cd 'F:\M548\avg_data\RPE_mean';
filename = append(file_name, "RPE_mean.mat");
save(filename, 'perm_stats','real_stat');

