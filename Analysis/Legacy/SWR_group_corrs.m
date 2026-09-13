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
cd 'F:\M533\M533_2024_08_21_recording3'; 
FP = load('M533_2024_08_21processed.mat'); % fiber data processed with my pipeline
load('M533_2024_08_21detectedSWRs.mat') % SWR intervals
load('M533_2024-08-21_track.mat') % for pseudo_outcomes % CHANGE THIS 
P = readtable('M533_2024_08_21-convertedDLC_resnet50_Linear TrackApr5shuffle1_100000.csv','PreserveVariableNames',true); % CHANGE THIS 

file_name = 'M533_2024_08_21'; 

%%
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

%%
% Step 1: Calculate time difference between SWRs
SWR_intervals = diff(SWR_time_mid); % Time difference between consecutive SWRs
figure;
histogram(SWR_intervals, 'BinWidth', 0.5); % Plot frequency distribution
xlim([0 10])
xlabel('Inter-SWR Interval (s)');
ylabel('Frequency');
title('Frequency Distribution of Inter-SWR Intervals');

%%

% Step 2: Define groups of SWRs based on a threshold (e.g., 2 seconds)
threshold = 4; % Adjust based on your data
group_idx = find(SWR_intervals < threshold); % Find indices of close SWRs

% Initialize variables to store grouped SWRs
SWR_groups = {}; % Cell array to store groups
group = 1;
SWR_groups{group} = SWR_time_mid(1); % First SWR in the first group

for i = 2:length(SWR_time_mid)
    if SWR_time_mid(i) - SWR_time_mid(i-1) < threshold
        SWR_groups{group} = [SWR_groups{group}, SWR_time_mid(i)]; % Add SWR to the group
    else
        group = group + 1; % Start a new group
        SWR_groups{group} = SWR_time_mid(i); % Add SWR to the new group
    end
end

%%

% Step 3: Perform linear regression between group size and fiber photometry signal
peak_signal = zeros(length(SWR_groups), 1);
group_size = zeros(length(SWR_groups), 1);

for i = 1:length(SWR_groups)
    group_size(i) = length(SWR_groups{i}); % Number of SWRs in the group
    % Find the middle SWR in the group
    middle_SWR_time = SWR_groups{i}(ceil(length(SWR_groups{i}) / 2)); % Middle SWR in the group
    fiber_idx = nearest_idx3(middle_SWR_time, time); % Find the closest fiber index
    window_idx = fiber_idx:fiber_idx + round(4 / mean(diff(time))); % 4-second window
    try
        peak_signal(i) = max(FP.zF_win_60s(window_idx)); % Get the peak signal in the window
    catch
        disp('out of bounds fiber signal after swr')
    end
end

% Plot group size vs peak signal
figure;
scatter(group_size, peak_signal, 'filled');
xlabel('Number of SWRs in Group');
ylabel('Peak Fiber Photometry Signal');
title('Linear Regression: SWR Groups vs Fiber Photometry Peak');
hold on 

% Perform linear regression
mdl = fitlm(group_size, peak_signal);
disp(mdl);

% Plot the regression line
x_range = linspace(min(group_size), max(group_size), 100); % Generate x values for plotting
y_fit = predict(mdl, x_range'); % Use the model to predict y values
plot(x_range, y_fit, 'r-', 'LineWidth', 2); % Plot the regression line

% Add the adjusted R^2 value to the plot
adj_R_squared = mdl.Rsquared.Adjusted; % Extract adjusted R^2
text(max(group_size) * 0.7, max(peak_signal) * 0.9, ...
    ['Adjusted R^2 = ', num2str(adj_R_squared, '%.3f')], ...
    'FontSize', 12, 'Color', 'k');

hold off;
