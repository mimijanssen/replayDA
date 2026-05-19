%% Mouse Averages for each session (Track days 1 - 6) . 
clear; clc;
cd 'D:\Reward_Omit'
Files=dir('*.*');

nDays = 6;
data_by_day = cell(1, nDays);

% Loading
for k=3:length(Files)
   FileNames=Files(k).name;
    
   % Extract day number
   day_num = regexp(FileNames, 'day(\d)', 'tokens');
   day_num = str2double(day_num{1}{1});
    
   % Load file
   S = load(FileNames);
    
   if isfield(S, 'avg_fiber_post') % I BE SWITCHING THIS OUT
       vec = S.avg_fiber_post(:)';  % force row vector
   else
       fn = fieldnames(S);
       vec = S.(fn{1}).avg_fiber_post;
   end
    data_by_day{day_num} = [data_by_day{day_num}; vec];
end


%% Averages and std per day
avg_by_day = cell(1, nDays);

for d = 1:nDays
    avg_by_day{d} = mean(data_by_day{d}, 1); % average across sessions
end

avg_by_day = cell(1, nDays);
std_by_day = cell(1, nDays);
sem_by_day = cell(1, nDays);


for d = 1:nDays
    avg_by_day{d} = mean(data_by_day{d}, 1);
    std_by_day{d} = std(data_by_day{d}, 0, 1);
    n = size(std_by_day{d},1);
    sem_by_day{d} = std(data_by_day{d},0,1) / sqrt(4);
end

%% time points
nPoints = length(avg_by_day{1,1});
t = linspace(-5, 5, nPoints);
%%
figure; hold on;

cmap1 = [ ...
    0.8 1.0 0.8;   % light green
    0.4 0.8 0.4;   % medium green
    0.0 0.5 0.0];  % dark green

for d = 1:3
    shaded_plot(t, avg_by_day{d}, sem_by_day{d}, cmap1(d,:));
end
xline(0, '--k', 'LineWidth', 1.5);
ylim([-0.6 0.4])

xlabel('Time (s)');
ylabel('DA Signal (z-score)');
title('Days 1–3');
legend({'','Day 1','','Day 2','','Day 3'});xlim([-5 5]);

figure; hold on;

cmap2 = [ ...
    0.7 1.0 0.9;
    0.3 0.7 0.6;
    0.0 0.4 0.3];

for d = 4:6
    shaded_plot(t, avg_by_day{d}, sem_by_day{d}, cmap2(d-3,:));
end

xline(0, '--k', 'LineWidth', 1.5);
ylim([-0.6 0.4])

xlabel('Time (s)');
ylabel('DA Signal (z-score)');
title('Days 4–6');
legend({'','Day 4','','Day 5','','Day 6'});
xlim([-5 5]);


%%  BAR PLOT 

post_idx = (t >= 0 & t <= 1);

% Extract single metric per session (0–2s mean)
metric_by_day = cell(1,nDays);

for d = 1:nDays
    data = data_by_day{d};
    metric_by_day{d} = max(data(:, post_idx), [], 2);
end

%% Compute stats
means = zeros(1,nDays);
sems  = zeros(1,nDays);

for d = 1:nDays
    x = metric_by_day{d};
    means(d) = mean(x);
    sems(d)  = std(x) / sqrt(numel(x));
end

%% 🔹 GREEN BAR PLOT
figure; hold on;

x = 1:6;

b = bar(x, means, 'FaceColor', 'flat');

% green gradient
b.CData = [
    0.85 1.00 0.85;
    0.70 0.95 0.70;
    0.55 0.90 0.55;
    0.40 0.80 0.40;
    0.25 0.65 0.25;
    0.10 0.45 0.10];

% error bars
errorbar(x, means, sems, ...
    'k', 'LineStyle', 'none', ...
    'LineWidth', 1.5, ...
    'CapSize', 10);

% Individual points + connecting lines
all_data = [];

for d = 1:nDays
    all_data = [all_data, metric_by_day{d}];
end

nMice = size(all_data,1);

for i = 1:nMice
    plot(1:6, all_data(i,:), '-o', ...
        'Color', [0.6 0.6 0.6], ...
        'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k');
end

% Labels
xticks(1:6);
xticklabels({'Day 1','Day 2','Day 3','Day 4','Day 5','Day 6'});

ylabel('Mean DA Response (0–1 s)');
title('SWR-DA Pre');

box off;

%
[h_14,p_14,ci_14,stats_14] = ttest(metric_by_day{1}, metric_by_day{4})
[h_36,p_36,ci_36,stats_36] = ttest(metric_by_day{3}, metric_by_day{6})
[h_16,p_16,ci_16,stats_16] = ttest(metric_by_day{1}, metric_by_day{6});

disp(['Day 1 vs 4 p = ', num2str(p_14)]);
disp(['Day 3 vs 6 p = ', num2str(p_36)]);
disp(['Day 1 vs 6 p = ', num2str(p_16)]);


%%
function shaded_plot(x, y, err, color)
    fill([x fliplr(x)], ...
         [y+err fliplr(y-err)], ...
         color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on;
    plot(x, y, 'Color', color, 'LineWidth', 2);
end
