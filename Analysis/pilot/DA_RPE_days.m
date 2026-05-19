%% Mouse Averages for each session (Track days 1 - 6) . 
clear; clc;
cd 'D:\Reward_Omit_RPE'
Files=dir('*.*');

nDays = 6;
rew_by_day = cell(1, nDays);
nonrew_by_day = cell(1, nDays);

% Loading
for k=3:length(Files)
    fname=Files(k).name;
    
    % extract day
    day_num = regexp(fname, 'day(\d)', 'tokens');
    day_num = str2double(day_num{1}{1});
    
    S = load(fname);
    
    % adjust if nested structure
    if isfield(S, 'avg_left')
        avg_left  = S.avg_left(:)';   % force row
        avg_right = S.avg_right(:)';
        med_side  = S.med_side;
    else
        fn = fieldnames(S);
        data = S.(fn{1});
        avg_left  = data.avg_left(:)';
        avg_right = data.avg_right(:)';
        med_side  = data.med_side;
    end
    
    % assign rewarded vs non-rewarded
    if strcmpi(med_side, 'left')
        rew = avg_left;
        nonrew = avg_right;
    else
        rew = avg_right;
        nonrew = avg_left;
    end
    
    % store (with length check)
    if isempty(rew_by_day{day_num})
        rew_by_day{day_num} = rew;
        nonrew_by_day{day_num} = nonrew;
    else
        if length(rew) ~= size(rew_by_day{day_num},2)
            warning('Skipping %s (length mismatch)', fname);
            continue;
        end
        
        rew_by_day{day_num} = [rew_by_day{day_num}; rew];
        nonrew_by_day{day_num} = [nonrew_by_day{day_num}; nonrew];
    end
end

%%
rew_mean = cell(1,nDays);
rew_std  = cell(1,nDays);

nonrew_mean = cell(1,nDays);
nonrew_std  = cell(1,nDays);

for d = [1 3 4 6]
    rew_mean{d} = mean(rew_by_day{d},1);
   % rew_std{d}  = std(rew_by_day{d},0,1);

    n = size(rew_by_day{d},1);
    rew_sem{d} = std(rew_by_day{d},0,1) / sqrt(n);

    
    nonrew_mean{d} = mean(nonrew_by_day{d},1);
   % nonrew_std{d}  = std(nonrew_by_day{d},0,1);
    n2 = size(nonrew_by_day{d},1);
    nonrew_sem{d} = std(nonrew_by_day{d},0,1) / sqrt(n2);
end

%% time points
t = linspace(-5, 5, length(rew_mean{1}));
%%
% rewarded (orange gradient)
rew_colors_13 = [ ...
    1.0 0.8 0.4;
    0.9 0.5 0.1];

rew_colors_46 = [ ...
    1.0 0.7 0.3;
    0.8 0.3 0.0];

% non-rewarded (purple gradient)
nonrew_colors_13 = [ ...
    0.8 0.6 1.0;
    0.5 0.2 0.8];

nonrew_colors_46 = [ ...
    0.7 0.5 0.9;
    0.4 0.0 0.6];

%%

figure; hold on;

days = [1 3];

for i = 1:2
    d = days(i);
    
    shaded_plot(t, rew_mean{d}, rew_sem{d}, rew_colors_13(i,:));
    shaded_plot(t, nonrew_mean{d}, nonrew_sem{d}, nonrew_colors_13(i,:));
end

xline(0,'--k','LineWidth',1.5);

xlabel('Time (s)');
ylabel('DA Signal (z-score)');
title('Days 1 & 3: Reward vs Omission');

legend({'','R D1','','O D1','','R D3','','O D3'});
xlim([-5 5]);

figure; hold on;

days = [4 6];

for i = 1:2
    d = days(i);
    
    shaded_plot(t, rew_mean{d}, rew_sem{d}, rew_colors_46(i,:));
    shaded_plot(t, nonrew_mean{d}, nonrew_sem{d}, nonrew_colors_46(i,:));
end

xline(0,'--k','LineWidth',1.5);

xlabel('Time (s)');
ylabel('DA Signal (z-score)');
title('Days 4 & 6: Reward vs Omission');

legend({'','R D4','','O D4','','R D6','','O D6'});
xlim([-5 5]);
%% BAR PLOTS AND STATS 
%post_idx = t >= 0;
idx = (t >= 0) & (t <= 1);

valid_days = [1 3 4 6];

%% 🔹 Extract peak per mouse
peak_rew = struct();

for d = valid_days
    data = nonrew_by_day{d}; %rew_by_day{d};
    
    % OPTIONAL: smooth (recommended)
    data_sm = smoothdata(data, 2, 'gaussian', 50);
    
    peaks = max(data_sm(:, idx), [], 2);
    %peaks = mean(data_sm(:, t>=0 & t<=2), 2);
    peak_rew.(sprintf('day%d', d)) = peaks;
end


%% 🔹 Collect data
p1 = peak_rew.day1;
p3 = peak_rew.day3;
p4 = peak_rew.day4;
p6 = peak_rew.day6;

means = [mean(p1), mean(p3), mean(p4), mean(p6)];

%%
sem = [ ...
    std(p1)/sqrt(numel(p1)), ...
    std(p3)/sqrt(numel(p3)), ...
    std(p4)/sqrt(numel(p4)), ...
    std(p6)/sqrt(numel(p6)) ];



%%
figure; hold on;

x = 1:4;

% bar plot
b = bar(x, means, 'FaceColor', 'flat');

b.CData = [
    1.0 0.8 0.4;
    1.0 0.6 0.2;
    0.9 0.4 0.1;
    0.7 0.2 0.0];

% SEM error bars
sem = [ ...
    std(p1)/sqrt(numel(p1)), ...
    std(p3)/sqrt(numel(p3)), ...
    std(p4)/sqrt(numel(p4)), ...
    std(p6)/sqrt(numel(p6)) ];

errorbar(x, means, sem, ...
    'k', 'LineStyle', 'none', ...
    'LineWidth', 1.5, ...
    'CapSize', 10);

% paired lines (mice trajectories)
all_peaks = [p1(:), p3(:), p4(:), p6(:)];
nMice = size(all_peaks, 1);

for i = 1:nMice
    plot(1:4, all_peaks(i,:), '-o', ...
        'Color', [0.6 0.6 0.6], ...
        'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k');
end


ymax = max([p1; p3; p4; p6]);

[p_43,p_43,ci_43,stats_43] = ttest(p4, p3);
[p_46,p_46,ci_46,stats_46] = ttest(p4, p6);

if p_43 < 0.05
    plot([2 3], [ymax ymax]*1.1, '-k');
    text(2.5, ymax*1.15, '*', 'HorizontalAlignment','center');
end

if p_46 < 0.05
    plot([3 4], [ymax ymax]*1.2, '-k');
    text(3.5, ymax*1.25, '*', 'HorizontalAlignment','center');
end

xticks(x);
xticklabels({'Day 1','Day 3','Day 4','Day 6'});

ylabel('DA Peak z-score (0–2s)');
title('Rewarded Peak Responses');
box off;

%% 🔹 Statistics

% Choose test type:
paired = true; % <-- set TRUE if same mice tracked across days

if paired
    [p_43,p_43,ci_43,stats_43] = ttest(p4, p3);
    [p_46,p_46,ci_46,stats_46] = ttest(p4, p6);
else
    [h_43,p_43] = ttest2(p4, p3);
    [p_46,~] = ttest2(p4, p6);
end

disp(['Day 4 vs Day 3 p = ', num2str(p_43)]);
disp(['Day 4 vs Day 6 p = ', num2str(p_46)]);

%% 🔹 Add significance stars
ymax = max([p1; p3; p4; p6]);

if p_43 < 0.05
    plot([2 3], [ymax ymax]*1.1, '-k');
    text(2.5, ymax*1.15, '*', 'HorizontalAlignment','center');
end

if p_46 < 0.05
    plot([3 4], [ymax ymax]*1.2, '-k');
    text(3.5, ymax*1.25, '*', 'HorizontalAlignment','center');
end

%% GRAND PETH AND STATS

% 🔹 AVERAGE WITHIN MOUSE → THEN ACROSS MICE

% --- identify mice ---
mouse_ids = {};
rew_mouse = containers.Map();
nonrew_mouse = containers.Map();

for k = 3:length(Files)
    fname = Files(k).name;
    
    % extract mouse ID (e.g., M646)
    mouse_id = regexp(fname, 'M\d+', 'match');
    mouse_id = mouse_id{1};
    
    % extract day
    day_num = regexp(fname, 'day(\d)', 'tokens');
    day_num = str2double(day_num{1}{1});
    
    if ~ismember(day_num, 1:6)
        continue;
    end
    
    S = load(fname);
    
    % handle structure
    if isfield(S, 'avg_left')
        avg_left  = S.avg_left(:)';
        avg_right = S.avg_right(:)';
        med_side  = S.med_side;
    else
        fn = fieldnames(S);
        data = S.(fn{1});
        avg_left  = data.avg_left(:)';
        avg_right = data.avg_right(:)';
        med_side  = data.med_side;
    end
    
    % assign rewarded / non-rewarded
    if strcmpi(med_side, 'left')
        rew = avg_left;
        nonrew = avg_right;
    else
        rew = avg_right;
        nonrew = avg_left;
    end
    
    % store per mouse
    if ~isKey(rew_mouse, mouse_id)
        rew_mouse(mouse_id) = rew;
        nonrew_mouse(mouse_id) = nonrew;
        mouse_ids{end+1} = mouse_id;
    else
        rew_mouse(mouse_id) = [rew_mouse(mouse_id); rew];
        nonrew_mouse(mouse_id) = [nonrew_mouse(mouse_id); nonrew];
    end
end

% 🔹 Average within each mouse (across days)
rew_avg_mouse = [];
nonrew_avg_mouse = [];

for i = 1:length(mouse_ids)
    m = mouse_ids{i};
    
    rew_data = rew_mouse(m);
    nonrew_data = nonrew_mouse(m);
    
    rew_avg_mouse(i,:) = mean(rew_data, 1);
    nonrew_avg_mouse(i,:) = mean(nonrew_data, 1);
end

%% 🔹 Average across mice + SEM
rew_mean = mean(rew_avg_mouse, 1);
rew_sem  = std(rew_avg_mouse, 0, 1) / sqrt(size(rew_avg_mouse,1));

nonrew_mean = mean(nonrew_avg_mouse, 1);
nonrew_sem  = std(nonrew_avg_mouse, 0, 1) / sqrt(size(nonrew_avg_mouse,1));

%% 🔹 Time vector
t = linspace(-5, 5, size(rew_mean,2));

%% 🔹 FINAL PETH PLOT (2 lines)
figure; hold on;

orange = [1.0 0.5 0.0];
purple = [0.5 0.0 0.7];

shaded_plot(t, rew_mean, rew_sem, orange);
shaded_plot(t, nonrew_mean, nonrew_sem, purple);

xline(0,'--k','LineWidth',1.5);

xlabel('Time (s)');
ylabel('DA Signal (z-score)');
title('Reward vs Omission');

xlim([-5 5]);
%
p_vals = zeros(1,length(t));

for i = 1:length(t)
    [~,p_vals(i)] = ttest(rew_avg_mouse(:,i), nonrew_avg_mouse(:,i));
end

sig = p_vals < 0.05;

y = min([rew_mean nonrew_mean]);

plot(t(sig), y*ones(sum(sig),1), '.', 'Color','k');
legend({'','Reward','','Omission'});

%% Bar plot stats 
idx = (t >= 0 & t <= 1);

%rew_sm  = smoothdata(rew_avg_mouse, 2, 'gaussian', 50);
%nonrew_sm = smoothdata(nonrew_avg_mouse, 2, 'gaussian', 50);

rew_metric  = max(rew_avg_mouse(:, idx), [], 2);
nonrew_metric = max(nonrew_avg_mouse(:, idx), [], 2);

[h,p,ci,stats] = ttest(rew_metric, nonrew_metric)

disp(['Reward vs Non-reward p = ', num2str(p)]);

figure; hold on;

means = [mean(rew_metric), mean(nonrew_metric)];
sems  = [std(rew_metric)/sqrt(numel(rew_metric)), ...
         std(nonrew_metric)/sqrt(numel(nonrew_metric))];

x = [1 2];

b = bar(x, means, 'FaceColor','flat');

b.CData = [
    1.0 0.5 0.0;   % orange (rewarded)
    0.5 0.0 0.7];  % purple (non-rewarded)

errorbar(x, means, sems, 'k', ...
    'LineStyle','none','LineWidth',1.5,'CapSize',10);

% paired lines (each mouse)
for i = 1:length(rew_metric)
    plot(x, [rew_metric(i), nonrew_metric(i)], '-o', ...
        'Color',[0.6 0.6 0.6], ...
        'MarkerFaceColor','k');
end

xticks(x);
xticklabels({'Rewarded','Omission'});

ylabel('Mean z-score DA (0–1 s)');
title('Reward vs Omission Arm');

box off;

ymax = max([rew_metric; nonrew_metric]);

if p < 0.05
    plot([1 2], [ymax ymax]*1.1, '-k');
    text(1.5, ymax*1.15, '*', 'HorizontalAlignment','center');
end

%%
function shaded_plot(x, y, err, color)
    fill([x fliplr(x)], ...
         [y+err fliplr(y-err)], ...
         color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on;
    plot(x, y, 'Color', color, 'LineWidth', 2);
end