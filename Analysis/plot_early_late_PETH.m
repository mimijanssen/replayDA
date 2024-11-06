%% Grand Average Over Mice : EARlY AND LATE SESSIONS 
addpath('C:\Users\mimia\Documents\Toolboxes\shadedErrorBar')

%% Mouse Average Plots
clear; clc;
cd 'F:\early_late_avg'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   sess.(['mouse',num2str(k-2)]) = load(FileNames);
end

%% LOAD THIS FOR TIME
cd 'F:\Mouse_avg'
Files=dir('*.*');
k = 6;
FileNames=Files(k).name;
time.(['mouse',num2str(k-2)]) = load(FileNames);


%% PETH AVG 

fiber_pre_early = [sess.mouse1.avg_fiber_pre_early; sess.mouse2.avg_fiber_pre_early; sess.mouse3.avg_fiber_pre_early; sess.mouse4.avg_fiber_pre_early; sess.mouse5.avg_fiber_pre_early;sess.mouse6.avg_fiber_pre_early;sess.mouse7.avg_fiber_pre_early;sess.mouse8.avg_fiber_pre_early]; %sess.sess4.avg_fiber_pre;sess.sess5.avg_fiber_pre;sess.sess6.avg_fiber_pre;sess.sess7.avg_fiber_pre;]; %sess.sess8.avg_fiber_pre
fiber_pre_late = [sess.mouse1.avg_fiber_pre_late; sess.mouse2.avg_fiber_pre_late; sess.mouse3.avg_fiber_pre_late; sess.mouse4.avg_fiber_pre_late; sess.mouse6.avg_fiber_pre_late;sess.mouse7.avg_fiber_pre_late;sess.mouse8.avg_fiber_pre_late]; %sess.sess4.avg_fiber_pre;sess.sess5.avg_fiber_pre;sess.sess6.avg_fiber_pre;sess.sess7.avg_fiber_pre;]; %sess.sess8.avg_fiber_pre
fiber_post_early = [sess.mouse1.avg_fiber_post_early; sess.mouse2.avg_fiber_post_early; sess.mouse3.avg_fiber_post_early; sess.mouse4.avg_fiber_post_early; sess.mouse5.avg_fiber_post_early;sess.mouse6.avg_fiber_post_early;sess.mouse7.avg_fiber_post_early;sess.mouse8.avg_fiber_post_early]; %sess.sess4.avg_fiber_post;sess.sess5.avg_fiber_post;sess.sess6.avg_fiber_post;sess.sess7.avg_fiber_post;]; %sess.sess8.avg_fiber_post
fiber_post_late = [sess.mouse1.avg_fiber_post_late; sess.mouse2.avg_fiber_post_late; sess.mouse3.avg_fiber_post_late; sess.mouse4.avg_fiber_post_late;sess.mouse6.avg_fiber_post_late;sess.mouse7.avg_fiber_post_late;sess.mouse8.avg_fiber_post_late]; %sess.sess4.avg_fiber_post;sess.sess5.avg_fiber_post;sess.sess6.avg_fiber_post;sess.sess7.avg_fiber_post;]; %sess.sess8.avg_fiber_post

avg_fiber_pre_early = mean(fiber_pre_early);
avg_fiber_pre_late = mean(fiber_pre_late);
avg_fiber_post_early = mean(fiber_post_early);
avg_fiber_post_late = mean(fiber_post_late);

sem_fiber_pre_early = std(fiber_pre_early)/sqrt(size(fiber_pre_early,1));
sem_fiber_pre_late = std(fiber_pre_late)/sqrt(size(fiber_pre_late,1));
sem_fiber_post_early = std(fiber_post_early)/sqrt(size(fiber_post_early,1));
sem_fiber_post_late = std(fiber_post_late)/sqrt(size(fiber_post_late,1));

%%
%dark_green = [78,178,101]./255;%[0,104,87]./255; 
%light_green = [144, 201, 135]./255;

%lighterp = [192,175,212]./255; % lilac
lightp =  [141,217,153]./255; %[167,241,121]./255; %[172,149,198]./225; % african violet
% cyan [28,144,153]./225; %
%darkp = [141, 109,176]./255; % amethyst
darkerp = [16,81,96]./255; %[101,71,133]./225; % ultra violet

%% 

figure(1)
plot([4, 4], [-0.5 0.5], '--k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on
shadedErrorBar(time.mouse4.sess.sess1.time(1,:),avg_fiber_pre_early,sem_fiber_pre_early,'lineProps',{'-','color',lightp,'MarkerFaceColor',lightp})
plot(time.mouse4.sess.sess1.time(1,:),avg_fiber_pre_early,'LineWidth',3,'Color',lightp)

shadedErrorBar(time.mouse4.sess.sess1.time(1,:),avg_fiber_pre_late,sem_fiber_pre_late,'lineProps',{'-','color',darkerp,'MarkerFaceColor',darkerp})
plot(time.mouse4.sess.sess1.time(1,:),avg_fiber_pre_late,'LineWidth',3,'Color',darkerp)

%xl = xline(4,'',{'SWR'});
%xl.LabelVerticalAlignment = 'top';
xlim([0 8])
xticks([0 4 8])
ylim([-0.1 0.2])
xticklabels({'-4','0','4'})
title('Pre-Track Rest [DA] after SWRs')
ylabel('Mean [DA] (z-score)')
xlabel('Time from SWR (s)')
legend('','early','','late','Location','northwest')
legend boxoff

set(gca,'fontsize', 18)
set(gcf, 'color', 'none');
set(gca, 'color', 'none');

set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Desktop')
exportgraphics(gcf, 'Pretrack_early_late.png', 'ContentType','vector');  % Export as PDF

hold off
%%
figure(2)
plot([4, 4], [-0.5 0.5], '--k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on
shadedErrorBar(time.mouse4.sess.sess1.time(1,:),avg_fiber_post_early,sem_fiber_post_early,'lineProps',{'-','color',lightp,'MarkerFaceColor',lightp})
plot(time.mouse4.sess.sess1.time(1,:),avg_fiber_post_early,'LineWidth',3,'Color',lightp)

shadedErrorBar(time.mouse4.sess.sess1.time(1,:),avg_fiber_post_late,sem_fiber_post_late,'lineProps',{'-','color',darkerp,'MarkerFaceColor',darkerp})
plot(time.mouse4.sess.sess1.time(1,:),avg_fiber_post_late,'LineWidth',3,'Color',darkerp)

%xl = xline(4,'',{'SWR'});
%xl.LabelVerticalAlignment = 'top';
xlim([0 8])
xticks([0 4 8])
ylim([-0.1 0.2])
xticklabels({'-4','0','4'})
title('Post-Track Rest [DA] after SWRs')
ylabel('Mean [DA] (z-score)')
xlabel('Time from SWR (s)')
legend('','early','','late','Location','northwest')
legend boxoff

set(gca,'fontsize', 18)
set(gcf, 'color', 'none');
set(gca, 'color', 'none');

set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Desktop')
exportgraphics(gcf, 'Posttrack_early_late.png', 'ContentType','vector');  % Export as PDF

hold off

%% Calculate dF 
n = 8; % number of mice 

% for mouse_fiber_post
% for 1:4000 (before swr)
x1 = 2001:1:4001;
% for 4001:8001 (after swr)
x2 = 4001:1:6001;
% area under the curve for each mouse ; row is mouse 

dF1_post_early = zeros(n,1);
dF2_post_early = zeros(n,1);
dF1_post_late = zeros(n,1);
dF2_post_late = zeros(n,1);

dF1_pre_early = zeros(n,1);
dF2_pre_early = zeros(n,1);
dF1_pre_late = zeros(n,1);
dF2_pre_late = zeros(n,1);


% dF for post 
% EARLY ~~~~~~~~~~~~~~~~~~~~~~~~
% before SWR
for i_post = 1:1:n
    dF1_post_early(i_post,:) = max(fiber_post_early(i_post,x1))-min(fiber_post_early(i_post,x1));
end
% after SWR
for i_post = 1:1:n
    dF2_post_early(i_post,:) = max(fiber_post_early(i_post,x2))-min(fiber_post_early(i_post,x2));
end
% LATE ~~~~~~~~~~~~~~~~~~~~~~~~~~
% before SWR
for i_post = 1:1:7
    dF1_post_late(i_post,:) = max(fiber_post_late(i_post,x1))-min(fiber_post_late(i_post,x1));
end
% after SWR
for i_post = 1:1:7
    dF2_post_late(i_post,:) = max(fiber_post_late(i_post,x2))-min(fiber_post_late(i_post,x2));
end

% dF for pre
for i_pre = 1:1:n
    dF1_pre_early(i_pre,:) = max(fiber_pre_early(i_pre,x1))-min(fiber_pre_early(i_pre,x1));
end
% after SWR
for i_pre = 1:1:n
    dF2_pre_early(i_pre,:) = max(fiber_pre_early(i_pre,x2))-min(fiber_pre_early(i_pre,x2));
end
% LATE ~~~~~~~~~~~~~~~~~~~~~~~~~~
% before SWR
for i_pre = 1:1:7
    dF1_pre_late(i_pre,:) = max(fiber_pre_late(i_pre,x1))-min(fiber_pre_late(i_pre,x1));
end
% after SWR
for i_pre = 1:1:7
    dF2_pre_late(i_pre,:) = max(fiber_pre_late(i_pre,x2))-min(fiber_pre_late(i_pre,x2));
end

%% OK PLOT THIS NOW
% colors for main boxplot
% early color 
earlyb = [148,219,157]./255;
% late color 
lateb = [30,148,176]./255;

%%
% Concatenate all data into a single column vector for boxchart input
combinedData = [dF1_pre_early(:); dF2_pre_early(:); dF1_pre_late(:); dF2_pre_late(:)];% Define group categories for each data point
categories = [ones(1, numel(dF1_pre_early)), 2 * ones(1, numel(dF2_pre_early)), ...
              3 * ones(1, numel(dF1_pre_late)), 4 * ones(1, numel(dF2_pre_late))];

% Create the boxchart figure
figure(4);
%h = boxchart(categories, combinedData);  % Use categories for x grouping
hold on;
%h.BoxFaceColor = lighterp;

boxchart(ones(size(dF1_pre_early)), dF1_pre_early, 'BoxFaceColor', earlyb);
boxchart(2 * ones(size(dF2_pre_early)), dF2_pre_early, 'BoxFaceColor', earlyb);
boxchart(3 * ones(size(dF1_pre_late)), dF1_pre_late, 'BoxFaceColor', lateb);
boxchart(4 * ones(size(dF2_pre_late)), dF2_pre_late, 'BoxFaceColor', lateb);


% Set x-axis positions for each group with consistent jitter
jitterAmount = 0.05;
x_jitter_early_before = 1 + randn(1, numel(dF1_pre_early)) * jitterAmount;
x_jitter_early_after = 2 + randn(1, numel(dF2_pre_early)) * jitterAmount;
x_jitter_late_before = 3 + randn(1, numel(dF1_pre_late)) * jitterAmount;
x_jitter_late_after = 4 + randn(1, numel(dF2_pre_late)) * jitterAmount;

% Scatter points for each group
scatter(x_jitter_early_before, dF1_pre_early, 60, 'MarkerFaceColor', lightp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
scatter(x_jitter_early_after, dF2_pre_early, 60, 'MarkerFaceColor', lightp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
scatter(x_jitter_late_before, dF1_pre_late, 60, 'MarkerFaceColor', darkerp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
scatter(x_jitter_late_after, dF2_pre_late, 60, 'MarkerFaceColor', darkerp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

% Draw connecting lines for each subject between before and after within early and late sessions
for i = 1:numel(dF1_pre_early)
    plot([x_jitter_early_before(i), x_jitter_early_after(i)], [dF1_pre_early(i), dF2_pre_early(i)], ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey lines for early
end

for i = 1:numel(dF1_pre_late)
    plot([x_jitter_late_before(i), x_jitter_late_after(i)], [dF1_pre_late(i), dF2_pre_late(i)], ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey lines for late
end

% Set x-tick labels and axis properties
xticks([1, 2, 3, 4]);
xticklabels({"Early Before", "Early After", "Late Before", "Late After"});
ylim([-0.05 0.32]);
ylabel("Mean \Delta [DA] (z-score)");
title("Pre-Track Rest");

% Set figure properties
set(gcf, 'color', 'none');
set(gca, 'color', 'none');
set(gca, 'fontsize', 18);

% Set renderer and export options
set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Desktop')
exportgraphics(gcf, 'dF_pre_early_late_box.png', 'ContentType', 'vector');  % Export as PDF

% Show the figure
hold off;

%% POST 
% Concatenate all data into a single column vector for boxchart input
combinedData = [dF1_post_early(:); dF2_post_early(:); dF1_post_late(:); dF2_post_late(:)];% Define group categories for each data point
categories = [ones(1, numel(dF1_post_early)), 2 * ones(1, numel(dF2_post_early)), ...
              3 * ones(1, numel(dF1_post_late)), 4 * ones(1, numel(dF2_post_late))];

% Create the boxchart figure
figure(5);
%h = boxchart(categories, combinedData);  % Use categories for x grouping
hold on;
%h.BoxFaceColor = lighterp;

boxchart(ones(size(dF1_post_early)), dF1_post_early, 'BoxFaceColor', earlyb);
boxchart(2 * ones(size(dF2_post_early)), dF2_post_early, 'BoxFaceColor', earlyb);
boxchart(3 * ones(size(dF1_post_late)), dF1_post_late, 'BoxFaceColor', lateb);
boxchart(4 * ones(size(dF2_post_late)), dF2_post_late, 'BoxFaceColor', lateb);


% Set x-axis positions for each group with consistent jitter
jitterAmount = 0.05;
x_jitter_early_before = 1 + randn(1, numel(dF1_post_early)) * jitterAmount;
x_jitter_early_after = 2 + randn(1, numel(dF2_post_early)) * jitterAmount;
x_jitter_late_before = 3 + randn(1, numel(dF1_post_late)) * jitterAmount;
x_jitter_late_after = 4 + randn(1, numel(dF2_post_late)) * jitterAmount;

% Scatter points for each group
scatter(x_jitter_early_before, dF1_post_early, 60, 'MarkerFaceColor', lightp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
scatter(x_jitter_early_after, dF2_post_early, 60, 'MarkerFaceColor', lightp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
scatter(x_jitter_late_before, dF1_post_late, 60, 'MarkerFaceColor', darkerp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
scatter(x_jitter_late_after, dF2_post_late, 60, 'MarkerFaceColor', darkerp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

% Draw connecting lines for each subject between before and after within early and late sessions
for i = 1:numel(dF1_post_early)
    plot([x_jitter_early_before(i), x_jitter_early_after(i)], [dF1_post_early(i), dF2_post_early(i)], ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey lines for early
end

for i = 1:numel(dF1_post_late)
    plot([x_jitter_late_before(i), x_jitter_late_after(i)], [dF1_post_late(i), dF2_post_late(i)], ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey lines for late
end

% Set x-tick labels and axis properties
xticks([1, 2, 3, 4]);
xticklabels({"Early Before", "Early After", "Late Before", "Late After"});
ylim([-0.05 0.32]);
ylabel("Mean \Delta [DA] (z-score)");
title("Post-Track Rest");

% Set figure properties
set(gcf, 'color', 'none');
set(gca, 'color', 'none');
set(gca, 'fontsize', 18);

% Set renderer and export options
set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Desktop')
exportgraphics(gcf, 'dF_post_early_late_box.png', 'ContentType', 'vector');  % Export as PDF

% Show the figure
hold off;

%% SAVE VARIABLES NEEDED FOR PLOTTING
% mouse_fiber_pre = sess_fiber_pre;
% mouse_fiber_post = sess_fiber_post; 
% 
% time_swr_da_plot = sess.mouse1.sess.sess1.time(1,:);
% time_rpe_plot = sess.mouse1.RPE.RPE1.t_shared;
% 
%  cd 'D:\Mouse_Avg\'
%  file_name = 'mouse_'; 
%  filename = append(file_name, "avg.mat");
%  save(filename, 'pre_count', 'post_count', 'mean_pre', 'mean_post', 'std_pre', 'std_post', 'mouse_high','mouse_med','mouse_low','mean_low','mean_med','mean_high','std_low','std_med','std_high','time_rpe_plot','avg_fiber_post','avg_fiber_pre','std_fiber_post','std_fiber_pre','time_swr_da_plot','mouse_fiber_post','mouse_fiber_pre');
