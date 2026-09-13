%% PETH AVG 

fiber_nrem = [sess.mouse1.avg_fiber_nrem; sess.mouse2.avg_fiber_nrem; sess.mouse3.avg_fiber_nrem; sess.mouse4.avg_fiber_nrem; sess.mouse5.avg_fiber_nrem;sess.mouse6.avg_fiber_nrem;sess.mouse7.avg_fiber_nrem;sess.mouse8.avg_fiber_nrem]; %sess.sess4.avg_fiber_nrem;sess.sess5.avg_fiber_nrem;sess.sess6.avg_fiber_nrem;sess.sess7.avg_fiber_nrem;]; %sess.sess8.avg_fiber_nrem
fiber_wake = [sess.mouse1.avg_fiber_wake; sess.mouse2.avg_fiber_wake; sess.mouse3.avg_fiber_wake; sess.mouse4.avg_fiber_wake; sess.mouse5.avg_fiber_wake;sess.mouse6.avg_fiber_wake;sess.mouse7.avg_fiber_wake;sess.mouse8.avg_fiber_wake]; %sess.sess4.avg_fiber_wake;sess.sess5.avg_fiber_wake;sess.sess6.avg_fiber_wake;sess.sess7.avg_fiber_wake;]; %sess.sess8.avg_fiber_wake

avg_fiber_nrem = mean(fiber_nrem);
avg_fiber_wake = mean(fiber_wake);

sem_fiber_nrem = std(fiber_nrem)/sqrt(size(fiber_nrem,1));
sem_fiber_wake = std(fiber_wake)/sqrt(size(fiber_wake,1));

%%
%dark_green = [78,178,101]./255;%[0,104,87]./255; 
%light_green = [144, 201, 135]./255;

%lighterp = [192,175,212]./255; % lilac
lightp =  [0,165,255]./255; %[163,207,95]./255; %[167,241,121]./255; %[172,149,198]./225; % african violet
% cyan [28,144,153]./225; %
%darkp = [141, 109,176]./255; % amethyst
darkerp = [0,78,100]./255 ; %[16,81,96]./255; %[101,71,133]./225; % ultra violet or 34 85 85

%% 

figure(1)
plot([4, 4], [-0.5 0.5], '--k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on
shadedErrorBar(time,avg_fiber_nrem,sem_fiber_nrem,'lineProps',{'-','color',darkerp,'MarkerFaceColor',darkerp})
plot(time,avg_fiber_nrem,'LineWidth',3,'Color',darkerp)

shadedErrorBar(time(1,:),avg_fiber_nrem,sem_fiber_nrem,'lineProps',{'-','color',lightp,'MarkerFaceColor',lightp})
plot(time,avg_fiber_nrem,'LineWidth',3,'Color',lightp)

%xl = xline(4,'',{'SWR'});
%xl.LabelVerticalAlignment = 'top';
xlim([0 8])
xticks([0 1 2 3 4 5 6 7 8])
ylim([-0.1 0.2])
xticklabels({'-4','','','','0','','','','4'})
title('NREM Rest [DA] after SWRs')
ylabel('Mean [DA] (z-score)')
xlabel('Time from SWR (s)')
legend('','','','','Location','northwest')
legend boxoff

set(gca,'fontsize', 18)
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');

set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Documents\ReplayDA Figures')
exportgraphics(gcf, 'nrem_ticks.png', 'ContentType','vector');  % Export as PDF

hold off
%%
figure(2)
plot([4, 4], [-0.5 0.5], '--k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on

shadedErrorBar(time,avg_fiber_wake,sem_fiber_wake,'lineProps',{'-','color',darkerp,'MarkerFaceColor',darkerp})
plot(time,avg_fiber_wake,'LineWidth',3,'Color',darkerp)

shadedErrorBar(time,avg_fiber_wake,sem_fiber_wake,'lineProps',{'-','color',lightp,'MarkerFaceColor',lightp})
plot(time,avg_fiber_wake,'LineWidth',3,'Color',lightp)

%xl = xline(4,'',{'SWR'});
%xl.LabelVerticalAlignment = 'top';
xlim([0 8])
xticks([0 1 2 3 4 5 6 7 8])
ylim([-0.1 0.2])
xticklabels({'-4','','','','0','','','','4'})
title('Quiet Wakefullness [DA] after SWRs')
ylabel('Mean [DA] (z-score)')
xlabel('Time from SWR (s)')
legend('','','','','Location','northwest')
legend boxoff

set(gca,'fontsize', 18)
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');

set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Documents\ReplayDA Figures')
exportgraphics(gcf, 'wake_ticks.eps', 'ContentType','vector');  % Export as PDF

hold off

%% Calculate dF 
n = 8; % number of mice 

% for mouse_fiber_wake
% for 1:4000 (before swr)
x1 = 2001:1:4001;
% for 4001:8001 (after swr)
x2 = 4001:1:6001;
% area under the curve for each mouse ; row is mouse 

dF1_wake = zeros(n,1);
dF2_wake = zeros(n,1);
%dF1_wake = zeros(7,1);
%dF2_wake = zeros(7,1);

dF1_nrem = zeros(n,1);
dF2_nrem = zeros(n,1);
%dF1_nrem = zeros(7,1);
%dF2_nrem = zeros(7,1);


% dF for wake 
% EARLY ~~~~~~~~~~~~~~~~~~~~~~~~
% before SWR
for i_wake = 1:1:n
    dF1_wake(i_wake,:) = max(fiber_wake(i_wake,x1))-min(fiber_wake(i_wake,x1));
end
% after SWR
for i_wake = 1:1:n
    dF2_wake(i_wake,:) = max(fiber_wake(i_wake,x2))-min(fiber_wake(i_wake,x2));
end

% dF for nrem
for i_nrem = 1:1:n
    dF1_nrem(i_nrem,:) = max(fiber_nrem(i_nrem,x1))-min(fiber_nrem(i_nrem,x1));
end
% after SWR
for i_nrem = 1:1:n
    dF2_nrem(i_nrem,:) = max(fiber_nrem(i_nrem,x2))-min(fiber_nrem(i_nrem,x2));
end

%% NREM AND WAKE BOXPLOT AND TTEST: 

% Concatenate all data into a single column vector for boxchart input
combinedData = [dF1_nrem(:); dF2_nrem(:); dF1_nrem(:); dF2_nrem(:)];% Define group categories for each data point
categories = [ones(1, numel(dF1_nrem)), 2 * ones(1, numel(dF2_nrem)), ...
              3 * ones(1, numel(dF1_nrem)), 4 * ones(1, numel(dF2_nrem))];

% Create the boxchart figure
figure(4);
%h = boxchart(categories, combinedData);  % Use categories for x grouping
hold on;
%h.BoxFaceColor = lighterp;

boxchart(ones(size(dF1_nrem)), dF1_nrem, 'BoxFaceColor', earlyb);
boxchart(2 * ones(size(dF2_nrem)), dF2_nrem, 'BoxFaceColor', earlyb);
boxchart(3 * ones(size(dF1_nrem)), dF1_nrem, 'BoxFaceColor', lateb);
boxchart(4 * ones(size(dF2_nrem)), dF2_nrem, 'BoxFaceColor', lateb);


% Set x-axis positions for each group with consistent jitter
jitterAmount = 0.05;
x_jitter_before = 1 + randn(1, numel(dF1_nrem)) * jitterAmount;
x_jitter_after = 2 + randn(1, numel(dF2_nrem)) * jitterAmount;
x_jitter_before = 3 + randn(1, numel(dF1_nrem)) * jitterAmount;
x_jitter_after = 4 + randn(1, numel(dF2_nrem)) * jitterAmount;

% Draw connecting lines for each subject between before and after within early and late sessions
for i = 1:numel(dF1_nrem)
    plot([x_jitter_before(i), x_jitter_after(i)], [dF1_nrem(i), dF2_nrem(i)], ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey lines for early
end

for i = 1:numel(dF1_nrem)
    plot([x_jitter_before(i), x_jitter_after(i)], [dF1_nrem(i), dF2_nrem(i)], ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey lines for late
end

% Scatter points for each group
scatter(x_jitter_before, dF1_nrem, 60, 'MarkerFaceColor', lightp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);
scatter(x_jitter_after, dF2_nrem, 60, 'MarkerFaceColor', lightp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);
scatter(x_jitter_before, dF1_nrem, 60, 'MarkerFaceColor', darkerp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);
scatter(x_jitter_after, dF2_nrem, 60, 'MarkerFaceColor', darkerp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);

% Set x-tick labels and axis properties
xticks([1, 2, 3, 4]);
xticklabels({"Early Before", "Early After", "Late Before", "Late After"});
ylim([-0.05 0.32]);
ylabel("Mean \Delta [DA] (z-score)");
title("nrem-Track Rest");

% Set figure properties
set(gcf, 'color', 'none');
set(gca, 'color', 'none');
set(gca, 'fontsize', 18);

% Set renderer and export options
set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Desktop')
exportgraphics(gcf, 'dF_nrem_box.eps', 'ContentType', 'vector');  % Export as PDF

% Show the figure
hold off;

%% wake 
% Concatenate all data into a single column vector for boxchart input
combinedData = [dF1_wake(:); dF2_wake(:); dF1_wake(:); dF2_wake(:)];% Define group categories for each data point
categories = [ones(1, numel(dF1_wake)), 2 * ones(1, numel(dF2_wake)), ...
              3 * ones(1, numel(dF1_wake)), 4 * ones(1, numel(dF2_wake))];

% Create the boxchart figure
figure(5);
%h = boxchart(categories, combinedData);  % Use categories for x grouping
hold on;
%h.BoxFaceColor = lighterp;

boxchart(ones(size(dF1_wake)), dF1_wake, 'BoxFaceColor', earlyb);
boxchart(2 * ones(size(dF2_wake)), dF2_wake, 'BoxFaceColor', earlyb);
boxchart(3 * ones(size(dF1_wake)), dF1_wake, 'BoxFaceColor', lateb);
boxchart(4 * ones(size(dF2_wake)), dF2_wake, 'BoxFaceColor', lateb);


% Set x-axis positions for each group with consistent jitter
jitterAmount = 0.05;
x_jitter_before = 1 + randn(1, numel(dF1_wake)) * jitterAmount;
x_jitter_after = 2 + randn(1, numel(dF2_wake)) * jitterAmount;
x_jitter_before = 3 + randn(1, numel(dF1_wake)) * jitterAmount;
x_jitter_after = 4 + randn(1, numel(dF2_wake)) * jitterAmount;

% Draw connecting lines for each subject between before and after within early and late sessions
for i = 1:numel(dF1_wake)
    plot([x_jitter_before(i), x_jitter_after(i)], [dF1_wake(i), dF2_wake(i)], ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey lines for early
end

for i = 1:numel(dF1_wake)
    plot([x_jitter_before(i), x_jitter_after(i)], [dF1_wake(i), dF2_wake(i)], ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey lines for late
end

% Scatter points for each group
scatter(x_jitter_before, dF1_wake, 60, 'MarkerFaceColor', lightp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);
scatter(x_jitter_after, dF2_wake, 60, 'MarkerFaceColor', lightp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);
scatter(x_jitter_before, dF1_wake, 60, 'MarkerFaceColor', darkerp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);
scatter(x_jitter_after, dF2_wake, 60, 'MarkerFaceColor', darkerp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);

% Set x-tick labels and axis properties
xticks([1, 2, 3, 4]);
xticklabels({"Early Before", "Early After", "Late Before", "Late After"});
ylim([-0.05 0.32]);
ylabel("Mean \Delta [DA] (z-score)");
title("wake-Track Rest");

% Set figure properties
set(gcf, 'color', 'none');
set(gca, 'color', 'none');
set(gca, 'fontsize', 18);

% Set renderer and export options
set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Desktop')
exportgraphics(gcf, 'dF_wake_box.eps', 'ContentType', 'vector');  % Export as PDF

% Show the figure
hold off;

%% T-Tests
% wake 
% early sessions
[h_wake_t,p_wake_t, ci_wake, stats_wake] = ttest(dF1_wake,dF2_wake);
% late
[h_wake_t,p_wake_t, ci_wake, stats_wake] = ttest(dF1_wake,dF2_wake);

% nrem
% early
[h_nrem_t,p_nrem_t, ci_nrem, stats_nrem] = ttest(dF1_nrem,dF2_nrem);
% late
[h_nrem_t,p_nrem_t, ci_nrem, stats_nrem] = ttest(dF1_nrem,dF2_nrem);


%% nrem Anova
% Combine response data into a single vector
response = [dF1_nrem(:); dF2_nrem(:); dF1_nrem(:); dF2_nrem(:);];

% Create a categorical variable for training (early vs. late sessions)
timing = [repmat("early", numel(dF1_nrem), 1); ...
          repmat("early", numel(dF2_nrem), 1); ...
          repmat("late", numel(dF1_nrem), 1); ...
          repmat("late", numel(dF2_nrem), 1)];

% Create a categorical variable for condition (before vs. after swr)
condition = [repmat("before", numel(dF1_nrem), 1); ...
             repmat("after", numel(dF2_nrem), 1); ...
             repmat("before", numel(dF1_nrem), 1); ...
             repmat("after", numel(dF2_nrem), 1)];

% Run the two-way ANOVA
[p, tbl, stats] = anovan(response, {timing, condition}, ...
                         'model', 'interaction', ...
                         'varnames', {'Training', 'Condition'});

% Display the results
disp(tbl);

%% wake Anova
response = [dF1_wake(:); dF2_wake(:); dF1_wake(:); dF2_wake(:);];

% Create a categorical variable for training (early vs. late sessions)
timing = [repmat("early", numel(dF1_wake), 1); ...
          repmat("early", numel(dF2_wake), 1); ...
          repmat("late", numel(dF1_wake), 1); ...
          repmat("late", numel(dF2_wake), 1)];

% Create a categorical variable for condition (before vs. after swr)
condition = [repmat("before", numel(dF1_wake), 1); ...
             repmat("after", numel(dF2_wake), 1); ...
             repmat("before", numel(dF1_wake), 1); ...
             repmat("after", numel(dF2_wake), 1)];

% Run the two-way ANOVA
[p_wake, tbl_wake, stats_wake] = anovan(response, {timing, condition}, ...
                         'model', 'interaction', ...
                         'varnames', {'Training', 'Condition'});

% Display the results
disp(tbl_wake);

%% Make sure it is two way: 
response = [dF1_wake(:); dF2_wake(:); dF1_wake(:); dF2_wake(:);];

% Create a categorical variable for training (early vs. late sessions)
timing = [repmat("early", numel(dF1_wake), 1); ...
          repmat("early", numel(dF2_wake), 1); ...
          repmat("late", numel(dF1_wake), 1); ...
          repmat("late", numel(dF2_wake), 1)];

% Create a categorical variable for condition (before vs. after swr)
condition = [repmat("before", numel(dF1_wake), 1); ...
             repmat("after", numel(dF2_wake), 1); ...
             repmat("before", numel(dF1_wake), 1); ...
             repmat("after", numel(dF2_wake), 1)];

% Run the two-way ANOVA
[p_wake, tbl_wake, stats_wake] = anovan(response, {timing, condition}, ...
                         'model', 2, ...
                         'varnames', {'Training', 'Condition'});

% Display the results
disp(tbl_wake);