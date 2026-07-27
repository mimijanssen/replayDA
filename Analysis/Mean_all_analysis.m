% AUC Analysis - All Data 
%% ttest 
clear ; clc
%% 
cd F:\Mouse_avg
load('mouse_avg.mat')
addpath('C:\Users\mimia\Documents\Toolboxes\shadedErrorBar')

%%
n = 8; % number of mice 

% for mouse_fiber_post
% for 1:4000 (before swr)
% 1.5 seconds 
x1 = 2001:1:4000;%1:1:4000;%2501:1:4001;
% for 4001:8001 (after swr)
x2 = 4001:6000;%4001:1:8000;%4501:1:6001;
% area under the curve for each mouse ; row is mouse 
A1_post = zeros(n,1);
for i_post = 1:1:n
    A1_post(i_post,:) = trapz(x1, mouse_fiber_post(i_post,x1));
end

A2_post = zeros(n,1);
for i_post = 1:1:n
    A2_post(i_post,:) = trapz(x2, mouse_fiber_post(i_post,x2));
end

%% AUC for Pre Track Rest
% area under the curve for each mouse ; row is mouse 
A1_pre = zeros(n,1);
for i_pre = 1:1:n
    A1_pre(i_pre,:) = trapz(x1, mouse_fiber_pre(i_pre,x1));
end

A2_pre = zeros(n,1);
for i_pre = 1:1:n
    A2_pre(i_pre,:) = trapz(x2, mouse_fiber_pre(i_pre,x2));
end

%% graphing 
% color 
low_c = [78,178,101]./255;%[0,104,87]./255; 
color2 = [144, 201, 135]./255;

% Calculate means and standard deviations
avg_A1_post = mean(A1_post);
avg_A2_post = mean(A2_post);

avg_A1_pre = mean(A1_pre);
avg_A2_pre = mean(A2_pre);

% Create the boxchart figure
figure(5);
h = boxchart([A1_post, A2_post]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;
% Convert XData to numeric for scatter plotting
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric

% Generate jittered x-coordinates for each condition
jitterAmount = 0.05;
x_jitter_A1 = xDataNumeric(1) + randn(1, numel(A1_post)) * jitterAmount;
x_jitter_A2 = xDataNumeric(2) + randn(1, numel(A2_post)) * jitterAmount;

% Draw light grey lines connecting each subject's data point in "Before" and "After" using the same jitter
for i = 1:numel(A1_post)
    plot([x_jitter_A1(i), x_jitter_A2(i)], [A1_post(i), A2_post(i)], ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey color, thinner line
end

% Scatter points with jitter for each condition
scatter(x_jitter_A1, A1_post, 60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
scatter(x_jitter_A2, A2_post, 60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

% Set x-tick labels and axis properties
xticklabels({"Before", "After"});
ylim([-100 400]);
ylabel("Mean [DA] AUC (z-score)");
title("Post-Track Rest AUC");

% Set figure background color to white
set(gcf, 'color', 'none');
set(gca, 'color', 'none');
set(gca,'fontsize', 18)

fontname("AvenirNext LT Pro Regular");

set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Desktop')

exportgraphics(gcf, 'AUC_post.eps', 'ContentType','vector');  % Export as PDF
hold off;

% Create the boxchart figure
figure(4);
h = boxchart([A1_pre, A2_pre]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;
% Add error bars
%errorbar(1:2, [avg_A1_post, avg_A2_post], [std_pre, std_post], 'LineStyle', 'none', 'Color', color2, 'LineWidth', 2);

% Convert XData to numeric for scatter plotting
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric

% Generate jittered x-coordinates for each condition
jitterAmount = 0.05;
x_jitter_A1 = xDataNumeric(1) + randn(1, numel(A1_pre)) * jitterAmount;
x_jitter_A2 = xDataNumeric(2) + randn(1, numel(A2_pre)) * jitterAmount;

% Draw light grey lines connecting each subject's data point in "Before" and "After" using the same jitter
for i = 1:numel(A1_pre)
    plot([x_jitter_A1(i), x_jitter_A2(i)], [A1_pre(i), A2_pre(i)], ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey color, thinner line
end

% Scatter points with jitter for each condition
scatter(x_jitter_A1, A1_pre, 60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
scatter(x_jitter_A2, A2_pre, 60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

% Set x-tick labels and axis properties
xticklabels({"Before", "After"});
ylim([-100 400]);
ylabel("Mean [DA] AUC (z-score)");
title("Pre-Track Rest AUC");

% Set figure background color to white
set(gcf, 'color', 'none');
set(gca, 'color', 'none');
set(gca,'fontsize', 18)

fontname("AvenirNext LT Pro Regular");

set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Desktop')

exportgraphics(gcf, 'AUC_pre.eps', 'ContentType','vector');  % Export as PDF

% Show the figure
hold off;

%% Paired t-test

[h_post_t,p_post_t, ci_post, stats_post] = ttest(A1_post,A2_post);
% p = 0.0142 
[h_pre_t,p_pre_t, ci_pre, stats_pre] = ttest(A1_pre,A2_pre);
% p = 0.2772 