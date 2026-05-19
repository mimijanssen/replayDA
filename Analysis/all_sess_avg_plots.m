%% plot all circ shifted lines and session averages per mice 
clear;
clc;
%% load data
cd 'D:\Mouse_avg'
load('M452_avg.mat')

% also pretty:
%addpath('C:\Users\mimia\Documents\Toolboxes\ScientificColourMaps8\batlow');
%load('batlow.mat');

addpath('C:\Users\mimia\Documents\Toolboxes\ScientificColourMaps8\navia');
load('navia.mat');

%% Pre Track Sleep Data Wrangling 
fields = fieldnames(sess);
numsess = length(fields);
t = time; % time is the same for all sessions

%%
figure(1);

% Add a gray dashed line at 4 seconds
plot([4, 4], [-0.5 0.5], '--k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);

hold on;

% Downsample t
t_down = downsample(t(1,:), 4);

% skip the top 20% ... too white
non_white_range = round(0.15*256):round(0.85*256);  % Skip the middle 20% of the colormap
new_colormap = navia(non_white_range, :);  % Create a new colormap excluding the middle part

% Apply the adjusted colormap
colormap(new_colormap);

% Select the colors for each session from the adjusted colormap
colors = new_colormap(round(linspace(1, length(non_white_range), numsess)), :); 

for i = 1:numsess
    field = fields{i};
    avg_fiber_pre = downsample([sess.(field).avg_fiber_pre], 4);  % Downsample to reduce file size
    % Plot with color corresponding to the session
    lineHandle = plot(t_down, avg_fiber_pre, 'LineWidth', 3, 'Color', colors(i, :));  % Apply color to each line
end

% Set x-ticks and labels centered around 0
xticks(0:1:8);
xticklabels({'-4', '', '-2', '', '0', '', '2', '', '4'});

yticks(-0.5:0.1:0.5);
yticklabels({'-0.5', '', '', '', '', '0', '', '', '', '', '0.5'});

ax = gca; % Get current axes handle
set(ax, 'Box', 'off');  % Remove the box around the plot
set(ax, 'XAxisLocation', 'bottom');  % Keep bottom axis
set(ax, 'YAxisLocation', 'left');    % Keep left axis

% Hide the top and right axis
ax.XColor = 'k'; % Set bottom axis color (black)
ax.YColor = 'k'; % Set left axis color (black)
ax.TickDir = 'out';  % Move the ticks to the outside

xlabel('Time from SWR (s)');
ylabel('Mean Signal (z-score)');
title('M452:Pre-Track Session');
xlim([0 8])
ylim([-0.5 0.5])
set(gcf, 'color', 'none');
set(gca, 'color', 'none');

% Apply font settings to the axes
set(gca, 'FontSize', 18);

% Add a color bar and label it
c = colorbar;
c.Label.String = 'Session Index';  % Customize this label if needed

% Set the color axis range from 1 to numsess
caxis([1 numsess]);

% Ensure vector output and export
set(gcf, 'renderer', 'painters');
exportgraphics(gcf, 'M452_v4_pre.eps');  % Export as PDF
fontname("AvenirNext LT Pro Regular");

hold off;

%%
figure(2);

fontname("AvenirNext LT Pro Regular");

% Add a gray dashed line at 4 seconds
plot([4, 4], [-0.5 0.5], '--k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);

hold on;

% Downsample t
t_down = downsample(t(1,:), 4);

% skip the top 20% ... too white
non_white_range = round(0.15*256):round(0.85*256);  % Skip the middle 20% of the colormap
new_colormap = navia(non_white_range, :);  % Create a new colormap excluding the middle part

% Apply the adjusted colormap
colormap(new_colormap);

% Select the colors for each session from the adjusted colormap
colors = new_colormap(round(linspace(1, length(non_white_range), numsess)), :); 

for i = 1:numsess
    field = fields{i};
    avg_fiber_post = downsample([sess.(field).avg_fiber_post], 4);  % Downsample to reduce file size
    % Plot with color corresponding to the session
    lineHandle = plot(t_down, avg_fiber_post, 'LineWidth', 3, 'Color', colors(i, :));  % Apply color to each line
    hold on
end

% Set x-ticks and labels centered around 0
xticks(0:1:8);
xticklabels({'-4', '', '-2', '', '0', '', '2', '', '4'});

yticks(-0.5:0.1:0.5);
yticklabels({'-0.5', '', '', '', '', '0', '', '', '', '', '0.5'});

ax = gca; % Get current axes handle
set(ax, 'Box', 'off');  % Remove the box around the plot
set(ax, 'XAxisLocation', 'bottom');  % Keep bottom axis
set(ax, 'YAxisLocation', 'left');    % Keep left axis

% Hide the top and right axis
ax.XColor = 'k'; % Set bottom axis color (black)
ax.YColor = 'k'; % Set left axis color (black)
ax.TickDir = 'out';  % Move the ticks to the outside

xlabel('Time from SWR (s)', 'FontName', 'AvenirNext LT Pro Regular');
ylabel('Mean Signal (z-score)', 'FontName', 'AvenirNext LT Pro Regular');
title('M547:Post-Track Session', 'FontName', 'AvenirNext LT Pro Regular');
xlim([0 8])
ylim([-0.5 0.5])
set(gcf, 'color', 'none');
set(gca, 'color', 'none');

% Apply font settings to the axes
set(gca, 'FontName', 'AvenirNext LT Pro', 'FontSize', 18);

% Add a color bar and label it
c = colorbar;
c.Label.String = 'Session Index';  % Customize this label if needed
set(c, 'FontName', 'AvenirNext LT Pro Regular');  % Apply font to color bar

% Set the color axis range from 1 to numsess
caxis([1 numsess]);

% Ensure vector output and export
set(gcf, 'renderer', 'painters');
exportgraphics(gcf, 'M452_v4_post.eps');  % Export as PDF

hold off;

