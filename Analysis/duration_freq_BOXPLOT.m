%% DURATION AND FREQUENCY PLOTS
cd('D:\Duration')

%%
Files=dir('*.*');

% load files for each mouse -- this is hard coded ew.

% M433
l = 2;
for k=3:9
   FileNames=Files(k).name;
   if k == 6 
       l = l - 1; % skip labeling the fourth session 
   end 
   M433.(['sess',num2str(k-l)]) = load(FileNames);
end

% M453 
l = 0;
for k=10:14
   FileNames=Files(k).name;
   if k == 10 % if the first file... skip to sess 2, if the 5th file, skip to next session (two skips.) 
       l = l + 1; % skip labeling the fourth session 
   elseif k == 14
       l = l + 2; 
   end 
   l = l + 1; 
   M453.(['sess',num2str(l)]) = load(FileNames);
end

% M460
l = 0;
for k=15:18
   FileNames=Files(k).name;
   if k == 16 % if the first file... skip to sess 2, if the 5th file, skip to next session (two skips.) 
       l = l + 3; % skip labeling the fourth session 
   end 
   l = l + 1; 
   M460.(['sess',num2str(l)]) = load(FileNames);
end

% M533
l = 0;
for k=19:25
   FileNames=Files(k).name;
   l = l + 1; 
   M533.(['sess',num2str(l)]) = load(FileNames);
end

% M534
l = 0;
for k=26:28
   FileNames=Files(k).name;
   if k == 27 % if the first file... skip to sess 2, if the 5th file, skip to next session (two skips.) 
       l = l + 1; % skip labeling the fourth session 
   end 
   l = l + 1; 
   M534.(['sess',num2str(l)]) = load(FileNames);
end

% M545
l = 3;
for k=29:32
   FileNames=Files(k).name;
   l = l + 1; 
   M545.(['sess',num2str(l)]) = load(FileNames);
end

% M547
l = 0;
for k=33:38
   FileNames=Files(k).name;
   l = l + 1; 
   M547.(['sess',num2str(l)]) = load(FileNames);
end

% M548
l = 0;
for k=39:45
   FileNames=Files(k).name;
   l = l + 1; 
   M548.(['sess',num2str(l)]) = load(FileNames);
end


%% SESSION PRE AND POST MATRIX 

% session 1 
sess1_pre_freq = [M433.sess1.freq(1), M460.sess1.freq(1), M533.sess1.freq(1), M534.sess1.freq(1), M547.sess1.freq(1), M548.sess1.freq(1)]; 
sess1_post_freq = [M433.sess1.freq(2), M460.sess1.freq(2), M533.sess1.freq(2), M534.sess1.freq(2), M547.sess1.freq(2), M548.sess1.freq(2)];
sess1_pre_dur = [M433.sess1.dur(1), M460.sess1.dur(1), M533.sess1.dur(1), M534.sess1.dur(1), M547.sess1.dur(1), M548.sess1.dur(1)];
sess1_post_dur = [M433.sess1.dur(2), M460.sess1.dur(2), M533.sess1.dur(2), M534.sess1.dur(2), M547.sess1.dur(2), M548.sess1.dur(2)];

% session 2 
sess2_pre_freq = [M433.sess2.freq(1), M453.sess2.freq(1), M533.sess2.freq(1), M547.sess2.freq(1), M548.sess2.freq(1)]; 
sess2_post_freq = [M433.sess2.freq(2), M453.sess2.freq(2), M533.sess2.freq(2), M547.sess2.freq(2), M548.sess2.freq(2)];
sess2_pre_dur = [M433.sess2.dur(1), M453.sess2.dur(1), M533.sess2.dur(1), M547.sess2.dur(1), M548.sess2.dur(1)];
sess2_post_dur = [M433.sess2.dur(2), M453.sess2.dur(2), M533.sess2.dur(2),  M547.sess2.dur(1), M548.sess2.dur(2)];

% session 3
sess3_pre_freq = [M433.sess3.freq(1), M453.sess3.freq(1), M533.sess3.freq(1), M534.sess3.freq(1), M547.sess3.freq(1), M548.sess3.freq(1)]; 
sess3_post_freq = [M433.sess3.freq(2), M453.sess3.freq(2), M533.sess3.freq(2), M534.sess3.freq(2), M547.sess3.freq(2), M548.sess3.freq(2)];
sess3_pre_dur = [M433.sess3.dur(1), M453.sess3.dur(1), M533.sess3.dur(1), M534.sess3.dur(1), M547.sess3.dur(1), M548.sess3.dur(1)];
sess3_post_dur = [M433.sess3.dur(2), M453.sess3.dur(2), M533.sess3.dur(2), M534.sess3.dur(1), M547.sess3.dur(2), M548.sess3.dur(2)];

% session 4 
sess4_pre_freq = [M453.sess4.freq(1), M533.sess4.freq(1), M534.sess4.freq(1), M545.sess4.freq(1),  M547.sess4.freq(1), M548.sess4.freq(1)]; 
sess4_post_freq = [M453.sess4.freq(2), M533.sess4.freq(2), M534.sess4.freq(2), M545.sess4.freq(2), M547.sess4.freq(2), M548.sess4.freq(2)];
sess4_pre_dur = [M453.sess4.dur(1), M533.sess4.dur(1), M534.sess4.dur(1), M545.sess4.dur(1), M547.sess4.dur(1), M548.sess4.dur(1)];
sess4_post_dur = [M453.sess4.dur(2), M533.sess4.dur(2), M534.sess4.dur(2), M545.sess4.dur(2), M547.sess4.dur(2), M548.sess4.dur(2)];

% session 5
sess5_pre_freq = [M433.sess5.freq(1), M453.sess5.freq(1), M460.sess5.freq(1), M533.sess5.freq(1), M545.sess5.freq(1), M547.sess5.freq(1), M548.sess5.freq(1)]; 
sess5_post_freq = [M433.sess5.freq(2), M453.sess5.freq(2), M460.sess5.freq(2), M533.sess5.freq(2), M545.sess5.freq(2),  M547.sess5.freq(2), M548.sess5.freq(2)];
sess5_pre_dur = [M433.sess5.dur(1), M453.sess5.dur(1), M460.sess5.dur(1), M533.sess5.dur(1), M545.sess5.dur(1),  M547.sess5.dur(1), M548.sess5.dur(1)];
sess5_post_dur = [M433.sess5.dur(2), M453.sess5.dur(2), M460.sess5.dur(2), M533.sess5.dur(2), M545.sess5.dur(2),  M547.sess5.dur(2), M548.sess5.dur(2)];

% session 6 
sess6_pre_freq = [M433.sess6.freq(1), M460.sess6.freq(1), M533.sess6.freq(1), M545.sess6.freq(1),  M547.sess6.freq(1), M548.sess6.freq(1)]; 
sess6_post_freq = [M433.sess6.freq(2), M460.sess6.freq(2), M533.sess6.freq(2), M545.sess6.freq(2), M547.sess6.freq(2), M548.sess6.freq(2)];
sess6_pre_dur = [M433.sess6.dur(1), M460.sess6.dur(1), M533.sess6.dur(1), M545.sess6.dur(1), M547.sess6.dur(1), M548.sess6.dur(1)];
sess6_post_dur = [M433.sess6.dur(2), M460.sess6.dur(2), M533.sess6.dur(2), M545.sess6.dur(2), M547.sess6.dur(2), M548.sess6.dur(2)];

% session 7
sess7_pre_freq = [M433.sess7.freq(1), M460.sess7.freq(1), M533.sess7.freq(1), M545.sess7.freq(1), M548.sess7.freq(1)]; 
sess7_post_freq = [M433.sess7.freq(2), M460.sess7.freq(2), M533.sess7.freq(2), M545.sess7.freq(2), M548.sess7.freq(2)];
sess7_pre_dur = [M433.sess7.dur(1), M533.sess7.dur(1), M545.sess7.dur(1), M548.sess7.dur(1)]; % rerun sess 7 pre dur for M460...
% I think I should either exclude M460 because lack of swrs... or just
% don't use the pre session... seems like only noise was being picke dup...M460.sess7.dur(1)
% 
sess7_post_dur = [M433.sess7.dur(2), M460.sess7.dur(2), M533.sess7.dur(2), M545.sess7.dur(2), M548.sess7.dur(2)];

% session 8 
sess8_pre_freq = [M433.sess8.freq(1), M453.sess8.freq(1)]; 
sess8_post_freq = [M433.sess8.freq(2), M453.sess8.freq(2)];
sess8_pre_dur = [M433.sess8.dur(1), M453.sess8.dur(1)];
sess8_post_dur = [M433.sess8.dur(2), M453.sess8.dur(2)];

%%
% Define colors
pre_color = [0.4, 0.6, 0.8]; % Example color for "pre" BoxFaceColor
post_color = [0.2, 0.4, 0.7]; % Example color for "post" BoxFaceColor
pre_scatter_color = [0.6, 0.8, 0.9]; % Example color for "pre" scatter points
post_scatter_color = [0.3, 0.5, 0.8]; % Example color for "post" scatter points

% Data preparation
freq_data = {
    sess1_pre_freq, sess1_post_freq;
    sess2_pre_freq, sess2_post_freq;
    sess3_pre_freq, sess3_post_freq;
    sess4_pre_freq, sess4_post_freq;
    sess5_pre_freq, sess5_post_freq;
    sess6_pre_freq, sess6_post_freq;
    sess7_pre_freq, sess7_post_freq;
    sess8_pre_freq, sess8_post_freq
};

dur_data = {
    sess1_pre_dur, sess1_post_dur;
    sess2_pre_dur, sess2_post_dur;
    sess3_pre_dur, sess3_post_dur;
    sess4_pre_dur, sess4_post_dur;
    sess5_pre_dur, sess5_post_dur;
    sess6_pre_dur, sess6_post_dur;
    sess7_pre_dur, sess7_post_dur;
    sess8_pre_dur, sess8_post_dur
};

% Number of sessions
num_sessions = size(freq_data, 1);

% Frequency Boxchart with scatter
figure;
hold on;
for i = 1:num_sessions
    % Combine pre and post data
    session_pre = freq_data{i, 1};
    session_post = freq_data{i, 2};
    combined_data = [session_pre, session_post];
    group_labels = [repmat(1, 1, numel(session_pre)), repmat(2, 1, numel(session_post))];
    
    % Plot boxcharts
    h1 = boxchart(group_labels(1:numel(session_pre)) + 2*(i-1), session_pre, 'BoxFaceColor', pre_color);
    h2 = boxchart(group_labels(numel(session_pre)+1:end) + 2*(i-1), session_post, 'BoxFaceColor', post_color);
    
    % Scatter points with jitter
    scatter(1 + 2*(i-1) + randn(size(session_pre)) * 0.05, session_pre, ...
        60, 'MarkerFaceColor', pre_scatter_color, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
    scatter(2 + 2*(i-1) + randn(size(session_post)) * 0.05, session_post, ...
        60, 'MarkerFaceColor', post_scatter_color, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
end
xlabel('Sessions');
ylabel('Frequency (SWRs per minute)');
xticks(1.5:2:2*num_sessions);
xticklabels(compose('Session %d', 1:num_sessions));
set(gca,'fontsize', 18)
legend('pre','post')
hold off;

% Duration Boxchart with scatter
figure;
hold on;
for i = 1:num_sessions
    % Combine pre and post data
    session_pre = dur_data{i, 1};
    session_post = dur_data{i, 2};
    combined_data = [session_pre, session_post];
    group_labels = [repmat(1, 1, numel(session_pre)), repmat(2, 1, numel(session_post))];
    
    % Plot boxcharts
    h1 = boxchart(group_labels(1:numel(session_pre)) + 2*(i-1), session_pre, 'BoxFaceColor', pre_color);
    h2 = boxchart(group_labels(numel(session_pre)+1:end) + 2*(i-1), session_post, 'BoxFaceColor', post_color);
    
    % Scatter points with jitter
    scatter(1 + 2*(i-1) + randn(size(session_pre)) * 0.05, session_pre, ...
        60, 'MarkerFaceColor', pre_scatter_color, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
    scatter(2 + 2*(i-1) + randn(size(session_post)) * 0.05, session_post, ...
        60, 'MarkerFaceColor', post_scatter_color, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
end
xlabel('Sessions');
ylabel('Duration (ms)');
xticks(1.5:2:2*num_sessions);
xticklabels(compose('Session %d', 1:num_sessions));
set(gca,'fontsize', 18)
legend('pre','post')
hold off;


%% Plot 
low_c = [82,137,199]./255;%[0,104,87]./255; 
color2 = [123,175,222]./255;

mean_pre = mean(pre_count);
mean_post = mean(post_count);
std_pre = std(pre_count);
std_post = std(post_count);

figure(4)
h = boxchart([pre_count, post_count]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;

% Scatter points with jitter for better visibility
scatter(repmat(xDataNumeric(1), 1, numel(pre_count)) + randn(1, numel(pre_count)) * 0.05, pre_count, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(2), 1, numel(post_count)) + randn(1, numel(post_count)) * 0.05, post_count, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

xticklabels(["Pre-Track" "Post-Track"])
ylabel("# of Detected SWRs")
xlabel("Rest Session")
title("SWR Count")

set(gcf, 'color', 'none');
set(gca, 'color', 'none');
set(gca,'fontsize', 18)

set(gcf, 'renderer', 'painters');
fontname("AvenirNext LT Pro Regular");

exportgraphics(gcf, 'SWR.eps', 'ContentType','vector');  % Export as PDF

hold off 