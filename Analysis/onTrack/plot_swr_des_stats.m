
cd('D:\trackcounts')

%%
Files=dir('*.*');

% load files for each mouse -- this is hard coded ew.

% M646
l = 2;
for k=3:10
   FileNames=Files(k).name
   M646.(['sess',num2str(k-l)]) = load(FileNames);
end
% check
M646

% M648
l = 0;
for k=11:18
   FileNames=Files(k).name
   l = l + 1; % skip labeling the fourth session 
   M648.(['sess',num2str(l)]) = load(FileNames) ;
end
M648

% M650
l = 0;
for k=19:26 %10:14
   FileNames=Files(k).name
   l = l + 1; % skip labeling the fourth session 
   M650.(['sess',num2str(l)]) = load(FileNames);
end
M650

% M654
l = 0;
for k=27:34 %15:18
   FileNames=Files(k).name
   l = l + 1; 
   M654.(['sess',num2str(l)]) = load(FileNames);
end
M654

%%
% Define colors
pre_color = [212, 162, 43]./255; %[0.4, 0.6, 0.8]; % Example color for "pre" BoxFaceColor
%post_color = [0.2, 0.4, 0.7]; % Example color for "post" BoxFaceColor
post_color = [44,123,116]./255; %[0.3, 0.5, 0.8];%[24,53,58]./255; %[0.5, 0.6, 0.9];
pre_scatter_color = [225, 190,106]./255; %[0.6, 0.8, 0.9]; % Example color for "pre" scatter points
post_scatter_color =  [64,176,166]./255; %[0.3, 0.5, 0.8]; % Example color for "post" scatter points

%% Descriptive stats 
% Average for M646
M646_freq.pre = mean([M646.sess1.freq; M646.sess2.freq;  M646.sess3.freq; M646.sess4.freq; M646.sess5.freq; M646.sess6.freq;  M646.sess7.freq;  M646.sess8.freq]); 
M646_dur.pre = mean([M646.sess1.dur; M646.sess2.dur;  M646.sess3.dur; M646.sess4.dur; M646.sess5.dur; M646.sess6.dur;  M646.sess7.dur;  M646.sess8.dur]); 
M646_swr_count.pre = mean([M646.sess1.swr_count; M646.sess2.swr_count;  M646.sess3.swr_count; M646.sess4.swr_count; M646.sess5.swr_count; M646.sess6.swr_count;  M646.sess7.swr_count;  M646.sess8.swr_count]); 

% Average for M648
M648_freq.pre = mean([M648.sess1.freq; M648.sess2.freq;  M648.sess3.freq; M648.sess4.freq; M648.sess5.freq; M648.sess6.freq;  M648.sess7.freq;  M648.sess8.freq]); 
M648_dur.pre = mean([M648.sess1.dur; M648.sess2.dur;  M648.sess3.dur; M648.sess4.dur; M648.sess5.dur; M648.sess6.dur;  M648.sess7.dur;  M648.sess8.dur]); 
M648_swr_count.pre = mean([M648.sess1.swr_count; M648.sess2.swr_count;  M648.sess3.swr_count; M648.sess4.swr_count; M648.sess5.swr_count; M648.sess6.swr_count;  M648.sess7.swr_count;  M648.sess8.swr_count]); 

% Average for M650
M650_freq.pre = mean([M650.sess1.freq; M650.sess2.freq;  M650.sess3.freq; M650.sess4.freq; M650.sess5.freq; M650.sess6.freq;  M650.sess7.freq;  M650.sess8.freq]); 
M650_dur.pre = mean([M650.sess1.dur; M650.sess2.dur;  M650.sess3.dur; M650.sess4.dur; M650.sess5.dur; M650.sess6.dur;  M650.sess7.dur;  M650.sess8.dur]); 
M650_swr_count.pre = mean([M650.sess1.swr_count; M650.sess2.swr_count;  M650.sess3.swr_count; M650.sess4.swr_count; M650.sess5.swr_count; M650.sess6.swr_count;  M650.sess7.swr_count;  M650.sess8.swr_count]); 

% Average for M654
M654_freq.pre = mean([M654.sess1.freq; M654.sess2.freq;  M654.sess3.freq; M654.sess4.freq; M654.sess5.freq; M654.sess6.freq;  M654.sess7.freq;  M654.sess8.freq]); 
M654_dur.pre = mean([M654.sess1.dur; M654.sess2.dur;  M654.sess3.dur; M654.sess4.dur; M654.sess5.dur; M654.sess6.dur;  M654.sess7.dur;  M654.sess8.dur]); 
M654_swr_count.pre = mean([M654.sess1.swr_count; M654.sess2.swr_count;  M654.sess3.swr_count; M654.sess4.swr_count; M654.sess5.swr_count; M654.sess6.swr_count;  M654.sess7.swr_count;  M654.sess8.swr_count]); 

%% mean and plot over mice

low_c = [82,137,199]./255;%[0,104,87]./255; 
color2 = [123,175,222]./255;

mean_freq = mean([M646_freq.pre; M648_freq.pre; M650_freq.pre; M654_freq.pre]);
mean_dur = mean([M646_dur.pre; M648_dur.pre; M650_dur.pre; M654_dur.pre]);
mean_swr_count = mean([M646_swr_count.pre; M648_swr_count.pre;  M650_swr_count.pre; M654_swr_count.pre]);

mean_freq_pre = mean_freq(1)
mean_freq_track = mean_freq(2)
mean_freq_post = mean_freq(3)

mean_dur_pre = mean_dur(1)
mean_dur_track = mean_dur(2)
mean_dur_post = mean_dur(3)

mean_swr_count_pre = mean_swr_count(1)
mean_swr_count_track = mean_swr_count(2)
mean_swr_count_post = mean_swr_count(3)

std_freq_pre = std([M646_freq.pre(1); M648_freq.pre(1); M650_freq.pre(1); M654_freq.pre(1)])
std_freq_track = std([M646_freq.pre(2); M648_freq.pre(2); M650_freq.pre(2); M654_freq.pre(2)])
std_freq_post= std([M646_freq.pre(3); M648_freq.pre(3); M650_freq.pre(3); M654_freq.pre(3)])

std_dur_pre = std([M646_dur.pre(1); M648_dur.pre(1); M650_dur.pre(1); M654_dur.pre(1)])
std_dur_track= std([M646_dur.pre(2); M648_dur.pre(2); M650_dur.pre(2); M654_dur.pre(2)])
std_dur_post= std([M646_dur.pre(3); M648_dur.pre(3); M650_dur.pre(3); M654_dur.pre(3)])

std_swr_count_pre = std([M646_swr_count.pre(1); M648_swr_count.pre(1); M650_swr_count.pre(1); M654_swr_count.pre(1)])
std_swr_count_track = std([M646_swr_count.pre(2); M648_swr_count.pre(2); M650_swr_count.pre(2); M654_swr_count.pre(2)])
std_swr_count_post= std([M646_swr_count.pre(3); M648_swr_count.pre(3); M650_swr_count.pre(3); M654_swr_count.pre(3)])

% PRE IS WAKE AND POST IS NREM !!! TOO LAZY TO CHANGE !!!
pre_freq = [M646_freq.pre(1); M648_freq.pre(1); M650_freq.pre(1); M654_freq.pre(1)];
track_freq = [M646_freq.pre(2); M648_freq.pre(2); M650_freq.pre(2); M654_freq.pre(2)];
post_freq = [M646_freq.pre(3); M648_freq.pre(3); M650_freq.pre(3); M654_freq.pre(3)];
pre_dur = [M646_dur.pre(1); M648_dur.pre(1); M650_dur.pre(1); M654_dur.pre(1)];
track_dur = [M646_dur.pre(2); M648_dur.pre(2); M650_dur.pre(2); M654_dur.pre(2)];
post_dur = [M646_dur.pre(3); M648_dur.pre(3); M650_dur.pre(3); M654_dur.pre(3)];
pre_swr_count = [M646_swr_count.pre(1); M648_swr_count.pre(1); M650_swr_count.pre(1); M654_swr_count.pre(1)];
track_swr_count = [M646_swr_count.pre(2); M648_swr_count.pre(2); M650_swr_count.pre(2); M654_swr_count.pre(2)];
post_swr_count = [M646_swr_count.pre(3); M648_swr_count.pre(3); M650_swr_count.pre(3); M654_swr_count.pre(3)];


%% PLOT FREQUENCY AND STATS
%cd ('C:\Users\mimia\Desktop\BOXCHARTS')
% mean freq
figure(3)
h = boxchart([pre_freq, track_freq, post_freq]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric
% Scatter points with jitter for better visibility
scatter(repmat(xDataNumeric(1), 1, numel(pre_freq)) + randn(1, numel(pre_freq)) * 0.05, pre_freq, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(2), 1, numel(track_freq)) + randn(1, numel(track_freq)) * 0.05, track_freq, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(3), 1, numel(post_freq)) + randn(1, numel(post_freq)) * 0.05, post_freq, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

xticklabels(["Pre" "Track" "Post"])
ylabel("Rate of SWRs (SWRs per minute)")
xlabel("Rest or Track Period")
set(gca,'fontsize', 18)
% p =    6.8861e-04     tstat: -4.3327

set(gcf, 'renderer', 'painters');%
%cd ('C:\Users\mimia\OneDrive\Desktop\SWR_supplement')
%exportgraphics(gcf, 'SWR_freq_mean.eps', 'ContentType','vector');  % Export as PDF
hold off 

% Arrange as matrix: rows = animals, columns = conditions
freq_matrix = [pre_freq, track_freq, post_freq]; % n_mice x 3
% Repeated measures ANOVA
t = array2table(freq_matrix, 'VariableNames', {'Pre','Track','Post'});
within = table({'Pre';'Track';'Post'}, 'VariableNames', {'Session'});
rm = fitrm(t, 'Pre,Track,Post ~ 1', 'WithinDesign', within);
ranova_results = ranova(rm);
disp(ranova_results)

% Post-hoc paired tests with Bonferroni correction
pairs = {'pre_freq', 'track_freq'; 'pre_freq', 'post_freq'; 'track_freq', 'post_freq'};
data  = {pre_freq, track_freq, post_freq};
names = {'Pre', 'Track', 'Post'};
alpha_corrected = 0.05 / 3; % Bonferroni

fprintf('\nPost-hoc paired t-tests (Bonferroni corrected alpha = %.4f):\n', alpha_corrected);
combos = nchoosek(1:3, 2);
for i = 1:size(combos, 1)
    a = combos(i,1); b = combos(i,2);
    [~, p, ~, s] = ttest(data{a}, data{b});
    sig = p < alpha_corrected;
    fprintf('%s vs %s: t(%d) = %.3f, p = %.4f %s\n', ...
        names{a}, names{b}, s.df, s.tstat, p, ...
        string(repmat('*', 1, sig)));
end

%% PLOT mean dur
figure(4)
h = boxchart([pre_dur track_dur, post_dur]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric
% Scatter points with jitter for better visibility
scatter(repmat(xDataNumeric(1), 1, numel(pre_dur)) + randn(1, numel(pre_dur)) * 0.05, pre_dur, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(2), 1, numel(track_dur)) + randn(1, numel(track_dur)) * 0.05, track_dur, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(3), 1, numel(post_dur)) + randn(1, numel(post_dur)) * 0.05, post_dur, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

xticklabels(["Pre" "Track" "Post"])
ylabel("SWR duration (ms)")
xlabel("Rest or Track Period")
set(gca,'fontsize', 18)

set(gcf, 'renderer', 'painters');
%cd ('C:\Users\mimia\OneDrive\Desktop\SWR_supplement')
%exportgraphics(gcf, 'SWR_dur_mean.eps', 'ContentType','vector');  % Export as PDF
hold off 

%[h_dur,p_dur,ci_dur,stats_dur] = ttest2(pre_dur, post_dur)

% Arrange as matrix: rows = animals, columns = conditions
dur_matrix = [pre_dur, track_dur, post_dur]; % n_mice x 3
% Repeated measures ANOVA
t = array2table(dur_matrix, 'VariableNames', {'Pre','Track','Post'});
within = table({'Pre';'Track';'Post'}, 'VariableNames', {'Session'});
rm = fitrm(t, 'Pre,Track,Post ~ 1', 'WithinDesign', within);
ranova_results = ranova(rm);
disp(ranova_results)

% Post-hoc paired tests with Bonferroni correction
pairs = {'pre_dur', 'track_dur'; 'pre_dur', 'post_dur'; 'track_dur', 'post_dur'};
data  = {pre_dur, track_dur, post_dur};
names = {'Pre', 'Track', 'Post'};
alpha_corrected = 0.05 / 3; % Bonferroni

fprintf('\nPost-hoc paired t-tests (Bonferroni corrected alpha = %.4f):\n', alpha_corrected);
combos = nchoosek(1:3, 2);
for i = 1:size(combos, 1)
    a = combos(i,1); b = combos(i,2);
    [~, p, ~, s] = ttest(data{a}, data{b});
    sig = p < alpha_corrected;
    fprintf('%s vs %s: t(%d) = %.3f, p = %.4f %s\n', ...
        names{a}, names{b}, s.df, s.tstat, p, ...
        string(repmat('*', 1, sig)));
end
%% mean swr_count
figure(5)
h = boxchart([pre_swr_count, track_swr_count, post_swr_count]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric
% Scatter points with jitter for better visibility
scatter(repmat(xDataNumeric(1), 1, numel(pre_swr_count)) + randn(1, numel(pre_swr_count)) * 0.05, pre_swr_count, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(2), 1, numel(track_swr_count)) + randn(1, numel(track_swr_count)) * 0.05, track_swr_count, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(3), 1, numel(post_swr_count)) + randn(1, numel(post_swr_count)) * 0.05, post_swr_count, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

xticklabels(["Pre" "Track" "Post"])
ylabel("# of SWRs")
xlabel("Rest or Track Period")
set(gca,'fontsize', 18)

set(gcf, 'renderer', 'painters');
hold off 

% Arrange as matrix: rows = animals, columns = conditions
count_matrix = [pre_swr_count, track_swr_count, post_swr_count]; % n_mice x 3
% Repeated measures ANOVA
t = array2table(count_matrix, 'VariableNames', {'Pre','Track','Post'});
within = table({'Pre';'Track';'Post'}, 'VariableNames', {'Session'});
rm = fitrm(t, 'Pre,Track,Post ~ 1', 'WithinDesign', within);
ranova_results = ranova(rm);
disp(ranova_results)

% Post-hoc paired tests with Bonferroni correction
pairs = {'pre_count', 'track_count'; 'pre_count', 'post_count'; 'track_count', 'post_count'};
data  = {pre_swr_count, track_swr_count, post_swr_count};
names = {'Pre', 'Track', 'Post'};
alpha_corrected = 0.05 / 3; % Bonferroni

fprintf('\nPost-hoc paired t-tests (Bonferroni corrected alpha = %.4f):\n', alpha_corrected);
combos = nchoosek(1:3, 2);
for i = 1:size(combos, 1)
    a = combos(i,1); b = combos(i,2);
    [~, p, ~, s] = ttest(data{a}, data{b});
    sig = p < alpha_corrected;
    fprintf('%s vs %s: t(%d) = %.3f, p = %.4f %s\n', ...
        names{a}, names{b}, s.df, s.tstat, p, ...
        string(repmat('*', 1, sig)));
end