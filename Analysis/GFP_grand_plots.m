%% ttest 
clear ; clc
%% 
cd D:\Mouse_avg
load('GFPmouse_avg.mat')
addpath('C:\Users\mimia\OneDrive\Documents\Toolboxes\raacampbell-shadedErrorBar')

%%
n = 3; % number of mice 

% for mouse_fiber_post
% for 1:4000 (before swr)
% 1.5 seconds 
x1 = 3200:1:6400; %2001:1:4000; %1:1:4000;%2501:1:4001;
% for 4001:8001 (after swr)
x2 = 6401:1:9600;  %4001:1:6000; %4001:1:8000;%4501:1:6001;
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
% A1_pre = zeros(n,1);
% for i_pre = 1:1:n
%     A1_pre(i_pre,:) = trapz(x1, mouse_fiber_pre(i_pre,x1));
% end
% 
% A2_pre = zeros(n,1);
% for i_pre = 1:1:n
%     A2_pre(i_pre,:) = trapz(x2, mouse_fiber_pre(i_pre,x2));
% end

%% dF for Pre and Post Track Rest 

% take the max and min in a 2 second period and subtract them... 

n = 3; % number of mice 

% for mouse_fiber_post
% for 1:4000 (before swr)
%x1 = 2001:1:4001;
% for 4001:8001 (after swr)
%x2 = 4001:1:6001;
% area under the curve for each mouse ; row is mouse 
t_before = linspace(-2,0,2000); % in seconds
t_after = linspace(0,2,2000); % in seconds


dF1_post = zeros(n,1);
dF2_post = zeros(n,1);
dF1_pre = zeros(n,1);
dF2_pre = zeros(n,1);

dF1_post_time = zeros(n,1);
dF2_post_time = zeros(n,1);
dF1_pre_time = zeros(n,1);
dF2_pre_time = zeros(n,1);

% dF for post 
% before SWR
for i_post = 1:1:n
    dF1_post(i_post,:) = max(mouse_fiber_post(i_post,x1));%-min(mouse_fiber_post(i_post,x1));
    dF1_post_time(i_post,:) = t_before(find(mouse_fiber_post(i_post,x1) == max(mouse_fiber_post(i_post,x1))));

end
% after SWR
for i_post = 1:1:n
    dF2_post(i_post,:) = max(mouse_fiber_post(i_post,x2));%-min(mouse_fiber_post(i_post,x2));
    dF2_post_time(i_post,:) = t_after(find(mouse_fiber_post(i_post,x2) == max(mouse_fiber_post(i_post,x2))));
end

% dF for pre
% before SWR
for i_pre = 1:1:n
    dF1_pre(i_pre,:) = max(mouse_fiber_pre(i_pre,x1));%-min(mouse_fiber_pre(i_pre,x1));
    dF1_pre_time(i_pre,:) = t_before(find(mouse_fiber_pre(i_pre,x1) == max(mouse_fiber_pre(i_pre,x1))));
end
% after SWR 
for i_pre = 1:1:n
    dF2_pre(i_pre,:) = max(mouse_fiber_pre(i_pre,x2));%-min(mouse_fiber_pre(i_pre,x2));
    dF2_pre_time(i_pre,:) = t_after(find(mouse_fiber_pre(i_pre,x2) == max(mouse_fiber_pre(i_pre,x2))));
end

% undid the - min for the significance test
% i used - min for calculating change in dF for the t-test. 

% edited to save time as well.

%%  POST TRACK REST dF
% color 
low_c = [78,178,101]./255;%[0,104,87]./255; 
color2 = [144, 201, 135]./255;

% Create the boxchart figure
figure(5);
h = boxchart([dF1_post, dF2_post]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;

% Add error bars
%errorbar(1:2, [avg_A1_post, avg_A2_post], [std_pre, std_post], 'LineStyle', 'none', 'Color', color2, 'LineWidth', 2);

% Convert XData to numeric for scatter plotting
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric

% Generate jittered x-coordinates for each condition
jitterAmount = 0.05;
x_jitter_A1 = xDataNumeric(1) + randn(1, numel(dF1_post)) * jitterAmount;
x_jitter_A2 = xDataNumeric(2) + randn(1, numel(dF2_post)) * jitterAmount;

% Draw light grey lines connecting each subject's data point in "Before" and "After" using the same jitter
for i = 1:numel(dF1_post)
    plot([x_jitter_A1(i), x_jitter_A2(i)], [dF1_post(i), dF2_post(i)], ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey color, thinner line
end

% Scatter points with jitter for each condition
scatter(x_jitter_A1, dF1_post, 60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
scatter(x_jitter_A2, dF2_post, 60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

% Set x-tick labels and axis properties
xticklabels({"Before", "After"});
ylabel("Mean \Delta [DA] (z-score)");
title("Post-Track Rest");
%ylim([0 0.3]);

% Set figure background color to white
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
set(gca,'fontsize', 18)

fontname("AvenirNext LT Pro Regular");

set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\OneDrive\Desktop\figures')

%exportgraphics(gcf, 'dF_post_ticks.png', 'ContentType','vector');  % Export as PDF

% Show the figure
hold off;


%% PRE REST dF 

% Create the boxchart figure
figure(6);
h = boxchart([dF1_pre, dF2_pre]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;

% Add error bars
%errorbar(1:2, [avg_A1_post, avg_A2_post], [std_pre, std_post], 'LineStyle', 'none', 'Color', color2, 'LineWidth', 2);

% Convert XData to numeric for scatter plotting
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric

% Generate jittered x-coordinates for each condition
jitterAmount = 0.05;
x_jitter_A1 = xDataNumeric(1) + randn(1, numel(dF1_pre)) * jitterAmount;
x_jitter_A2 = xDataNumeric(2) + randn(1, numel(dF2_pre)) * jitterAmount;

% Draw light grey lines connecting each subject's data point in "Before" and "After" using the same jitter
for i = 1:numel(dF1_pre)
    plot([x_jitter_A1(i), x_jitter_A2(i)], [dF1_pre(i), dF2_pre(i)], ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey color, thinner line
end

% Scatter points with jitter for each condition
scatter(x_jitter_A1, dF1_pre, 60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
scatter(x_jitter_A2, dF2_pre, 60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

% Set x-tick labels and axis properties
xticklabels({"Before", "After"});
ylabel("Mean \Delta [DA] (z-score)");
title("Pre-Track Rest");
%ylim([0 0.3]);

% Set figure background color to white
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
set(gca,'fontsize', 18)

%fontname("AvenirNext LT Pro Regular");

set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Desktop')

%exportgraphics(gcf, 'dF_pre.eps', 'ContentType','vector');  % Export as PDF

% Show the figure
hold off;
%% Finding grand mean (average over mice) for the swr post peak
disp(dF2_post) % post sessions

disp(dF2_pre) % pre sessions

mouse_grand_dF = mean([dF2_post, dF2_pre],2) ; % same number of n so this is ok.
disp(mouse_grand_dF)

grand_dF = mean(mouse_grand_dF);
disp('grand dF')
disp(grand_dF)

% grandstd 
disp(std(mouse_grand_dF))

% time mean
mouse_grand_time = mean([dF2_post_time, dF2_pre_time],2) ; % same number of n so this is ok.
grand_time = mean(mouse_grand_time);
disp(mouse_grand_time)
disp(grand_time)

grand_time_std = std(mouse_grand_time);
disp(grand_time_std)

%
disp('pre session mean and std ')
pre_mean = mean(dF2_pre);
disp(pre_mean)
disp(std(dF2_pre))

pre_time_mean = mean(dF2_pre_time);
disp(pre_time_mean)
disp(std(dF2_pre_time))

disp('post')
post_mean = mean(dF2_post);
disp(post_mean)
disp(std(dF2_post))

post_time_mean = mean(dF2_post_time);
disp(post_time_mean)
disp(std(dF2_post_time))

%% paired t test 
[h_post_t,p_post_t, ci_post, stats_post] = ttest(dF1_post,dF2_post);
% p = 0.0142 
[h_pre_t,p_pre_t, ci_pre, stats_pre] = ttest(dF1_pre,dF2_pre);

%% Test for normality : 
% otherwise you would need to use a Wilcoxin test 

[h,p] = kstest(dF1_post)


%% signrank tests 
[p,h,stats] = signrank(dF1_post, dF2_post); 
[p,h,stats] = signrank(dF1_pre,dF2_pre); 

%% chat POST TRACK REST AUC 
% color 
low_c = [78,178,101]./255;%[0,104,87]./255; 
color2 = [144, 201, 135]./255;

% Calculate means and standard deviations
avg_A1_post = mean(A1_post);
avg_A2_post = mean(A2_post);
std_pre = std(A1_post);
std_post = std(A2_post);

% Create the boxchart figure
figure(5);
h = boxchart([A1_post, A2_post]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;


% Add error bars
%errorbar(1:2, [avg_A1_post, avg_A2_post], [std_pre, std_post], 'LineStyle', 'none', 'Color', color2, 'LineWidth', 2);

% Convert XData to numeric for scatter plotting
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric

% Scatter points with jitter for better visibility
scatter(repmat(xDataNumeric(1), 1, numel(A1_post)) + randn(1, numel(A1_post)) * 0.05, A1_post, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(2), 1, numel(A2_post)) + randn(1, numel(A2_post)) * 0.05, A2_post, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

% Set x-tick labels and axis properties
xticklabels({"Before", "After"});
ylabel("Mean [DA] AUC (z-score)");
title("Post-Track Rest AUC");
%ylim([0.03 0.2]);

% Set figure background color to white
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
set(gca,'fontsize', 18)

fontname("AvenirNext LT Pro Regular");

set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Desktop')

%exportgraphics(gcf, 'AUC_post.eps', 'ContentType','vector');  % Export as PDF

% Show the figure
hold off;

%% PRE TRACK REST

% Calculate means and standard deviations
avg_A1_pre = mean(A1_pre);
avg_A2_pre = mean(A2_pre);
std_A1_pre = std(A1_pre);
std_A2_pre = std(A2_pre);

% Create the boxchart figure
figure(4);
h = boxchart([A1_pre, A2_pre]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;
% Add error bars
%errorbar(1:2, [avg_A1_post, avg_A2_post], [std_pre, std_post], 'LineStyle', 'none', 'Color', color2, 'LineWidth', 2);

% Convert XData to numeric for scatter plotting
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric

% Scatter points with jitter for better visibility
scatter(repmat(xDataNumeric(1), 1, numel(A1_pre)) + randn(1, numel(A1_pre)) * 0.05, A1_pre, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(2), 1, numel(A2_pre)) + randn(1, numel(A2_pre)) * 0.05, A2_pre, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

% Set x-tick labels and axis properties
xticklabels({"Before", "After"});
%ylim([-100 400]);
ylabel("Mean [DA] AUC (z-score)");
title("Pre-Track Rest AUC");

% Set figure background color to white
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
set(gca,'fontsize', 18)

fontname("AvenirNext LT Pro Regular");

set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Desktop')

%exportgraphics(gcf, 'AUC_pre.eps', 'ContentType','vector');  % Export as PDF

% Show the figure
hold off;


%% Save the data to a .mat file
cd 'D:\Mouse_Avg\';
file_name = 'mouse_'; 
filename = append(file_name, "auc.mat");
save(filename,'A1_post','A1_pre','A2_pre','A2_post');

%% SWR COUNT PLOT 
% for M433 with 8 sessions
%pre_count = [RPE.RPE1.swr_count(1); RPE.RPE2.swr_count(1); RPE.RPE3.swr_count(1); RPE.RPE4.swr_count(1); RPE.RPE5.swr_count(1); RPE.RPE6.swr_count(1); RPE.RPE7.swr_count(1); RPE.RPE8.swr_count(1)];  
%post_count = [RPE.RPE1.swr_count(2); RPE.RPE2.swr_count(2); RPE.RPE3.swr_count(2); RPE.RPE4.swr_count(2); RPE.RPE5.swr_count(2); RPE.RPE6.swr_count(2); RPE.RPE7.swr_count(2); RPE.RPE8.swr_count(2)];  
% pre_count = [sess.mouse1.mean_pre(1); sess.mouse2.mean_pre(1); sess.mouse3.mean_pre(1); sess.mouse4.mean_pre(1); sess.mouse5.mean_pre(1) ];%RPE.RPE4.swr_count(1); RPE.RPE5.swr_count(1); RPE.RPE6.swr_count(1); RPE.RPE7.swr_count(1)];  
% post_count = [sess.mouse1.mean_post(1); sess.mouse2.mean_post(1); sess.mouse3.mean_post(1); sess.mouse4.mean_post(1); sess.mouse5.mean_post(1)]; %RPE.RPE4.swr_count(2); RPE.RPE5.swr_count(2); RPE.RPE6.swr_count(2); RPE.RPE7.swr_count(2)];  

% color 
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
%title("SWR Count")

%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
set(gca,'fontsize', 18)

set(gcf, 'renderer', 'painters');
fontname("AvenirNext LT Pro Regular");

%exportgraphics(gcf, 'SWR.eps', 'ContentType','vector');  % Export as PDF

hold off 
%% RPE PLOT 
% overall RPE average
% each row is a session. 
%SEM_amp = std([M1_mean_amp;M2_mean_amp],0,1)/sqrt(length([M1_mean_amp;M2_mean_amp]));

% high 
sess_high = mouse_high; %RPE.RPE4.avg_high; RPE.RPE5.avg_high; RPE.RPE6.avg_high; RPE.RPE7.avg_high; ];  %RPE.RPE8.avg_high
mean_high = mean(sess_high);
sem_high = std(sess_high)/sqrt(size(sess_high,1)); % this is actually SEM 
% medium
sess_med = mouse_med; %RPE.RPE4.avg_med; RPE.RPE5.avg_med; RPE.RPE6.avg_med; RPE.RPE7.avg_med;];  % RPE.RPE8.avg_med
mean_med = mean(sess_med);
sem_med = std(sess_med)/sqrt(size(sess_med,1));
% low
sess_low = mouse_low; %RPE.RPE4.avg_low; RPE.RPE5.avg_low; RPE.RPE6.avg_low; RPE.RPE7.avg_low; ];  %RPE.RPE8.avg_low
mean_low = mean(sess_low);
sem_low = std(sess_low)/sqrt(size(sess_low,1));

med_c = [78,178,101]./255;%[104,187,225]./255; % rgb(167, 199, 231) rgb(255, 165, 0)  blue: rgb(104,187,227)
low_c = [157,130,187]./255;
%low_c = [82,137,199]./255;%[0,104,87]./255; 
high_c = [255, 165,0]./255; % rgb(204, 204, 255) periwinkle

mean_high2 = decimate(mean_high, 10);
time_rpe_plot2 = decimate(time_rpe_plot,10);
sem_high2 = decimate(sem_high,10);
mean_med2 = decimate(mean_med, 10);
sem_med2 = decimate(sem_med,10);
mean_low2 = decimate(mean_low, 10);
sem_low2 = decimate(sem_low,10);

%fig = figure('units','inch','position', [0, 0, 6, 4]);
figure(1)
shadedErrorBar(time_rpe_plot2,mean_high2,sem_high2,'lineprops',{'-','color',high_c,'MarkerFaceColor',high_c});
hold on
plot(time_rpe_plot2,mean_high2,'LineWidth',3,'Color',high_c)
hold on
shadedErrorBar(time_rpe_plot2,mean_med2,sem_med2,'lineprops',{'-','color',med_c,'MarkerFaceColor',med_c});
hold on
plot(time_rpe_plot2,mean_med2,'LineWidth',3,'Color',med_c)
shadedErrorBar(time_rpe_plot2,mean_low2,sem_low2,'lineprops',{'-','color',low_c,'MarkerFaceColor',low_c});
plot(time_rpe_plot2,mean_low2,'LineWidth',3,'Color',low_c)

set(gca,'fontsize', 18)
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');

set(gcf, 'renderer', 'painters');

xlabel('Time from photobeam break (s)')
xlim([0 16])
ylim([-1 1.5])
xticks([0 2 4 6 8 10 12 14 16])
xticklabels({'-8','','','','0','','','','8'})
legend('high','','medium','','low','Location','Best');
legend boxoff
ylabel('Mean [DA] (z-score)')
title('Fiber Signal RPE')
%fontname("AvenirNext LT Pro Regular");
cd ('C:\Users\mimia\OneDrive\Desktop\updated figures')
exportgraphics(gcf, 'RPE_GFP_resave.eps', 'ContentType','vector');  % Export as PDF

hold off

%% RPE T-Test Paired
% are the peaks different from before
%peak1_high_entry = zeros(5,1);

% 1.5 seconds 
%x1 = 3200:1:6400; %2001:1:4000; %1:1:4000;%2501:1:4001;
% for 4001:8001 (after swr)
%x2 = 6401:1:9600;  %4001:1:6000; %4001:1:8000;%4501:1:6001;

peak1_high = zeros(3,1);
transient1 = 3200:1:6400; %4001:1:8000; % changed this so that I only take the peak or trough for 4 seconds not 8. 

for i = 1:1:3
    max_val = max(sess_high(i,transient1));
    peak1_high(i,:) = max_val;
end

% second half of high reward rpe 
peak2_high = zeros(3,1);
transient2 = 6401:1:9600; %8001:1:12000;

for i = 1:1:3
    max_val = max(sess_high(i,transient2));
    peak2_high(i,:) = max_val;
end

[h,p_high] = ttest(peak1_high, peak2_high);
high_mean = mean(peak2_high);
high_std = std(peak2_high);

%[p,h,stats] = signrank(peak1_high, peak2_high);
% p value = 0.0625 

% medium 
peak1_med = zeros(3,1);
peak2_med = zeros(3,1);
for i = 1:1:3
    max_val = max(sess_med(i,transient1));
    peak1_med(i,:) = max_val;
end
for i = 1:1:3
    max_val = max(sess_med(i,transient2));
    peak2_med(i,:) = max_val;
end
%[p,h,stats] = signrank(peak1_med, peak2_med);
[h,p_med] = ttest(peak1_med, peak2_med);
med_mean = mean(peak2_med);
med_std = std(peak2_med);

% p value = 0.3125
peak1_low = zeros(3,1);
peak2_low = zeros(3,1);
for i = 1:1:3
    min_val = min(sess_low(i,transient1));
    peak1_low(i,:) = min_val;
end
for i = 1:1:3
    min_val = min(sess_low(i,transient2));
    peak2_low(i,:) = min_val;
end

[h,p_low] = ttest(peak1_low, peak2_low);
low_mean = mean(peak2_low);
low_std = std(peak2_low);
%[p,h,stats] = signrank(peak1_low, peak2_low);
% p = 0.1250 

%%% against each other? 
[h,p_against,ci,stats] = ttest2(peak2_low, peak2_high)
%[h,p_against,ci,stats] = ttest2(peak2_low, peak2_med)
%[h,p_against,ci,stats] = ttest2(peak2_med, peak2_high)


rpe_peaks = [peak2_low, peak2_med, peak2_high];
x = {'omission','medium','high'};

boxchart(rpe_peaks);


figure(2)
boxplot(rpe_peaks,x);
hold on
scatter([1,2,3],rpe_peaks)
ylabel('Mean [DA] Peak Values (z-score)')
xlabel('Reward Type')


% are the peaks different from each other??
% pearson correlations
% linear... need to do this laters... 


% need to test if this is normal... so you would need to run chi2gof on all
% rpe values... 
[h_c_high,p_c_high,stats_c_high] = chi2gof(peak2_high)
[h_c_low,p_c_low,stats_c_low] = chi2gof(peak2_low)


%%
grand_dF = mean(mouse_grand_dF);
grand_rpe_high = mean(peak2_high);

percentage_diff = grand_dF/grand_rpe_high; 
%% t test for AUC

%save(filename,'A1_post','A1_pre','A2_pre','A2_post');


[p_post,h_post,stats_post] = signrank(A1_post, A2_post);
% p = 0.0391
[p_pre,h_pre,stats_pre] = signrank(A1_pre, A2_pre);
% p = 0.2500


% [h,p,ci,stats]
[h_post_t,p_post_t, ci_post, stats_post] = ttest(A1_post,A2_post);
% p = 0.0142 
[h_pre_t,p_pre_t, ci_pre, stats_pre] = ttest(A1_pre,A2_pre);
% p = 0.2772 



% paired t test 
% [h,p,ci,stats]
[h_post_t,p_post_t, ci_post, stats_post] = ttest(dF1_post,dF2_post);
% p = 0.0142 
[h_pre_t,p_pre_t, ci_pre, stats_pre] = ttest(dF1_pre,dF2_pre);
% p = 0.2772 


%% PRE TRACK PETH 
% Pre-Track Rest PETH 
% Average of Session circshifted signals 
% Average of Session DA signals 
% 
% %sess_circ_pre = [sess.sess1.circ_avg_pre;sess.sess2.circ_avg_pre;sess.sess3.circ_avg_pre;]; %sess.sess4.circ_avg_pre;sess.sess5.circ_avg_pre;sess.sess6.circ_avg_pre;sess.sess7.circ_avg_pre;]; %sess.sess8.circ_avg_pre
% sess_fiber_pre = [sess.mouse1.avg_fiber_pre; sess.mouse2.avg_fiber_pre; sess.mouse3.avg_fiber_pre; sess.mouse4.avg_fiber_pre; sess.mouse5.avg_fiber_pre]; %sess.sess4.avg_fiber_pre;sess.sess5.avg_fiber_pre;sess.sess6.avg_fiber_pre;sess.sess7.avg_fiber_pre;]; %sess.sess8.avg_fiber_pre
% %sess_circ_post = [sess.sess1.circ_avg_post;sess.sess2.circ_avg_post;sess.sess3.circ_avg_post;]; %sess.sess4.circ_avg_post;sess.sess5.circ_avg_post;sess.sess6.circ_avg_post;sess.sess7.circ_avg_post;]; %sess.sess8.circ_avg_post
% sess_fiber_post = [sess.mouse1.avg_fiber_post; sess.mouse2.avg_fiber_post; sess.mouse3.avg_fiber_post; sess.mouse4.avg_fiber_post; sess.mouse5.avg_fiber_post]; %sess.sess4.avg_fiber_post;sess.sess5.avg_fiber_post;sess.sess6.avg_fiber_post;sess.sess7.avg_fiber_post;]; %sess.sess8.avg_fiber_post
% 
% %circ_avg_fiber_post = mean(sess_circ_post);
% %circ_avg_fiber_pre = mean(sess_circ_pre);
% %circ_std_fiber_post = 2*std(sess_circ_post);
% %circ_std_fiber_pre = 2*std(sess_circ_pre);
% avg_fiber_pre = mean(sess_fiber_pre);
% avg_fiber_post = mean(sess_fiber_post);
sem_fiber_pre = std(mouse_fiber_pre)/sqrt(size(mouse_fiber_pre,1)); % just to make sure we are calculating SEM
sem_fiber_post = std(mouse_fiber_post)/sqrt(size(mouse_fiber_post,1));
% % took the average of the circ shifted signal... is that legit?
dark_green = [78,178,101]./255;%[0,104,87]./255; 
light_green = [144, 201, 135]./255;

figure(2)
%shadedErrorBar(sess.sess1.time(1,:),circ_avg_fiber_pre,circ_std_fiber_pre,'lineProps','-k','transparent',1)
%hold on
%plot(sess.sess1.time(1,:),circ_avg_fiber_pre,'LineWidth',3,'Color','k')
% plot average on top with larger line
%hold on
plot([4, 4], [-0.5 0.5], '--k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on
shadedErrorBar(time_swr_da_plot,avg_fiber_pre,sem_fiber_pre,'lineProps',{'-','color',light_green,'MarkerFaceColor',light_green})
%shadedErrorBar(time_rpe_plot,mean_low,sem_low,'lineprops',{'-','color',low_c,'MarkerFaceColor',low_c});

plot(time_swr_da_plot,avg_fiber_pre,'LineWidth',3,'Color',dark_green)
%xl = xline(4,'',{'SWR'});
%xl.LabelVerticalAlignment = 'top';
xlim([0 8])
xticks([0 1 2 3 4 5 6 7 8])
ylim([-0.1 0.2])
xticklabels({'-4','','','','0','','','','4'})
title('Pre-Track Rest [DA] after SWRs')
ylabel('Mean [DA] (z-score)')
xlabel('Time from SWR (s)')
legend('','signal','Location','northwest')
legend boxoff

set(gca,'fontsize', 18)
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');

set(gcf, 'renderer', 'painters');
%fontname("AvenirNext LT Pro Regular");
cd ('C:\Users\mimia\OneDrive\Desktop\updated figures')
exportgraphics(gcf, 'Pretrack_PETH_GFP.png', 'ContentType','vector');  % Export as PDF

hold off

%% POST TRACK

dark_green = [78,178,101]./255;%[0,104,87]./255; 
light_green = [144, 201, 135]./255;

figure(2)

plot([4, 4], [-0.5 0.5], '--k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on
shadedErrorBar(time_swr_da_plot,avg_fiber_post,std_fiber_post,'lineProps',{'-','color',light_green,'MarkerFaceColor',light_green})
%shadedErrorBar(time_rpe_plot,mean_low,sem_low,'lineprops',{'-','color',low_c,'MarkerFaceColor',low_c});

plot(time_swr_da_plot,avg_fiber_post,'LineWidth',3,'Color',dark_green)
%xl = xline(4,'',{'SWR'});
%xl.LabelVerticalAlignment = 'top';
xlim([0 8])
xticks([0 1 2 3 4 5 6 7 8])
ylim([-0.1 0.2])
xticklabels({'-4','','','','0','','','','4'})
title('Post-Track Rest [DA] after SWRs')
ylabel('Mean [DA] (z-score)')
xlabel('Time from SWR (s)')
legend('','signal','Location','northwest')
legend boxoff

set(gca,'fontsize', 18)
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');

set(gcf, 'renderer', 'painters');
fontname("AvenirNext LT Pro Regular");
cd ('C:\Users\mimia\OneDrive\Desktop\updated figures')
exportgraphics(gcf, 'Posttrack_PETH_ticks_GFP.png', 'ContentType','vector');  % Export as PDF

hold off



%% PEAK ANALYSIS FOR PRE AND POST 

% % pre
% peak1_pre = zeros(5,1);
% trans1 = 1:1:4000;
% % second half 
% peak2_pre = zeros(5,1);
% trans2 = 4001:1:8001;
% 
% for i = 1:1:5 
%     max_val = max(sess_fiber_pre(i,trans1));
%     peak1_pre(i,:) = max_val;
% end
% for i = 1:1:5 
%     max_val = max(sess_fiber_pre(i,trans2));
%     peak2_pre(i,:) = max_val;
% end
% 
% % post 
% peak1_post = zeros(5,1);
% trans1 = 1:1:4000;
% % second half 
% peak2_post = zeros(5,1);
% trans2 = 4001:1:8001;
% 
% for i = 1:1:5 
%     max_val = max(sess_fiber_post(i,trans1));
%     peak1_post(i,:) = max_val;
% end
% for i = 1:1:5 
%     max_val = max(sess_fiber_post(i,trans2));
%     peak2_post(i,:) = max_val;
% end
% 
% [h,p_pre] = ttest(peak1_pre, peak2_pre); 
% [h,p_post] = ttest(peak1_post, peak2_post); 
% 
% swr_da_peaks_pre = [peak1_pre, peak2_pre];
% swr_da_peaks_post = [peak1_post, peak2_post];
% swr_da_x = {'Before SWR', 'After SWR'};
% 
% figure(3)
% boxplot(swr_da_peaks_pre,swr_da_x);
% hold on
% scatter([1,2],swr_da_peaks_pre)
% ylabel('Mean [DA] Peak Values (z-score)')
% title('Pre-Track Rest')
% 
% figure(4)
% boxplot(swr_da_peaks_post,swr_da_x);
% hold on
% scatter([1,2],swr_da_peaks_post)
% ylabel('Mean [DA] Peak Values (z-score)')
% title('Post-Track Rest')


%% SAVE VARIABLES NEEDED FOR PLOTTING
% VARIABLES 

% pre_count = [sess.mouse1.mean_pre(1); sess.mouse2.mean_pre(1); sess.mouse3.mean_pre(1); sess.mouse4.mean_pre(1) ];%RPE.RPE4.swr_count(1); RPE.RPE5.swr_count(1); RPE.RPE6.swr_count(1); RPE.RPE7.swr_count(1)];  
% post_count = [sess.mouse1.mean_post(1); sess.mouse2.mean_post(1); sess.mouse3.mean_post(1); sess.mouse4.mean_post(1) ]; %RPE.RPE4.swr_count(2); RPE.RPE5.swr_count(2); RPE.RPE6.swr_count(2); RPE.RPE7.swr_count(2)];  
% mean_pre = mean(pre_count);
% mean_post = mean(post_count);
% std_pre = std(pre_count);
% std_post = std(post_count);

mouse_high = sess_high ; %RPE.RPE4.avg_high; RPE.RPE5.avg_high; RPE.RPE6.avg_high; RPE.RPE7.avg_high; ];  %RPE.RPE8.avg_high
% mean_high = mean(sess_high);
% std_high = std(sess_high)/sqrt(size(sess_high,1));
% % medium
mouse_med = sess_med;
% sess_med = [sess.mouse1.mean_med; sess.mouse2.mean_med; sess.mouse3.mean_med; sess.mouse4.mean_med]; %RPE.RPE4.avg_med; RPE.RPE5.avg_med; RPE.RPE6.avg_med; RPE.RPE7.avg_med;];  % RPE.RPE8.avg_med
% mean_med = mean(sess_med);
% std_med = std(sess_med)/sqrt(size(sess_med,1));
% % low
mouse_low = sess_low;
% sess_low = [sess.mouse1.mean_low; sess.mouse2.mean_low; sess.mouse3.mean_low; sess.mouse4.mean_low;]; %RPE.RPE4.avg_low; RPE.RPE5.avg_low; RPE.RPE6.avg_low; RPE.RPE7.avg_low; ];  %RPE.RPE8.avg_low
% mean_low = mean(sess_low);
% std_low = std(sess_low)/sqrt(size(sess_low,1));
% 
mouse_fiber_pre = sess_fiber_pre;
mouse_fiber_post = sess_fiber_post; 
% sess_fiber_pre = [sess.mouse1.avg_fiber_pre; sess.mouse2.avg_fiber_pre; sess.mouse3.avg_fiber_pre; sess.mouse4.avg_fiber_pre]; %sess.sess4.avg_fiber_pre;sess.sess5.avg_fiber_pre;sess.sess6.avg_fiber_pre;sess.sess7.avg_fiber_pre;]; %sess.sess8.avg_fiber_pre
% sess_fiber_post = [sess.mouse1.avg_fiber_post; sess.mouse2.avg_fiber_post; sess.mouse3.avg_fiber_post; sess.mouse4.avg_fiber_post]; %sess.sess4.avg_fiber_post;sess.sess5.avg_fiber_post;sess.sess6.avg_fiber_post;sess.sess7.avg_fiber_post;]; %sess.sess8.avg_fiber_post
% 
% avg_fiber_pre = mean(sess_fiber_pre);
% avg_fiber_post = mean(sess_fiber_post);
% std_fiber_pre = std(sess_fiber_pre)/sqrt(size(sess_fiber_pre,1));
% std_fiber_post = std(sess_fiber_post)/sqrt(size(sess_fiber_post,1));

time_swr_da_plot = sess.mouse1.sess.sess1.time(1,:);
time_rpe_plot = sess.mouse1.RPE.RPE1.t_shared;

cd 'D:\Mouse_Avg\'
file_name = 'mouse_'; 
filename = append(file_name, "avg.mat");
save(filename, 'pre_count', 'post_count', 'mean_pre', 'mean_post', 'std_pre', 'std_post', 'mouse_high','mouse_med','mouse_low','mean_low','mean_med','mean_high','std_low','std_med','std_high','time_rpe_plot','avg_fiber_post','avg_fiber_pre','std_fiber_post','std_fiber_pre','time_swr_da_plot','mouse_fiber_post','mouse_fiber_pre');


