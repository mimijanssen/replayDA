%% Grand Average Over Mice 

%% Mouse Average Plots
clear; clc;
cd 'D:\GFPavg'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   sess.(['mouse',num2str(k-2)]) = load(FileNames);
end


%% averaging an average because each set has the same number of sessions - for grand mean 
fiber_1 = mean([sess.mouse1.avg_fiber_pre ; sess.mouse1.avg_fiber_post]);
fiber_2 = mean([sess.mouse2.avg_fiber_pre ; sess.mouse2.avg_fiber_post]);
fiber_3 = mean([sess.mouse3.avg_fiber_pre ; sess.mouse3.avg_fiber_post]);

mouse_fiber = [fiber_1; fiber_2; fiber_3;];

avg_mouse = mean(mouse_fiber);
sem_mouse = std(mouse_fiber)/sqrt(size(mouse_fiber,1));
% took the average of the circ shifted signal... is that legit?

%% colors

dark_green = [78,178,101]./255;%[0,104,87]./255; 
light_green = [144, 201, 135]./255;



%%
figure(10)
%shadedErrorBar(sess.sess1.time(1,:),circ_avg_fiber_pre,circ_std_fiber_pre,'lineProps','-k','transparent',1)
%hold on
%plot(sess.sess1.time(1,:),circ_avg_fiber_pre,'LineWidth',3,'Color','k')
% plot average on top with larger line
%hold on
shadedErrorBar(sess.mouse1.sess.sess1.time(1,:),avg_mouse,sem_mouse,'lineProps','-g','transparent',1)
hold on
plot(sess.mouse1.sess.sess1.time(1,:),avg_mouse,'LineWidth',3,'Color',low_c)
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
xlim([0 8])
xticks([0 4 8])
xticklabels({'-4','0','4'})
title('Pre Track Rest PETH')
ylabel('Avg Signal (detrended & z-scored)')
xlabel('Time from SWR (s)')
legend('','signal','Location','northwest')
legend boxoff
% 
% set(gcf,'Color',[1,1,1])
% set(gca,'fontsize', 24)
% shg
hold off

%exportgraphics(gcf,'pretrackrest.eps','BackgroundColor','none','ContentType','vector');

%% BETTER PLOTTING? ~~~
figure(2)
plot([4, 4], [-0.5 0.5], '--k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on
shadedErrorBar(sess.mouse1.sess.sess1.time(1,:),avg_mouse,sem_mouse,'lineProps',{'-','color',light_green,'MarkerFaceColor',light_green})
%shadedErrorBar(time_rpe_plot,mean_low,sem_low,'lineprops',{'-','color',low_c,'MarkerFaceColor',low_c});

plot(sess.mouse1.sess.sess1.time(1,:),avg_mouse,'LineWidth',3,'Color',dark_green)
%xl = xline(4,'',{'SWR'});
%xl.LabelVerticalAlignment = 'top';
xlim([0 8])
xticks([0 1 2 3 4 5 6 7 8])
ylim([-0.1 0.2])
xticklabels({'-4','','','','0','','','','4'})
title('[DA] after SWRs')
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
exportgraphics(gcf, 'GrandGFPmouse.png', 'ContentType','vector');  % Export as PDF

hold off

%%
figure(3)
%shadedErrorBar(sess.mouse1.sess.sess1.time(1,:),circ_avg_fiber_post,circ_std_fiber_post,'lineProps','-k','transparent',1)
%hold on
%plot(sess.sess1.time(1,:),circ_avg_fiber_post,'LineWidth',3,'Color','k')
%hold on
shadedErrorBar(sess.mouse1.sess.sess1.time(1,:),avg_fiber_post,std_fiber_post,'lineProps','-g','transparent',1)
hold on
plot(sess.mouse1.sess.sess1.time(1,:),avg_fiber_post,'LineWidth',3,'Color',low_c)
hold on
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
xlim([0 8])
xticks([0 4 8])
xticklabels({'-4','0','4'})
title('Post Track Rest PETH')
ylabel('Avg Signal (detrended & z-scored)')
xlabel('Time from SWR (s)')
legend('','signal','Location','northwest')
legend boxoff 

% set(gcf,'Color',[1,1,1])
% set(gca,'fontsize', 24)
% shg
hold off

exportgraphics(gcf,'posttrackrest.eps','BackgroundColor','none','ContentType','vector');

%% PRE TRACK PETH  zoom
figure(4)
%shadedErrorBar(sess.sess1.time(1,:),circ_avg_fiber_pre,circ_std_fiber_pre,'lineProps','-k','transparent',1)
%hold on
%plot(sess.sess1.time(1,:),circ_avg_fiber_pre,'LineWidth',3,'Color','k')
% plot average on top with larger line
%hold on
shadedErrorBar(sess.mouse1.sess.sess1.time(1,:),avg_fiber_pre,std_fiber_pre,'lineProps','-g','transparent',1)
hold on
plot(sess.mouse1.sess.sess1.time(1,:),avg_fiber_pre,'LineWidth',3,'Color',low_c)
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
hold off
xlim([3.75 4.25])
xticks([3.75 4 4.25])
xticklabels({'-0.25','0','0.25'})
title('Pre Track Rest PETH')
ylabel('Averaged Signal (zdF)')
xlabel('Time from SWR (s)')

set(gcf,'Color',[1,1,1])
set(gca,'fontsize', 24)
shg
hold off

figure(5)
%shadedErrorBar(sess.mouse1.sess.sess1.time(1,:),circ_avg_fiber_post,circ_std_fiber_post,'lineProps','-k','transparent',1)
%hold on
%plot(sess.sess1.time(1,:),circ_avg_fiber_post,'LineWidth',3,'Color','k')
%hold on
shadedErrorBar(sess.mouse1.sess.sess1.time(1,:),avg_fiber_post,std_fiber_post,'lineProps','-g','transparent',1)
hold on
plot(sess.mouse1.sess.sess1.time(1,:),avg_fiber_post,'LineWidth',3,'Color',low_c)
hold on
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
xlim([3.75 4.25])
xticks([3.75 4 4.25])
xticklabels({'-0.25','0','0.25'})
title('Post Track Rest PETH')
ylabel('Averaged Signal (zdF)')
xlabel('Time from SWR (s)')

set(gcf,'Color',[1,1,1])
set(gca,'fontsize', 24)
shg
hold off

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
