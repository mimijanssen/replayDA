%% Mouse Average Plots
clear; clc;
cd 'D:\TrackAvg'
Files=dir('*.*');
for k=3:length(Files)-2
   FileNames=Files(k).name;
   sess.(['mouse',num2str(k-2)]) = load(FileNames);
end

%%
figure(1)
plot([5, 5], [-0.5 0.5], '--k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on
shadedErrorBar(sess.mouse1.sess.sess1.time(1,:),avg_fiber_pre,std_fiber_pre,'lineProps',{'-','color',low_c,'MarkerFaceColor',low_c})
plot(sess.mouse1.sess.sess1.time(1,:),avg_fiber_pre,'LineWidth',3,'Color',low_c)
xlim([1 9])
xticks([1 2 3 4 5 6 7 8 9])
ylim([-0.3 0.3])
xticklabels({'-4','','','','0','','','','4'})
title('Pre-Track Rest')
ylabel('Mean [DA] (z-score)')
xlabel('Time from SWR (s)')
legend('','signal','Location','northwest')
legend boxoff
set(gca,'fontsize', 18)
set(gcf, 'color', 'none');
set(gca, 'color', 'none');
set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Desktop\Track')
exportgraphics(gcf, 'pre_track_track.eps', 'ContentType','vector');  % Export as PDF
hold off

%%
figure(2)
plot([5, 5], [-0.5 0.5], '--k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on
shadedErrorBar(sess.mouse1.sess.sess1.time(1,:),avg_fiber_track,std_fiber_track,'lineProps',{'-','color',low_c,'MarkerFaceColor',low_c})
%shadedErrorBar(time_rpe_plot,mean_low,sem_low,'lineprops',{'-','color',low_c,'MarkerFaceColor',low_c});

plot(sess.mouse1.sess.sess1.time(1,:),avg_fiber_track,'LineWidth',3,'Color',low_c)
%xl = xline(4,'',{'SWR'});
%xl.LabelVerticalAlignment = 'top';
xlim([1 9])
xticks([1 2 3 4 5 6 7 8 9])
ylim([-0.3 0.3])
xticklabels({'-4','','','','0','','','','4'})
title('Track')
ylabel('Mean [DA] (z-score)')
xlabel('Time from SWR (s)')
legend('','signal','Location','northwest')
legend boxoff

set(gca,'fontsize', 18)
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');

set(gcf, 'renderer', 'painters');
%fontname("AvenirNext LT Pro Regular");

cd ('C:\Users\mimia\Desktop\Track')
exportgraphics(gcf, 'track.eps', 'ContentType','vector');  % Export as PDF

hold off

%%
figure(3)
plot([5, 5], [-0.5 0.5], '--k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on
shadedErrorBar(sess.mouse1.sess.sess1.time(1,:),avg_fiber_post,std_fiber_post,'lineProps',{'-','color',low_c,'MarkerFaceColor',low_c})
plot(sess.mouse1.sess.sess1.time(1,:),avg_fiber_post,'LineWidth',3,'Color',low_c)
xlim([1 9])
xticks([1 2 3 4 5 6 7 8 9])
ylim([-0.3 0.3])
xticklabels({'-4','','','','0','','','','4'})
title('Post-Track Rest')
ylabel('Mean [DA] (z-score)')
xlabel('Time from SWR (s)')
legend('','signal','Location','northwest')
legend boxoff
set(gca,'fontsize', 18)
set(gcf, 'color', 'none');
set(gca, 'color', 'none');
set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\Desktop\Track')
exportgraphics(gcf, 'post_track_track.eps', 'ContentType','vector');  % Export as PDF
hold off
