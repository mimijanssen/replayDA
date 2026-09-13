%% Mouse Average Plots for Sleep/Wake 
clear; clc;
cd 'D:\M545\sleep'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   sess.(['sess',num2str(k-2)]) = load(FileNames);
end

%% COLOR
med_c = [104,187,225]./255; % rgb(167, 199, 231) rgb(255, 165, 0)  blue: rgb(104,187,227)
low_c = [0,104,87]./255; 
high_c = [255, 165,0]./255; % rgb(204, 204, 255) periwinkle

%% NREM AND WAKE PETH 
% nrem-Track Rest PETH 
% Average of Session circshifted signals 
% Average of Session DA signals 

sess_circ_nrem = [sess.sess1.circ_avg_nrem;sess.sess2.circ_avg_nrem;sess.sess3.circ_avg_nrem;sess.sess4.circ_avg_nrem;];%sess.sess5.circ_avg_nrem;sess.sess6.circ_avg_nrem;];%sess.sess7.circ_avg_nrem;]; %sess.sess8.circ_avg_nrem
sess_fiber_nrem = [sess.sess1.avg_fiber_nrem;sess.sess2.avg_fiber_nrem;sess.sess3.avg_fiber_nrem; sess.sess4.avg_fiber_nrem;];%sess.sess5.avg_fiber_nrem;sess.sess6.avg_fiber_nrem;];%sess.sess7.avg_fiber_nrem;]; %sess.sess8.avg_fiber_nrem
sess_circ_wake = [sess.sess1.circ_avg_wake;sess.sess2.circ_avg_wake;sess.sess3.circ_avg_wake;sess.sess4.circ_avg_wake;];%sess.sess5.circ_avg_wake;sess.sess6.circ_avg_wake;];%sess.sess7.circ_avg_wake;]; %sess.sess8.circ_avg_wake
sess_fiber_wake = [sess.sess1.avg_fiber_wake;sess.sess2.avg_fiber_wake;sess.sess3.avg_fiber_wake;sess.sess4.avg_fiber_wake;];%sess.sess5.avg_fiber_wake;sess.sess6.avg_fiber_wake;];%sess.sess7.avg_fiber_wake;]; %sess.sess8.avg_fiber_wake

circ_avg_fiber_wake = mean(sess_circ_wake);
circ_avg_fiber_nrem = mean(sess_circ_nrem);
circ_std_fiber_wake = 2*std(sess_circ_wake);
circ_std_fiber_nrem = 2*std(sess_circ_nrem);
avg_fiber_nrem = mean(sess_fiber_nrem);
avg_fiber_wake = mean(sess_fiber_wake);
std_fiber_nrem = std(sess_fiber_nrem);
std_fiber_wake = std(sess_fiber_wake);
% took the average of the circ shifted signal... is that legit?

time = linspace(0,8,8001);
figure(2)
%shadedErrorBar(sess.sess1.time(1,:),circ_avg_fiber_nrem,circ_std_fiber_nrem,'lineProps','-k','transparent',1)
%hold on
%plot(sess.sess1.time(1,:),circ_avg_fiber_nrem,'LineWidth',3,'Color','k')
% plot average on top with larger line
%hold on
shadedErrorBar(time,avg_fiber_nrem,std_fiber_nrem,'lineProps','-g','transparent',1)
hold on
plot(time,avg_fiber_nrem,'LineWidth',3,'Color',low_c)
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
xlim([0 8])
xticks([0 4 8])
xticklabels({'-4','0','4'})
title('M545: NREM Rest PETH')
ylabel('Averaged Signal (zdF)')
xlabel('Time from SWR (s)')
legend('','signal','','signal','Location','northwest')
legend boxoff

set(gcf,'Color',[1,1,1])
shg
hold off

figure(3)
shadedErrorBar(time,avg_fiber_wake,std_fiber_wake,'lineProps','-g','transparent',1)
hold on
plot(time,avg_fiber_wake,'LineWidth',3,'Color',low_c)
hold on
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
xlim([0 8])
xticks([0 4 8])
xticklabels({'-4','0','4'})
title('M545: Quiet Wakefulness PETH')
ylabel('Averaged Signal (zdF)')
xlabel('Time from SWR (s)')
legend('','signal','','signal','Location','northwest')
legend boxoff 

set(gcf,'Color',[1,1,1])
shg
hold off

%% save variables 
cd 'D:\M545'
file_name = 'M545'; 
filename = append(file_name, "avg_sleep.mat");
save(filename, 'avg_fiber_wake','avg_fiber_nrem','std_fiber_wake','std_fiber_nrem','sess');