% SAVE SESS AVG EARLY AND LATE 
%% Mouse Average Plots
clear; clc;
cd 'F:\M545\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   sess.(['sess',num2str(k-2)]) = load(FileNames);
end

%% EXTRACT EARLY AND LATE PETH 
% M433 
% sess_fiber_pre_early = [sess.sess1.avg_fiber_pre;sess.sess2.avg_fiber_pre;sess.sess3.avg_fiber_pre;];% sess.sess7.avg_fiber_pre;]; %sess.sess8.avg_fiber_pre
% sess_fiber_pre_late = [sess.sess4.avg_fiber_pre;sess.sess5.avg_fiber_pre;sess.sess6.avg_fiber_pre;sess.sess7.avg_fiber_pre;];% sess.sess7.avg_fiber_pre;]; %sess.sess8.avg_fiber_pre
% sess_fiber_post_early = [sess.sess1.avg_fiber_post;sess.sess2.avg_fiber_post;sess.sess3.avg_fiber_post;];%sess.sess7.avg_fiber_post;]; %sess.sess8.avg_fiber_post
% sess_fiber_post_late = [sess.sess4.avg_fiber_post;sess.sess5.avg_fiber_post;sess.sess6.avg_fiber_post;sess.sess7.avg_fiber_post;];%sess.sess7.avg_fiber_post;]; %sess.sess8.avg_fiber_post

% M453 
% sess_fiber_pre_early = [sess.sess1.avg_fiber_pre;sess.sess2.avg_fiber_pre;sess.sess3.avg_fiber_pre;];
% sess_fiber_pre_late = [sess.sess4.avg_fiber_pre;sess.sess5.avg_fiber_pre;];
% sess_fiber_post_early = [sess.sess1.avg_fiber_post;sess.sess2.avg_fiber_post;sess.sess3.avg_fiber_post;];
% sess_fiber_post_late = [sess.sess4.avg_fiber_post;sess.sess5.avg_fiber_post;];

% M460 or M545
sess_fiber_pre_early = [sess.sess1.avg_fiber_pre;];
sess_fiber_pre_late = [sess.sess2.avg_fiber_pre;sess.sess3.avg_fiber_pre;sess.sess4.avg_fiber_pre;];
sess_fiber_post_early = [sess.sess1.avg_fiber_post;];
sess_fiber_post_late = [sess.sess2.avg_fiber_post;sess.sess3.avg_fiber_post;sess.sess4.avg_fiber_post];

% M533 or M548
% sess_fiber_pre_early = [sess.sess1.avg_fiber_pre;sess.sess2.avg_fiber_pre;sess.sess3.avg_fiber_pre;sess.sess4.avg_fiber_pre;];% sess.sess7.avg_fiber_pre;]; %sess.sess8.avg_fiber_pre
% sess_fiber_pre_late = [sess.sess5.avg_fiber_pre;sess.sess6.avg_fiber_pre;sess.sess7.avg_fiber_pre;];% sess.sess7.avg_fiber_pre;]; %sess.sess8.avg_fiber_pre
% sess_fiber_post_early = [sess.sess1.avg_fiber_post;sess.sess2.avg_fiber_post;sess.sess3.avg_fiber_post;sess.sess4.avg_fiber_post;];%sess.sess7.avg_fiber_post;]; %sess.sess8.avg_fiber_post
% sess_fiber_post_late = [sess.sess5.avg_fiber_post;sess.sess6.avg_fiber_post;sess.sess7.avg_fiber_post;];%sess.sess7.avg_fiber_post;]; %sess.sess8.avg_fiber_post

% M534 
%sess_fiber_pre_early = [sess.sess1.avg_fiber_pre;sess.sess2.avg_fiber_pre;sess.sess3.avg_fiber_pre;];
% sess_fiber_pre_late = [sess.sess4.avg_fiber_pre;sess.sess5.avg_fiber_pre;];
%sess_fiber_post_early = [sess.sess1.avg_fiber_post;sess.sess2.avg_fiber_post;sess.sess3.avg_fiber_post;];
% sess_fiber_post_late = [sess.sess4.avg_fiber_post;sess.sess5.avg_fiber_post;];

% M547
% sess_fiber_pre_early = [sess.sess1.avg_fiber_pre;sess.sess2.avg_fiber_pre;sess.sess3.avg_fiber_pre;sess.sess4.avg_fiber_pre;];% sess.sess7.avg_fiber_pre;]; %sess.sess8.avg_fiber_pre
% sess_fiber_pre_late = [sess.sess5.avg_fiber_pre;sess.sess6.avg_fiber_pre;];% sess.sess7.avg_fiber_pre;]; %sess.sess8.avg_fiber_pre
% sess_fiber_post_early = [sess.sess1.avg_fiber_post;sess.sess2.avg_fiber_post;sess.sess3.avg_fiber_post;sess.sess4.avg_fiber_post;];%sess.sess7.avg_fiber_post;]; %sess.sess8.avg_fiber_post
% sess_fiber_post_late = [sess.sess5.avg_fiber_post;sess.sess6.avg_fiber_post;];%sess.sess7.avg_fiber_post;]; %sess.sess8.avg_fiber_post


%% average
avg_fiber_pre_early = mean(sess_fiber_pre_early,1);
avg_fiber_pre_late = mean(sess_fiber_pre_late,1);
avg_fiber_post_early = mean(sess_fiber_post_early,1);
avg_fiber_post_late = mean(sess_fiber_post_late,1);

%% 
cd 'F:\M545'
file_name = 'M545_'; 
filename = append(file_name, "avg_split.mat");
save(filename, 'avg_fiber_post_early','avg_fiber_post_late','avg_fiber_pre_early','avg_fiber_pre_early');