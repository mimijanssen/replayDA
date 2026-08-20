
%% DURATION AND FREQUENCY PLOTS
cd('F:\Duration')

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

 %M452 8-13. 
 l = 0;
 for k=10:15
    FileNames=Files(k).name;
    if k == 12 % if the first file... skip to sess 2, if the 5th file, skip to next session (two skips.) 
        l = l + 1; % skip labeling the fourth session 
    end 
    l = l + 1; 
    M452.(['sess',num2str(l)]) = load(FileNames);
 end

l = 0;
for k=16:20
   FileNames=Files(k).name;
   if k == 16 % if the first file... skip to sess 2, if the 5th file, skip to next session (two skips.) 
       l = l + 1; % skip labeling the fourth session 
   elseif k == 20
       l = l + 2; 
   end 
   l = l + 1; 
   M453.(['sess',num2str(l)]) = load(FileNames);
end

% M460
l = 0;
for k=21:24
   FileNames=Files(k).name;
   if k == 22 % if the first file... skip to sess 2, if the 5th file, skip to next session (two skips.) 
       l = l + 3; % skip labeling the fourth session 
   end 
   l = l + 1; 
   M460.(['sess',num2str(l)]) = load(FileNames);
end

% M533
l = 0;
for k=25:31
   FileNames=Files(k).name;
   l = l + 1; 
   M533.(['sess',num2str(l)]) = load(FileNames);
end

% M534
% l = 0;
% for k=26:28
%    FileNames=Files(k).name;
%    if k == 27 % if the first file... skip to sess 2, if the 5th file, skip to next session (two skips.) 
%        l = l + 1; % skip labeling the fourth session 
%    end 
%    l = l + 1; 
%    M534.(['sess',num2str(l)]) = load(FileNames);
% end

% M545
l = 3;
for k=32:35
   FileNames=Files(k).name;
   l = l + 1; 
   M545.(['sess',num2str(l)]) = load(FileNames);
end

% M547
l = 0;
for k=36:41
   FileNames=Files(k).name;
   l = l + 1; 
   M547.(['sess',num2str(l)]) = load(FileNames);
end

% M548
l = 0;
for k=42:48
   FileNames=Files(k).name;
   l = l + 1; 
   M548.(['sess',num2str(l)]) = load(FileNames);
end

%%
% session 1 
sess1_pre_freq = [M433.sess1.freq(1), M452.sess1.freq(1), M460.sess1.freq(1), M533.sess1.freq(1),  M547.sess1.freq(1), M548.sess1.freq(1)]; 
sess1_post_freq = [M433.sess1.freq(2), M452.sess1.freq(2), M460.sess1.freq(2), M533.sess1.freq(2),  M547.sess1.freq(2), M548.sess1.freq(2)];
sess1_pre_swr_count = [M433.sess1.swr_count(1), M452.sess1.swr_count(1), M460.sess1.swr_count(1), M533.sess1.swr_count(1),  M547.sess1.swr_count(1), M548.sess1.swr_count(1)];
sess1_post_swr_count = [M433.sess1.swr_count(2), M452.sess1.swr_count(2), M460.sess1.swr_count(2), M533.sess1.swr_count(2),  M547.sess1.swr_count(2), M548.sess1.swr_count(2)];

% session 2 
sess2_pre_freq = [M433.sess2.freq(1),  M452.sess2.freq(1), M453.sess2.freq(1), M533.sess2.freq(1), M547.sess2.freq(1), M548.sess2.freq(1)]; 
sess2_post_freq = [M433.sess2.freq(2), M452.sess2.freq(2), M453.sess2.freq(2), M533.sess2.freq(2), M547.sess2.freq(2), M548.sess2.freq(2)];
sess2_pre_swr_count = [M433.sess2.swr_count(1), M452.sess2.swr_count(1), M453.sess2.swr_count(1), M533.sess2.swr_count(1), M547.sess2.swr_count(1), M548.sess2.swr_count(1)];
sess2_post_swr_count = [M433.sess2.swr_count(2), M452.sess2.swr_count(2), M453.sess2.swr_count(2), M533.sess2.swr_count(2),  M547.sess2.swr_count(1), M548.sess2.swr_count(2)];

% session 3
sess3_pre_freq = [M433.sess3.freq(1), M453.sess3.freq(1), M533.sess3.freq(1), M547.sess3.freq(1), M548.sess3.freq(1)]; 
sess3_post_freq = [M433.sess3.freq(2), M453.sess3.freq(2), M533.sess3.freq(2), M547.sess3.freq(2), M548.sess3.freq(2)];
sess3_pre_swr_count = [M433.sess3.swr_count(1), M453.sess3.swr_count(1), M533.sess3.swr_count(1), M547.sess3.swr_count(1), M548.sess3.swr_count(1)];
sess3_post_swr_count = [M433.sess3.swr_count(2), M453.sess3.swr_count(2), M533.sess3.swr_count(2), M547.sess3.swr_count(2), M548.sess3.swr_count(2)];

% session 4 
sess4_pre_freq = [M453.sess4.freq(1), M452.sess4.freq(1), M533.sess4.freq(1), M545.sess4.freq(1),  M547.sess4.freq(1), M548.sess4.freq(1)]; 
sess4_post_freq = [M453.sess4.freq(2), M452.sess4.freq(2), M533.sess4.freq(2), M545.sess4.freq(2), M547.sess4.freq(2), M548.sess4.freq(2)];
sess4_pre_swr_count = [M453.sess4.swr_count(1), M452.sess4.swr_count(1), M533.sess4.swr_count(1), M545.sess4.swr_count(1), M547.sess4.swr_count(1), M548.sess4.swr_count(1)];
sess4_post_swr_count = [M453.sess4.swr_count(2), M452.sess4.swr_count(2), M533.sess4.swr_count(2), M545.sess4.swr_count(2), M547.sess4.swr_count(2), M548.sess4.swr_count(2)];

% session 5
sess5_pre_freq = [M433.sess5.freq(1), M452.sess5.freq(1), M453.sess5.freq(1), M460.sess5.freq(1), M533.sess5.freq(1), M545.sess5.freq(1), M547.sess5.freq(1), M548.sess5.freq(1)]; 
sess5_post_freq = [M433.sess5.freq(2), M452.sess5.freq(2), M453.sess5.freq(2), M460.sess5.freq(2), M533.sess5.freq(2), M545.sess5.freq(2),  M547.sess5.freq(2), M548.sess5.freq(2)];
sess5_pre_swr_count = [M433.sess5.swr_count(1), M452.sess5.swr_count(1), M453.sess5.swr_count(1), M460.sess5.swr_count(1), M533.sess5.swr_count(1), M545.sess5.swr_count(1),  M547.sess5.swr_count(1), M548.sess5.swr_count(1)];
sess5_post_swr_count = [M433.sess5.swr_count(2), M452.sess5.swr_count(2), M453.sess5.swr_count(2), M460.sess5.swr_count(2), M533.sess5.swr_count(2), M545.sess5.swr_count(2),  M547.sess5.swr_count(2), M548.sess5.swr_count(2)];

% session 6 
sess6_pre_freq = [M433.sess6.freq(1), M452.sess6.freq(1), M460.sess6.freq(1), M533.sess6.freq(1), M545.sess6.freq(1),  M547.sess6.freq(1), M548.sess6.freq(1)]; 
sess6_post_freq = [M433.sess6.freq(2), M452.sess6.freq(2), M460.sess6.freq(2), M533.sess6.freq(2), M545.sess6.freq(2), M547.sess6.freq(2), M548.sess6.freq(2)];
sess6_pre_swr_count = [M433.sess6.swr_count(1), M452.sess6.swr_count(1), M460.sess6.swr_count(1), M533.sess6.swr_count(1), M545.sess6.swr_count(1), M547.sess6.swr_count(1), M548.sess6.swr_count(1)];
sess6_post_swr_count = [M433.sess6.swr_count(2),M452.sess6.swr_count(2), M460.sess6.swr_count(2), M533.sess6.swr_count(2), M545.sess6.swr_count(2), M547.sess6.swr_count(2), M548.sess6.swr_count(2)];

% session 7
sess7_pre_freq = [M433.sess7.freq(1), M452.sess7.freq(1), M460.sess7.freq(1), M533.sess7.freq(1), M545.sess7.freq(1), M548.sess7.freq(1)]; 
sess7_post_freq = [M433.sess7.freq(2), M452.sess7.freq(2), M460.sess7.freq(2), M533.sess7.freq(2), M545.sess7.freq(2), M548.sess7.freq(2)];
sess7_pre_swr_count = [M433.sess7.swr_count(1), M452.sess7.swr_count(1), M460.sess7.swr_count(1), M533.sess7.swr_count(1), M545.sess7.swr_count(1), M548.sess7.swr_count(1)]; % rerun sess 7 pre swr_count for M460...
% I think I should either exclude M460 because lack of swrs... or just
% don't use the pre session... seems like only noise was being picke dup...M460.sess7.swr_count(1)
sess7_post_swr_count = [M433.sess7.swr_count(2), M452.sess7.swr_count(2), M460.sess7.swr_count(2), M533.sess7.swr_count(2), M545.sess7.swr_count(2), M548.sess7.swr_count(2)];

% session 8 
sess8_pre_freq = [M433.sess8.freq(1), M453.sess8.freq(1)]; 
sess8_post_freq = [M433.sess8.freq(2), M453.sess8.freq(2)];
sess8_pre_swr_count = [M433.sess8.swr_count(1), M453.sess8.swr_count(1)];
sess8_post_swr_count = [M433.sess8.swr_count(2), M453.sess8.swr_count(2)];



%% 
all_sess_count = [sess1_pre_swr_count, sess1_post_swr_count, sess2_pre_swr_count, sess2_post_swr_count, sess3_pre_swr_count, sess3_post_swr_count, sess4_pre_swr_count, sess4_post_swr_count, sess5_pre_swr_count, sess5_post_swr_count, sess6_pre_swr_count, sess6_post_swr_count, sess7_pre_swr_count, sess7_post_swr_count, sess8_pre_swr_count, sess8_post_swr_count]

%%
mean(all_sess_count)
std(all_sess_count)
