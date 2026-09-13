%% DURATION AND FREQUENCY PLOTS
% AND COUNTS 

cd('F:\sleepcounts')

%%
Files=dir('*.*');

% load files for each mouse -- this is hard coded ew.

% M433
l = 2;
for k=3:9
   FileNames=Files(k).name
   if k == 6 
       l = l - 1; % skip labeling the fourth session 
   end 
   M433.(['sess',num2str(k-l)]) = load(FileNames);
end
% check
M433

% M452
l = 0;
for k=10:15
   FileNames=Files(k).name
   if k == 12 
       l = l + 1; % skip labeling the fourth session 
   end
   l = l + 1; 
   M452.(['sess',num2str(l)]) = load(FileNames) ;
end
M452


% M453 
l = 0;
for k=16:20 %10:14
   FileNames=Files(k).name
   if k == 16 % if the first file... skip to sess 2, if the 5th file, skip to next session (two skips.) 
       l = l + 1; % skip labeling the fourth session 
   elseif k == 20
       l = l + 2; 
   end 
   l = l + 1; 
   M453.(['sess',num2str(l)]) = load(FileNames);
end
M453

% M460
l = 0;
for k=21:24 %15:18
   FileNames=Files(k).name
   if k == 22 % if the first file... skip to sess 2, if the 5th file, skip to next session (two skips.) 
       l = l + 3; % skip labeling the fourth session 
   end 
   l = l + 1; 
   M460.(['sess',num2str(l)]) = load(FileNames);
end
M460

% M533
l = 0;
for k=25:31%19:25
   FileNames=Files(k).name;
   l = l + 1; 
   M533.(['sess',num2str(l)]) = load(FileNames);
end
M533

% M534
%l = 0;
%for k=32:34%26:28
%   FileNames=Files(k).name
%   if k == 33 % if the first file... skip to sess 2, if the 5th file, skip to next session (two skips.) 
%       l = l + 1; % skip labeling the fourth session 
%   end 
%   l = l + 1; 
%   M534.(['sess',num2str(l)]) = load(FileNames);
%end

% M545
l = 3;
for k=32:35%29:32
   FileNames=Files(k).name
   l = l + 1; 
   M545.(['sess',num2str(l)]) = load(FileNames);
end
M545

% M547
l = 0;
for k=36:41 %33:38
   FileNames=Files(k).name
   l = l + 1; 
   M547.(['sess',num2str(l)]) = load(FileNames);
end
M547

% M548
l = 0;
for k=42:48
   FileNames=Files(k).name
   l = l + 1; 
   M548.(['sess',num2str(l)]) = load(FileNames);
end
M548


%% SESSION PRE AND POST MATRIX 
% NEED TO make one for count! with same structure . 

% session 1 
sess1_pre_freq = [M433.sess1.freq(1), M452.sess1.freq(1), M460.sess1.freq(1), M533.sess1.freq(1), M547.sess1.freq(1), M548.sess1.freq(1)]; 
sess1_post_freq = [M433.sess1.freq(2), M452.sess1.freq(2), M460.sess1.freq(2), M533.sess1.freq(2), M547.sess1.freq(2), M548.sess1.freq(2)];
sess1_pre_dur = [M433.sess1.dur(1), M452.sess1.dur(1), M460.sess1.dur(1), M533.sess1.dur(1), M547.sess1.dur(1), M548.sess1.dur(1)];
sess1_post_dur = [M433.sess1.dur(2), M452.sess1.dur(2), M460.sess1.dur(2), M533.sess1.dur(2), M547.sess1.dur(2), M548.sess1.dur(2)];
sess1_pre_swr_count = [M433.sess1.swr_count(1), M452.sess1.swr_count(1), M460.sess1.swr_count(1), M533.sess1.swr_count(1), M547.sess1.swr_count(1), M548.sess1.swr_count(1)];
sess1_post_swr_count = [M433.sess1.swr_count(2), M452.sess1.swr_count(2), M460.sess1.swr_count(2), M533.sess1.swr_count(2), M547.sess1.swr_count(2), M548.sess1.swr_count(2)];
sess1_pre_theta = [M433.sess1.theta(1), M452.sess1.theta(1), M460.sess1.theta(1), M533.sess1.theta(1), M547.sess1.theta(1), M548.sess1.theta(1)];
sess1_post_theta = [M433.sess1.theta(2), M452.sess1.theta(2), M460.sess1.theta(2), M533.sess1.theta(2), M547.sess1.theta(2), M548.sess1.theta(2)];
sess1_pre_delta = [M433.sess1.delta(1), M452.sess1.delta(1), M460.sess1.delta(1), M533.sess1.delta(1), M547.sess1.delta(1), M548.sess1.delta(1)];
sess1_post_delta = [M433.sess1.delta(2), M452.sess1.delta(2), M460.sess1.delta(2), M533.sess1.delta(2), M547.sess1.delta(2), M548.sess1.delta(2)];
sess1_pre_swr = [M433.sess1.swr(1), M452.sess1.swr(1), M460.sess1.swr(1), M533.sess1.swr(1), M547.sess1.swr(1), M548.sess1.swr(1)];
sess1_post_swr = [M433.sess1.swr(2), M452.sess1.swr(2), M460.sess1.swr(2), M533.sess1.swr(2), M547.sess1.swr(2), M548.sess1.swr(2)];

% session 2 
sess2_pre_freq = [M433.sess2.freq(1), M452.sess2.freq(1), M453.sess2.freq(1), M533.sess2.freq(1), M547.sess2.freq(1), M548.sess2.freq(1)]; 
sess2_post_freq = [M433.sess2.freq(2), M452.sess2.freq(2), M453.sess2.freq(2), M533.sess2.freq(2), M547.sess2.freq(2), M548.sess2.freq(2)];
sess2_pre_dur = [M433.sess2.dur(1), M452.sess2.dur(1), M453.sess2.dur(1), M533.sess2.dur(1), M547.sess2.dur(1), M548.sess2.dur(1)];
sess2_post_dur = [M433.sess2.dur(2), M452.sess2.dur(2), M453.sess2.dur(2), M533.sess2.dur(2),  M547.sess2.dur(1), M548.sess2.dur(2)];
sess2_pre_swr_count = [M433.sess2.swr_count(1), M452.sess2.swr_count(1), M453.sess2.swr_count(1), M533.sess2.swr_count(1), M547.sess2.swr_count(1), M548.sess2.swr_count(1)];
sess2_post_swr_count = [M433.sess2.swr_count(2), M452.sess2.swr_count(2), M453.sess2.swr_count(2), M533.sess2.swr_count(2),  M547.sess2.swr_count(1), M548.sess2.swr_count(2)];
sess2_pre_theta = [M433.sess2.theta(1), M452.sess2.theta(1), M453.sess2.theta(1), M533.sess2.theta(1), M547.sess2.theta(1), M548.sess2.theta(1)];
sess2_post_theta = [M433.sess2.theta(2), M452.sess2.theta(2), M453.sess2.theta(2), M533.sess2.theta(2),  M547.sess2.theta(1), M548.sess2.theta(2)];
sess2_pre_delta = [M433.sess2.delta(1), M452.sess2.delta(1), M453.sess2.delta(1), M533.sess2.delta(1), M547.sess2.delta(1), M548.sess2.delta(1)];
sess2_post_delta = [M433.sess2.delta(2), M452.sess2.delta(2), M453.sess2.delta(2), M533.sess2.delta(2),  M547.sess2.delta(1), M548.sess2.delta(2)];
sess2_pre_swr = [M433.sess2.swr(1), M452.sess2.swr(1), M453.sess2.swr(1), M533.sess2.swr(1), M547.sess2.swr(1), M548.sess2.swr(1)];
sess2_post_swr = [M433.sess2.swr(2), M452.sess2.swr(2), M453.sess2.swr(2), M533.sess2.swr(2),  M547.sess2.swr(1), M548.sess2.swr(2)];

% session 3
sess3_pre_freq = [M433.sess3.freq(1), M453.sess3.freq(1), M533.sess3.freq(1), M547.sess3.freq(1), M548.sess3.freq(1)]; 
sess3_post_freq = [M433.sess3.freq(2), M453.sess3.freq(2), M533.sess3.freq(2), M547.sess3.freq(2), M548.sess3.freq(2)];
sess3_pre_dur = [M433.sess3.dur(1), M453.sess3.dur(1), M533.sess3.dur(1), M547.sess3.dur(1), M548.sess3.dur(1)];
sess3_post_dur = [M433.sess3.dur(2), M453.sess3.dur(2), M533.sess3.dur(2), M547.sess3.dur(2), M548.sess3.dur(2)];
sess3_pre_swr_count = [M433.sess3.swr_count(1), M453.sess3.swr_count(1), M533.sess3.swr_count(1), M547.sess3.swr_count(1), M548.sess3.swr_count(1)];
sess3_post_swr_count = [M433.sess3.swr_count(2), M453.sess3.swr_count(2), M533.sess3.swr_count(2), M547.sess3.swr_count(2), M548.sess3.swr_count(2)];
sess3_pre_theta = [M433.sess3.theta(1), M453.sess3.theta(1), M533.sess3.theta(1), M547.sess3.theta(1), M548.sess3.theta(1)];
sess3_post_theta = [M433.sess3.theta(2), M453.sess3.theta(2), M533.sess3.theta(2), M547.sess3.theta(2), M548.sess3.theta(2)];
sess3_pre_delta = [M433.sess3.delta(1), M453.sess3.delta(1), M533.sess3.delta(1), M547.sess3.delta(1), M548.sess3.delta(1)];
sess3_post_delta = [M433.sess3.delta(2), M453.sess3.delta(2), M533.sess3.delta(2), M547.sess3.delta(2), M548.sess3.delta(2)];
sess3_pre_swr = [M433.sess3.swr(1), M453.sess3.swr(1), M533.sess3.swr(1), M547.sess3.swr(1), M548.sess3.swr(1)];
sess3_post_swr = [M433.sess3.swr(2), M453.sess3.swr(2), M533.sess3.swr(2), M547.sess3.swr(2), M548.sess3.swr(2)];

% session 4 
sess4_pre_freq = [M453.sess4.freq(1), M452.sess4.freq(1), M533.sess4.freq(1), M545.sess4.freq(1),  M547.sess4.freq(1), M548.sess4.freq(1)]; 
sess4_post_freq = [M453.sess4.freq(2), M452.sess4.freq(2), M533.sess4.freq(2), M545.sess4.freq(2), M547.sess4.freq(2), M548.sess4.freq(2)];
sess4_pre_dur = [M453.sess4.dur(1), M452.sess4.dur(1), M533.sess4.dur(1), M545.sess4.dur(1), M547.sess4.dur(1), M548.sess4.dur(1)];
sess4_post_dur = [M453.sess4.dur(2), M452.sess4.dur(2), M533.sess4.dur(2), M545.sess4.dur(2), M547.sess4.dur(2), M548.sess4.dur(2)];
sess4_pre_swr_count = [M453.sess4.swr_count(1), M452.sess4.swr_count(1), M533.sess4.swr_count(1), M545.sess4.swr_count(1), M547.sess4.swr_count(1), M548.sess4.swr_count(1)];
sess4_post_swr_count = [M453.sess4.swr_count(2), M452.sess4.swr_count(2), M533.sess4.swr_count(2), M545.sess4.swr_count(2), M547.sess4.swr_count(2), M548.sess4.swr_count(2)];
sess4_pre_theta = [M453.sess4.theta(1), M452.sess4.theta(1), M533.sess4.theta(1), M545.sess4.theta(1), M547.sess4.theta(1), M548.sess4.theta(1)];
sess4_post_theta = [M453.sess4.theta(2), M452.sess4.theta(2), M533.sess4.theta(2), M545.sess4.theta(2), M547.sess4.theta(2), M548.sess4.theta(2)];
sess4_pre_delta = [M453.sess4.delta(1), M452.sess4.delta(1), M533.sess4.delta(1), M545.sess4.delta(1), M547.sess4.delta(1), M548.sess4.delta(1)];
sess4_post_delta = [M453.sess4.delta(2), M452.sess4.delta(2), M533.sess4.delta(2), M545.sess4.delta(2), M547.sess4.delta(2), M548.sess4.delta(2)];
sess4_pre_swr = [M453.sess4.swr(1), M452.sess4.swr(1), M533.sess4.swr(1), M545.sess4.swr(1), M547.sess4.swr(1), M548.sess4.swr(1)];
sess4_post_swr = [M453.sess4.swr(2), M452.sess4.swr(2), M533.sess4.swr(2), M545.sess4.swr(2), M547.sess4.swr(2), M548.sess4.swr(2)];

% session 5
sess5_pre_freq = [M433.sess5.freq(1), M452.sess5.freq(1), M453.sess5.freq(1), M460.sess5.freq(1), M533.sess5.freq(1), M545.sess5.freq(1), M547.sess5.freq(1), M548.sess5.freq(1)]; 
sess5_post_freq = [M433.sess5.freq(2), M452.sess5.freq(2), M453.sess5.freq(2), M460.sess5.freq(2), M533.sess5.freq(2), M545.sess5.freq(2),  M547.sess5.freq(2), M548.sess5.freq(2)];
sess5_pre_dur = [M433.sess5.dur(1), M452.sess5.dur(1), M453.sess5.dur(1), M460.sess5.dur(1), M533.sess5.dur(1), M545.sess5.dur(1),  M547.sess5.dur(1), M548.sess5.dur(1)];
sess5_post_dur = [M433.sess5.dur(2), M452.sess5.dur(2), M453.sess5.dur(2), M460.sess5.dur(2), M533.sess5.dur(2), M545.sess5.dur(2),  M547.sess5.dur(2), M548.sess5.dur(2)];
sess5_pre_swr_count = [M433.sess5.swr_count(1), M452.sess5.swr_count(1), M453.sess5.swr_count(1), M460.sess5.swr_count(1), M533.sess5.swr_count(1), M545.sess5.swr_count(1),  M547.sess5.swr_count(1), M548.sess5.swr_count(1)];
sess5_post_swr_count = [M433.sess5.swr_count(2), M452.sess5.swr_count(2), M453.sess5.swr_count(2), M460.sess5.swr_count(2), M533.sess5.swr_count(2), M545.sess5.swr_count(2),  M547.sess5.swr_count(2), M548.sess5.swr_count(2)];
sess5_pre_theta = [M433.sess5.theta(1), M452.sess5.theta(1), M453.sess5.theta(1), M460.sess5.theta(1), M533.sess5.theta(1), M545.sess5.theta(1),  M547.sess5.theta(1), M548.sess5.theta(1)];
sess5_post_theta = [M433.sess5.theta(2), M452.sess5.theta(2), M453.sess5.theta(2), M460.sess5.theta(2), M533.sess5.theta(2), M545.sess5.theta(2),  M547.sess5.theta(2), M548.sess5.theta(2)];
sess5_pre_delta = [M433.sess5.delta(1), M452.sess5.delta(1), M453.sess5.delta(1), M460.sess5.delta(1), M533.sess5.delta(1), M545.sess5.delta(1),  M547.sess5.delta(1), M548.sess5.delta(1)];
sess5_post_delta = [M433.sess5.delta(2), M452.sess5.delta(2), M453.sess5.delta(2), M460.sess5.delta(2), M533.sess5.delta(2), M545.sess5.delta(2),  M547.sess5.delta(2), M548.sess5.delta(2)];
sess5_pre_swr = [M433.sess5.swr(1), M452.sess5.swr(1), M453.sess5.swr(1), M460.sess5.swr(1), M533.sess5.swr(1), M545.sess5.swr(1),  M547.sess5.swr(1), M548.sess5.swr(1)];
sess5_post_swr = [M433.sess5.swr(2), M452.sess5.swr(2), M453.sess5.swr(2), M460.sess5.swr(2), M533.sess5.swr(2), M545.sess5.swr(2),  M547.sess5.swr(2), M548.sess5.swr(2)];

% session 6 
sess6_pre_freq = [M433.sess6.freq(1), M452.sess6.freq(1), M460.sess6.freq(1), M533.sess6.freq(1), M545.sess6.freq(1),  M547.sess6.freq(1), M548.sess6.freq(1)]; 
sess6_post_freq = [M433.sess6.freq(2), M452.sess6.freq(2), M460.sess6.freq(2), M533.sess6.freq(2), M545.sess6.freq(2), M547.sess6.freq(2), M548.sess6.freq(2)];
sess6_pre_dur = [M433.sess6.dur(1), M452.sess6.dur(1), M460.sess6.dur(1), M533.sess6.dur(1), M545.sess6.dur(1), M547.sess6.dur(1), M548.sess6.dur(1)];
sess6_post_dur = [M433.sess6.dur(2), M452.sess6.dur(2), M460.sess6.dur(2), M533.sess6.dur(2), M545.sess6.dur(2), M547.sess6.dur(2), M548.sess6.dur(2)];
sess6_pre_swr_count = [M433.sess6.swr_count(1), M452.sess6.swr_count(1), M460.sess6.swr_count(1), M533.sess6.swr_count(1), M545.sess6.swr_count(1), M547.sess6.swr_count(1), M548.sess6.swr_count(1)];
sess6_post_swr_count = [M433.sess6.swr_count(2), M452.sess6.swr_count(2), M460.sess6.swr_count(2), M533.sess6.swr_count(2), M545.sess6.swr_count(2), M547.sess6.swr_count(2), M548.sess6.swr_count(2)];
sess6_pre_theta = [M433.sess6.theta(1), M452.sess6.theta(1), M460.sess6.theta(1), M533.sess6.theta(1), M545.sess6.theta(1), M547.sess6.theta(1), M548.sess6.theta(1)];
sess6_post_theta = [M433.sess6.theta(2), M452.sess6.theta(2), M460.sess6.theta(2), M533.sess6.theta(2), M545.sess6.theta(2), M547.sess6.theta(2), M548.sess6.theta(2)];
sess6_pre_delta = [M433.sess6.delta(1), M452.sess6.delta(1), M460.sess6.delta(1), M533.sess6.delta(1), M545.sess6.delta(1), M547.sess6.delta(1), M548.sess6.delta(1)];
sess6_post_delta = [M433.sess6.delta(2), M452.sess6.delta(2), M460.sess6.delta(2), M533.sess6.delta(2), M545.sess6.delta(2), M547.sess6.delta(2), M548.sess6.delta(2)];
sess6_pre_swr = [M433.sess6.swr(1), M452.sess6.swr(1), M460.sess6.swr(1), M533.sess6.swr(1), M545.sess6.swr(1), M547.sess6.swr(1), M548.sess6.swr(1)];
sess6_post_swr = [M433.sess6.swr(2), M452.sess6.swr(2), M460.sess6.swr(2), M533.sess6.swr(2), M545.sess6.swr(2), M547.sess6.swr(2), M548.sess6.swr(2)];

% session 7
sess7_pre_freq = [M433.sess7.freq(1), M452.sess7.freq(1), M460.sess7.freq(1), M533.sess7.freq(1), M545.sess7.freq(1), M548.sess7.freq(1)]; 
sess7_post_freq = [M433.sess7.freq(2), M452.sess7.freq(2), M460.sess7.freq(2), M533.sess7.freq(2), M545.sess7.freq(2), M548.sess7.freq(2)];
sess7_pre_dur = [M433.sess7.dur(1), M452.sess7.dur(1), M460.sess7.dur(1), M533.sess7.dur(1), M545.sess7.dur(1), M548.sess7.dur(1)]; 
sess7_post_dur = [M433.sess7.dur(2), M452.sess7.dur(2), M460.sess7.dur(2), M533.sess7.dur(2), M545.sess7.dur(2), M548.sess7.dur(2)];
sess7_pre_swr_count = [M433.sess7.swr_count(1), M452.sess7.swr_count(1), M460.sess7.swr_count(1), M533.sess7.swr_count(1), M545.sess7.swr_count(1), M548.sess7.swr_count(1)]; 
sess7_post_swr_count = [M433.sess7.swr_count(2), M452.sess7.swr_count(2), M460.sess7.swr_count(2), M533.sess7.swr_count(2), M545.sess7.swr_count(2), M548.sess7.swr_count(2)];
sess7_pre_theta = [M433.sess7.theta(1), M452.sess7.theta(1), M460.sess7.theta(1), M533.sess7.theta(1), M545.sess7.theta(1), M548.sess7.theta(1)]; 
sess7_post_theta = [M433.sess7.theta(2), M452.sess7.theta(2), M460.sess7.theta(2), M533.sess7.theta(2), M545.sess7.theta(2), M548.sess7.theta(2)];
sess7_pre_delta = [M433.sess7.delta(1), M452.sess7.delta(1), M460.sess7.delta(1), M533.sess7.delta(1), M545.sess7.delta(1), M548.sess7.delta(1)]; 
sess7_post_delta = [M433.sess7.delta(2), M452.sess7.delta(2), M460.sess7.delta(2), M533.sess7.delta(2), M545.sess7.delta(2), M548.sess7.delta(2)];
sess7_pre_swr = [M433.sess7.swr(1), M452.sess7.swr(1), M460.sess7.swr(1), M533.sess7.swr(1), M545.sess7.swr(1), M548.sess7.swr(1)]; 
sess7_post_swr = [M433.sess7.swr(2), M452.sess7.swr(2), M460.sess7.swr(2), M533.sess7.swr(2), M545.sess7.swr(2), M548.sess7.swr(2)];

% session 8 
sess8_pre_freq = [M433.sess8.freq(1), M453.sess8.freq(1)]; 
sess8_post_freq = [M433.sess8.freq(2), M453.sess8.freq(2)];
sess8_pre_dur = [M433.sess8.dur(1), M453.sess8.dur(1)];
sess8_post_dur = [M433.sess8.dur(2), M453.sess8.dur(2)];
sess8_pre_swr_count = [M433.sess8.swr_count(1), M453.sess8.swr_count(1)];
sess8_post_swr_count = [M433.sess8.swr_count(2), M453.sess8.swr_count(2)];
sess8_pre_theta = [M433.sess8.theta(1), M453.sess8.theta(1)];
sess8_post_theta = [M433.sess8.theta(2), M453.sess8.theta(2)];
sess8_pre_delta = [M433.sess8.delta(1), M453.sess8.delta(1)];
sess8_post_delta = [M433.sess8.delta(2), M453.sess8.delta(2)];
sess8_pre_swr = [M433.sess8.swr(1), M453.sess8.swr(1)];
sess8_post_swr = [M433.sess8.swr(2), M453.sess8.swr(2)];

%%
%cd('Desktop') 
% Define colors
pre_color = [212, 162, 43]./255; %[0.4, 0.6, 0.8]; % Example color for "pre" BoxFaceColor
%post_color = [0.2, 0.4, 0.7]; % Example color for "post" BoxFaceColor
post_color = [44,123,116]./255; %[0.3, 0.5, 0.8];%[24,53,58]./255; %[0.5, 0.6, 0.9];
pre_scatter_color = [225, 190,106]./255; %[0.6, 0.8, 0.9]; % Example color for "pre" scatter points
post_scatter_color =  [64,176,166]./255; %[0.3, 0.5, 0.8]; % Example color for "post" scatter points

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

swr_count_data = {
    sess1_pre_swr_count, sess1_post_swr_count;
    sess2_pre_swr_count, sess2_post_swr_count;
    sess3_pre_swr_count, sess3_post_swr_count;
    sess4_pre_swr_count, sess4_post_swr_count;
    sess5_pre_swr_count, sess5_post_swr_count;
    sess6_pre_swr_count, sess6_post_swr_count;
    sess7_pre_swr_count, sess7_post_swr_count;
    sess8_pre_swr_count, sess8_post_swr_count
};

delta_data = {
    sess1_pre_delta, sess1_post_delta;
    sess2_pre_delta, sess2_post_delta;
    sess3_pre_delta, sess3_post_delta;
    sess4_pre_delta, sess4_post_delta;
    sess5_pre_delta, sess5_post_delta;
    sess6_pre_delta, sess6_post_delta;
    sess7_pre_delta, sess7_post_delta;
    sess8_pre_delta, sess8_post_delta
};

theta_data = {
    sess1_pre_theta, sess1_post_theta;
    sess2_pre_theta, sess2_post_theta;
    sess3_pre_theta, sess3_post_theta;
    sess4_pre_theta, sess4_post_theta;
    sess5_pre_theta, sess5_post_theta;
    sess6_pre_theta, sess6_post_theta;
    sess7_pre_theta, sess7_post_theta;
    sess8_pre_theta, sess8_post_theta
};

swr_data = {
    sess1_pre_swr, sess1_post_swr;
    sess2_pre_swr, sess2_post_swr;
    sess3_pre_swr, sess3_post_swr;
    sess4_pre_swr, sess4_post_swr;
    sess5_pre_swr, sess5_post_swr;
    sess6_pre_swr, sess6_post_swr;
    sess7_pre_swr, sess7_post_swr;
    sess8_pre_swr, sess8_post_swr
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
ylabel('SWR Rate (SWRs per minute)');
xticks(1.5:2:2*num_sessions);
xticklabels(compose('%d', 1:num_sessions));
set(gca,'fontsize', 18)
legend('wake','nrem','Location','best')
set(gcf,'Color',[1,1,1])
set(gca,'fontsize', 24)
shg
hold off

%exportgraphics(gcf,'SWRfreq_sleep.png','BackgroundColor','none','ContentType','vector');

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
ylabel('SWR Duration (ms)');
xticks(1.5:2:2*num_sessions);
xticklabels(compose('%d', 1:num_sessions));
set(gca,'fontsize', 18)
legend('wake','nrem','Location','best')

set(gcf,'Color',[1,1,1])
set(gca,'fontsize', 24)
shg
hold off

%exportgraphics(gcf,'SWRduration.png','BackgroundColor','none','ContentType','vector');

% swr_countation Boxchart with scatter
figure;
hold on;
for i = 1:num_sessions
    % Combine pre and post data
    session_pre = swr_count_data{i, 1};
    session_post = swr_count_data{i, 2};
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
ylabel('# of SWRs');
xticks(1.5:2:2*num_sessions);
xticklabels(compose('%d', 1:num_sessions));
set(gca,'fontsize', 18)
legend('wake','nrem','Location','best')

set(gcf,'Color',[1,1,1])
set(gca,'fontsize', 24)
shg
hold off



% thetaation Boxchart with scatter
figure ;
hold on;
for i = 1:num_sessions
    % Combine pre and post data
    session_pre = theta_data{i, 1};
    session_post = theta_data{i, 2};
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
ylabel('Theta Power');
xticks(1.5:2:2*num_sessions);
xticklabels(compose('%d', 1:num_sessions));
set(gca,'fontsize', 18)
legend('wake','nrem','Location','best')


% deltaation Boxchart with scatter
figure;
hold on;
for i = 1:num_sessions
    % Combine pre and post data
    session_pre = theta_data{i,1}./delta_data{i, 1};
    session_post = theta_data{i,2}./delta_data{i, 2};
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
ylabel('Theta/Delta Power');
ylim([-10, 10])
xticks(1.5:2:2*num_sessions);
xticklabels(compose('%d', 1:num_sessions));
set(gca,'fontsize', 18)
legend('wake','nrem','Location','best')

%exportgraphics(gcf,'SWRcount.png','BackgroundColor','none','ContentType','vector');


% deltaation Boxchart with scatter
figure;
hold on;
for i = 1:num_sessions
    % Combine pre and post data
    session_pre = delta_data{i, 1};
    session_post = delta_data{i, 2};
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
ylabel('Delta Power (ms)');
xticks(1.5:2:2*num_sessions);
xticklabels(compose('%d', 1:num_sessions));
set(gca,'fontsize', 18)
legend('wake','nrem','Location','best')


figure;
hold on;
for i = 1:num_sessions
    % Combine pre and post data
    session_pre = swr_data{i, 1};
    session_post = swr_data{i, 2};
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
ylabel('SWR Band Power (ms)');
xticks(1.5:2:2*num_sessions);
xticklabels(compose('%d', 1:num_sessions));
set(gca,'fontsize', 18)
legend('wake','nrem','Location','best')


%% Kruskal Wallis test 
% to determine if the frequency or duration is different in any of the
% sessions

% KW test is nonparametric. You can use ANOVA or t-tests if you suspect
% data is normal. To test normality use a KS test. 

% curate data format so that each column is a new session... 
% ~ Pre FREQ ~
freq_data_col = [sess1_pre_freq',
    sess2_pre_freq', 
    sess3_pre_freq', 
    sess4_pre_freq',
    sess5_pre_freq',
    sess6_pre_freq', 
    sess7_pre_freq',
    sess8_pre_freq', 
];

h = kstest(freq_data_col) % data does not come from a normal distribution 

% grouping variable 
freq_data_group = [ones(length(sess1_pre_freq),1); ones(length(sess2_pre_freq),1)*2; ones(length(sess3_pre_freq),1)*3; ones(length(sess4_pre_freq),1)*4; ones(length(sess5_pre_freq),1)*5; ones(length(sess6_pre_freq),1)*6; ones(length(sess7_pre_freq),1)*7; ones(length(sess8_pre_freq),1)*8];

% Anova
disp('Pre SWR rate ANOVA')
[p,t,stats] = anova1(freq_data_col, freq_data_group);

% Kruskalwallis Test 
disp('Pre SWR rate Kruskalwallis test')
[p,tbl,stats] = kruskalwallis(freq_data_col, freq_data_group,'off')
% The returned value of p indicates that the test rejects the null hypothesis at the 1% significance level. 

% Dunn's test using multicompare
disp('Pre Dunn test')
[results,m,h,gnames] = multcompare(stats,"CriticalValueType","dunn-sidak");

% no groups have a mean significantly different from each other. 
tbl = array2table(results,"VariableNames", ...
    ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

% ~ Pre Dur ~ 
dur_data_col = [sess1_pre_dur',
    sess2_pre_dur', 
    sess3_pre_dur', 
    sess4_pre_dur',
    sess5_pre_dur',
    sess6_pre_dur', 
    sess7_pre_dur',
    sess8_pre_dur', 
];

% grouping variable 
dur_data_group = [ones(length(sess1_pre_dur),1); ones(length(sess2_pre_dur),1)*2; ones(length(sess3_pre_dur),1)*3; ones(length(sess4_pre_dur),1)*4; ones(length(sess5_pre_dur),1)*5; ones(length(sess6_pre_dur),1)*6; ones(length(sess7_pre_dur),1)*7; ones(length(sess8_pre_dur),1)*8];

[p,t,stats] = anova1(dur_data_col, dur_data_group);

disp('Pre SWR duration Kruskalwallis test')
[p,tbl,stats] = kruskalwallis(dur_data_col, dur_data_group,'off')

disp('Pre Dunn test for duration')
[results,m,h,gnames] = multcompare(stats,"CriticalValueType","dunn-sidak");
% no groups have mean ranks significantly different from each other. 
% no groups have a mean significantly different from each other. 
tbl = array2table(results,"VariableNames", ...
    ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

% ~ Pre swr_count ~ 
swr_count_data_col = [sess1_pre_swr_count',
    sess2_pre_swr_count', 
    sess3_pre_swr_count', 
    sess4_pre_swr_count',
    sess5_pre_swr_count',
    sess6_pre_swr_count', 
    sess7_pre_swr_count',
    sess8_pre_swr_count', 
];

% grouping variable 
swr_count_data_group = [ones(length(sess1_pre_swr_count),1); ones(length(sess2_pre_swr_count),1)*2; ones(length(sess3_pre_swr_count),1)*3; ones(length(sess4_pre_swr_count),1)*4; ones(length(sess5_pre_swr_count),1)*5; ones(length(sess6_pre_swr_count),1)*6; ones(length(sess7_pre_swr_count),1)*7; ones(length(sess8_pre_swr_count),1)*8];

[p,t,stats] = anova1(swr_count_data_col, swr_count_data_group);

disp('Pre SWR swr_countation Kruskalwallis test')
[p,tbl,stats] = kruskalwallis(swr_count_data_col, swr_count_data_group,'off')

disp('Pre Dunn test for swr_countation')
[results,m,h,gnames] = multcompare(stats,"CriticalValueType","dunn-sidak")
% no groups have mean ranks significantly different from each other. 
tbl = array2table(results,"VariableNames", ...
    ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

% ~ Post FREQ ~
freq_data_col_post = [sess1_post_freq',
    sess2_post_freq', 
    sess3_post_freq', 
    sess4_post_freq',
    sess5_post_freq',
    sess6_post_freq', 
    sess7_post_freq',
    sess8_post_freq', 
];

% grouping variable 
freq_data_group_post = [ones(length(sess1_post_freq),1); ones(length(sess2_post_freq),1)*2; ones(length(sess3_post_freq),1)*3; ones(length(sess4_post_freq),1)*4; ones(length(sess5_post_freq),1)*5; ones(length(sess6_post_freq),1)*6; ones(length(sess7_post_freq),1)*7; ones(length(sess8_post_freq),1)*8];

[p,t,stats] = anova1(freq_data_col_post, freq_data_group_post);

disp('Post SWR rate Kruskalwallis test')
[p,tbl,stats] = kruskalwallis(freq_data_col_post, freq_data_group_post,'off')

disp('Post Dunn test freq')
[results,m,h,gnames] = multcompare(stats,"CriticalValueType","dunn-sidak")

% no groups have mean ranks significantly different from each other. 
tbl = array2table(results,"VariableNames", ...
    ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

% ~ Post Dur ~ 
dur_data_col_post = [sess1_post_dur',
    sess2_post_dur', 
    sess3_post_dur', 
    sess4_post_dur',
    sess5_post_dur',
    sess6_post_dur', 
    sess7_post_dur',
    sess8_post_dur', 
];

% grouping variable 
dur_data_group_post = [ones(length(sess1_post_dur),1); ones(length(sess2_post_dur),1)*2; ones(length(sess3_post_dur),1)*3; ones(length(sess4_post_dur),1)*4; ones(length(sess5_post_dur),1)*5; ones(length(sess6_post_dur),1)*6; ones(length(sess7_post_dur),1)*7; ones(length(sess8_post_dur),1)*8];

[p,t,stats] = anova1(dur_data_col_post, dur_data_group_post);

disp('Post SWR duration Kruskalwallis test')
[p,tbl,stats] = kruskalwallis(dur_data_col_post, dur_data_group_post,'off')

% if anything is significant, run a post hoc Dunn test. 
disp('Post Dunn test duration')
[results,m,h,gnames] = multcompare(stats,"CriticalValueType","dunn-sidak")
% no groups have mean ranks significantly different from each other. 
tbl = array2table(results,"VariableNames", ...
    ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

% ~ Post swr_count ~ 
swr_count_data_col_post = [sess1_post_swr_count',
    sess2_post_swr_count', 
    sess3_post_swr_count', 
    sess4_post_swr_count',
    sess5_post_swr_count',
    sess6_post_swr_count', 
    sess7_post_swr_count',
    sess8_post_swr_count', 
];

% grouping variable 
swr_count_data_group_post = [ones(length(sess1_post_swr_count),1); ones(length(sess2_post_swr_count),1)*2; ones(length(sess3_post_swr_count),1)*3; ones(length(sess4_post_swr_count),1)*4; ones(length(sess5_post_swr_count),1)*5; ones(length(sess6_post_swr_count),1)*6; ones(length(sess7_post_swr_count),1)*7; ones(length(sess8_post_swr_count),1)*8];

[p,t,stats] = anova1(swr_count_data_col_post, swr_count_data_group_post);

disp('Post SWR swr_countation Kruskalwallis test')
[p,tbl,stats] = kruskalwallis(swr_count_data_col_post, swr_count_data_group_post,'off')

% if anything is significant, run a post hoc Dunn test. 
disp('Post Dunn test swr_countation')
[results,m,h,gnames] = multcompare(stats,"CriticalValueType","dunn-sidak")
% no groups have mean ranks significantly different from each other. 

tbl = array2table(results,"VariableNames", ...
    ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

%% Descriptive stats 
% Average for M433 
M433_freq.pre = mean([M433.sess1.freq; M433.sess2.freq;  M433.sess3.freq; M433.sess5.freq; M433.sess6.freq;  M433.sess7.freq;  M433.sess8.freq]); 
M433_dur.pre = mean([M433.sess1.dur; M433.sess2.dur;  M433.sess3.dur; M433.sess5.dur; M433.sess6.dur;  M433.sess7.dur;  M433.sess8.dur]); 
M433_swr_count.pre = mean([M433.sess1.swr_count; M433.sess2.swr_count;  M433.sess3.swr_count; M433.sess5.swr_count; M433.sess6.swr_count;  M433.sess7.swr_count;  M433.sess8.swr_count]); 
M433_theta.pre = mean([M433.sess1.theta; M433.sess2.theta;  M433.sess3.theta; M433.sess5.theta; M433.sess6.theta;  M433.sess7.theta;  M433.sess8.theta]); 
M433_delta.pre = mean([M433.sess1.delta; M433.sess2.delta;  M433.sess3.delta; M433.sess5.delta; M433.sess6.delta;  M433.sess7.delta;  M433.sess8.delta]); 
M433_swr.pre = mean([M433.sess1.swr; M433.sess2.swr;  M433.sess3.swr; M433.sess5.swr; M433.sess6.swr;  M433.sess7.swr;  M433.sess8.swr]); 
M433_theta_delta.pre = mean([M433.sess1.theta; M433.sess2.theta;  M433.sess3.theta; M433.sess5.theta; M433.sess6.theta;  M433.sess7.theta;  M433.sess8.theta]./[M433.sess1.delta; M433.sess2.delta;  M433.sess3.delta; M433.sess5.delta; M433.sess6.delta;  M433.sess7.delta;  M433.sess8.delta]);


M452_freq.pre = mean([M452.sess1.freq; M452.sess2.freq;  M452.sess4.freq; M452.sess5.freq; M452.sess6.freq;  M452.sess7.freq;]); 
M452_dur.pre = mean([M452.sess1.dur; M452.sess2.dur;  M452.sess4.dur; M452.sess5.dur; M452.sess6.dur;  M452.sess7.dur;]); 
M452_swr_count.pre = mean([M452.sess1.swr_count; M452.sess2.swr_count;  M452.sess4.swr_count; M452.sess5.swr_count; M452.sess6.swr_count;  M452.sess7.swr_count;]); 
M452_theta.pre = mean([M452.sess1.theta; M452.sess2.theta;  M452.sess4.theta; M452.sess5.theta; M452.sess6.theta;  M452.sess7.theta;]); 
M452_delta.pre = mean([M452.sess1.delta; M452.sess2.delta;  M452.sess4.delta; M452.sess5.delta; M452.sess6.delta;  M452.sess7.delta;]); 
M452_swr.pre = mean([M452.sess1.swr; M452.sess2.swr;  M452.sess4.swr; M452.sess5.swr; M452.sess6.swr;  M452.sess7.swr;]); 
M452_theta_delta.pre = mean([M452.sess1.theta; M452.sess2.theta;  M452.sess4.theta; M452.sess5.theta; M452.sess6.theta;  M452.sess7.theta;]./[M452.sess1.delta; M452.sess2.delta;  M452.sess4.delta; M452.sess5.delta; M452.sess6.delta;  M452.sess7.delta;]);


M453_freq.pre = mean([M453.sess2.freq;  M453.sess3.freq;  M453.sess4.freq; M453.sess5.freq; M453.sess8.freq]); 
M453_dur.pre = mean([M453.sess2.dur;  M453.sess3.dur;  M453.sess4.dur; M453.sess5.dur; M453.sess8.dur]); 
M453_swr_count.pre = mean([M453.sess2.swr_count;  M453.sess3.swr_count;  M453.sess4.swr_count; M453.sess5.swr_count; M453.sess8.swr_count]); 
M453_theta.pre = mean([M453.sess2.theta;  M453.sess3.theta;  M453.sess4.theta; M453.sess5.theta; M453.sess8.theta]); 
M453_delta.pre = mean([M453.sess2.delta;  M453.sess3.delta;  M453.sess4.delta; M453.sess5.delta; M453.sess8.delta]); 
M453_swr.pre = mean([M453.sess2.swr;  M453.sess3.swr;  M453.sess4.swr; M453.sess5.swr; M453.sess8.swr]); 
M453_theta_delta.pre = mean([M453.sess2.theta;  M453.sess3.theta;  M453.sess4.theta; M453.sess5.theta; M453.sess8.theta]./[M453.sess2.delta;  M453.sess3.delta;  M453.sess4.delta; M453.sess5.delta; M453.sess8.delta]);

M460_freq.pre = mean([M460.sess1.freq;  M460.sess5.freq;  M460.sess6.freq; M460.sess7.freq]); 
M460_dur.pre = mean([M460.sess1.dur;  M460.sess5.dur;  M460.sess6.dur; M460.sess7.dur]); 
M460_swr_count.pre = mean([M460.sess1.swr_count;  M460.sess5.swr_count;  M460.sess6.swr_count; M460.sess7.swr_count]); 
M460_theta.pre = mean([M460.sess1.theta;  M460.sess5.theta;  M460.sess6.theta; M460.sess7.theta]); 
M460_delta.pre = mean([M460.sess1.delta;  M460.sess5.delta;  M460.sess6.delta; M460.sess7.delta]); 
M460_swr.pre = mean([M460.sess1.swr;  M460.sess5.swr;  M460.sess6.swr; M460.sess7.swr]); 
M460_theta_delta.pre = mean([M460.sess1.theta;  M460.sess5.theta;  M460.sess6.theta; M460.sess7.theta]./[M460.sess1.delta;  M460.sess5.delta;  M460.sess6.delta; M460.sess7.delta]);

M533_freq.pre = mean([M533.sess1.freq; M533.sess2.freq;  M533.sess3.freq;  M533.sess4.freq; M533.sess5.freq; M533.sess6.freq;  M533.sess7.freq]); 
M533_dur.pre = mean([M533.sess1.dur; M533.sess2.dur;  M533.sess3.dur;  M533.sess4.dur; M533.sess5.dur; M533.sess6.dur;  M533.sess7.dur]); 
M533_swr_count.pre = mean([M533.sess1.swr_count; M533.sess2.swr_count;  M533.sess3.swr_count;  M533.sess4.swr_count; M533.sess5.swr_count; M533.sess6.swr_count;  M533.sess7.swr_count]); 
M533_theta.pre = mean([M533.sess1.theta; M533.sess2.theta;  M533.sess3.theta;  M533.sess4.theta; M533.sess5.theta; M533.sess6.theta;  M533.sess7.theta]); 
M533_delta.pre = mean([M533.sess1.delta; M533.sess2.delta;  M533.sess3.delta;  M533.sess4.delta; M533.sess5.delta; M533.sess6.delta;  M533.sess7.delta]); 
M533_swr.pre = mean([M533.sess1.swr; M533.sess2.swr;  M533.sess3.swr;  M533.sess4.swr; M533.sess5.swr; M533.sess6.swr;  M533.sess7.swr]); 
M533_theta_delta.pre = mean([M533.sess1.theta; M533.sess2.theta;  M533.sess3.theta;  M533.sess4.theta; M533.sess5.theta; M533.sess6.theta;  M533.sess7.theta]./[M533.sess1.delta; M533.sess2.delta;  M533.sess3.delta;  M533.sess4.delta; M533.sess5.delta; M533.sess6.delta;  M533.sess7.delta]);

%M534_freq.pre = mean([M534.sess1.freq;  M534.sess3.freq;  M534.sess4.freq]); 
%M534_dur.pre = mean([M534.sess1.dur;  M534.sess3.dur;  M534.sess4.dur]); 

M545_freq.pre = mean([M545.sess4.freq;  M545.sess5.freq;  M545.sess6.freq; M545.sess7.freq]); 
M545_dur.pre = mean([M545.sess4.dur;  M545.sess5.dur;  M545.sess6.dur; M545.sess7.dur]); 
M545_swr_count.pre = mean([M545.sess4.swr_count;  M545.sess5.swr_count;  M545.sess6.swr_count; M545.sess7.swr_count]); 
M545_theta.pre = mean([M545.sess4.theta;  M545.sess5.theta;  M545.sess6.theta; M545.sess7.theta]); 
M545_delta.pre = mean([M545.sess4.delta;  M545.sess5.delta;  M545.sess6.delta; M545.sess7.delta]); 
M545_swr.pre = mean([M545.sess4.swr;  M545.sess5.swr;  M545.sess6.swr; M545.sess7.swr]); 
M545_theta_delta.pre = mean([M545.sess4.theta;  M545.sess5.theta;  M545.sess6.theta; M545.sess7.theta]./[M545.sess4.swr;  M545.sess5.swr;  M545.sess6.swr; M545.sess7.swr]);

M547_freq.pre = mean([M547.sess1.freq; M547.sess2.freq;  M547.sess3.freq;  M547.sess4.freq; M547.sess5.freq; M547.sess6.freq]); 
M547_dur.pre = mean([M547.sess1.dur; M547.sess2.dur;  M547.sess3.dur;  M547.sess4.dur; M547.sess5.dur; M547.sess6.dur]); 
M547_swr_count.pre = mean([M547.sess1.swr_count; M547.sess2.swr_count;  M547.sess3.swr_count;  M547.sess4.swr_count; M547.sess5.swr_count; M547.sess6.swr_count]); 
M547_theta.pre = mean([M547.sess1.theta; M547.sess2.theta;  M547.sess3.theta;  M547.sess4.theta; M547.sess5.theta; M547.sess6.theta]); 
M547_delta.pre = mean([M547.sess1.delta; M547.sess2.delta;  M547.sess3.delta;  M547.sess4.delta; M547.sess5.delta; M547.sess6.delta]); 
M547_swr.pre = mean([M547.sess1.swr; M547.sess2.swr;  M547.sess3.swr;  M547.sess4.swr; M547.sess5.swr; M547.sess6.swr]); 
M547_theta_delta.pre = mean([M547.sess1.theta; M547.sess2.theta;  M547.sess3.theta;  M547.sess4.theta; M547.sess5.theta; M547.sess6.theta]./[M547.sess1.delta; M547.sess2.delta;  M547.sess3.delta;  M547.sess4.delta; M547.sess5.delta; M547.sess6.delta]); 

M548_freq.pre = mean([M548.sess1.freq; M548.sess2.freq;  M548.sess3.freq;  M548.sess4.freq; M548.sess5.freq; M548.sess6.freq;  M548.sess7.freq]); 
M548_dur.pre = mean([M548.sess1.dur; M548.sess2.dur;  M548.sess3.dur;  M548.sess4.dur; M548.sess5.dur; M548.sess6.dur;  M548.sess7.dur]); 
M548_swr_count.pre = mean([M548.sess1.swr_count; M548.sess2.swr_count;  M548.sess3.swr_count;  M548.sess4.swr_count; M548.sess5.swr_count; M548.sess6.swr_count;  M548.sess7.swr_count]); 
M548_theta.pre = mean([M548.sess1.theta; M548.sess2.theta;  M548.sess3.theta;  M548.sess4.theta; M548.sess5.theta; M548.sess6.theta;  M548.sess7.theta]); 
M548_delta.pre = mean([M548.sess1.delta; M548.sess2.delta;  M548.sess3.delta;  M548.sess4.delta; M548.sess5.delta; M548.sess6.delta;  M548.sess7.delta]); 
M548_swr.pre = mean([M548.sess1.swr; M548.sess2.swr;  M548.sess3.swr;  M548.sess4.swr; M548.sess5.swr; M548.sess6.swr;  M548.sess7.swr]); 
M548_theta_delta.pre = mean([M548.sess1.theta; M548.sess2.theta;  M548.sess3.theta;  M548.sess4.theta; M548.sess5.theta; M548.sess6.theta;  M548.sess7.theta]./[M548.sess1.delta; M548.sess2.delta;  M548.sess3.delta;  M548.sess4.delta; M548.sess5.delta; M548.sess6.delta;  M548.sess7.delta]);

%% mean and plot over mice

low_c = [82,137,199]./255;%[0,104,87]./255; 
color2 = [123,175,222]./255;

mean_freq = mean([M433_freq.pre; M452_freq.pre; M453_freq.pre; M460_freq.pre; M533_freq.pre; M545_freq.pre; M547_freq.pre; M548_freq.pre]);
mean_dur = mean([M433_dur.pre; M452_dur.pre; M453_dur.pre; M460_dur.pre; M533_dur.pre; M545_dur.pre; M547_dur.pre; M548_dur.pre]);
mean_swr_count = mean([M433_swr_count.pre; M452_swr_count.pre;  M453_swr_count.pre; M460_swr_count.pre; M533_swr_count.pre; M545_swr_count.pre; M547_swr_count.pre; M548_swr_count.pre]);
mean_theta = mean([M433_theta.pre; M452_theta.pre;  M453_theta.pre; M460_theta.pre; M533_theta.pre; M545_theta.pre; M547_theta.pre; M548_theta.pre]);
mean_delta = mean([M433_delta.pre; M452_delta.pre;  M453_delta.pre; M460_delta.pre; M533_delta.pre; M545_delta.pre; M547_delta.pre; M548_delta.pre]);
mean_theta_delta = mean([M433_theta_delta.pre; M452_theta_delta.pre;  M453_theta_delta.pre; M460_theta_delta.pre; M533_theta_delta.pre; M545_theta_delta.pre; M547_theta_delta.pre; M548_theta_delta.pre]);
mean_swr = mean([M433_swr.pre; M452_swr.pre;  M453_swr.pre; M460_swr.pre; M533_swr.pre; M545_swr.pre; M547_swr.pre; M548_swr.pre]);

mean_freq_pre = mean_freq(1)
mean_freq_post = mean_freq(2)
mean_dur_pre = mean_dur(1)
mean_dur_post = mean_dur(2)
mean_swr_count_pre = mean_swr_count(1)
mean_swr_count_post = mean_swr_count(2)
mean_theta_pre = mean_theta(1)
mean_theta_post = mean_theta(2)
mean_delta_pre = mean_delta(1)
mean_delta_post = mean_delta(2)
mean_swr_pre = mean_swr(1)
mean_swr_post = mean_swr(2)
mean_theta_delta_pre = mean_theta_delta(1)
mean_theta_delta_post = mean_theta_delta(2)

std_freq_pre = std([M433_freq.pre(1); M452_freq.pre(1); M453_freq.pre(1); M460_freq.pre(1); M533_freq.pre(1); M545_freq.pre(1); M547_freq.pre(1); M548_freq.pre(1)])
std_freq_post= std([M433_freq.pre(2); M452_freq.pre(2); M453_freq.pre(2); M460_freq.pre(2); M533_freq.pre(2); M545_freq.pre(2); M547_freq.pre(2); M548_freq.pre(2)])
std_dur_pre = std([M433_dur.pre(1); M452_dur.pre(1); M453_dur.pre(1); M460_dur.pre(1); M533_dur.pre(1); M545_dur.pre(1); M547_dur.pre(1); M548_dur.pre(1)])
std_dur_post= std([M433_dur.pre(2); M452_dur.pre(2); M453_dur.pre(2); M460_dur.pre(2); M533_dur.pre(2); M545_dur.pre(2); M547_dur.pre(2); M548_dur.pre(2)])
std_swr_count_pre = std([M433_swr_count.pre(1); M452_swr_count.pre(1); M453_swr_count.pre(1); M460_swr_count.pre(1); M533_swr_count.pre(1); M545_swr_count.pre(1); M547_swr_count.pre(1); M548_swr_count.pre(1)])
std_swr_count_post= std([M433_swr_count.pre(2); M452_swr_count.pre(2); M453_swr_count.pre(2); M460_swr_count.pre(2); M533_swr_count.pre(2); M545_swr_count.pre(2); M547_swr_count.pre(2); M548_swr_count.pre(2)])
std_theta_pre = std([M433_theta.pre(1); M452_theta.pre(1); M453_theta.pre(1); M460_theta.pre(1); M533_theta.pre(1); M545_theta.pre(1); M547_theta.pre(1); M548_theta.pre(1)])
std_theta_post= std([M433_theta.pre(2); M452_theta.pre(2); M453_theta.pre(2); M460_theta.pre(2); M533_theta.pre(2); M545_theta.pre(2); M547_theta.pre(2); M548_theta.pre(2)])
std_delta_pre = std([M433_delta.pre(1); M452_delta.pre(1); M453_delta.pre(1); M460_delta.pre(1); M533_delta.pre(1); M545_delta.pre(1); M547_delta.pre(1); M548_delta.pre(1)])
std_delta_post= std([M433_delta.pre(2); M452_delta.pre(2); M453_delta.pre(2); M460_delta.pre(2); M533_delta.pre(2); M545_delta.pre(2); M547_delta.pre(2); M548_delta.pre(2)])
std_swr_pre = std([M433_swr.pre(1); M452_swr.pre(1); M453_swr.pre(1); M460_swr.pre(1); M533_swr.pre(1); M545_swr.pre(1); M547_swr.pre(1); M548_swr.pre(1)])
std_swr_post= std([M433_swr.pre(2); M452_swr.pre(2); M453_swr.pre(2); M460_swr.pre(2); M533_swr.pre(2); M545_swr.pre(2); M547_swr.pre(2); M548_swr.pre(2)])
std_theta_delta_pre = std([M433_theta_delta.pre(1); M452_theta_delta.pre(1); M453_theta_delta.pre(1); M460_theta_delta.pre(1); M533_theta_delta.pre(1); M545_theta_delta.pre(1); M547_theta_delta.pre(1); M548_theta_delta.pre(1)])
std_theta_delta_post= std([M433_theta_delta.pre(2); M452_theta_delta.pre(2); M453_theta_delta.pre(2); M460_theta_delta.pre(2); M533_theta_delta.pre(2); M545_theta_delta.pre(2); M547_theta_delta.pre(2); M548_theta_delta.pre(2)])

% PRE IS WAKE AND POST IS NREM !!! TOO LAZY TO CHANGE !!!
pre_freq = [M433_freq.pre(1); M452_freq.pre(1); M453_freq.pre(1); M460_freq.pre(1); M533_freq.pre(1); M545_freq.pre(1); M547_freq.pre(1); M548_freq.pre(1)];
post_freq = [M433_freq.pre(2); M452_freq.pre(2); M453_freq.pre(2); M460_freq.pre(2); M533_freq.pre(2); M545_freq.pre(2); M547_freq.pre(2); M548_freq.pre(2)];
pre_dur = [M433_dur.pre(1); M452_dur.pre(1); M453_dur.pre(1); M460_dur.pre(1); M533_dur.pre(1); M545_dur.pre(1); M547_dur.pre(1); M548_dur.pre(1)];
post_dur = [M433_dur.pre(2); M452_dur.pre(2); M453_dur.pre(2); M460_dur.pre(2); M533_dur.pre(2); M545_dur.pre(2); M547_dur.pre(2); M548_dur.pre(2)];
pre_swr_count = [M433_swr_count.pre(1); M452_swr_count.pre(1); M453_swr_count.pre(1); M460_swr_count.pre(1); M533_swr_count.pre(1); M545_swr_count.pre(1); M547_swr_count.pre(1); M548_swr_count.pre(1)];
post_swr_count = [M433_swr_count.pre(2); M452_swr_count.pre(2); M453_swr_count.pre(2); M460_swr_count.pre(2); M533_swr_count.pre(2); M545_swr_count.pre(2); M547_swr_count.pre(2); M548_swr_count.pre(2)];
pre_theta = [M433_theta.pre(1); M452_theta.pre(1); M453_theta.pre(1); M460_theta.pre(1); M533_theta.pre(1); M545_theta.pre(1); M547_theta.pre(1); M548_theta.pre(1)];
post_theta = [M433_theta.pre(2); M452_theta.pre(2); M453_theta.pre(2); M460_theta.pre(2); M533_theta.pre(2); M545_theta.pre(2); M547_theta.pre(2); M548_theta.pre(2)];
pre_delta = [M433_delta.pre(1); M452_delta.pre(1); M453_delta.pre(1); M460_delta.pre(1); M533_delta.pre(1); M545_delta.pre(1); M547_delta.pre(1); M548_delta.pre(1)];
post_delta = [M433_delta.pre(2); M452_delta.pre(2); M453_delta.pre(2); M460_delta.pre(2); M533_delta.pre(2); M545_delta.pre(2); M547_delta.pre(2); M548_delta.pre(2)];
pre_swr = [M433_swr.pre(1); M452_swr.pre(1); M453_swr.pre(1); M460_swr.pre(1); M533_swr.pre(1); M545_swr.pre(1); M547_swr.pre(1); M548_swr.pre(1)];
post_swr = [M433_swr.pre(2); M452_swr.pre(2); M453_swr.pre(2); M460_swr.pre(2); M533_swr.pre(2); M545_swr.pre(2); M547_swr.pre(2); M548_swr.pre(2)];
pre_theta_delta = [M433_theta_delta.pre(1); M452_theta_delta.pre(1); M453_theta_delta.pre(1); M460_theta_delta.pre(1); M533_theta_delta.pre(1); M545_theta_delta.pre(1); M547_theta_delta.pre(1); M548_theta_delta.pre(1)];
post_theta_delta = [M433_theta_delta.pre(2); M452_theta_delta.pre(2); M453_theta_delta.pre(2); M460_theta_delta.pre(2); M533_theta_delta.pre(2); M545_theta_delta.pre(2); M547_theta_delta.pre(2); M548_theta_delta.pre(2)];

%cd ('C:\Users\mimia\Desktop\BOXCHARTS')
% mean freq
figure(3)
h = boxchart([pre_freq, post_freq]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric
% Scatter points with jitter for better visibility
scatter(repmat(xDataNumeric(1), 1, numel(pre_freq)) + randn(1, numel(pre_freq)) * 0.05, pre_freq, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(2), 1, numel(post_freq)) + randn(1, numel(post_freq)) * 0.05, post_freq, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

xticklabels(["WAKE" "NREM"])
ylabel("Rate of SWRs (SWRs per minute)")
xlabel("Sleep-Wake State")
set(gca,'fontsize', 18)
% p =    6.8861e-04     tstat: -4.3327





set(gcf, 'renderer', 'painters');%
%cd ('C:\Users\mimia\OneDrive\Desktop\SWR_supplement')
%exportgraphics(gcf, 'SWR_freq_mean.eps', 'ContentType','vector');  % Export as PDF
hold off 

[h_freq,p_freq,ci_freq,stats_freq] = ttest(pre_freq, post_freq)

% mean dur
figure(4)
h = boxchart([pre_dur, post_dur]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric
% Scatter points with jitter for better visibility
scatter(repmat(xDataNumeric(1), 1, numel(pre_dur)) + randn(1, numel(pre_dur)) * 0.05, pre_dur, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(2), 1, numel(post_dur)) + randn(1, numel(post_dur)) * 0.05, post_dur, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

xticklabels(["WAKE" "NREM"])
ylabel("SWR duration (ms)")
xlabel("Sleep-Wake State")
set(gca,'fontsize', 18)

set(gcf, 'renderer', 'painters');
%cd ('C:\Users\mimia\OneDrive\Desktop\SWR_supplement')
%exportgraphics(gcf, 'SWR_dur_mean.eps', 'ContentType','vector');  % Export as PDF
hold off 

[h_dur,p_dur,ci_dur,stats_dur] = ttest(pre_dur, post_dur)


% mean swr_count
figure(5)
h = boxchart([pre_swr_count, post_swr_count]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric
% Scatter points with jitter for better visibility
scatter(repmat(xDataNumeric(1), 1, numel(pre_swr_count)) + randn(1, numel(pre_swr_count)) * 0.05, pre_swr_count, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(2), 1, numel(post_swr_count)) + randn(1, numel(post_swr_count)) * 0.05, post_swr_count, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

xticklabels(["WAKE" "NREM"])
ylabel("# of SWRs")
xlabel("Sleep-Wake State")
set(gca,'fontsize', 18)

set(gcf, 'renderer', 'painters');
%cd ('C:\Users\mimia\OneDrive\Desktop\SWR_supplement')
exportgraphics(gcf, 'SWR_count_mean.eps', 'ContentType','vector');  % Export as PDF
hold off 

[h_swr_count,p_swr_count,ci_swr_count,stats_swr_count] = ttest2(pre_swr_count, post_swr_count)


figure(6)
h = boxchart([pre_theta, post_theta]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric
% Scatter points with jitter for better visibility
scatter(repmat(xDataNumeric(1), 1, numel(pre_theta)) + randn(1, numel(pre_theta)) * 0.05, pre_theta, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(2), 1, numel(post_theta)) + randn(1, numel(post_theta)) * 0.05, post_theta, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

xticklabels(["WAKE" "NREM"])
ylabel("Theta Power")
xlabel("Sleep-Wake State")
set(gca,'fontsize', 18)

set(gcf, 'renderer', 'painters');
%cd ('C:\Users\mimia\OneDrive\Desktop\SWR_supplement')
exportgraphics(gcf, 'theta_mean.eps', 'ContentType','vector');  % Export as PDF
hold off 

[h_theta,p_theta,ci_theta,stats_theta] = ttest(pre_theta, post_theta)
% p = 00.48 ; tstat = 3.34 SIGNIFICANTLY DIFFERENT~!



% mean delta
figure(7)
h = boxchart([pre_delta, post_delta]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric
% Scatter points with jitter for better visibility
scatter(repmat(xDataNumeric(1), 1, numel(pre_delta)) + randn(1, numel(pre_delta)) * 0.05, pre_delta, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(2), 1, numel(post_delta)) + randn(1, numel(post_delta)) * 0.05, post_delta, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

xticklabels(["WAKE" "NREM"])
ylabel("Delta Power")
xlabel("Sleep-Wake State")
set(gca,'fontsize', 18)

set(gcf, 'renderer', 'painters');
%cd ('C:\Users\mimia\OneDrive\Desktop\SWR_supplement')
exportgraphics(gcf, 'delta_mean.eps', 'ContentType','vector');  % Export as PDF
hold off 

[h_delta,p_delta,ci_delta,stats_delta] = ttest(pre_delta, post_delta)


figure(8)
h = boxchart([pre_swr, post_swr]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric
% Scatter points with jitter for better visibility
scatter(repmat(xDataNumeric(1), 1, numel(pre_swr)) + randn(1, numel(pre_swr)) * 0.05, pre_swr, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(2), 1, numel(post_swr)) + randn(1, numel(post_swr)) * 0.05, post_swr, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

xticklabels(["WAKE" "NREM"])
ylabel("SWR Power")
xlabel("Sleep-Wake State")
set(gca,'fontsize', 18)

set(gcf, 'renderer', 'painters');
%cd ('C:\Users\mimia\OneDrive\Desktop\SWR_supplement')
exportgraphics(gcf, 'swr__power_mean.eps', 'ContentType','vector');  % Export as PDF
hold off 

[h_swr,p_swr,ci_swr,stats_swr] = ttest(pre_swr, post_swr)


% mean theta_delta
figure(9)
h = boxchart([pre_theta_delta, post_theta_delta]);  % Combine the data for the boxchart
hold on;
h.BoxFaceColor = color2;
xDataNumeric = double(h.XData);  % Convert categorical XData to numeric
% Scatter points with jitter for better visibility
scatter(repmat(xDataNumeric(1), 1, numel(pre_theta_delta)) + randn(1, numel(pre_theta_delta)) * 0.05, pre_theta_delta, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

scatter(repmat(xDataNumeric(2), 1, numel(post_theta_delta)) + randn(1, numel(post_theta_delta)) * 0.05, post_theta_delta, ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

xticklabels(["WAKE" "NREM"])
ylabel("Theta/Delta")
xlabel("Sleep-Wake State")
set(gca,'fontsize', 18)

set(gcf, 'renderer', 'painters');
%cd ('C:\Users\mimia\OneDrive\Desktop\theta_delta_supplement')
exportgraphics(gcf, 'theta_delta_mean.eps', 'ContentType','vector');  % Export as PDF
hold off 

[h_theta_delta,p_theta_delta,ci_theta_delta,stats_theta_delta] = ttest2(pre_theta_delta, post_theta_delta)