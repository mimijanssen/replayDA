%% Load Each File 
% save each variable into a new matrix (sess_all_events)

clear all; clc;
D = 'C:\Data\M437\processed_data';
cd C:\Data\M437\processed_data;

S = dir(fullfile(D,'*'));

for i = 3:1:8
    load(S(i).name)
    sess_all_evt{i} = evtt;
    sess_all_t{i} = t;
    sess_all_zdF{i} = zdF; %zdF; 
end

%% Extract all events ordered 
all_events = [];

for i = 3:1:8 % for all of the sessions
    load(S(i).name)
    all_events = [all_events ;sort([evtt{1,1}, evtt{1,2}, evtt{1,3}])];
end

% initialize values for zdF matrix 
window = 8/0.001; 


%% Extract only time points we care about. 
high_sess = [];
med_sess = [];
low_sess = [];

for i = 3:1:8 % for all of the sessions
    load(S(i).name)
    high_sess = [high_sess ;evtt{1,1}];
    med_sess = [med_sess ;evtt{1,2}];
    low_sess = [low_sess ;evtt{1,3}];
end

% initialize values for zdF matrix 
window = 8/0.001; 


%% zdF Matrix 
count_h = 0 ;
count_m = 0 ;
count_l = 0 ; 

high = zeros(length(high_sess),window*2); %rows = trial, columns = signal 
med = zeros(length(med_sess),window*2);
low = zeros(length(low_sess),window*2);

clear high_sess med_sess low_sess

high_sess_avg = [];
med_sess_avg = [];
low_sess_avg = [];
count = 0;

for i = 3:1:8 
    count = count + 1;
    high_sess = [];
    med_sess = [];
    low_sess = [];
    load(S(i).name)
    
    for iter = 1:1:length(evtt{1,1})
        try
        count_h = count_h +1 ; 
        time_h_ind = evtt{1,1}(iter); % each time index 
        indxpb = find(abs(t-time_h_ind) < 0.0005); % find the time of the event. 
        init_trial = (indxpb - (window)); % 8 seconds; remember 0.0002 is the time in seconds for each voltage point
        end_trial = (indxpb + (window));
        high(count_h,:) = F_detrend(init_trial:end_trial-1); % zdF 
        high_sess(count_h,:) = F_detrend(init_trial:end_trial-1); 
        catch
            disp('you cut this')
            continue
        end
    end
    high_sess_avg(count,:) = mean(high_sess);

    for iter = 1:1:length(evtt{1,2})
        try 
        count_m = count_m +1 ; 
        time_m_ind = evtt{1,2}(iter); % each time index 
        indxpb = find(abs(t-time_m_ind) < 0.0005); % find the time of the event. 
        init_trial = (indxpb - (window)); % 8 seconds; remember 0.0002 is the time in seconds for each voltage point
        end_trial = (indxpb + (window));
        med(count_m,:) = F_detrend(init_trial:end_trial-1); 
        med_sess(count_m,:) = F_detrend(init_trial:end_trial-1); 
        catch
            disp('you cut this')
            continue
        end
    end
    med_sess_avg(count,:) = mean(med_sess);

    for iter = 1:1:length(evtt{1,3})
        count_l = count_l +1 ; 
        time_l_ind = evtt{1,3}(iter); % each time index 
        indxpb = find(abs(t-time_l_ind) < 0.0005); % find the time of the event. 
        init_trial = (indxpb - (window)); % 8 seconds; remember 0.0002 is the time in seconds for each voltage point
        end_trial = (indxpb + (window));
        low(count_l,:) = F_detrend(init_trial:end_trial-1); 
        low_sess(count_l,:) = F_detrend(init_trial:end_trial-1); 
        t_shared = t(init_trial:end_trial-1)-t(init_trial); % double the window
    end
    low_sess_avg(count,:) = mean(low_sess);

end


%% zdF Matrix for all events 
count = 0; 
count1 = 0; 
count2 = 0; 
count3 = 0; 
count4 = 0; 

q1 =  zeros(length(all_events)/4,window*2); 
q2 =  zeros(length(all_events)/4,window*2); 
q3 =  zeros(length(all_events)/4,window*2); 
q4 =  zeros(length(all_events)/4,window*2); 

q1_sess_avg = [];
q2_sess_avg = [];
q3_sess_avg = [];
q4_sess_avg = [];

for i = 3:1:8 
    count = count + 1;
    q1_sess = [];
    q2_sess = [];
    q3_sess = [];
    q4_sess = [];

    load(S(i).name)
    
    for iter = 1:1:length(all_events)/4 % for 1:15 events 
        try
        count1 = count1 + 1; 
        time_ind = all_events(count,iter); % each time 
        indxpb = find(abs(t-time_ind) < 0.0005); % find the time of the event. 
        init_trial = (indxpb - (window)); % 8 seconds; remember 0.0002 is the time in seconds for each voltage point
        end_trial = (indxpb + (window));
        q1_sess(count1,:) = zdF(init_trial:end_trial-1); % zdF 
        catch
            disp('you cut this')
            continue
        end
    end
    q1_sess_avg(count,:) = mean(q1_sess);
    
    for iter = 16:1:30 % for 1:15 events 
        try
        count2 = count2 +1 ; 
        time_ind = all_events(count,iter); % each time 
        indxpb = find(abs(t-time_ind) < 0.0005); % find the time of the event. 
        init_trial = (indxpb - (window)); % 8 seconds; remember 0.0002 is the time in seconds for each voltage point
        end_trial = (indxpb + (window));
        q2_sess(count2,:) = zdF(init_trial:end_trial-1); % zdF 
        catch
            disp('you cut this')
            continue
        end
    end
    q2_sess_avg(count,:) = mean(q2_sess);

    for iter = 31:1:45 % for 1:15 events 
        try
        count3 = count3 +1 ; 
        time_ind = all_events(count,iter); % each time 
        indxpb = find(abs(t-time_ind) < 0.0005); % find the time of the event. 
        init_trial = (indxpb - (window)); % 8 seconds; remember 0.0002 is the time in seconds for each voltage point
        end_trial = (indxpb + (window));
        q3_sess(count3,:) = zdF(init_trial:end_trial-1); % zdF 
        catch
            disp('you cut this')
            continue
        end
    end
    q3_sess_avg(count,:) = mean(q3_sess);

    for iter = 46:1:60 % for 1:15 events 
        try
        count4 = count4 +1 ; 
        time_ind = all_events(count,iter); % each time 
        indxpb = find(abs(t-time_ind) < 0.0005); % find the time of the event. 
        init_trial = (indxpb - (window)); % 8 seconds; remember 0.0002 is the time in seconds for each voltage point
        end_trial = (indxpb + (window));
        q4_sess(count4,:) = zdF(init_trial:end_trial-1); % zdF 
        t_shared = t(init_trial:end_trial-1)-t(init_trial); % double the window
        catch
            disp('you cut this')
            continue
        end
    end
    q4_sess_avg(count,:) = mean(q4_sess);

end
%% Scrambling the session averages 
all_sessions = [high_sess_avg; med_sess_avg; low_sess_avg]; %18 rows
r = randperm(18); 
random_sess1 = [all_sessions(r(1:6),:)] ; % first 6 
random_sess2 = [all_sessions(r(7:12),:)] ; 
random_sess3 = [all_sessions(r(13:18),:)] ; 

%% Plot Averaged Values 
addpath(genpath('C:\Users\mimia\Documents\GitHub\shadedErrorBar'));

% Color
med_c = [104,187,225]./255; % rgb(167, 199, 231) rgb(255, 165, 0)  blue: rgb(104,187,227)
low_c = [0,104,87]./255; 
high_c = [255, 165,0]./255; % rgb(204, 204, 255) periwinkle

all_colors = [{high_c},{med_c},{low_c}];

% Parameters 
window = 8/0.001; 

% average by columns
avg_high = mean(high_sess_avg); 
avg_med = mean(med_sess_avg);
avg_low = mean(low_sess_avg);

% std by columns
std_high = std(high_sess_avg)/sqrt(size(high_sess_avg,1)); 
std_med = std(med_sess_avg)/sqrt(size(med_sess_avg,1));
std_low = std(low_sess_avg)/sqrt(size(low_sess_avg,1));

figure(1)
shadedErrorBar(t_shared,avg_high,std_high,'lineprops',{'-','color',high_c,'MarkerFaceColor',high_c});
hold on
shadedErrorBar(t_shared,avg_med,std_med,'lineprops',{'-','color',med_c,'MarkerFaceColor',med_c});
shadedErrorBar(t_shared,avg_low,std_low,'lineprops',{'-','color',low_c,'MarkerFaceColor',low_c});

xlabel('Time from Photobeam Break (s)')
xlim([0 16])
xticks([0 8 16])
xticklabels({'-8','0','8'})
legend('high','medium','low');
legend boxoff
ylabel('Mean Signal (dF z-scored)')
title('Averaged Sessions')

set(gcf,'Color',[1,1,1])
shg

orient(gcf,'landscape')

figure(2)
shadedErrorBar(t_shared,avg_high,std_high,'lineprops',{'-','color',high_c,'MarkerFaceColor',high_c});
hold on
shadedErrorBar(t_shared,avg_med,std_med,'lineprops',{'-','color',med_c,'MarkerFaceColor',med_c});
shadedErrorBar(t_shared,avg_low,std_low,'lineprops',{'-','color',low_c,'MarkerFaceColor',low_c});

xlabel('Time from Photobeam Break (s)')
xlim([8 16])
xticks([8 12 16])
xticklabels({'0','4','8'})
legend('high','medium','low');
legend boxoff
ylabel('Mean Signal (dF z-scored)')
title('Averaged Sessions')

set(gcf,'Color',[1,1,1])
shg

orient(gcf,'landscape')


%% Plot Random Sessions
addpath(genpath('C:\Users\mimia\Documents\GitHub\shadedErrorBar'));

% Color
med_c = [104,187,225]./255; % rgb(167, 199, 231) rgb(255, 165, 0)  blue: rgb(104,187,227)
low_c = [0,104,87]./255; 
high_c = [255, 165,0]./255; % rgb(204, 204, 255) periwinkle

all_colors = [{high_c},{med_c},{low_c}];

% Parameters 
window = 8/0.001; 

% average by columns
avg_high = mean(random_sess1); 
avg_med = mean(random_sess2);
avg_low = mean(random_sess3);

% std by columns
std_high = std(random_sess1)/sqrt(size(random_sess1,1)); 
std_med = std(random_sess2)/sqrt(size(random_sess2,1));
std_low = std(random_sess3)/sqrt(size(random_sess3,1));

figure(1)
shadedErrorBar(t_shared,avg_high,std_high,'lineprops',{'-','color',high_c,'MarkerFaceColor',high_c});
hold on
shadedErrorBar(t_shared,avg_med,std_med,'lineprops',{'-','color',med_c,'MarkerFaceColor',med_c});
shadedErrorBar(t_shared,avg_low,std_low,'lineprops',{'-','color',low_c,'MarkerFaceColor',low_c});

xlabel('Time from Photobeam Break (s)')
xlim([0 16])
xticks([0 8 16])
xticklabels({'-8','0','8'})
legend('randsess1','randsess2','randsess3');
legend boxoff
ylabel('Mean Signal (dF z-scored)')
title('Averaged Sessions')

set(gcf,'Color',[1,1,1])
shg

orient(gcf,'landscape')

figure(2)
shadedErrorBar(t_shared,avg_high,std_high,'lineprops',{'-','color',high_c,'MarkerFaceColor',high_c});
hold on
shadedErrorBar(t_shared,avg_med,std_med,'lineprops',{'-','color',med_c,'MarkerFaceColor',med_c});
shadedErrorBar(t_shared,avg_low,std_low,'lineprops',{'-','color',low_c,'MarkerFaceColor',low_c});

xlabel('Time from Photobeam Break (s)')
xlim([8 16])
xticks([8 12 16])
xticklabels({'0','4','8'})
legend('randsess1','randsess2','randsess3');
legend boxoff
ylabel('Mean Signal (dF z-scored)')
title('Averaged Sessions')

set(gcf,'Color',[1,1,1])
shg

orient(gcf,'landscape')



%% Plot Quarters
addpath(genpath('C:\Users\mimia\Documents\GitHub\shadedErrorBar'));

% Color
q1_c = [104,187,225]./255; % rgb(167, 199, 231) rgb(255, 165, 0)  blue: rgb(104,187,227)
q2_c = [0,104,87]./255; 
q3_c = [255, 165,0]./255; % rgb(204, 204, 255) periwinkle
q4_c = [219, 112,147]./255; % rgb(204, 204, 255) periwinklergb(219,112,147)

all_colors = [{q1_c},{q2_c},{q3_c},{q4_c}];

% Parameters 
window = 8/0.001; 

% average by columns
avg_q1 = mean(q1_sess_avg); 
avg_q2 = mean(q2_sess_avg);
avg_q3 = mean(q3_sess_avg);
avg_q4 = mean(q4_sess_avg);


% std by columns
std_q1 = std(q1_sess_avg)/sqrt(size(q1_sess_avg,1)); 
std_q2 = std(q2_sess_avg)/sqrt(size(q2_sess_avg,1));
std_q3 = std(q3_sess_avg)/sqrt(size(q3_sess_avg,1));
std_q4 = std(q4_sess_avg)/sqrt(size(q4_sess_avg,1));

figure(1)
shadedErrorBar(t_shared,avg_q1,std_q1,'lineprops',{'-','color',q1_c,'MarkerFaceColor',q1_c});
hold on
shadedErrorBar(t_shared,avg_q2,std_q2,'lineprops',{'-','color',q2_c,'MarkerFaceColor',q2_c});
shadedErrorBar(t_shared,avg_q3,std_q3,'lineprops',{'-','color',q3_c,'MarkerFaceColor',q3_c});
shadedErrorBar(t_shared,avg_q4,std_q4,'lineprops',{'-','color',q4_c,'MarkerFaceColor',q4_c});


xlabel('Time from Photobeam Break (s)')
xlim([0 16])
xticks([0 8 16])
xticklabels({'-8','0','8'})
legend('q1','q2','q3','q4');
legend boxoff
ylabel('Mean Signal (detrended)')
title('Averaged Sessions')

set(gcf,'Color',[1,1,1])
shg

orient(gcf,'landscape')

figure(2)
shadedErrorBar(t_shared,avg_q1,std_q1,'lineprops',{'-','color',q1_c,'MarkerFaceColor',q1_c});
hold on
shadedErrorBar(t_shared,avg_q2,std_q2,'lineprops',{'-','color',q2_c,'MarkerFaceColor',q2_c});
shadedErrorBar(t_shared,avg_q3,std_q3,'lineprops',{'-','color',q3_c,'MarkerFaceColor',q3_c});
shadedErrorBar(t_shared,avg_q4,std_q4,'lineprops',{'-','color',q4_c,'MarkerFaceColor',q4_c});

xlabel('Time from Photobeam Break (s)')
xlim([8 16])
xticks([8 12 16])
xticklabels({'0','4','8'})
legend('q1','q2','q3','q4');
legend boxoff
ylabel('Mean Signal (detrended)')
title('Averaged Sessions')

set(gcf,'Color',[1,1,1])
shg

orient(gcf,'landscape')