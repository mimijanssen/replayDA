%% Load
clear; clc;
cd 'C:\Users\mimia\Desktop\linearfit_constrained';
D = 'C:\Users\mimia\Desktop\linearfit_constrained';
S = dir(fullfile(D,'*processed*'));
N = length(S); 
% Consider saving sampling rate as well in the future 

% Initializing variables and timestamps 
inject_start = zeros(1,8);
inject_end = zeros(1,8); % to find latest time after injection 
ped_start = zeros(1,8); 
ped2_end = zeros(1,8); % to find earliest time of pedestal 2 ends

all_dF = []; % can optimize this more later
all_FP_z = [];
all_zdF = [];
all_shared_t = []; 

% Finding timestamps 
for iit = 1:N
    load(S(iit).name);
    inject_start(iit) = evtt(2); % consider changing this to automatically find based on evtlabel
    inject_end(iit) = evtt(3);
    ped_start(iit) = evtt(1);
    ped2_end(iit) = evtt(4);
end

% I know the sampling rate is 1000 Hz --> 1000 samples per second, so a
% sample every 0.001 s

% Find the session and time we want to align to 
align_time = max(inject_end); % max time 
align_session = find(inject_end == align_time); % session with max time stamp 
load(S(align_session).name)
align_index = find(align_time == t);

% Find xlim 
% max pedestal start time index
max_ped_start = max( ped_start);
max_ped_start_ind = find(ped_start == max_ped_start);
load(S(max_ped_start_ind).name)
xlim_init = find(max_ped_start == t);
% MORE IMPORTANT TO KNOW THE SESSION.. I think I need to look for first Nan
% in that session... 

% min pedestal 2 end time index
% this will also be used to truncate the data. 
min_ped2_end = min( ped2_end);
min_ped2_end_ind = find(ped2_end == min_ped2_end);
load(S(min_ped2_end_ind).name)
xlim_end = find(min_ped2_end == t);

% Saving variables 
% Save each files variables in a new row in a array, but only the
% shared time stamps 
% fiber data is shifted to latest time after injection 

for iiv = 1:N 
    clear vars dF evtlabel evtt FP_z zdF
    load(S(iiv).name); 

    % inject_end time index
    index = find(inject_end(iiv) == t); %finds index of inject time 
    
    % find difference from that index and align time 
    index_diff = align_index(1) - index; % will always be positive 
    if index_diff < 0 % mental check
        disp('error');
    end

    % append NaN to the front of the data based on index difference. 
    % create a row with number of elements of Nans needed
    buffer = NaN(1,index_diff);
    
    % Shift over all event codes to correspond with this buffer now... 
    additional_time = 0.001*length(buffer); % this is the number of indicies to shift the data, where each index is 0.001 s apart 
    % add additional time to ped start variable 
    ped_start(iiv) = ped_start(iiv)+additional_time; % updated pd_start with new buffer 
    
    if iiv == 6 
        add_time_buffer6 = additional_time;
    end

    % append to data in a column
    one_dF = [buffer'; dF];
    one_FP_z = [buffer'; FP_z];
    one_zdF = [buffer'; zdF];
    one_shared_t = [buffer'; (t + additional_time)]; % add time to where data is shifted?

    % it is possible I loose some points this way becuase if I shift over
    % the data, it might now be longer than the session I'm aligning to. 
    % IS aligned session end earlier than min ped2 session end? 
    load(S(align_session).name)
    align_ped2_index = find(ped2_end(align_session) == t); 
    
    if xlim_end > align_ped2_index % if the index is longer than aligned end index then use the aligned end index else, use xlim_end
        % use align_ped2_index 
        one_dF = one_dF(1:align_ped2_index,1)';
        one_FP_z = one_FP_z(1:align_ped2_index,1)';
        one_zdF = one_zdF(1:align_ped2_index,1)';
        one_shared_t = one_shared_t(1:align_ped2_index,1)';
    else
        % use xlim_end
        one_dF = one_dF(1:xlim_end,1)';
        one_FP_z = one_FP_z(1:xlim_end,1)';
        one_zdF = one_zdF(1:xlim_end,1)';
        one_shared_t = one_shared_t(1:xlim_end,1)';
    end

    % save data
    all_dF = [all_dF; one_dF];
    all_FP_z = [all_FP_z; one_FP_z];
    all_zdF = [all_zdF; one_zdF];
    all_shared_t = [all_shared_t; one_shared_t]; 

end

% truncate start time for all data 
max_ped_start2 = max(ped_start); % time of max pedestal 1 time -- 106.11
align_ped_start2 = ped_start(align_session); % 4.99
% now we know that align is longer ... 
% we actually want max start time... becuase it will be shorter when we
% truncate. 

% load max start index 
load(S(max_ped_start_ind).name)
xlim_init = find(abs(all_shared_t(1,:)-(max_ped_start2+add_time_buffer6)) < 0.001); % index for new time 

% find difference from that index and align time 
index_diff = align_index - xlim_init; % will always be positive 
if index_diff < 0 % mental check
   disp('error');
end
%%
%xlim_init = xlim_init + index_diff; % the new shifted index 

% 1,634,127 --> 828,041
all_dF = all_dF(:,xlim_init(1):end);
all_FP_z = all_FP_z(:,xlim_init(1):end);
all_zdF = all_zdF(:,xlim_init(1):end);
all_shared_t = all_shared_t(:,xlim_init(1):end); 


%% Plot data
addpath(genpath('C:\Users\mimia\Documents\GitHub\shadedErrorBar'));
% Amphetamine 
color_amp = [114 147 203]./255;
% Over M383 (2 & 4)
M1_mean_amp = mean(all_zdF([2,4],:));
%SEM1_amp = std(all_zdF([2,4],:),0,1)/sqrt(length(all_zdF([2,4],:)));
% Over M411 (5 & 7)
M2_mean_amp = mean(all_zdF([5,7],:));
%SEM2_amp = std(all_zdF([5,7],:),0,1)/sqrt(length(all_zdF([5,7],:)));
% Avg_Amp
%avg_amp = mean([M1_mean_amp;M2_mean_amp]);
%SEM_amp = std([M1_mean_amp;M2_mean_amp],0,1)/sqrt(length([M1_mean_amp;M2_mean_amp]));

% Saline 
color_sal = [128 133 133]./255; 
%loyolagreen = 1/255*[0,104,87];

% Over M383 (1 & 3)
M1_mean_sal = mean(all_zdF([1,3],:));
%SEM1_sal = std(all_zdF([1,3],:),0,1)/sqrt(length(all_zdF([1,3],:)));
% Over M411 (6 & 8)
M2_mean_sal = mean(all_zdF([6,8],:));
%SEM2_sal = std(all_zdF([6,8],:),0,1)/sqrt(length(all_zdF([6,8],:)));
% Avg_Sal
avg_sal = mean([M1_mean_sal;M2_mean_sal]);
SEM_sal = std([M1_mean_sal;M2_mean_sal],0,1)/sqrt(length([M1_mean_sal;M2_mean_sal]));

c = [0.8 0.7 0.8];
%L = patch([linspace(time(ind_drug),time(ind_admin),40) fliplr(linspace(time(ind_drug),time(ind_admin),40))],[(linspace(y(1),y(1),40)) fliplr((linspace(y(2),y(2),40)))],c,'EdgeColor','none');
%uistack(L,'bottom')
% legend('Injection','Location','southeast')
% legend('boxoff')
% align_index is when injection ends 
% find when injection starts 
% Find the session and time we want to align to

% overlay transparent background for each mouse... 
% This is just using session 1
load(S(align_session).name)
align_index_start = find(abs(inject_start(align_session) -all_shared_t(1,:)) < 0.001); % find the time that matches new time index 
align_index = find(abs(inject_end(align_session) - all_shared_t(1,:)) < 0.001);
time2 = all_shared_t(1,:);
center = (time2(align_index_start) + time2(align_index(1)))/2;

figure(1)
clf
%line_B = shadedErrorBar(ages,meanB,CIB,'lineProps',{'-b','MarkerFaceColor','b','LineWidth', 3});
%shadedErrorBar(all_shared_t(1,:),avg_sal,SEM_sal,'lineProps',{'-b','MarkerFaceColor',color_sal,'LineWidth',3}); %'-g','transparent',1)
shadedErrorBar(all_shared_t(1,:),avg_sal,SEM_sal,'lineProps','-k','transparent',1)
hold on
plot(all_shared_t(1,:),avg_sal)
hold on
shadedErrorBar(all_shared_t(1,:),avg_amp,SEM_amp,'lineProps','-b','transparent',1)
hold on
plot(all_shared_t(1,:),avg_amp)
hold on 
y = ylim;
L = patch([linspace(time2(align_index_start(1)),time2(align_index(1)),40) fliplr(linspace(time2(align_index_start(1)),time2(align_index(1)),40))],[(linspace(y(1),y(1),40)) fliplr((linspace(y(2),y(2),40)))],c,'EdgeColor','none');
uistack(L,'bottom')
hold off
title('Mean FP Trace')
ylabel('Signal (Z-scored)')
xlabel('Time from Injection (min)')
xlim([(center-600) (center+900)]) % change x lim based on below 
xticks([(center-600) (center-300) (center) center+300 center+600 center+900])
xticklabels({'-10','-5','0','5','10','15'})
legend('','Saline','','Amphetamine','','Location','southeast')
legend('boxoff')

set(gcf,'Color',[1,1,1])
shg

orient(gcf,'landscape')

% try changing aligment to be in the center of the injection? 
