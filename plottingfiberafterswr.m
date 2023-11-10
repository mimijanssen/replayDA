%% Plot all detected SWRS? 
%% Load data SWR information:
% start and end interval for SWR
load('2023-09-19_M433_recording1detectedSWRs.mat')

% CSC used 
csc_name = [];
csc_name.fc = {'CSC15.ncs'};
csc = LoadCSC(csc_name);

%% Load fiber info 
%z-scorded df/f (zdf)
load('M433-2023-09-19processed.mat')

%% Plot Fiber Photometry data after SWR events 

% initialize LFP 
lfp_time = csc.tvec - csc.tvec(1);
lfp = csc.data;


% initialize SWR interval
SWR_start = evt.tstart - csc.tvec(1);
SWR_end = evt.tend - csc.tvec(1);
SWR_iv = [SWR_start SWR_end]; 

%% plot first example of a SWR
% find the closest time inedex in lfp from the SWR interval 
SWR_ind_start = nearest_idx3(SWR_iv(:,1),lfp_time); %find(abs(lfp-SWR_iv(1,1)) < 0.0005); % only for one but can we extend this to everything??
SWR_ind_end = nearest_idx3(SWR_iv(:,2),lfp_time); %find(abs(lfp-SWR_iv(1,2)) < 0.0005); % only for one but can we extend this to everything??

%%
%example plots !
figure (1)
plot((SWR_ind_start(1):SWR_ind_end(1))/1000, lfp(SWR_ind_start(1):SWR_ind_end(1)));
ylabel('lfp (V)'); xlabel('Time (s)');
figure (2)
plot((SWR_ind_start(10):SWR_ind_end(10))/1000, lfp(SWR_ind_start(10):SWR_ind_end(10)));
ylabel('lfp (V)'); xlabel('Time (s)');

%% Pick midpoint of SWR and plot fiber data 
% 1) 
SWR_ind_mid = (SWR_ind_start + SWR_ind_end)/2;  %middle index 

% parameters 
% 1 second before and after SWR event 
% single example

% seconds you want plotted: 
% 1000 Hz ... plot two seconds which is 1000 samples/s * 2 seconds = 2000
% samples in 2 seconds 
% samples 

% since you are using two different time series: 1) Find the index for the
% SWR, 2) then find the time for that, 3) then the time the corresponds in fiber
% data, 4) then the fiber index (same as time index).

% 2) & 3) & 4)
% using lfpTime 
% 1052 detected SWRs
SWR_time_mid = zeros(length(SWR_ind_mid),1); 
SWR_fiber_ind = zeros(length(SWR_ind_mid),1); 
for i = 1:1:size(SWR_ind_mid)
    SWR_time_mid(i) = lfp_time(round(SWR_ind_mid(i)));
    SWR_fiber_ind(i) = nearest_idx3(SWR_time_mid(i),time);
end

% event number 
enum = 1100 ; 
%enum100 is great

% reset time 
timeset = time((SWR_fiber_ind(enum)-10000):(SWR_fiber_ind(enum)+10000)); 

% SWR is 100 milliseconds so it won't show up with the fiber timescale...
% make a subplot 
figure(3) 
%yyaxis left
subplot(2,1,1)
plot((time((SWR_fiber_ind(enum)-10000):(SWR_fiber_ind(enum)+10000))-timeset(1)), (zdF((SWR_fiber_ind(enum)-10000):(SWR_fiber_ind(enum)+10000))), 'Color', [0 0.5 0])
hold on
xl = xline(4,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
xlabel('Time (s)')
xlim([0 8])
xticks([0 4 8])
xticklabels({'-4','0','4'})
ylabel('Signal (dF z-scored)')
title('Fiber Signal Centered on Detected SWR #900')
subplot(2,1,2)
plot((lfp_time(round(SWR_ind_start(enum)):round(SWR_ind_end(enum)))/1000), lfp(round(SWR_ind_start(enum)):round(SWR_ind_end(enum))));
ylabel('lfp (V)'); xlabel('Time (s)');

%% PETH of fiber signal following SWR 
% for all events 
figure(4) 


%%
figure(1)
for ptime = 1:1: length(photobeam_times)
    % tolerance 
    indxpb = find(abs(time-photobeam_times(ptime)) < 0.0001); % why does this keep chaning???
    init_trial = (indxpb - (8/0.0002)); % 8 seconds; remember 0.0002 is the time in seconds for each voltage point 
    end_trial = (indxpb + (8/0.0002));
    %if (max(left_end_peth.data*32) - min(left_end_peth.data*32)) < 0.15
    if abs(max(FP_detrended(init_trial:end_trial))) + abs(min(FP_detrended(init_trial:end_trial))) < 0.4

    plot((time(init_trial:end_trial))-time(init_trial), FP_detrended(init_trial:end_trial), 'Color', [jet(ptime,:)])
    else 
    disp('removed artifact')
    end
    hold on
colormap(jet);
cbh = colorbar;
cbh.Ticks = linspace (0,1,5); % 9 - need to change this based on length
cbh.TickLabels = num2cell(0:5:20); % 30; 40 - need to change this based on length
cbh.Label.String = 'Trial';
xlabel('time (s)')
xlim([0 16])
xticks([0 8 16])
xticklabels({'-8','0','8'})
ylabel('F')
title('F PET with 16s window')
end





