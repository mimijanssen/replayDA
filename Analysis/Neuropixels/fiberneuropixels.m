%% PLOTTING FIBER NEUROPIXELS

load ('all_events.mat')
load('dorsal_SWR.mat')
load('ventral_SWR.mat')
load('XA7_synced.mat')

addpath('C:\Users\mimia\Documents\Toolboxes\shadedErrorBar')


%% fiber restructure
% rename variables and store in a tsd struct
FP = tsd;
FP.data = data; %voltage multiplier
FP.tvec = tvec; 
FP.cfg.hdr{1}.SamplingFrequency = fs;
%% plot
sessionTitle = 'CW_';
last_time = length(FP.tvec); %the value of the last time point is how many seconds the recording was
timerange1 = 10*FP.cfg.hdr{1}.SamplingFrequency; 
timerange2= 100*FP.cfg.hdr{1}.SamplingFrequency; 
timerange3= length(FP.tvec); 
time_ranges = [timerange1, timerange2, timerange3];

figure(1)
for t_i = 1:length(time_ranges)
    t_range = 1:time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(FP.tvec(t_range), FP.data(t_range), 'Color', [0 0.5 0])
    %title([num2str(time_ranges(t_i)),' samples'], 'Interpreter','none')
    ylabel('Fiber Signal (V)'); xlabel('Time (s)');
end


%% filtered

% Median filter
FP_denoised = medfilt1(FP.data);

figure(2)
plot(FP.tvec,FP.data);
hold on
plot(FP.tvec,FP_denoised);
title('Median Filtered FP Signal')
ylabel('Signal (V)')
xlabel('Time (s)')

% Butterworth filter
fc = 20; % cut off frequency
fs = FP.cfg.hdr{1}.SamplingFrequency; % sampling rate 
[b,a] = butter(8,fc/(fs/2)); % 8th order filter
FP.FP_denoised = filtfilt(b, a, FP_denoised);

figure(3)
plot(FP.tvec,FP.data);
hold on
plot(FP.tvec,FP.FP_denoised)
hold off
title('20 Hz Butterworth Filtered FP Signal over Raw Data')
ylabel('Signal (V)')
xlabel('Time (s)')
legend('raw','butterworth','Location','northeast')

% pwelch to check if butterworth worked: 
% pwelch of raw data
wsize = FP.cfg.hdr{1}.SamplingFrequency * 4; % sampling frequency * 4s is the amount of samples (the window size needed)
[Pxx,F] = pwelch(FP.data,hanning(wsize),wsize/2,[],FP.cfg.hdr{1}.SamplingFrequency);
figure(4)
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 100]);
title('PSD for Raw Fiber')

[Pxx_butter,F_butter] =  pwelch(FP.FP_denoised,hanning(wsize),wsize/2,[],FP.cfg.hdr{1}.SamplingFrequency);
figure(5)
plot(F_butter,10*log10(Pxx_butter),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 100]);
title('PSD for Denoised Fiber')
%% Local Detrend: Windowed Detrend using Locdetrend
% another way to detrend the signal is a simple linear detrend. 
addpath('C:\Users\mimia\Documents\GitHub\vandermeerlab-replay-da\chronux-master\spectral_analysis\helper')
addpath('C:\Users\mimia\Documents\GitHub\vandermeerlab-replay-da\chronux-master\spectral_analysis\continuous')

% 1 min detrend, 1 sample stepsize  
FP_detrend_60s = locdetrend(FP.FP_denoised,FP.cfg.hdr{1}.SamplingFrequency,[60 1]);
% 10 s detrend, 1 sample stepsize 
FP_detrend_10s = locdetrend(FP.FP_denoised,FP.cfg.hdr{1}.SamplingFrequency,[10 1]); 

% 1 min detrend, almost no overlap
% if your step size is the size of your window, then you will have no
% overlap
FP_detrend_60s_no = locdetrend(FP.FP_denoised,FP.cfg.hdr{1}.SamplingFrequency,[61 60]);
% There is slight overlap on the edges since I had a edge artifact when I
% did 60, 60. 

FP.detrend_60s = FP_detrend_60s; 
FP.detrend_10s = FP_detrend_10s;
FP.detrend_60s_no = FP_detrend_60s; 

% detrended and filtered signal (using 10s window for detrend)
figure(8)
plot(FP.tvec, FP.FP_denoised)
hold on
plot(FP.tvec,FP_detrend_60s_no)
title('Detrended and Filtered FP Signal (60s)')
ylabel('Signal (V)')
xlabel('Time (s)')

%% Normalization for Windowed Detrend (locdetrend)
% z-scored locdetrended signal
%zF_win_60s_no = (FP_detrend_60s_no - mean(FP_detrend_60s_no))./std(FP_detrend_60s_no);
%zF_win_60s = (FP_detrend_60s - mean(FP_detrend_60s))./std(FP_detrend_60s); % you don't need to apply the . becase it should find the standard deviation for everything... try just doing zscore... and see if you get a different output
%zF_win_10s = (FP_detrend_10s - mean(FP_detrend_10s))./std(FP_detrend_10s);
t = FP.tvec; % added

zF_win_60s_no = zscore(FP_detrend_60s_no);
zF_win_60s = zscore(FP_detrend_60s);
zF_win_10s = zscore(FP_detrend_10s);

FP.zF_win_60s = zF_win_60s; 
FP.zF_win_10s = zF_win_10s;
FP.zF_win_60s_no = zF_win_60s_no; 

% filtered, detrended, and normalized signal (using 60s window for detrend)
figure(9)
subplot(2,1,1)
plot(t,zF_win_60s)
title('Filtered, Detrended, and Normalized FP Signal (60s)')
ylabel('Signal z-scored (V)')
xlabel('Time (s)')
subplot(2,1,2)
plot(t,zF_win_60s_no)
title('Filtered, Detrended, Normalized FP Signal (60s, almost no overlap)')
ylabel('Signal z-scored (V)')
xlabel('Time (s)')

%%
%filename = append(file_name, "processed.mat");
%save(filename, '-struct','FP')

%% 
seconds = 10; % I want to save 5 seconds now. 55555 but am i saving that on the x? 
samples = (seconds*FP.cfg.hdr{1,1}.SamplingFrequency)/2;

% find lfp time index and fiber signal index 
SWR_fiber_ind_post = zeros(length(dorsal_swr.iv.tstart),1); 
for i = 1:1:length(dorsal_swr.iv.tstart)
    SWR_fiber_ind_post(i) = nearest_idx3(dorsal_swr.iv.tstart(i),tvec);
end

SWR_fiber_ind_pre = zeros(length(ventral_swr.iv.tstart),1); 
for i = 1:1:length(ventral_swr.iv.tstart)
    SWR_fiber_ind_pre(i) = nearest_idx3(ventral_swr.iv.tstart(i),tvec);
end
%% dorsal
zdF_extract_post = zeros(length(dorsal_swr.iv.tstart)-1, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 

for ievt = 1:1:length(dorsal_swr.iv.tstart)-1 % -2 for M654 track 1 ; -3 for M600 recording 1; -3 for M646 recording 7; and M654 recording 6 -4 for M646 recording 4; -1 for everyone else. !!! start from 3 for recording 2 M654
    zdF_extract_post(ievt,:) = FP.zF_win_60s((SWR_fiber_ind_post(ievt)-samples):(SWR_fiber_ind_post(ievt)+samples));
end

avg_fiber_dorsal = nanmean(zdF_extract_post);
std_fiber_dorsal = 2*std(zdF_extract_post);

X = FP.zF_win_60s; 
N = 1000; % number of circshifts 
K=randi([1 length(FP.zF_win_60s)],1, N); % pick a random number between 1 and number of samples ... 100 times 
events_num = length(dorsal_swr.iv.tstart)-1;

circ_zdF_extract = zeros(events_num, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
avg_circ_zdF_extract_post = zeros(N,seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for iter_circ = 1:1:N 
    Y = circshift(X,K(iter_circ)); 
    for ievt = 1:1:length(dorsal_swr.iv.tstart)-1 % ALSO CHANGE HERE!!!!!
    % find time for x axis
       circ_zdF_extract(ievt,:) = (Y((SWR_fiber_ind_post(ievt)-samples):(SWR_fiber_ind_post(ievt)+samples))); %resets every time
    end
    avg_circ_zdF_extract_post(iter_circ,:) = nanmean(circ_zdF_extract); % this is a mean over all the trials .... 
end

circ_avg_fiber_dorsal = nanmean(avg_circ_zdF_extract_post);
circ_std_fiber_dorsal = 2*std(avg_circ_zdF_extract_post);


%% ventral
zdF_extract_pre = zeros(length(ventral_swr.iv.tstart)-1, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 

for ievt = 1:1:length(ventral_swr.iv.tstart)-1 % -2 for M654 track 1 ; -3 for M600 recording 1; -3 for M646 recording 7; and M654 recording 6 -4 for M646 recording 4; -1 for everyone else. !!! start from 3 for recording 2 M654
    zdF_extract_pre(ievt,:) = FP.zF_win_60s((SWR_fiber_ind_pre(ievt)-samples):(SWR_fiber_ind_pre(ievt)+samples));
end

avg_fiber_ventral = nanmean(zdF_extract_pre);
std_fiber_ventral = 2*std(zdF_extract_pre);

X = FP.zF_win_60s; 
N = 1000; % number of circshifts 
K=randi([1 length(FP.zF_win_60s)],1, N); % pick a random number between 1 and number of samples ... 100 times 
events_num = length(ventral_swr.iv.tstart)-1;

circ_zdF_extract = zeros(events_num, seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1); 
avg_circ_zdF_extract_pre = zeros(N,seconds* FP.cfg.hdr{1,1}.SamplingFrequency+1);

for iter_circ = 1:1:N 
    Y = circshift(X,K(iter_circ)); 
    for ievt = 1:1:length(ventral_swr.iv.tstart)-1 % ALSO CHANGE HERE!!!!!
    % find time for x axis
       circ_zdF_extract(ievt,:) = (Y((SWR_fiber_ind_pre(ievt)-samples):(SWR_fiber_ind_pre(ievt)+samples))); %resets every time
    end
    avg_circ_zdF_extract_pre(iter_circ,:) = nanmean(circ_zdF_extract); % this is a mean over all the trials .... 
end

circ_avg_fiber_ventral = nanmean(avg_circ_zdF_extract_pre);
circ_std_fiber_ventral = 2*std(avg_circ_zdF_extract_pre);

time_extract_pre = linspace(-5,5,10001);

%%
figure (1)
shadedErrorBar(time_extract_pre(1,:),circ_avg_fiber_dorsal,circ_std_fiber_dorsal,'lineProps','-k','transparent',1) % subtracted the circ mean here
hold on
plot(time_extract_pre(1,:),circ_avg_fiber_dorsal,'LineWidth',2,'Color','k') %subtracted the circ mean here
hold on
plot(time_extract_pre(1,:),avg_fiber_dorsal,'LineWidth',2)
hold on
xl = xline(5,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
ylim([-1.5 1.5])
xticks([-5 0 5])
xticklabels({'-5','0','5'})
title('Dorsal SWR PETH','FontSize', 20)
ylabel('Mean [DA] (z-score)','FontSize', 16)
xlabel('Time from SWR (s)','FontSize', 16)
legend('','shuffle','signal','Location','northwest')
legend boxoff 

%% 
figure (2)
shadedErrorBar(time_extract_pre(1,:),circ_avg_fiber_ventral,circ_std_fiber_ventral,'lineProps','-k','transparent',1) % subtracted the circ mean here
hold on
plot(time_extract_pre(1,:),circ_avg_fiber_ventral,'LineWidth',2,'Color','k') %subtracted the circ mean here
hold on
plot(time_extract_pre(1,:),avg_fiber_ventral,'LineWidth',2)
hold on
xl = xline(5,'-',{'SWR'});
xl.LabelVerticalAlignment = 'top';
%hold off
ylim([-1.5 1.5])
xticks([-5 0 5])
xticklabels({'-5','0','5'})
title('Ventral SWR PETH','FontSize', 20)
ylabel('Mean [DA] (z-score)','FontSize', 16)
xlabel('Time from SWR (s)','FontSize', 16)
legend('','shuffle','signal','Location','northwest')
legend boxoff 


%% 

win = [-2 2]; % peri-event window in seconds
fs = 1/median(diff(tvec)); % sampling rate
win_idx = round(win * fs); % convert to indices
t_axis = linspace(win(1), win(2), diff(win_idx)+1);

figure;

for e = 1:length(evt.t)

    event_times = evt.t{1,e};
    segments = [];

    for i = 1:length(event_times)

        % Find closest index in signal
        idx = nearest_idx3(event_times(i), tvec);

        % Define window
        idx_win = (idx + win_idx(1)):(idx + win_idx(2));

        % Skip if out of bounds
        if idx_win(1) < 1 || idx_win(end) > length(FP.zF_win_60s)
            continue
        end

        % Extract segment
        segments(i,:) = FP.zF_win_60s(idx_win);

    end

    % Remove empty rows (in case some skipped)
    segments = segments(~all(segments==0,2),:);

    % Compute mean and SEM
    mean_trace = mean(segments,1,'omitnan');
    sem_trace = std(segments,[],1,'omitnan') / sqrt(size(segments,1));

    % Plot
    subplot(3,3,e)
    hold on
    fill([t_axis fliplr(t_axis)], ...
         [mean_trace+sem_trace fliplr(mean_trace-sem_trace)], ...
         [0.8 0.8 1], 'EdgeColor','none');

    plot(t_axis, mean_trace, 'b', 'LineWidth', 2);
    xline(0,'k--');

    title(['Event ' num2str(e) ' (n=' num2str(size(segments,1)) ')'])
    xlabel('Time (s)')
    ylabel('zF')
end
