% SWR peaks vs. Frequency Bands and Movement Speed 

% 1) find SWRs 
% 2) find SWR-DA peaks within one second
% 3) fit a spline
% 4) run a spectrogram
% 5) save movement speed
% 6) run a correlation between signals 

%%  Load Data
clear; clc;
rng(pi)
cd 'F:\M533\M533_2024_08_20_recording2'; 
file_name = 'M533_2024_08_20'; 

FP_file=dir('*processed*');
FP = load(FP_file.name); % preprocessed data. 

SWR_file = dir('*detectedSWRs*');
load(SWR_file.name) % SWR intervals.

Track_file = dir('*_track.mat*'); 
load(Track_file.name); % for pseudo_outcomes

DLC_file = dir('*convertedDLC*'); 
P = readtable(DLC_file.name,'PreserveVariableNames',true);

% Display a message if everything loaded correctly 
% meaning the right number of files were found : 
if length(FP_file) == 1 && length(SWR_file)==1 && length(Track_file) ==1 && length(DLC_file)==1
    disp('loading looks good')
else
    disp('check files')
end

%% 1) Find SWRs 
% extract events (times of sleep sessions)
LoadExpKeys
cfg_evt = [];
evt2 = LoadEvents(cfg_evt);

% load raw fiber data 
%cfg_fiber.fc = {'CSC30.ncs'};
%raw_fiber = LoadCSC(cfg_fiber);
%raw_fiber_time = raw_fiber.tvec - raw_fiber.tvec(1); 

% extract LFP 
csc_name = [];
csc_name.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(csc_name); % csc with good ripples
% initialize LFP 
lfp_time = csc.tvec- csc.tvec(1); % lfp time 
lfp = csc.data; 

%time = FP.tvec; % fiber time - processed

% time of post recording 
post = ExpKeys.postrecord(1) - csc.tvec(1); % time of post sleep period, initialized  

% initialize SWR interval
SWR_start = evt.tstart- csc.tvec(1);
SWR_end = evt.tend- csc.tvec(1);
SWR_iv = [SWR_start SWR_end];
SWR_ind_start = nearest_idx3(SWR_iv(:,1),lfp_time); %find(abs(lfp-SWR_iv(1,1)) < 0.0005); % only for one but can we extend this to everything??
SWR_ind_end = nearest_idx3(SWR_iv(:,2),lfp_time); %find(abs(lfp-SWR_iv(1,2)) < 0.0005); % only for one but can we extend this to everything??

SWR_ind_mid = (SWR_ind_start + SWR_ind_end)/2;  %middle timepoint for each SWR

% find that corresponding SWR time index closest to that 
SWR_ind_start_post = nearest_idx3(post,SWR_iv(:,1)); %find(abs(lfp-SWR_iv(1,1)) < 0.0005); % only for one but can we extend this to everything??
SWR_ind_end_post = nearest_idx3(post,SWR_iv(:,2)); %find(abs(lfp-SWR_iv(1,2)) < 0.0005); % only for one but can we extend this to everything??
SWR_ind_mid_post = (SWR_ind_start_post + SWR_ind_end_post)/2;  %middle index 

% keep all SWR after that time -- should be 654
post_SWR_ind = SWR_ind_mid(SWR_ind_mid > SWR_ind_mid(round(SWR_ind_mid_post))); %SWR index that is greater than first post time point

% Pre track rest ^.^ 
pre_SWR_ind = SWR_ind_mid(SWR_ind_mid < SWR_ind_mid(round(SWR_ind_mid_post))); %all middle SWR timepoints that are less then the 

% for fiber after swr
SWR_time_mid = zeros(length(SWR_ind_mid),1); 
SWR_fiber_ind = zeros(length(SWR_ind_mid),1); 
for i = 1:1:size(SWR_ind_mid,1)
    SWR_time_mid(i) = lfp_time(round(SWR_ind_mid(i))); % lfp time in (s) for a swr
    SWR_fiber_ind(i) = nearest_idx3(SWR_time_mid(i),FP.tvec); % find the corresponidng time in fibFPer and saves the index.
end

pre_count = length(pre_SWR_ind);
post_count = length(post_SWR_ind);

prerecord_init = ExpKeys.prerecord(1)-csc.tvec(1);
prerecord_end = ExpKeys.prerecord(2)-csc.tvec(1);

postrecord_init = ExpKeys.postrecord(1)-csc.tvec(1);
postrecord_end = ExpKeys.postrecord(2)-csc.tvec(1);

%% initialize matrix 
% each swr has it's own row 
matrix_sess = array2table(zeros(length(post_SWR_ind),5),'VariableNames',{'swrID','TwosBeforePeak','TwosAfterPeak','TimeAfterPeak','peak_time'});

matrix_sess.('swrID')(:) = linspace(0,1,length(matrix_sess.('swrID')))';

%% FIND SWR-DA peaks within One Second. 
% dF from 2 seconds 
% for post session only. 

seconds = 2; 
samples = (seconds*FP.cfg.hdr{1,1}.SamplingFrequency)/2; % samples I will take before and after swrs

x1 = 1:1:1000;%1001:1:3000; %1:1:2000;% % two seconds before for preproc data
x2 = 1001:1:2001; %3001:1:5000; %2001:1:4000;%  % two seconds after for preproc data

for i = 1:height(matrix_sess) % iterate through each swr. 
    swr_time = lfp_time(round(post_SWR_ind(i))); % swr lfp time  
    fiber_index = nearest_idx3(swr_time, FP.tvec); % fiber index closest to middle swr_time 
    signal = FP.zF_win_60s(fiber_index-samples:fiber_index+samples); 
    matrix_sess.('OnesBeforePeak')(i) = max(signal(x1)); %-min(signal(x1)); % maybe the average signal might be better than the lowest signal?? 
    [matrix_sess.('OnesAfterPeak')(i),I] = max(signal(x2)); %-min(signal(x2)); 
    matrix_sess.('TimeAfterPeak')(i) = FP.tvec(I); % time of the peak post swr OK THIS IS WRONG... I NEED TO TAKE IT OF SIGNAL and then 
    matrix_sess.('peak_time')(i) = FP.tvec(fiber_index+I); % time of the peak post swr
end

%% Spectrogram
% extract behavioral events
LoadExpKeys
cfg_evt = [];
evt2 = LoadEvents(cfg_evt);

% extract LFP with SWRs
csc_name = [];
csc_name.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(csc_name);
% initialize LFP 
lfp_time = csc.tvec- csc.tvec(1);
lfp = csc.data;
FS = csc.cfg.hdr{1}.SamplingFrequency; % set FP_data.acq.Fs to sampling frequency rate (5000 points per second) 
% zscore signal
zlfp = zscore_tsd(csc);


% time of post recording 
post = nearest_idx3(ExpKeys.postrecord(1) - csc.tvec(1),FP.tvec); % time of post sleep period, initialized  
post_end = nearest_idx3(ExpKeys.postrecord(2) - csc.tvec(1),FP.tvec); % time of post sleep period, initialized  

% initialize SWR interval
SWR_start = evt.tstart- csc.tvec(1);
SWR_end = evt.tend- csc.tvec(1);
SWR_iv = [SWR_start SWR_end];
SWR_ind_start = nearest_idx3(SWR_iv(:,1),lfp_time);
SWR_ind_end = nearest_idx3(SWR_iv(:,2),lfp_time); 
SWR_ind_mid = (SWR_ind_start + SWR_ind_end)/2;  %middle timepoint for each SWR

% find that corresponding SWR time index closest to that 
SWR_ind_start_post = nearest_idx3(post,SWR_iv(:,1));
SWR_ind_end_post = nearest_idx3(post,SWR_iv(:,2)); 

% I need to make a new FP structure where data is zdF win so that it can be
% restricted correctly 
FP_win = FP;
FP_win.data = FP.zdF_win; 
FP_restrict_win = restrict(FP_win, FP.tvec(post), FP.tvec(post_end)); % does not change the tvec data, just the data. '
FP_restrict_win.tvec = FP_restrict_win.tvec - FP_restrict_win.tvec(1);

% restrict to only post session 
CSC_restrict = restrict(zlfp, ExpKeys.postrecord(1), ExpKeys.postrecord(2));
FP_restrict = restrict(FP, FP.tvec(post), FP.tvec(post_end)); % does not change the tvec data, just the data. '
FP_restrict = zscore_tsd(FP_restrict);

win_size = 5 * FS; % 10-second moving window
overlap =round(win_size * 0.95); % 90% overlap
freq_range = 1:250; % Frequency range

% Define custom frequency ranges
low_freq_range = linspace(1, 30, 250);   % More bins for low frequencies (1-30 Hz) 900 events. so need at least 450?
mid_freq_range = linspace(30, 100, 100); % Fewer bins for mid frequencies (30-100 Hz) 5000
high_freq_range = linspace(100, 250, 10); % Fewest bins for high frequencies (100-250 Hz)

% Combine all frequency ranges into a single frequency vector
freq_bins = [low_freq_range, mid_freq_range, high_freq_range];

%freq_bins = logspace(log10(1), log10(250), 100); % Log-spaced bins

[S, F, T, P] = spectrogram(CSC_restrict.data, hanning(win_size), overlap, [1:250], FS);
% instead of 1:250 you can use freq_bins 

%% mean power band over spline times 

n_peaks = length(matrix_sess.TimeAfterPeak);
peak_photometry = matrix_sess.OnesAfterPeak;

% Define frequency bands
delta_band = [2 5];
theta_band = [6 10];
beta_band = [12 35];
low_gamma_band = [35 70];
high_gamma_band = [70 100];
swr_band   = [140 250];

% Initialize arrays
delta_power = nan(n_peaks, 1);
theta_power = nan(n_peaks, 1);
beta_power = nan(n_peaks,1);
swr_power   = nan(n_peaks, 1);

for i = 1:n_peaks
    t_peak = post_time_sec(i);

    % Find index in T closest to this time
    [~, idx_time] = min(abs(T - t_peak));

    % Get power in each band at this time point
    delta_idx = F >= delta_band(1) & F <= delta_band(2);
    theta_idx = F >= theta_band(1) & F <= theta_band(2);
    beta_idx
    swr_idx   = F >= swr_band(1)   & F <= swr_band(2);

    % Take average over frequency band, at this time point
    delta_power(i) = mean(S(delta_idx, idx_time), 'omitnan');
    theta_power(i) = mean(S(theta_idx, idx_time), 'omitnan');
    swr_power(i)   = mean(S(swr_idx, idx_time), 'omitnan');
end
%% max lags
max_lag = 2;  % Allow at least ±6 seconds to capture meaningful shifts
time_step = mean(diff(T));  % Spectrogram time step (e.g., 3s)
max_lag_samples = round(max_lag / time_step);

%% xcorr
%[xc_delta, lags] = xcorr(delta_power - mean(delta_power), FP_interp - mean(FP_interp), max_lag_samples, 'coeff');


%%
figure;

subplot(1,3,1)
[r, p] = xcorr(real(delta_power)-mean(real(delta_power)), peak_photometry - mean(peak_photometry),max_lags_samples,'coeff');
plot(p, r, 'g'); title('2-5 Hz LFP Power'); xlabel('Fiber Lag (s)'); ylabel('Cross-Correlation (R-value)');

scatter(real(delta_power), peak_photometry, 'filled');
xlabel('Delta Power (2-5 Hz)');
ylabel('Peak Fiber Value');
title(sprintf('Delta: r=%.2f, p=%.3f', r, p));

subplot(1,3,2)
[r, p] = xcorr(real(theta_power), peak_photometry, 'rows', 'complete');
scatter(real(theta_power), peak_photometry, 'filled');
xlabel('Theta Power (6-10 Hz)');
title(sprintf('Theta: r=%.2f, p=%.3f', r, p));

subplot(1,3,3)
[r, p] = xcorr(real(swr_power), peak_photometry, 'rows', 'complete');
scatter(real(swr_power), peak_photometry, 'filled');
xlabel('SWR Power (140-250 Hz)');
title(sprintf('SWR: r=%.2f, p=%.3f', r, p));


%% 
% time initialized 
post_time_sec = matrix_sess.peak_time- matrix_sess.peak_time(1);
spline_x = linspace(0.01,post_time_sec(end),length(post_time_sec)); 

% remove any non unique points?? 
[~,ia,ic] = unique(peak_photometry, 'rows','stable');          % Unique Elements
v = accumarray(ic, 1);                  % Tally Occurrences Of Rows
B = peak_photometry(ia(v==1),:) ;                      % Keep Rows That Only Appear Once
B_time_spline = spline_x(ia(v==1,:));
post_time_sec_B = post_time_sec(ia(v==1,:));


s = interp1(post_time_sec_B, B,B_time_spline,'spline');

figure (3)
plot(post_time_sec, peak_photometry) 
%xq2 = 0:0.01:15;
%s = spline(x,y,xq2);
hold on
scatter(post_time_sec , peak_photometry, '*')
plot(B_time_spline,s,'--')
title('Peak SWR-DA Values (1 sec window)')
xlabel(['time (s)']);
ylabel(['SWR-[DA] Peaks (z-score)']);
