% code to see if there is any predictive value in the LFP signal
% frequencies and the dopamine signal. 

%% load data 
clear; clc; 
%rng(pi);

cd 'F:\M548\M548_2024_08_25_recording1';
FP = load('M548_2024_08_25processed.mat');

% extract SWR intervals
load('M548_2024_08_25detectedSWRs.mat')

addpath('C:\Users\mimia\Documents\Toolboxes\shadedErrorBar')

%% extract events
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

% Fiber time
time = FP.tvec; % fiber time (previously initialized)
FP.data = FP.data' ;

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

% restrict to only post session 
CSC_restrict = restrict(zlfp, ExpKeys.postrecord(1), ExpKeys.postrecord(2));
FP_restrict = restrict(FP, FP.tvec(post), FP.tvec(post_end)); % does not change the tvec data, just the data. '
FP_restrict = zscore_tsd(FP_restrict);
% I need to make a new FP structure where data is zdF win so that it can be
% restricted correctly 
FP_win = FP;
FP_win.data = FP.zdF_win; 
FP_restrict_win = restrict(FP_win, FP.tvec(post), FP.tvec(post_end)); % does not change the tvec data, just the data. '
FP_restrict_win.tvec = FP_restrict_win.tvec - FP_restrict_win.tvec(1);

%% fourier transform to power spectral densities 
% 2 second moving window 
% for frequencies 1 - 300 

%% from chat 
%% Define Parameters
win_size = 10 * FS; % 10-second moving window
overlap =round(win_size * 0.95); % 90% overlap
freq_range = 1:250; % Frequency range

% Define custom frequency ranges
low_freq_range = linspace(1, 30, 250);   % More bins for low frequencies (1-30 Hz) 900 events. so need at least 450?
mid_freq_range = linspace(30, 100, 100); % Fewer bins for mid frequencies (30-100 Hz) 5000
high_freq_range = linspace(100, 250, 10); % Fewest bins for high frequencies (100-250 Hz)

% Combine all frequency ranges into a single frequency vector
freq_bins = [low_freq_range, mid_freq_range, high_freq_range];

%freq_bins = logspace(log10(1), log10(250), 100); % Log-spaced bins

[S, F, T, P] = spectrogram(CSC_restrict.data, hanning(win_size), overlap, freq_bins, FS);

% Normalize power across frequencies (to mitigate bias towards lower frequencies)
%P_normalized = P ./ max(P, [], 1);  % Normalize each time point by its max value
%P_log = 10 * log10(P_normalized);   % Log-transform to dB


% High-pass filter (optional, if low-frequency drift is a problem)
%[b, a] = butter(4, 1 / (FS / 2), 'high');  % 0.5 Hz high-pass filter
%filtered_data = filtfilt(b, a, CSC_restrict.data);
% Recompute spectrogram for filtered data
%[S_filtered, F_filtered, T_filtered, P_filtered] = spectrogram(filtered_data, hanning(win_size), overlap, freq_bins, FS);
%P_filtered_log = 10 * log10(P_filtered);

% Compute spectrogram with a Hanning window
%[S, F, T, P] = spectrogram(CSC_restrict.data, hanning(win_size), overlap, freq_range, FS);
% where: 
% S is short-time Fourier transform returned as a matrix. Time increases
% across teh columns of s and frequency increases down the rows starting
% from zero. 
% F is cyclical frequencies returned as a vector. f has a length equal to
% the number of rows of s
% T is the midpoint of each time segment 
% P is the power spectral density 

%col = [0.09*10^-10 3*10^-10]; % some arbitrary range for the color scaling
                
% Log-transform power for normalization
P_log = 10 * log10(P);  

% Align time bins with fiber photometry signal
%T_aligned = T + CSC_restrict.tvec(1); % Adjust spectrogram time to match actual recording time
FP_interp = interp1(FP_restrict_win.tvec, FP_restrict_win.data, T, 'linear', 'extrap'); % Interpolate FP signal

% spectrogram function post rest
%[S_post,F_post,T_post,P_post] = spectrogram(zlfpr_post.data,hanning(7500),3750,1:20,csc.cfg.hdr{1,1}.SamplingFrequency);

figure(2);
imagesc(T, F, P_log);  % Log-transformed spectrogram (normalized)
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
%colormap jet;
colorbar;
%caxis([prctile(P_filtered_log(:), 5), prctile(P_filtered_log(:), 95)]); % Better contrast
title('Spectrogram of LFP Signal');

%% Regression: Predict FP from LFP Power
X = P'; % Predictor: Power across frequencies
y = FP_interp; % Response: Fiber photometry signal
mdl = fitlm(X, y); % Linear regression model
figure(3);
plot(mdl)
title('linear regression between freq power & fiber signal')
xlabel('log transform power across frequencies')
ylabel('fiber voltage')

%% Plot results
figure;
subplot(2,1,1);
imagesc(T, F, P_log);
axis xy; colorbar;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Spectrogram of LFP');

subplot(2,1,2);
plot(T, FP_interp, 'r');
xlabel('Time (s)'); ylabel('Fiber Signal (V)');
title('Fiber Signal');

disp(mdl)

%%
%% Plot Spectrogram
figure;
subplot(2,1,1);
imagesc(t_aligned, F, P_log);
axis xy; colorbar;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Spectrogram with 2s Hanning Window');

subplot(2,1,2);
plot(t_aligned, FP_interp, 'r');
xlabel('Time (s)'); ylabel('Dopamine Signal');
title('Interpolated Dopamine Signal');

%% Regression: Predict FP from LFP Power
X = P_log'; % Predictor: Power across frequencies
y = FP_interp; % Response: Fiber photometry signal
mdl = fitlm(X, y); % Linear regression model

disp(mdl)

        
%% windowed prediction
% Initialize an array to store predicted fiber signal values
predicted_FP = zeros(size(FP_interp)); 

% Loop through each time point and make predictions using the model
for t = 1:length(T)
    % Extract the frequency power at time point t (all frequency bins for that time point)
    X_t = P(:, t)';  % Transpose to make it a row vector of power values at time t
    
    % Predict the fiber signal at time t using the regression model
    predicted_FP(t) = predict(mdl, X_t);  % Use the linear regression model for prediction
end

% Plot the results
figure;
subplot(2,1,1);
imagesc(T, F, P_log);  % Log-transformed spectrogram (normalized)
axis xy; colorbar;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Spectrogram of LFP');

subplot(2,1,2);
plot(T, FP_interp, 'r'); hold on;
plot(T, predicted_FP, 'b--');  % Plot the predicted fiber signal
xlabel('Time (s)'); ylabel('Fiber Signal (V)');
title('Fiber Signal and Predicted Fiber Signal');
legend('Actual Fiber Signal', 'Predicted Fiber Signal');

%%
% Initialize an array to store the correlation values
corr_values = zeros(length(F), length(T));  % Dimensions: frequencies x time points

% Loop through each frequency and calculate the correlation with the fiber signal
for f = 1:length(F)
    for t = 1:length(T)
        % Extract the frequency power at time point t (for frequency f)
        power_at_t = P(f, t);  % Power at frequency f and time t

        % Compute correlation between fiber signal and power at frequency f
        corr_values(f, t) = corr(power_at_t, FP_interp(t), 'Type', 'Pearson');  % Pearson correlation
    end
end

% Plot correlation matrix over time for each frequency
figure;
imagesc(T, F, corr_values);  % Plot correlation over time and frequency
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
caxis([-1, 1]);  % Set correlation range from -1 to 1
colorbar;
colormap('jet');  % Use 'jet' colormap for better contrast
title('Correlation between LFP Power and Fiber Signal');


%% 
% Define the frequency ranges and bin sizes
low_freq_range = linspace(1, 30, 250);   % More bins for low frequencies (1-30 Hz)
mid_freq_range = linspace(30, 100, 100); % Fewer bins for mid frequencies (30-100 Hz)
high_freq_range = linspace(100, 250, 10); % Fewest bins for high frequencies (100-250 Hz)

% Combine all frequency ranges into a single frequency vector
freq_bins = [low_freq_range, mid_freq_range, high_freq_range];

% Create labels for each frequency bin (Low, Mid, High)
freq_labels = cell(length(freq_bins), 1);
freq_labels(1:length(low_freq_range)) = {'Low'};
freq_labels(length(low_freq_range)+1:length(low_freq_range)+length(mid_freq_range)) = {'Mid'};
freq_labels(length(low_freq_range)+length(mid_freq_range)+1:end) = {'High'};

% P is the power matrix [n_freqs x n_timepoints]
% Transpose P to have time points as rows and frequencies as columns
P_transposed = P';

% Compute the correlation matrix between frequency bands
corr_matrix = corrcoef(P_transposed); 

% Plot the correlation matrix as a heatmap
figure;
imagesc(corr_matrix);
colorbar;
colormap jet;
title('Correlation between Frequency Bands');
xlabel('Frequency Bin');
ylabel('Frequency Bin');
axis tight;

%%
% Make sure P is in the correct shape: [n_freqs x n_timepoints]
% P_log is already the log-transformed power spectrogram
P_log = 10 * log10(P);  % Log-transform the power for normalization

% Align time points for fiber photometry signal
T_aligned = T;  % Assuming T is already aligned to the right time points

% Interpolate fiber photometry data to match the spectrogram time points
FP_interp = interp1(FP_restrict_win.tvec, FP_restrict_win.data, T_aligned, 'linear', 'extrap');

% Run the correlation for each frequency band with the fiber signal
% We will calculate the correlation for each time window across the frequencies
corr_values = nan(size(P, 1), 1);  % Placeholder for correlation values (one for each frequency bin)

% Compute correlation for each frequency band
for i = 1:size(P_log, 1)  % Loop over frequency bins
    corr_values(i) = corr(P(i, :)', FP_interp');  % Correlate frequency power with fiber signal
end

%plot(corrcoef(corr_values))
% Plot the correlation heatmap
figure;
imagesc(corr_values);  % Display correlation values
colorbar;
colormap jet;
xlabel('Frequency Bins');
ylabel('Correlation with Fiber Signal');
title('Correlation between Frequency Power and Fiber Photometry Signal');

%% start with a few frequency bands 
% fiber data oscillates at 0.
% delta 2-5 
% theta 6-10 
% swr 140-250
% and compute cross correlations (time axis from -0.5 to 0.5 offsets for
% each) 
% Define frequency bands

% Extract frequency bands
delta_idx = F >= 2 & F <= 5;
theta_idx = F >= 6 & F <= 10;
swr_idx = F >= 140 & F <= 250;

% Compute mean power over each band (aligned with T)
% this gives you the mean power over T
delta_power = mean(P_log(delta_idx, :), 1);
theta_power = mean(P_log(theta_idx, :), 1);
swr_power = mean(P_log(swr_idx, :), 1);

% Define cross-correlation time window (-0.5 to 0.5 sec)
max_lag = 2;  % Allow at least ±6 seconds to capture meaningful shifts
time_step = mean(diff(T));  % Spectrogram time step (e.g., 3s)
max_lag_samples = round(max_lag / time_step);
% to do -0.5 to 0.5 seconds, i need to increase my temporal resolution. 
% try 10*fs next time with an overlap of 95%

%% make signal that oscillates at 1 Hz and interp the same as fiber data 
sim = sin(1:length(FP_restrict_win.tvec));
FP_sim = interp1(FP_restrict_win.tvec, sim, T, 'linear', 'extrap'); % Interpolate FP signal

% %% xcorr 2 way 
% XC_delta = xcorr2(delta_power, FP_interp);
% XC_theta = xcorr2(theta_power, FP_interp);
% XC_swr   = xcorr2(swr_power, FP_interp);
% % lags 
% max_lag = 6;  % Choose a reasonable lag window (e.g., ±6s)
% time_step = mean(diff(T));  % Spectrogram time step (e.g., 3s)
% max_lag_samples = round(max_lag / time_step);
% 
% lags = (-max_lag_samples:max_lag_samples) * time_step;  % Convert to seconds
% % confused 
% 
% %% xcorr2 plot
% figure;
% subplot(3,1,1);
% plot(lags, XC_delta);
% xlabel('Lag (s)'); ylabel('XCorr');
% title('Delta Band Cross-Correlation');
% 
% subplot(3,1,2);
% plot(lags, XC_theta);
% xlabel('Lag (s)'); ylabel('XCorr');
% title('Theta Band Cross-Correlation');
% 
% subplot(3,1,3);
% plot(lags, XC_swr);
% xlabel('Lag (s)'); ylabel('XCorr');
% title('SWR Band Cross-Correlation');


%% xcorr way
% Compute cross-correlations
[xc_delta, lags] = xcorr(delta_power - mean(delta_power), FP_interp - mean(FP_interp), max_lag_samples, 'coeff');
[xc_theta, ~] = xcorr(theta_power - mean(theta_power), FP_interp - mean(FP_interp), max_lag_samples, 'coeff');
[xc_swr, ~] = xcorr(swr_power - mean(swr_power), FP_interp - mean(FP_interp), max_lag_samples, 'coeff');
[xc_sine, ~] = xcorr(FP_sim - mean(FP_sim), FP_interp - mean(FP_interp), max_lag_samples, 'coeff');

% Convert lags to time
lag_time = lags * mean(diff(T));  % Use spectrogram time step

% Plot results
figure (2);
subplot(4,1,1);
plot(lag_time, xc_delta, 'g'); title('Delta Band x Fiber Cross-Correlation'); xlabel('Lag (s)'); ylabel('Correlation');
ylim([0 0.6]);
subplot(4,1,2);
plot(lag_time, xc_theta, 'b'); title('Theta Band x Fiber Cross-Correlation'); xlabel('Lag (s)'); ylabel('Correlation');
ylim([0 0.6]);
subplot(4,1,3);
plot(lag_time, xc_swr, 'r'); title('SWR Band x Fiber Cross-Correlation'); xlabel('Lag (s)'); ylabel('Correlation');
ylim([0 0.6]);
subplot(4,1,4);
plot(lag_time, xc_sine, 'k'); title('Sine x Fiber Cross-Correlation'); xlabel('Lag (s)'); ylabel('Correlation');
ylim([0 0.6]);