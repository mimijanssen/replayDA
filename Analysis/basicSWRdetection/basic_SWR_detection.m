%% basic SWR detection
% 1) restrict to periods of rest (pre- and post- track)
% 2) filter based on ripple band (120-250 Hz)
% 3) Hilbert transform and calculate mean power 
% 4) Keep anything with power above 1 z-score 
% 5) Optional: Keep anything above 1 z-score from a control channel with no
% ripples. 
% optional : take the difference from a control channel and save ripples
% that are unique to the ripple channel (avoids movement and other
% artifacts)
clear; clc; 

cd('F:\M548\M548_2024_08_31_recording7')
file_name = 'M548_2024_08_31';
folderPath = ('C:\Users\mimia\Documents\ReplayDA Figures\M548\recording7\basic');

%% Load LFP and metadata
LoadExpKeys;
LoadMetadata;

cfg = [];
cfg.fc = ExpKeys.goodSWR(1);
SWR_lfp_raw = LoadCSC(cfg);

cfg_temp = [];
cfg_temp.getRatings = 0;
cfg_temp.load_questionable_cells = 1;
cfg_temp.verbose = 1;
cfg_temp.fc = ExpKeys.goodSpikes;
S = LoadNST(cfg_temp);

%% Restrict to sleep epochs only
SWR_lfp_raw = restrict(SWR_lfp_raw, ...
    [ExpKeys.prerecord(1) ExpKeys.postrecord(1)], ...
    [ExpKeys.prerecord(2) ExpKeys.postrecord(2)]);
S = restrict(S, ...
    [ExpKeys.prerecord(1) ExpKeys.postrecord(1)], ...
    [ExpKeys.prerecord(2) ExpKeys.postrecord(2)]);

%% PSD check
Fs = SWR_lfp_raw.cfg.hdr{1,1}.SamplingFrequency;
wsize = round(Fs * 0.5);
[Pxx, F] = pwelch(SWR_lfp_raw.data, hanning(wsize), wsize/2, [], Fs);

figure;
plot(F, 10*log10(Pxx), 'k');
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([50 350]);
title('PSD - confirm ripple band power');
xline(120, 'b--'); xline(250, 'b--');

%% Bandpass filter: 120-250 Hz (consistent with SWR literature)
cfg = [];
cfg.f = [120 250];  
cfg.display_filter = 0;
SWR_lfp_filtered = FilterLFP(cfg, SWR_lfp_raw);

%% Smooth the power envelope before z-scoring
% Key step: raw Hilbert power is noisy; smoothing reveals true ripple bursts
cfg = [];
cfg.output = 'power';
SWR_power = LFPpower(cfg, SWR_lfp_filtered);

% Gaussian smoothing (~10 ms window is standard for SWR detection)
smooth_win_sec = 0.01;  % 10 ms
smooth_win_samples = round(smooth_win_sec * Fs);
SWR_power_smooth = SWR_power;
SWR_power_smooth.data = smoothdata(SWR_power.data, 'gaussian', smooth_win_samples);

% Z-score the smoothed envelope
SWR_power_z = zscore_tsd(SWR_power_smooth);

% %% Optional: control channel artifact rejection
% cfg = [];
% cfg.fc = ExpKeys.goodRef(1);
% SWR_lfp_control = LoadCSC(cfg);
% SWR_lfp_control = restrict(SWR_lfp_control, ...
%     [ExpKeys.prerecord(1) ExpKeys.postrecord(1)], ...
%     [ExpKeys.prerecord(2) ExpKeys.postrecord(2)]);
% 
% cfg = [];
% cfg.output = 'power';
% control_power = LFPpower(cfg, SWR_lfp_control);
% control_power_smooth = control_power;
% control_power_smooth.data = smoothdata(control_power.data, 'gaussian', smooth_win_samples);
% 
% % Subtract control power to suppress movement/common-mode artifacts
% diff_power = SWR_power_smooth;
% diff_power.data = SWR_power_smooth.data - control_power_smooth.data;
% diff_power_z = zscore_tsd(diff_power);

%% Detect candidate events (low threshold) then filter by peak (high threshold)
% Using two-threshold approach: detect at 2 SD, require peak >= 3 SD
% Adjust thresholds based on your data — noisier recordings need higher values

detection_z  = 2;   % event boundary threshold
peak_z       = 3;   % minimum peak z-score within each candidate event

cfg = [];
cfg.method    = 'zscore';
cfg.threshold = detection_z;
cfg.operation = '>';
cfg.merge_thr = 0.05;  % merge events within 50 ms (catches doublets)
cfg.minlen    = 0.01; % minimum 20 ms duration (at least 2-3 cycles at 120 Hz)

SWR_evt = TSDtoIV(cfg, SWR_power_z);  % use diff_power_z if using control channel

%% Add peak z-score per event and filter out weak detections
cfg = [];
cfg.method = 'max';
cfg.label  = 'peakZ';
SWR_evt = AddTSDtoIV(cfg, SWR_evt, SWR_power_z);

cfg = [];
cfg.method = 'mean';
cfg.label  = 'meanZ';
SWR_evt = AddTSDtoIV(cfg, SWR_evt, SWR_power_z);

% Keep only events whose peak exceeds the higher threshold
keep_idx = SWR_evt.usr.peakZ >= peak_z;
SWR_evt.tstart = SWR_evt.tstart(keep_idx);
SWR_evt.tend   = SWR_evt.tend(keep_idx);
SWR_evt.usr.peakZ  = SWR_evt.usr.peakZ(keep_idx);
SWR_evt.usr.meanZ  = SWR_evt.usr.meanZ(keep_idx);

fprintf('Detected %d SWR events\n', length(SWR_evt.tstart));

%% Visualize detections
cfg_mr = [];
cfg_mr.lfpMax     = Inf;
cfg_mr.lfpHeight  = 15;
cfg_mr.lfpSpacing = 5;
cfg_mr.lfpColor   = 'k';
cfg_mr.lfp(1) = SWR_lfp_raw;
cfg_mr.lfp(2) = SWR_power_z;
cfg_mr.evt    = SWR_evt;   % <-- was SWR_iv (undefined) in your original

MultiRaster(cfg_mr, S);
xlim([7005.5 7008.5]);

%% Plot ripple-aligned LFP traces
cfg_plot = [];
cfg_plot.display = 'iv';
cfg_plot.mode    = 'center';
cfg_plot.fgcol   = 'k';
PlotTSDfromIV(cfg_plot, SWR_evt, SWR_lfp_raw);

%% Save
filename = fullfile(pwd, append(file_name, '_detectedSWRs_basic.mat'));
save(filename, '-struct', 'SWR_evt');

%% Save all open figures
figList = findobj(allchild(0), 'flat', 'Type', 'figure');

cd (folderPath) 

for i = 1:length(figList)
    figHandle = figList(i);
    figNum = num2str(get(figHandle, 'Number'));
    saveas(figHandle, ['Figure_' figNum '.png']);
end