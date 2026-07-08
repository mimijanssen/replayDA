%% Multi-session PETH - Dorsal and Ventral SWRs
clear; clc;

addpath('C:\Users\mimia\Documents\Toolboxes\shadedErrorBar')

%% Settings
base_dir    = 'D:\FiberSWRDataZip'; % change to your path
seconds     = 10;  % window: +/- 5 seconds around SWR
n_samples   = seconds * 1000 + 1; % 10001 samples assuming 1kHz -- adjust if needed
tvec_plot   = linspace(-seconds/2, seconds/2, n_samples);
N_circ      = 1000; % number of circshifts for shuffle

%% Find all session folders
session_dirs = dir(fullfile(base_dir, 'M*'));
session_dirs = session_dirs([session_dirs.isdir]); % keep only folders
fprintf('Found %d session folders\n', length(session_dirs));

%% Initialize accumulators across sessions
% We will collect one mean trace per session, then SEM across sessions
dorsal_session_means  = [];
ventral_session_means = [];
dorsal_circ_means     = [];
ventral_circ_means    = [];
dorsal_circ_stds     = [];
ventral_circ_stds    = [];

%% Loop over sessions
addpath('C:\Users\mimia\Documents\GitHub\vandermeerlab-replay-da\chronux-master\spectral_analysis\helper')
addpath('C:\Users\mimia\Documents\GitHub\vandermeerlab-replay-da\chronux-master\spectral_analysis\continuous')

for s = 1:length(session_dirs)
    sess_path = fullfile(base_dir, session_dirs(s).name);
    fprintf('\nProcessing session: %s\n', session_dirs(s).name);

    % Check which SWR files exist
    has_dorsal  = exist(fullfile(sess_path, 'dorsal_imec0_SWR.mat'),  'file');
    has_ventral = exist(fullfile(sess_path, 'ventral_imec0_SWR.mat'), 'file');
    has_fiber   = exist(fullfile(sess_path, 'XA7_synced.mat'),        'file');

    if ~has_fiber
        fprintf('  No fiber file found, skipping\n');
        continue
    end

    % Load fiber signal
    fiber_data = load(fullfile(sess_path, 'XA7_synced.mat'));
    data = fiber_data.data;
    tvec = fiber_data.tvec;
    fs   = fiber_data.fs;

    samples = round((seconds/2) * fs); % samples before/after SWR

    % Preprocess fiber signal (matching your existing pipeline)
    FP        = struct();
    FP.data   = data;
    FP.tvec   = tvec;
    FP.fs     = fs;

    % Median filter
    FP_denoised = medfilt1(FP.data);

    % Butterworth filter
    fc      = 20;
    [b, a]  = butter(8, fc/(fs/2));
    FP_filt = filtfilt(b, a, FP_denoised);

    % Local detrend (60s window)
    FP_detrend = locdetrend(FP_filt, fs, [60 1]);

    % Z-score
    zF = zscore(FP_detrend);
    % edit: tstart is actually tcenter
    % Process dorsal SWRs
    if has_dorsal
        dorsal_data = load(fullfile(sess_path, 'dorsal_imec0_SWR.mat'));
        % adjust field name if yours differs
        if isfield(dorsal_data, 'dorsal_swr')
            swr_center = (dorsal_data.dorsal_swr.iv.tstart + dorsal_data.dorsal_swr.iv.tend)./2;
        elseif isfield(dorsal_data, 'evt')
            swr_center = (dorsal_data.evt.tstart + dorsal_data.evt.tend)./2;
        else
            fields = fieldnames(dorsal_data);
            swr_center = (dorsal_data.(fields{1}).iv.tstart + dorsal_data.(fields{1}).iv.tend)./2;
        end

        [d_mean, d_circ, d_std] = extract_peth(swr_center, tvec, zF, fs, samples, n_samples, N_circ);
        dorsal_session_means = [dorsal_session_means; d_mean];
        dorsal_circ_means    = [dorsal_circ_means;    d_circ];
        dorsal_circ_stds    = [dorsal_circ_stds;    d_std];
        fprintf('  Dorsal: %d SWRs\n', length(swr_center));
    end

    % Process ventral SWRs
    if has_ventral
        ventral_data = load(fullfile(sess_path, 'ventral_imec0_SWR.mat'));
        if isfield(ventral_data, 'ventral_swr')
            swr_center = (ventral_data.ventral_swr.iv.tstart + ventral_data.ventral_swr.iv.tend)./2;
        elseif isfield(ventral_data, 'evt')
            swr_center = (ventral_data.evt.tstart + ventral_data.evt.tend)./2;
        else
            fields = fieldnames(ventral_data);
            swr_center = (ventral_data.(fields{1}).iv.tstart + ventral_data.(fields{1}).iv.tend)./2;
        end

        [v_mean, v_circ, v_stds] = extract_peth(swr_center, tvec, zF, fs, samples, n_samples, N_circ);
        ventral_session_means = [ventral_session_means; v_mean];
        ventral_circ_means    = [ventral_circ_means;    v_circ];
        ventral_circ_stds    = [ventral_circ_stds;    v_stds];
        fprintf('  Ventral: %d SWRs\n', length(swr_center));
    end

end

%% Grand average and SEM across sessions
fprintf('\n--- Summary ---\n');
fprintf('Dorsal:  %d sessions\n', size(dorsal_session_means, 1));
fprintf('Ventral: %d sessions\n', size(ventral_session_means, 1));

% Dorsal
dorsal_grand_mean  = nanmean(dorsal_session_means,  1);
dorsal_grand_sem   = nanstd(dorsal_session_means,  0, 1) / sqrt(size(dorsal_session_means,  1));
dorsal_circ_mean   = nanmean(dorsal_circ_means,     1);
dorsal_circ_sem    = nanstd(dorsal_circ_means,      0, 1) / sqrt(size(dorsal_circ_means,    1));

% Ventral
ventral_grand_mean = nanmean(ventral_session_means, 1);
ventral_grand_sem  = nanstd(ventral_session_means, 0, 1) %/ sqrt(size(ventral_session_means, 1));
ventral_circ_mean  = nanmean(ventral_circ_means,    1);
ventral_circ_sem   = nanstd(ventral_circ_means,     0, 1); %/ sqrt(size(ventral_circ_means,   1));

%% Plot Dorsal
figure(1); clf;

% Shuffle
%shadedErrorBar(tvec_plot, dorsal_circ_mean, dorsal_circ_sem, ...
%    'lineProps', {'Color', [0.5 0.5 0.5], 'LineWidth', 2}, 'transparent', 1);
xline(0, 'k--', 'LineWidth', 1.5);
yline(0, 'k:',  'LineWidth', 0.8);

hold on;

% Signal
shadedErrorBar(tvec_plot, dorsal_grand_mean, dorsal_grand_sem, ...
    'lineProps', {'Color', [0 0.45 0.74], 'LineWidth', 2}, 'transparent', 1);

xlim([-seconds/2 seconds/2]);
ylim([-0.6 0.6]);
xlabel('Time from SWR (s)', 'FontSize', 16);
ylabel('Mean [DA] (z-score)',  'FontSize', 16);
title(sprintf('Dorsal SWR PETH (n = %d sessions)', size(dorsal_session_means,1)), 'FontSize', 20);
%legend('', 'Signal', 'Location', 'northwest');
%legend boxoff;
set(gca, 'FontSize', 14); box off;
%set(gcf, 'color', 'none'); set(gca, 'color', 'none');
set(gcf, 'renderer', 'painters');
%fontname("AvenirNext LT Pro Regular");

%% Dorsal zoom 

figure(2); clf;

% Shuffle
%shadedErrorBar(tvec_plot, dorsal_circ_mean, dorsal_circ_sem, ...
%    'lineProps', {'Color', [0.5 0.5 0.5], 'LineWidth', 2}, 'transparent', 1);
xline(0, 'k--', 'LineWidth', 1.5);
yline(0, 'k:',  'LineWidth', 0.8);

hold on;

% Signal
shadedErrorBar(tvec_plot, dorsal_grand_mean, dorsal_grand_sem, ...
    'lineProps', {'Color', [0 0.45 0.74], 'LineWidth', 2}, 'transparent', 1);

xlim([-seconds/2 seconds/2]);
ylim([-0.2 0.6]);
xlim([-0.1 0.5]);
xlabel('Time from SWR (s)', 'FontSize', 16);
ylabel('Mean [DA] (z-score)',  'FontSize', 16);
title(sprintf('Dorsal SWR PETH (n = %d sessions)', size(dorsal_session_means,1)), 'FontSize', 20);
%legend('', 'Signal', 'Location', 'northwest');
%legend boxoff;
set(gca, 'FontSize', 14); box off;
%set(gcf, 'color', 'none'); set(gca, 'color', 'none');
set(gcf, 'renderer', 'painters');

%% Plot Ventral
figure(3); clf;
hold on;

shadedErrorBar(tvec_plot, ventral_circ_mean, ventral_circ_stds, ...
    'lineProps', {'Color', [0.5 0.5 0.5], 'LineWidth', 2}, 'transparent', 1);

shadedErrorBar(tvec_plot, ventral_grand_mean, ventral_grand_sem, ...
    'lineProps', {'Color', [0.85 0.33 0.1], 'LineWidth', 2}, 'transparent', 1);

xline(0, 'k--', 'LineWidth', 1.5);
yline(0, 'k:',  'LineWidth', 0.8);
xlim([-seconds/2 seconds/2]);
ylim([-1.5 1.5]);
xlabel('Time from SWR (s)', 'FontSize', 16);
ylabel('Mean [DA] (z-score)',  'FontSize', 16);
title(sprintf('Ventral SWR PETH (n = %d sessions)', size(ventral_session_means,1)), 'FontSize', 20);
legend( 'Shuffle','Signal','', '', 'Location', 'northwest');
legend boxoff;
set(gca, 'FontSize', 14); box off;
%set(gcf, 'color', 'none'); set(gca, 'color', 'none');
set(gcf, 'renderer', 'painters');
%fontname("AvenirNext LT Pro Regular");

%% Plot individual session PETHs for Dorsal with shuffle

n_dorsal_sess = size(dorsal_session_means, 1);
n_cols = 3;
n_rows = ceil(n_dorsal_sess / n_cols);

figure(4); clf;
sgtitle('Dorsal SWR PETH - Individual Sessions', 'FontSize', 16);

% track which sessions had dorsal SWRs (need to re-loop to get session names)
dorsal_sess_names = {};
for s = 1:length(session_dirs)
    sess_path  = fullfile(base_dir, session_dirs(s).name);
    has_dorsal = exist(fullfile(sess_path, 'dorsal_imec0_SWR.mat'), 'file');
    has_fiber  = exist(fullfile(sess_path, 'XA7_synced.mat'),       'file');
    if has_dorsal && has_fiber
        dorsal_sess_names{end+1} = session_dirs(s).name;
    end
end

for s = 1:n_dorsal_sess
    subplot(n_rows, n_cols, s);
    hold on;

    sess_mean = dorsal_session_means(s, :);
    circ_mean = dorsal_circ_means(s, :);
    circ_std = dorsal_circ_stds(s, :);


    % Shuffle trace
    plot(tvec_plot, circ_mean, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
        'DisplayName', 'Shuffle');

    shadedErrorBar(tvec_plot, circ_mean, circ_std, ...
    'lineProps', {'Color', [0.5 0.5 0.5], 'LineWidth', 2}, 'transparent', 1);

    % Signal trace
    plot(tvec_plot, sess_mean, 'Color', [0 0.45 0.74], 'LineWidth', 2, ...
        'DisplayName', 'Signal');

    xline(0, 'k--', 'LineWidth', 1.2);
    yline(0, 'k:',  'LineWidth', 0.8);

    xlim([-seconds/2 seconds/2]);
    ylim([-1 1.5]);
    xlabel('Time from SWR (s)', 'FontSize', 11);
    ylabel('[DA] (z-score)',     'FontSize', 11);

    if s <= length(dorsal_sess_names)
        title(dorsal_sess_names{s}, 'FontSize', 10, 'Interpreter', 'none');
    else
        title(sprintf('Session %d', s), 'FontSize', 10);
    end

    if s == 1
        legend('','Shuffle','Signal','','Location', 'northwest', 'FontSize', 8);
        legend boxoff;
    end

    box off;
    set(gca, 'FontSize', 10);
end

%set(gcf, 'color', 'none');
set(gcf, 'renderer', 'painters');
%fontname("AvenirNext LT Pro Regular");

%% Helper: extract PETH for a given set of SWR times
function [sess_mean, sess_circ_mean, sess_circ_std] = extract_peth(swr_center, tvec, zF, fs, samples, n_samples, N_circ)

% Find fiber indices for each SWR
    swr_idx = zeros(length(swr_center), 1);
    for i = 1:length(swr_center)
        swr_idx(i) = nearest_idx3(swr_center(i), tvec);
    end

    % Extract segments
    segments = NaN(length(swr_center), n_samples);
    for ievt = 1:length(swr_center)
        idx_start = swr_idx(ievt) - samples;
        idx_end   = swr_idx(ievt) + samples;
        if idx_start < 1 || idx_end > length(zF)
            continue
        end
        seg = zF(idx_start:idx_end);
        if length(seg) == n_samples
            segments(ievt,:) = seg;
        end
     end

     % Session mean (ignoring NaN rows)
     sess_mean = nanmean(segments, 1);

     % Circshift shuffle
     X        = zF;
     K        = randi([1 length(zF)], 1, N_circ);
     circ_mat = NaN(N_circ, n_samples);

     for iter = 1:N_circ
         Y        = circshift(X, K(iter));

         circ_seg = NaN(length(swr_center), n_samples);
         for ievt = 1:length(swr_center)
             idx_start = swr_idx(ievt) - samples;
             idx_end   = swr_idx(ievt) + samples;
             if idx_start < 1 || idx_end > length(Y)
                  continue
              end
              seg = Y(idx_start:idx_end);
              if length(seg) == n_samples
                  circ_seg(ievt,:) = seg;
              end
         end
         circ_mat(iter,:) = nanmean(circ_seg, 1);
     end
     sess_circ_mean = nanmean(circ_mat, 1);
     sess_circ_std  = 2*std(circ_mat,0,1);
end

    