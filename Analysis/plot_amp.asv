% Initialize
clear; clc;
addpath(genpath('C:\Users\mimia\Documents\Toolboxes\shadedErrorBar'));

% Define folder paths
salineFolder = 'F:\Saline';
amphetamineFolder = 'F:\Amphetamine';

% Get list of mouse folders
salineMice = dir(fullfile(salineFolder, 'M*'));
amphetamineMice = dir(fullfile(amphetamineFolder, 'M*'));

% Choose a reference file for alignment
referenceMouse = 'M433'; % Update with your choice
referenceFile = 'M433_2023_10_18processed'; % Update with your choice

% Load reference data
refData = load(fullfile(salineFolder, referenceMouse, referenceFile));
refKeys = (fullfile(salineFolder, referenceMouse, strrep(referenceFile, 'processed', ''))); % load
cd (refKeys)
LoadExpKeys(); 

% Extract reference injection time
cfg_evt = [];
cfg_evt.eventList = ExpKeys.eventList;
cfg_evt.eventLabel = ExpKeys.eventLabel;
evt = LoadEvents(cfg_evt);
refInjectTime = evt.t{1, 7}; 
refAlignIdx = nearest_idx3(refInjectTime, refData.tvec);


% Initialize storage and determine maximum length across all time series
maxLength = 0;

% First pass: Determine the maximum length for both saline and amphetamine groups
for folderType = {'salineFolder', 'amphetamineFolder'}
    folderPath = eval(folderType{1});
    miceList = dir(folderPath);
    for i = 1:length(miceList)
        if miceList(i).isdir && ~startsWith(miceList(i).name, '.')
            mouseFolder = fullfile(folderPath, miceList(i).name);
            dataFiles = dir(fullfile(mouseFolder, '*processed*'));
            for j = 1:length(dataFiles)
                data = load(fullfile(mouseFolder, dataFiles(j).name));
                maxLength = max(maxLength, length(data.F_zscored_base));
            end
        end
    end
end

% Preallocate storage matrices
alignedSaline = NaN(maxLength, length(salineMice));
alignedAmphetamine = NaN(maxLength, length(amphetamineMice));
tAligned = NaN(maxLength, length(salineMice) + length(amphetamineMice));

%%
% --- Process saline group ---
for i = 1:length(salineMice)
    mouseFolder = fullfile(salineFolder, salineMice(i).name);
    dataFiles = dir(fullfile(mouseFolder, '*processed*'));
    
    for j = 1:length(dataFiles)
        % Load data and keys
        data = load(fullfile(mouseFolder, dataFiles(j).name));
        keys = (fullfile(mouseFolder, strrep(dataFiles(j).name, 'processed.mat', '')));
        cd (keys);
        cfg_evt.eventList = ExpKeys.eventList;
        cfg_evt.eventLabel = ExpKeys.eventLabel;
        evt = LoadEvents(cfg_evt);
        
        % Extract injection time
        injectTime = evt.t{1, 7};
        alignIdx = nearest_idx3(injectTime, data.tvec);
        
        % Align data to reference
        shift = refAlignIdx - alignIdx;
        alignedFP = circshift(data.F_zscored_base, shift);
        alignedTime = circshift(data.tvec, shift);
        
        % Pad with NaNs and store in matrices
        paddedFP = NaN(maxLength, 1);
        paddedFP(1:length(alignedFP)) = alignedFP;
        alignedSaline(:, i) = paddedFP';
        
        paddedTime = NaN(maxLength, 1);
        paddedTime(1:length(alignedTime)) = alignedTime;
        tAligned(:, i) = paddedTime';
    end
end

% --- Process amphetamine group ---
for i = 1:length(amphetamineMice)
    mouseFolder = fullfile(amphetamineFolder, amphetamineMice(i).name);
    dataFiles = dir(fullfile(mouseFolder, '*processed*'));
    
    for j = 1:length(dataFiles)
        % Load data and keys
        data = load(fullfile(mouseFolder, dataFiles(j).name));
        keys = (fullfile(mouseFolder, strrep(dataFiles(j).name, 'processed.mat', '')));
        cd (keys);
        cfg_evt.eventList = ExpKeys.eventList;
        cfg_evt.eventLabel = ExpKeys.eventLabel;
        evt = LoadEvents(cfg_evt);
        
        % Extract injection time
        injectTime = evt.t{1, 7};
        alignIdx = nearest_idx3(injectTime, data.tvec);
        
        % Align data to reference
        shift = refAlignIdx - alignIdx;
        alignedFP = circshift(data.F_zscored_base, shift);
        alignedTime = circshift(data.tvec, shift);
        
        % Pad with NaNs and store in matrices
        paddedFP = NaN(maxLength, 1);
        paddedFP(1:length(alignedFP)) = alignedFP;
        alignedAmphetamine(:, i) = paddedFP;
        
        paddedTime = NaN(maxLength, 1);
        paddedTime(1:length(alignedTime)) = alignedTime;
        tAligned(:, length(salineMice) + i) = paddedTime;
    end
end


%% FIRST PASS AT PLOTTING 

% Number of mice in each group
numMice = 8;

% Calculate averages and SEM for each group
avgSaline = nanmean(alignedSaline, 2);
semSaline = nanstd(alignedSaline, 0, 2) / sqrt(numMice);

avgAmphetamine = nanmean(alignedAmphetamine, 2);
semAmphetamine = nanstd(alignedAmphetamine, 0, 2) / sqrt(numMice);

% Time vector (using the first column of tAligned as the reference)
timeVector = tAligned(:, 1); % Assuming time alignment is consistent

% Plotting
figure;
hold on;

% Plot SEM as shaded area (saline group)
%fill([timeVector; flipud(timeVector)], ...
%     [avgSaline - semSaline; flipud(avgSaline + semSaline)], ...
%     [0.5, 0.5, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Overlay average trace (saline group)
plot(timeVector, avgSaline, 'b', 'LineWidth', 2, 'DisplayName', 'Saline (Avg)');

% Plot SEM as shaded area (amphetamine group)
%fill([timeVector; flipud(timeVector)], ...
%     [avgAmphetamine - semAmphetamine; flipud(avgAmphetamine + semAmphetamine)], ...
%     [1, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Overlay average trace (amphetamine group)
plot(timeVector, avgAmphetamine, 'r', 'LineWidth', 2, 'DisplayName', 'Amphetamine (Avg)');

% Add dashed line at reference injection time
xline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Injection Time');

% Customize plot
xlabel('Time from Injection (s)');
ylabel('Signal (Z-scored)');
title('Aligned Fiber Photometry Traces');
legend('Location', 'best');
set(gca, 'FontSize', 12);
grid on;

% Adjust figure aesthetics
set(gcf, 'Color', 'w');
%xlim([-500, 1500]); % Adjust based on your data range
hold off;



%% Aligning everything to the earliest injection time. 
% Number of mice in each group
numMice = 8;

% Calculate averages and SEM for each group
avgSaline = nanmean(alignedSaline, 2);
semSaline = nanstd(alignedSaline, 0, 2) / sqrt(numMice);

avgAmphetamine = nanmean(alignedAmphetamine, 2);
semAmphetamine = nanstd(alignedAmphetamine, 0, 2) / sqrt(numMice);

% Time vector (using the first column of tAligned as the reference)
timeVector = tAligned(:, 1); % Assuming time alignment is consistent

% Injection start and end times for each mouse (example structure)
% Replace these with your actual injection start and end times (in seconds)
% Each row corresponds to a mouse, with [start_time, end_time]
injectionWindowsSaline = [-10, 50; -8, 52; -11, 48; -12, 55; -9, 50; -10, 51; -10, 49; -11, 52];
injectionWindowsAmphetamine = [-15, 45; -13, 47; -14, 48; -16, 46; -12, 50; -15, 44; -13, 49; -14, 46];

% Plotting
figure;
hold on;

% Transparent background for each saline mouse
for i = 1:numMice
    patch([injectionWindowsSaline(i, 1), injectionWindowsSaline(i, 2), ...
           injectionWindowsSaline(i, 2), injectionWindowsSaline(i, 1)], ...
          [min(ylim), min(ylim), max(ylim), max(ylim)], ...
          [0.5, 0.5, 1], 'FaceAlpha', 0.05 * i, 'EdgeColor', 'none');
end

% Transparent background for each amphetamine mouse
for i = 1:numMice
    patch([injectionWindowsAmphetamine(i, 1), injectionWindowsAmphetamine(i, 2), ...
           injectionWindowsAmphetamine(i, 2), injectionWindowsAmphetamine(i, 1)], ...
          [min(ylim), min(ylim), max(ylim), max(ylim)], ...
          [1, 0.5, 0.5], 'FaceAlpha', 0.05 * i, 'EdgeColor', 'none');
end

% Plot SEM as shaded area (saline group)
fill([timeVector; flipud(timeVector)], ...
     [avgSaline - semSaline; flipud(avgSaline + semSaline)], ...
     [0.5, 0.5, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Overlay average trace (saline group)
plot(timeVector, avgSaline, 'b', 'LineWidth', 2, 'DisplayName', 'Saline (Avg)');

% Plot SEM as shaded area (amphetamine group)
fill([timeVector; flipud(timeVector)], ...
     [avgAmphetamine - semAmphetamine; flipud(avgAmphetamine + semAmphetamine)], ...
     [1, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Overlay average trace (amphetamine group)
plot(timeVector, avgAmphetamine, 'r', 'LineWidth', 2, 'DisplayName', 'Amphetamine (Avg)');

% Add dashed line at reference injection time
xline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Injection Time');

% Customize plot
xlabel('Time from Injection (s)');
ylabel('Signal (Z-scored)');
title('Aligned Fiber Photometry Traces');
legend('Location', 'best');
set(gca, 'FontSize', 12);
grid on;

% Adjust figure aesthetics
set(gcf, 'Color', 'w');
xlim([-500, 1500]); % Adjust based on your data range
hold off;


%%
% Align saline data
for i = 1 %1:length(salineMice)
    mouseFolder = fullfile(salineFolder, salineMice(i).name);
    dataFiles = dir(fullfile(mouseFolder, '*processed*'));
    
    for j = 1:length(dataFiles) % iterate through sessions
        data = load(fullfile(mouseFolder, dataFiles(j).name));
        keys = (fullfile(mouseFolder, strrep(dataFiles(j).name, 'processed.mat', '')));
        cd (keys)
        cfg_evt.eventList = ExpKeys.eventList;
        cfg_evt.eventLabel = ExpKeys.eventLabel;
        evt = LoadEvents(cfg_evt);
        
        injectTime = evt.t{1, 7};
        alignIdx = nearest_idx3(injectTime, data.tvec); % find the time point that matches the injection index
        
        shift = refAlignIdx - alignIdx; % shift based on time difference 
        alignedFP = circshift(data.F_zscored_base, shift);
        alignedSaline = [alignedSaline; alignedFP]; % Collect aligned data
    end
end

% Repeat for amphetamine data
for i = 1:length(amphetamineMice)
    mouseFolder = fullfile(amphetamineFolder, amphetamineMice(i).name);
    dataFiles = dir(fullfile(mouseFolder, '*processed*'));
    
    for j = 1:length(dataFiles)
        % Load data and keys
        data = load(fullfile(mouseFolder, dataFiles(j).name));
        keys = (fullfile(mouseFolder, strrep(dataFiles(j).name, 'processed.mat', '')));
        cd (keys)
        cfg_evt.eventList = ExpKeys.eventList;
        cfg_evt.eventLabel = ExpKeys.eventLabel;
        evt = LoadEvents(cfg_evt);
        
        % Extract injection time
        injectTime = evt.t{1, 7};
        alignIdx = nearest_idx3(injectTime, data.tvec);
        
        % Align data to reference
        shift = refAlignIdx - alignIdx;
        alignedFP = circshift(data.F_zscored_base, shift);
        alignedAmphetamine = [alignedAmphetamine; alignedFP]; % Collect aligned data
    end
end

% Compute average and SEM
meanSaline = nanmean(alignedSaline, 1);
semSaline = nanstd(alignedSaline, 0, 1) / sqrt(size(alignedSaline, 1));
meanAmphetamine = nanmean(alignedAmphetamine, 1);
semAmphetamine = nanstd(alignedAmphetamine, 0, 1) / sqrt(size(alignedAmphetamine, 1));

% Plot the data
figure;
hold on;

% Saline
shadedErrorBar(refData.tvec, meanSaline, semSaline, 'lineProps', {'-k', 'LineWidth', 2});

% Amphetamine
shadedErrorBar(refData.tvec, meanAmphetamine, semAmphetamine, 'lineProps', {'-b', 'LineWidth', 2});

% Injection window shading
injectStart = refData.tvec(refAlignIdx);
injectEnd = refData.tvec(refAlignIdx + 500); % Example: adjust based on injection duration
patch([injectStart, injectEnd, injectEnd, injectStart], ...
      [min(meanSaline), min(meanSaline), max(meanAmphetamine), max(meanAmphetamine)], ...
      [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Dashed line for injection onset
xline(injectStart, '--k', 'LineWidth', 1);

% Labels and legend
xlabel('Time from Injection (s)');
ylabel('Z-scored Signal');
legend({'Saline', 'Amphetamine'}, 'Location', 'best');
title('Average Fiber Photometry Traces');
set(gca, 'FontSize', 12);
hold off;
