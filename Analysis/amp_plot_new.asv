%% Initialize
clear; clc;
addpath('C:\Users\mimia\OneDrive\Documents\Toolboxes\raacampbell-shadedErrorBar'); % path to shaded error bar 

% Define folder paths
salineFolder = 'D:\Saline';
amphetamineFolder = 'D:\Amphetamine';

% Get list of mouse folders
salineMice = dir(fullfile(salineFolder, 'M*'));
amphetamineMice = dir(fullfile(amphetamineFolder, 'M*'));

%% 1) save all session data (-10 min to 14 min from injection time) 
% --> and find longest injection time 

% initialize variables
% 25 min = 1500 s ; 1500 s * 1000 samples/s = 1500000 samples
% double check what I expect ^
salineMice = dir(fullfile(salineFolder, 'M*'));
ampMice = dir(fullfile(amphetamineFolder, 'M*'));
referenceMouse = 'M433'; % Update with your choice
referenceFile = 'M433_2023_10_18processed'; % Update with your choice
refData = load(fullfile(salineFolder, referenceMouse, referenceFile));
minutes = 21; 
samples = minutes*60*refData.cfg.hdr{1,1}.SamplingFrequency;

samples_before = 10*60*refData.cfg.hdr{1,1}.SamplingFrequency;
samples_after = 11*60*refData.cfg.hdr{1,1}.SamplingFrequency;

Saline_sess = {}; 
Amp_sess = {}; 
Saline_time = {};
Amp_time = {};
Sal_inj_time = {}; 
Amp_inj_time = {}; 

% find injection time : avg between drug admin and admin done 
% --- Process saline group ---
for i = 1:length(salineMice) % iterate through the mice
    mouseFolder = fullfile(salineFolder, salineMice(i).name);
    dataFiles = dir(fullfile(mouseFolder, '*processed*')); % list of sessions for that mouse
    mouseName = salineMice(i).name;
    numSessions = length(dataFiles);
    sessionData = zeros(numSessions, samples); % assuming all sessions have 1500000 samples
    sess_salineTime = zeros(numSessions, samples);
    sess_injTime = zeros(numSessions,1);

    for j = 1:numSessions % now open each of those sessions 

        data = load(fullfile(mouseFolder, dataFiles(j).name)); % data.tvec; data.F_zscored_base for detrended and z_scored data
        keys = (fullfile(mouseFolder, strrep(dataFiles(j).name, 'processed.mat', '')));
        cd (keys);
        LoadExpKeys();
        cfg_evt = [];
        cfg_evt.eventList = ExpKeys.eventList;
        cfg_evt.eventLabel = ExpKeys.eventLabel;
        evt = LoadEvents(cfg_evt); % 7 is injection start ; 8 is injection end

        % ok, I did not save the original tvec. so now I need to extract
        % these injection times and then Load the original data and
        % subtract the original tvec. Then I will take the average. 
        inj_init = evt.t{1,7};
        inj_end = evt.t{1,8}; 

        raw_data_path = fullfile(strrep(fileparts(mouseFolder), 'Saline', salineMice(i).name), strrep(dataFiles(j).name, 'processed.mat', ''));
        cd(raw_data_path)
        cfg.fc = {'CSC30.ncs'};
        csc_photo = LoadCSC(cfg);

        inj_init = inj_init - csc_photo.tvec(1); 
        inj_end = inj_end - csc_photo.tvec(1); 
        
        inj_time = (inj_init + inj_end)/2; 
        index_time = nearest_idx3(inj_time, data.tvec); % values, look-up values

        init_bound = index_time - samples_before; 
        end_bound = index_time + samples_after; 
        %check = end_bound - init_bound % should be number of samples 
        
        % take -10 to 15 min of that and put it in a structure 
        sessionData(j, :) = data.zdF_base(init_bound:end_bound-1); %F_zscored_base
        sess_salineTime(j,:) = data.tvec(init_bound:end_bound-1); 

        % save injection time
        sess_injTime(j,:) = inj_time; 
    end
    
    Saline_sess.(mouseName) = sessionData;
    Saline_time.(mouseName) = sess_salineTime;
    Sal_inj_time.(mouseName) = sess_injTime; 

end

% --- Process amphetamine group ---
for i = 1:length(ampMice) % iterate through the mice
    mouseFolder = fullfile(amphetamineFolder, ampMice(i).name);
    dataFiles = dir(fullfile(mouseFolder, '*processed*')); % list of sessions for that mouse
    mouseName = ampMice(i).name;
    numSessions = length(dataFiles);
    sessionData = zeros(numSessions, samples); % assuming all sessions have 1500000 samples
    sess_ampTime = zeros(numSessions, samples);
    sess_injTime = zeros(numSessions,1);

    for j = 1:numSessions % now open each of those sessions 

        data = load(fullfile(mouseFolder, dataFiles(j).name)); % data.tvec; data.F_zscored_base for detrended and z_scored data
        keys = (fullfile(mouseFolder, strrep(dataFiles(j).name, 'processed.mat', '')));
        cd (keys);
        LoadExpKeys();
        cfg_evt = [];
        cfg_evt.eventList = ExpKeys.eventList;
        cfg_evt.eventLabel = ExpKeys.eventLabel;
        evt = LoadEvents(cfg_evt); % 7 is injection start ; 8 is injection end

        % ok, I did not save the original tvec. so now I need to extract
        % these injection times and then Load the original data and
        % subtract the original tvec. Then I will take the average. 
        inj_init = evt.t{1,7}(1);
        inj_end = evt.t{1,8}(1); 

        raw_data_path = fullfile(strrep(fileparts(mouseFolder), 'Amphetamine', ampMice(i).name), strrep(dataFiles(j).name, 'processed.mat', ''));
        cd(raw_data_path)
        cfg.fc = {'CSC30.ncs'};
        csc_photo = LoadCSC(cfg);

        inj_init = inj_init - csc_photo.tvec(1); 
        inj_end = inj_end - csc_photo.tvec(1); 
        
        inj_time = (inj_init + inj_end)/2; 
        index_time = nearest_idx3(inj_time, data.tvec); % values, look-up values

        init_bound = index_time - samples_before; 
        end_bound = index_time + samples_after; 
        %check = end_bound - init_bound % should be number of samples 
        
        % take -10 to 15 min of that and put it in a structure 
        sessionData(j, :) = data.zdF_base(init_bound:end_bound-1); 
        sess_ampTime(j,:) = data.tvec(init_bound:end_bound-1); 

        % save injection time
        sess_injTime(j,:) = inj_time; 
    end
    
    Amp_sess.(mouseName) = sessionData;
    Amp_time.(mouseName) = sess_ampTime;
    Amp_inj_time.(mouseName) = sess_injTime; 

end

% find max injection time for both saline and amphetamine group 

% find max mouse for saline sessions: 
[M,I] = max([max(Sal_inj_time.M433),max(Sal_inj_time.M453),max(Sal_inj_time.M460),max(Sal_inj_time.M533),max(Sal_inj_time.M534),max(Sal_inj_time.M545),max(Sal_inj_time.M547),max(Sal_inj_time.M548)]); 
% 808.3070 ; mouse 6 
[M2,I2] = max([max(Amp_inj_time.M433),max(Amp_inj_time.M453),max(Amp_inj_time.M460),max(Amp_inj_time.M533),max(Amp_inj_time.M534),max(Amp_inj_time.M545),max(Amp_inj_time.M547),max(Amp_inj_time.M548)]); 
% 786 ; mouse 6 
% Saline session of M545 is the longest injection time 

% find session: 
[M3,I3] = max(Sal_inj_time.M545); 
% session 1. 

%% 2) save all average mouse data (-10 to 14 min from max injection time)
% average each field and save in an overall matrix 

% Get the names of the fields in the structure
fields_saline = fieldnames(Saline_sess);
fields_amp = fieldnames(Amp_sess); 

% Initialize an empty matrix to store the averages
avg_sal = [];
avg_amp = []; 

% Loop through each field
for iter = 1:length(fields_saline) % iterate through each mouse 
    % Get the data for the current field
    data_saline_sess = Saline_sess.(fields_saline{iter});
    % time_saline_sess = Saline_time.(fields_saline{iter});  
    % inj_saline_sess = Sal_inj_time.(fields_saline{iter});
    % Calculate the average of the current field
    field_average = mean(data_saline_sess,1);
    % Append the average to the matrix as a new row
    avg_sal = [avg_sal; field_average];
end

for iter = 1:length(fields_amp)
    data_amp_sess = Amp_sess.(fields_amp{iter});
    field_average = mean(data_amp_sess,1);
    avg_amp = [avg_amp; field_average];
end


%% 3) find mean and SEM for plotting 
% average and SEM 
meanSal = mean(avg_sal,1); 
meanAmp = mean(avg_amp,1); 
SEMSal = std(avg_sal)/sqrt(size(avg_sal, 1));
SEMAmp = std(avg_amp)/sqrt(size(avg_amp, 1)); 

%% 4) plotting 
time_for_plotting = Saline_time.M433(1,:) - Saline_time.M433(1,1); % intialized. 
% This means that the injection variable should be + samples before + 1
% this is for 1,000 Hz. How would we plot minutes? can relabel the x axis. 
time_for_plotting = time_for_plotting/(60); % in minutes now? 
injection_time_index = 600001;

figure (1);
hold on;

% Saline
shadedErrorBar(time_for_plotting, meanSal, SEMSal, 'lineProps', {'-k', 'LineWidth', 2});

% Amphetamine
shadedErrorBar(time_for_plotting, meanAmp, SEMAmp, 'lineProps', {'-b', 'LineWidth', 2});

% Injection window shading
%injectStart = refData.tvec(refAlignIdx);
%injectEnd = refData.tvec(refAlignIdx + 500); % Example: adjust based on injection duration
%patch([injectStart, injectEnd, injectEnd, injectStart], ...
    %  [min(meanSaline), min(meanSaline), max(meanAmphetamine), max(meanAmphetamine)], ...
    %  [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Dashed line for injection onset
% xline(injection_time_index, '--k', 'LineWidth', 1);
% 
% 
% x_ticks = get(gca, 'XTick');  % Get current x-tick positions
% x_ticks_in_minutes = (x_ticks - injection_time_index) / 60*1000;  % Convert to minutes
% 
% % Set the new x-tick positions and labels
% xticks(x_ticks);  % Set x-tick positions
% xticklabels(arrayfun(@(x) sprintf('%d', x), x_ticks_in_minutes, 'UniformOutput', false));  % Set x-tick labels

% Labels and legend
xlabel('Time from Injection (s)');
ylabel('Z-scored Signal');
legend({'Saline', 'Amphetamine'}, 'Location', 'best');
title('Average Fiber Photometry Traces');
set(gca, 'FontSize', 12);
hold off;
