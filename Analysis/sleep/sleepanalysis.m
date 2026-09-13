 restoredefaultpath; clear classes; % start with a clean slate
 cd('C:\Users\mimia\Documents\GitHub\AccuSleep_X-main'); % or, wherever your code is located -- NOTE \shared subfolder!
 p = genpath(pwd); % create list of all folders from here
 addpath(p);

%% Save data as matlab file 
cd 'F:\M460\M460-2024-01-22_recording7'; 
LoadExpKeys
csc_name = [];
csc_name.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(csc_name); % csc with good ripples
EEG = csc.data; 

% SAVE EEG 
save('EEG.mat', 'EEG');

% legnth baased resampling 
N_target = length(EEG);
P = readtable('M460-2024-01-22-VT1-convertedDLC_resnet50_Linear TrackApr5shuffle1_100000.csv','PreserveVariableNames',true); % CHANGE THIS 
midbody = [P(:,1), P(:,17) P(:,18), P(:,19)];
midbody = table2array(midbody); 
old_axis = linspace(0, 1, length(midbody));
new_axis = linspace(0, 1, N_target);

EMG = interp1(old_axis, midbody(:,2), new_axis, 'linear');

% SAVE EEG 
save('EMG.mat', 'EMG');

% if you have 10 second epochs. how many bins do you need to represent your
% data 

% SAVE labels as an empty thing? 


fs = 1500;                 % sampling rate (Hz)
epoch_length = 2;         % seconds
samples_per_epoch = fs * epoch_length;

N_samples = length(EEG);
n_epochs = floor(N_samples / samples_per_epoch);  % number of full epochs

% Create label vector (change values later if needed)
labels = ones(1,n_epochs);   % column vector is usually safer

save('labels.mat', 'labels');

% Save Position as EMG?? BUT IS IT THE SAME SIZE.. 


%%
save('AS_config_mimi.mat','cfg_colors','cfg_names','')
