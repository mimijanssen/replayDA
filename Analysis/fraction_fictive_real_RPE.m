%% FRACTION OF SWR peaks to RPE peaks... 

% ~ saving: ~ 
% RPE strength (difference between omission signal and high-value signal) 
% SWR-DA strength (peak signal - 1 sd of circ shuffle) 
% motivation (average speed on track) 

% OTHER NOTE: MAKE SURE Excitation light is in exp keys 
% Excitation light (power uW) -- should be in exp keys. 
%  other things to do : list of sig swr index: pos and negative 
clear; clc;
cd 'D:\Mouse_avg'
load ('colors.mat')
% edited colors.mat so that M433 is a unique color and then resaved it...
% as colors.mat in the harddrive...

%% LOAD ALL MICE INFO 

% ~~~~~~~~~~~~~~ SWR-DA ~~~~~~~~~~~~~~
cd 'D:\M433\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M433sd.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M453\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M453sd.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M460\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M460sd.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M533\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M533sd.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M534\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M534sd.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M545\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M545sd.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M547\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M547sd.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M548\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M548sd.(['sess',num2str(k-2)]) = load(FileNames);
end

%% SAVE SWR-DA Strength ? 4 seconds....
SWR_DA_strength = {}; 

% List of structure names
structure_names = {'M433sd',  'M453sd', 'M460sd','M533sd', 'M534sd','M545sd','M547sd','M548sd'}; % Add all your structure names here

% Initialize an empty structure to hold the SWR-DA strengths for each dataset
SWR_DA_strength_all = struct();

% Iterate through each structure
for s = 1:length(structure_names)
    curr_structure_name = structure_names{s}; % Get current structure name as a string
    curr_structure = eval(curr_structure_name); % Get the structure itself using eval
    
    % Initialize SWR_DA_strength for the current structure
    SWR_DA_strength = {};
    num_sessions = length(fieldnames(curr_structure)); % Number of sessions (assuming each field is a session)
    
    SWR_DA_strength.pre = zeros(1, num_sessions); % Initialize pre values
    SWR_DA_strength.post = zeros(1, num_sessions); % Initialize post values
    
    % Get session names (assuming they are sess1, sess2, ..., sessN)
    session_names = fieldnames(curr_structure);
    
    % Iterate through each session in the current structure
    for i = 1:num_sessions
        session = session_names{i}; % Get current session name
        
        first_half = floor(length(curr_structure.(session).avg_fiber_pre)/2);

        % Pre-condition SWR-DA strength
        SWR_DA_strength.pre(1,i) = max(curr_structure.(session).avg_fiber_pre(first_half+1:end-(first_half/2)));
        
        % Post-condition SWR-DA strength
        SWR_DA_strength.post(1,i) = max(curr_structure.(session).avg_fiber_post(first_half+1:end-(first_half/2)));
    end
    
    % Save the results in SWR_DA_strength_all under the current structure name
    SWR_DA_strength_all.(curr_structure_name) = SWR_DA_strength;
end

% SWR_DA_strenght_all has all of the SWR_DA strengths. 

%% ~~~~~~~~~~~~~~ RPE ~~~~~~~~~~~~~~

cd 'D:\M433\avg_data\RPE_ttest'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M433rpe.(['RPE',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M453\avg_data\RPE_ttest'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M453rpe.(['RPE',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M460\avg_data\RPE_ttest'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M460rpe.(['RPE',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M533\avg_data\RPE_ttest'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M533rpe.(['RPE',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M534\avg_data\RPE_ttest'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M534rpe.(['RPE',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M545\avg_data\RPE_ttest'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M545rpe.(['RPE',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M547\avg_data\RPE_ttest'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M547rpe.(['RPE',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M548\avg_data\RPE_ttest'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M548rpe.(['RPE',num2str(k-2)]) = load(FileNames);
end

%% RPE STRENGTH
RPE_strength_dF = {}; 
RPE_strength_all = struct();

% List of structure names
structure_names = {'M433rpe',  'M453rpe','M460rpe', 'M533rpe', 'M534rpe','M545rpe','M547rpe','M548rpe'}; % Add all your structure names here

for s = 1:length(structure_names)
    curr_structure_name = structure_names{s}; % Get current structure name as a string
    curr_structure = eval(curr_structure_name); % Get the structure itself using eval

    % Initialize SWR_DA_strength for the current structure
    RPE_strength_dF = {};
    num_sessions = length(fieldnames(curr_structure)); % Number of sessions (assuming each field is a session)
    
    RPE_strength_dF_all = zeros(1, num_sessions); % Initialize RPE values
    
    % Get session names (assuming they are sess1, sess2, ..., sessN)
    session_names = fieldnames(curr_structure);
    
    % Iterate through each session in the current structure
    for i = 1:num_sessions
        session = session_names{i}; % Get current session name
        %[h,p]= ttest(max(curr_structure.(session).avg_high),min(curr_structure.(session).avg_low));

        % Pre-condition SWR-DA strength
        RPE_strength_dF_all(1,i) = mean(curr_structure.(session).dF_high); % was previously taking the max...idk why but it doesn't change anything since there is one value. 
        % this is for 2 seconds following photobeam break!
    end
    
    % Save the results in SWR_DA_strength_all under the current structure name
    RPE_strength_all.(curr_structure_name) = RPE_strength_dF_all;
end

% RPE_strength_all has all of the RPE high - low. 

%% CORRELATION OF RPE STRENGTH AND SWR_DA 
% PRE TASK REST
structure_names = {'M433rpe', 'M453rpe', 'M460rpe','M533rpe', 'M534rpe','M545rpe','M547rpe','M548rpe'}; % Add all your structure names here
structure_names_swrda = {'M433sd',  'M453sd','M460sd', 'M533sd', 'M534sd','M545sd','M547sd','M548sd'}; % Add all your structure names here

% Initialize arrays to hold all pre values
SWR_DA_pre_all = [];
SWR_DA_post_all = [];
RPE_all = [];
group_labels = []; % To store group labels for coloring the plot

% Iterate through each structure
for s = 1:length(structure_names_swrda)
    curr_structure_name = structure_names_swrda{s}; % Get current structure name as a string
    curr_structure_name_rpe = structure_names{s}; % Get current structure name as a string

    % Extract the SWR_DA_strength and RPE_strength for this structure
    SWR_DA_strength = SWR_DA_strength_all.(curr_structure_name); % Assuming you've already computed this
    RPE_strength = RPE_strength_all.(curr_structure_name_rpe); % Assuming this is computed similarly
    
    % Append the pre values to the arrays
    SWR_DA_pre_all = [SWR_DA_pre_all, SWR_DA_strength.pre]; % Append SWR_DA pre values
    SWR_DA_post_all = [SWR_DA_post_all, SWR_DA_strength.post]; % Append SWR_DA pre values

    RPE_all = [RPE_all, RPE_strength]; % Append RPE pre values

    % Append the group labels for coloring
    group_labels = [group_labels, repmat(s, 1, length(SWR_DA_strength.pre))];
end

%% Pre-Task Rest Scatter Plot - Linear Model 
figure (1);
hold on;
for s = 1:length(structure_names)
    % Get the logical indices for the current group
    idx = group_labels == s;
    
    % Check if there are any points for this group
    if sum(idx) > 0
        % Plot each structure's values with different colors and make points semi-transparent
        scatter(RPE_all(idx),SWR_DA_pre_all(idx), ...
            100, colors(s,:), 'filled', 'DisplayName', structure_names{s}, ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7); % Transparency level of 0.6
    end
end

% Optionally, fit a regression line to all data
mdl = fitlm(RPE_all,SWR_DA_pre_all);
h = plot(mdl); % Plot the fit line

delete(h(1))
legend('hide')

adj_R_squared = mdl.Rsquared.Adjusted; % Extract adjusted R^2
text(max(SWR_DA_pre_all) * 0.5, max(RPE_all) * 0.8, ...
    ['Adjusted R^2 = ', num2str(adj_R_squared, '%.3f')], ...
    'FontSize', 16, 'Color', 'k');

set(gca,'fontsize', 16)
%set(gcf, 'color','none');
%set(gca,'color','none');
set(gcf, 'renderer','painters');
%fontname("AvenirNext LT Pro Regular");

%ylim([-1.25 1.75]);
%xlim([0 20]);

ylabel('SWR-DA Average Transient');
xlabel('DA Value Average Transient for High Volume Trials');
title('Pre SWR-DA Transient vs Value Transient');
hold off;
% 
%cd 'C:\Users\mimia\OneDrive\Desktop\Figures_Correlations'
%exportgraphics(gcf,'ValuevsSWRDA_Pre_dF.png','ContentType','vector'); 

%% Post-Task Rest Scatter Plot - Linear Model 
figure (2);
hold on;

adj_R_squared = mdl.Rsquared.Adjusted; % Extract adjusted R^2
text(max(SWR_DA_post_all) * 0.5, max(RPE_all) * 0.8, ...
    ['Adjusted R^2 = ', num2str(adj_R_squared, '%.3f')], ...
    'FontSize', 16, 'Color', 'k');
legend('hide')

for s = 1:length(structure_names)
    % Get the logical indices for the current group
    idx = group_labels == s;
    
    % Check if there are any points for this group
    if sum(idx) > 0
        % Plot each structure's values with different colors and make points semi-transparent
        scatter(RPE_all(idx),SWR_DA_post_all(idx), ...
            100, colors(s,:), 'filled', 'DisplayName', structure_names{s}, ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7); % Transparency level of 0.6
    end
end

mdl = fitlm( RPE_all,SWR_DA_post_all);
h = plot(mdl); % Plot the fit line

delete(h(1))
legend('hide')

set(gca,'fontsize', 16)
%set(gcf, 'color','none');
%set(gca,'color','none');
set(gcf, 'renderer','painters');
%fontname("AvenirNext LT Pro Regular");

%xlim([0 20]);
%ylim([-1.25 1.75]);

% Title add labels
%legend('hide')
ylabel('SWR-DA Average Transient');
xlabel('DA Value Average Transient for High Volume Trials');
title('Post SWR-DA Transient vs Value Transient');
hold off;

%% 
mean(RPE_all)
mean(SWR_DA_post_all)
mean(SWR_DA_pre_all)
% and std

 mean(SWR_DA_pre_all)*100 / mean(RPE_all)
 mean(SWR_DA_post_all)*100 / mean(RPE_all)