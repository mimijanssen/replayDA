%% EXTRA SESSION DATA to save: 

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

%% SAVE SWR-DA Strength 
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
        
        % Pre-condition SWR-DA strength
        SWR_DA_strength.pre(1,i) = max(curr_structure.(session).avg_fiber_pre)/max(curr_structure.(session).circ_std_pre);
        
        % Post-condition SWR-DA strength
        SWR_DA_strength.post(1,i) = max(curr_structure.(session).avg_fiber_post)/max(curr_structure.(session).circ_std_post);
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
        RPE_strength_dF_all(1,i) = max(curr_structure.(session).dF_tstats.tstat);%dF_tstats.tstat); %p;
 
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
        scatter(SWR_DA_pre_all(idx), RPE_all(idx), ...
            100, colors(s,:), 'filled', 'DisplayName', structure_names{s}, ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7); % Transparency level of 0.6
    end
end

% Optionally, fit a regression line to all data
mdl = fitlm(SWR_DA_pre_all, RPE_all);
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

xlim([-1 4]);
%ylim([-0.01 0.03]);

xlabel('SWR-DA Strength (Pre-Task Rest)');
ylabel('DA Value Strength');
title('Pre SWR-DA Strength vs Value Strength');
hold off;
% 
 cd 'C:\Users\mimia\OneDrive\Desktop\figures'
 exportgraphics(gcf,'ValuevsSWRDA_Pre_dF.eps','ContentType','vector'); 

%% Post-Task Rest Scatter Plot - Linear Model 
figure (2);
hold on;

mdl = fitlm(SWR_DA_post_all, RPE_all);
h = plot(mdl); % Plot the fit line

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
        scatter(SWR_DA_post_all(idx), RPE_all(idx), ...
            100, colors(s,:), 'filled', 'DisplayName', structure_names{s}, ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7); % Transparency level of 0.6
    end
end

delete(h(1))
legend('hide')

set(gca,'fontsize', 16)
%set(gcf, 'color','none');
%set(gca,'color','none');
set(gcf, 'renderer','painters');
fontname("AvenirNext LT Pro Regular");

xlim([-1 4]);
%ylim([-0.01 0.03]);

% Title add labels
%legend('hide')
xlabel('SWR-DA Strength (Post-Task Rest)');
ylabel('DA Value Strength');
title('Post SWR-DA Strength vs Value Strength');
hold off;

% cd 'C:\Users\mimia\OneDrive\Desktop\figures'
% exportgraphics(gcf,'ValuevsSWRDA_Post_dF.eps','ContentType','vector'); 

%% Making Matrix 
% subject, session, pre/post, SWR-DA, dF-value RPE 
% rpe_swrda_tbl =  

%% Linear Mixed Effects Model of Data

%% ~~~~~~~~~~~~~~ Speed ~~~~~~~~~~~~~~

cd 'D:\M433\avg_data\pos'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M433pos.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M453\avg_data\pos'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M453pos.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M460\avg_data\pos'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M460pos.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M533\avg_data\pos'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M533pos.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M534\avg_data\pos'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M534pos.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M545\avg_data\pos'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M545pos.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M547\avg_data\pos'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M547pos.(['sess',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M548\avg_data\pos'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M548pos.(['sess',num2str(k-2)]) = load(FileNames);
end

%% SAVE SPEED
% List of structure names
structure_names = {'M433pos', 'M453pos','M460pos', 'M533pos', 'M534pos','M545pos','M547pos','M548pos'}; % Add all your structure names here

% Initialize an empty structure to hold the SWR-DA strengths for each dataset
speed_all = struct();

% Iterate through each structure
for s = 1:length(structure_names)
    curr_structure_name = structure_names{s}; % Get current structure name as a string
    curr_structure = eval(curr_structure_name); % Get the structure itself using eval
    
    % Initialize SWR_DA_strength for the current structure
    track_speed = {};
    num_sessions = length(fieldnames(curr_structure)); % Number of sessions (assuming each field is a session)
    
    track_speed = zeros(1, num_sessions); % Initialize post values
    
    % Get session names (assuming they are sess1, sess2, ..., sessN)
    session_names = fieldnames(curr_structure);
    
    % Iterate through each session in the current structure
    for i = 1:num_sessions
        session = session_names{i}; % Get current session name
        
        % Pre-condition SWR-DA strength
        track_speed(1,i) = curr_structure.(session).speed_track;
    end
    
    % Save the results in SWR_DA_strength_all under the current structure name
    speed_all.(curr_structure_name) = track_speed;
end

%% SWR_DA vs Speed

structure_names = {'M433pos', 'M453pos','M460pos', 'M533pos', 'M534pos','M545pos','M547pos','M548pos'}; % Add all your structure names here
structure_names_swrda = {'M433sd',  'M453sd', 'M460sd','M533sd', 'M534sd','M545sd','M547sd','M548sd'}; % Add all your structure names here

% Initialize arrays to hold all pre values
SWR_DA_pre_all = [];
SWR_DA_post_all = [];
speed_all2 = [];
group_labels = []; % To store group labels for coloring the plot

% Iterate through each structure
for s = 1:length(structure_names_swrda)
    curr_structure_name = structure_names_swrda{s}; % Get current structure name as a string
    curr_structure_name_pos = structure_names{s}; % Get current structure name as a string

    % Extract the SWR_DA_strength and RPE_strength for this structure
    SWR_DA_strength = SWR_DA_strength_all.(curr_structure_name); % Assuming you've already computed this
    track_speed = speed_all.(curr_structure_name_pos); % Assuming this is computed similarly
    
    % Append the pre values to the arrays
    SWR_DA_pre_all = [SWR_DA_pre_all, SWR_DA_strength.pre]; % Append SWR_DA pre values
    SWR_DA_post_all = [SWR_DA_post_all, SWR_DA_strength.post]; % Append SWR_DA pre values
    speed_all2 = [speed_all2,track_speed]; % Append RPE pre values

    % Append the group labels for coloring
    group_labels = [group_labels, repmat(s, 1, length(SWR_DA_strength.pre))];
end
%%
% Create the scatter plot
fg3 = figure(3);
hold on;

% Optionally, fit a regression line to all data
mdl = fitlm(speed_all2,SWR_DA_pre_all);
h = plot(mdl); % Plot the fit line

for s = 1:length(structure_names)
    % Get the logical indices for the current group
    idx = group_labels == s;
    if sum(idx) > 0
        scat = scatter(speed_all2(idx), SWR_DA_pre_all(idx),...
            100, colors(s,:), 'filled', 'DisplayName', structure_names{s}, ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7); % Transparency level of 0.6
    end
end

xlabel('Track Speed (pixels/s)');
ylabel('SWR-DA Strength');
title('"Motivation" vs. SWR-DA (Pre)');

delete(h(1))
legend('hide')
adj_R_squared = mdl.Rsquared.Adjusted; % Extract adjusted R^2
text( max(speed_all2) * 0.6,max(SWR_DA_pre_all) * 0.8, ...
    ['Adjusted R^2 = ', num2str(adj_R_squared, '%.3f')], ...
    'FontSize', 16, 'Color', 'k');
set(gca,'fontsize', 16)
%set(gcf, 'color','none');
%set(gca,'color','none');
set(gcf, 'renderer','painters');
%fontname("AvenirNext LT Pro Regular");
ylim([-1 4]);
%legend('Fit','CI','M433','M453','M460','M533','M534','M545','M547','M548')
hold off;

%fg3.WindowState = 'maximized';

cd 'C:\Users\mimia\OneDrive\Desktop\figures'
exportgraphics(gcf,'MotivationvsSWRDA_Pre.eps','ContentType','vector'); 

%%
% POST TASK RESTTTTTT
fg4 = figure(4);
hold on;

mdl = fitlm(speed_all2,SWR_DA_post_all);
h = plot(mdl); % Plot the fit line

for s = 1:length(structure_names)
    % Get the logical indices for the current group
    idx = group_labels == s;
    % Check if there are any points for this group
    if sum(idx) > 0
        % Plot each structure's values with different colors and make points semi-transparent
        scatter(speed_all2(idx), SWR_DA_post_all(idx),...
            100, colors(s,:), 'filled', 'DisplayName', structure_names{s}, ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7); % Transparency level of 0.6
    end
end

xlabel('Track Speed (pixels/s)');
ylabel('SWR-DA Strength');
title('"Motivation" vs. SWR-DA (Post)');

delete(h(1))
legend('hide')

adj_R_squared = mdl.Rsquared.Adjusted; % Extract adjusted R^2
text( max(speed_all2) * 0.6,max(SWR_DA_post_all) * 0.8, ...
    ['Adjusted R^2 = ', num2str(adj_R_squared, '%.3f')], ...
    'FontSize', 16, 'Color', 'k');
set(gca,'fontsize', 16)
%set(gcf, 'color','none');
%set(gca,'color','none');
set(gcf, 'renderer','painters');
%fontname("AvenirNext LT Pro Regular");
ylim([-1 4]);
hold off;

cd 'C:\Users\mimia\OneDrive\Desktop\figures'
exportgraphics(gcf,'MotivationvsSWRDA_Post.eps','ContentType','vector'); 


%% ~~~~~~~~~~~~~~ Fiber Power and Voltage Difference ~~~~~~~~~~~~~~
% 43 entries... 'M433','M453','M460','M533','M534','M545','M547','M548'
power = [50,50,50,50,50,50,50,160,160,160,160.2,160, 50,50,50.2,50,50,50,50,50.2,50,50.3,50.1,100,50.2,50.1,70.6,85.7,NaN,96.6,100,100.5,100,100.3,100,100,50,50,50.1,50.3,50,50.3,50]; 
voltdiff = [0.575,0.665,0.905,NaN,0.725,0.61,0.645,0.235,0.215,0.345,0.245,0.1, 0.27,0.29,0.135,0.16,0.51,0.44,0.615,1.105,0.71,1.22,0.69,0.31,0.355,0.345,0.535,0.77,NaN,0.8,0.225,0.165,0.28,0.335,0.30,0.36,0.6,0.74,0.84,0.67,0.77,0.71,0.65]; 

%% Plot POWER

fg5 = figure(5);
hold on;

mdl = fitlm(power,SWR_DA_post_all);
h = plot(mdl); % Plot the fit line

for s = 1:length(structure_names)
    % Get the logical indices for the current group
    idx = group_labels == s;
    % Check if there are any points for this group
    if sum(idx) > 0
        % Plot each structure's values with different colors and make points semi-transparent
        scatter(power(idx), SWR_DA_post_all(idx),...
            100, colors(s,:), 'filled', 'DisplayName', structure_names{s}, ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7); % Transparency level of 0.6
    end
end

xlabel('Fiber Power (uW)');
ylabel('SWR-DA Strength');
title('SWR-DA (Post) vs. Power');

delete(h(1))
legend('hide')

adj_R_squared = mdl.Rsquared.Adjusted; % Extract adjusted R^2
text( max(power) * 0.6,max(SWR_DA_post_all) * 0.8, ...
    ['Adjusted R^2 = ', num2str(adj_R_squared, '%.3f')], ...
    'FontSize', 16, 'Color', 'k');
set(gca,'fontsize', 16)
%set(gcf, 'color','none');
%set(gca,'color','none');
set(gcf, 'renderer','painters');
fontname("AvenirNext LT Pro Regular");
ylim([-1 4]);
hold off;

cd 'C:\Users\mimia\Desktop\PrePrint Figures'
exportgraphics(gcf,'PowervsSWRDA_Post.png','ContentType','vector'); 

%% POWER PRE
fg6 = figure(6);
hold on;

mdl = fitlm(power,SWR_DA_pre_all);
h = plot(mdl); % Plot the fit line

for s = 1:length(structure_names)
    % Get the logical indices for the current group
    idx = group_labels == s;
    % Check if there are any points for this group
    if sum(idx) > 0
        % Plot each structure's values with different colors and make points semi-transparent
        scatter(power(idx), SWR_DA_pre_all(idx),...
            100, colors(s,:), 'filled', 'DisplayName', structure_names{s}, ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7); % Transparency level of 0.6
    end
end

xlabel('Fiber Power (uW)');
ylabel('SWR-DA Strength');
title('SWR-DA (Pre) vs. Power');

delete(h(1))
legend('hide')

adj_R_squared = mdl.Rsquared.Adjusted; % Extract adjusted R^2
text( max(power) * 0.6,max(SWR_DA_pre_all) * 0.8, ...
    ['Adjusted R^2 = ', num2str(adj_R_squared, '%.3f')], ...
    'FontSize', 16, 'Color', 'k');
set(gca,'fontsize', 16)
%set(gcf, 'color','none');
%set(gca,'color','none');
set(gcf, 'renderer','painters');
fontname("AvenirNext LT Pro Regular");
ylim([-1 4]);
hold off;


cd 'C:\Users\mimia\Desktop\PrePrint Figures'
exportgraphics(gcf,'PowervsSWRDA_Pre.png','ContentType','vector'); 

%% POST VOLT DIFF

fg7 = figure(7);
hold on;

mdl = fitlm(voltdiff,SWR_DA_post_all);
h = plot(mdl); % Plot the fit line

for s = 1:length(structure_names)
    % Get the logical indices for the current group
    idx = group_labels == s;
    % Check if there are any points for this group
    if sum(idx) > 0
        % Plot each structure's values with different colors and make points semi-transparent
        scatter(voltdiff(idx), SWR_DA_post_all(idx),...
            100, colors(s,:), 'filled', 'DisplayName', structure_names{s}, ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7); % Transparency level of 0.6
    end
end

xlabel('Voltage (V)');
ylabel('SWR-DA Strength');
title('SWR-DA (Post) vs. Voltage Diff');

delete(h(1))
legend('hide')

adj_R_squared = mdl.Rsquared.Adjusted; % Extract adjusted R^2
text( max(voltdiff) * 0.6,max(SWR_DA_post_all) * 0.8, ...
    ['Adjusted R^2 = ', num2str(adj_R_squared, '%.3f')], ...
    'FontSize', 16, 'Color', 'k');
set(gca,'fontsize', 16)
%set(gcf, 'color','none');
%set(gca,'color','none');
set(gcf, 'renderer','painters');
fontname("AvenirNext LT Pro Regular");
ylim([-1 4]);
hold off;

cd 'C:\Users\mimia\Desktop\PrePrint Figures'
exportgraphics(gcf,'VoltvsSWRDA_Post.png','ContentType','vector'); 

%% PRE VOLT DIFF

fg8 = figure(8);
hold on;

mdl = fitlm(voltdiff,SWR_DA_pre_all);
h = plot(mdl); % Plot the fit line

for s = 1:length(structure_names)
    % Get the logical indices for the current group
    idx = group_labels == s;
    % Check if there are any points for this group
    if sum(idx) > 0
        % Plot each structure's values with different colors and make points semi-transparent
        scatter(voltdiff(idx), SWR_DA_pre_all(idx),...
            100, colors(s,:), 'filled', 'DisplayName', structure_names{s}, ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7); % Transparency level of 0.6
    end
end

xlabel('Voltage (V)');
ylabel('SWR-DA Strength');
title('SWR-DA (Pre) vs. Voltage Diff');

delete(h(1))
legend('hide')

adj_R_squared = mdl.Rsquared.Adjusted; % Extract adjusted R^2
text( max(voltdiff) * 0.6,max(SWR_DA_pre_all) * 0.8, ...
    ['Adjusted R^2 = ', num2str(adj_R_squared, '%.3f')], ...
    'FontSize', 16, 'Color', 'k');
set(gca,'fontsize', 16)
%set(gcf, 'color','none');
%set(gca,'color','none');
set(gcf, 'renderer','painters');
fontname("AvenirNext LT Pro Regular");
ylim([-1 4]);
hold off;

cd 'C:\Users\mimia\Desktop\PrePrint Figures'
exportgraphics(gcf,'VoltvsSWRDA_Pre.png','ContentType','vector'); 
