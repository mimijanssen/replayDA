%% Mixed Linear Effects Model: Value Strenght - SWR_DA 
% subject, session, pre/post, SWR-DA, dF-value RPE 

clear; clc;
cd 'C:\Users\mimia\Documents\GitHub\replayDA\Analysis'
load ('colors.mat')

cd 'F:\Mouse_avg'

matrix_valswr = zeros(92,5);

%%  Populate Matrix with Mouse Names 

% ~~~~~~~~~~~~~~ SWR-DA ~~~~~~~~~~~~~~
cd 'F:\M433\avg_data\avg_data'
Files=dir('*.*');
count_mouse = 1; 
for k=3:length(Files)
   FileNames=Files(k).name;
   M433sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 1; % M433 is mouse 1 
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 1; % M433 is mouse 1 
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

cd 'F:\M452\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M452sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 2; % M452 is mouse 2 
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 2; % M452 is mouse 2
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end


cd 'F:\M453\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M453sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 3; % M453 is mouse 3
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 3; % M453 is mouse 3
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

cd 'F:\M460\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M460sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 4; % M460 is mouse 4
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 4; % M460 is mouse 4
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

cd 'F:\M533\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M533sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 5; % M453 is mouse 5
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 5; % M453 is mouse 5
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

cd 'F:\M545\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M545sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 6; % M453 is mouse 6
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 6; % M453 is mouse 6
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

cd 'F:\M547\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M547sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 7; % M453 is mouse 7
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 7; % M453 is mouse 7
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

cd 'F:\M548\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M548sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 8; % M453 is mouse 8
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 8; % M453 is mouse 8
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

% session needs to be hard coded. 
sessions = [1,1,2,2,3,3,5,5,6,6,7,7,8,8,1,1,2,2,4,4,5,5,6,6,7,7,2,2,3,3,4,4,5,5,8,8,1,1,5,5,6,6,7,7,1,1,2,2,3,3,4,4,5,5,6,6,7,7,4,4,5,5,6,6,7,7,1,1,2,2,3,3,4,4,5,5,6,6,1,1,2,2,3,3,4,4,5,5,6,6,7,7];

matrix_valswr(:,2) = sessions';

%% SAVE SWR-DA Strength 
SWR_DA_strength = {}; 
count_mouse = 1; 
% List of structure names
structure_names = {'M433sd', 'M452sd','M453sd', 'M460sd','M533sd','M545sd','M547sd','M548sd'}; % Add all your structure names here

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

        one_sec = first_half + 1000; % ok I hardcoded this
        mu_pre = mean(curr_structure.(session).circ_avg_pre(first_half+1:one_sec));
        mu_post = mean(curr_structure.(session).circ_avg_post(first_half+1:one_sec));
        std_pre = mean(curr_structure.(session).circ_std_pre(first_half+1:one_sec)); 
        std_post = mean(curr_structure.(session).circ_std_post(first_half+1:one_sec)); 


        SWR_DA_strength.pre(1,i) = max((curr_structure.(session).avg_fiber_pre(first_half+1:one_sec)-mu_pre)/std_pre);  % currently this is dividing by 2 sd. So I want to divide by one
        matrix_valswr(count_mouse,4) = max((curr_structure.(session).avg_fiber_pre(first_half+1:one_sec)-mu_pre)/std_pre); % pre is 1
        count_mouse = count_mouse + 1;
        % did i z-score first for correlation plots?
        % can't compare the size of these two peaks because I'm z-scoring
        % based on two different means and stds. consider changing this in
        % the future. 

        % Post-condition SWR-DA strength
        SWR_DA_strength.post(1,i) = max((curr_structure.(session).avg_fiber_post(first_half+1:one_sec)-mu_post)/std_post);
        matrix_valswr(count_mouse,4) =max((curr_structure.(session).avg_fiber_post(first_half+1:one_sec)-mu_post)/std_post);
        count_mouse = count_mouse + 1;

    end
    
    % Save the results in SWR_DA_strength_all under the current structure name
    SWR_DA_strength_all.(curr_structure_name) = SWR_DA_strength;
end

%% ~~~~~~~~~~~~~~ RPE ~~~~~~~~~~~~~~

cd 'D:\M433\ratio'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M433rat.(['ratio',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M452\ratio'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M452rat.(['ratio',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M453\ratio'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M453rat.(['ratio',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M460\ratio'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M460rat.(['ratio',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M533\ratio'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M533rat.(['ratio',num2str(k-2)]) = load(FileNames);
end

% cd 'F:\M534\avg_data\RPE_ttest'
% Files=dir('*.*');
% for k=3:length(Files)
%    FileNames=Files(k).name;
%    M534rpe.(['RPE',num2str(k-2)]) = load(FileNames);
% end

cd 'D:\M545\ratio'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M545rat.(['ratio',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M547\ratio'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M547rat.(['ratio',num2str(k-2)]) = load(FileNames);
end

cd 'D:\M548\ratio'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M548rat.(['ratio',num2str(k-2)]) = load(FileNames);
end


%% RPE STRENGTH
Ratio_strength_dF = {}; 
Ratio_strength_all = struct();
Ratio_strength_dF_all = struct();
count_mouse = 1; 

% CHANGE SUB TO RATIO DEPENDING ON WHAT YOU ARE TESTING! 

% List of structure names
structure_names = {'M433rat', 'M452rat','M453rat','M460rat', 'M533rat','M545rat','M547rat','M548rat'}; % Add all your structure names here

for s = 1:length(structure_names)
    curr_structure_name = structure_names{s}; % Get current structure name as a string
    curr_structure = eval(curr_structure_name); % Get the structure itself using eval

    % Initialize SWR_DA_strength for the current structure
    Ratio_strength_dF_all.pre = zeros(1, num_sessions); % Initialize pre values
    Ratio_strength_dF_all.post = zeros(1, num_sessions); % Initialize pre values

    num_sessions = length(fieldnames(curr_structure)); % Number of sessions (assuming each field is a session)
        
    % Get session names (assuming they are sess1, sess2, ..., sessN)
    session_names = fieldnames(curr_structure);
    
    % Iterate through each session in the current structure
    for i = 1:num_sessions
        session = session_names{i}; % Get current session name
        Ratio_strength_dF_all.pre(1,i) = curr_structure.(session).wake_nrem_sub.pre;%dF_tstats.tstat); %p;
        matrix_valswr(count_mouse,5) = curr_structure.(session).wake_nrem_sub.pre; 
        count_mouse = count_mouse + 1;
        Ratio_strength_dF_all.post(1,i) = curr_structure.(session).wake_nrem_sub.post;%dF_tstats.tstat); %p;

        matrix_valswr(count_mouse,5) = curr_structure.(session).wake_nrem_sub.post; 
        count_mouse = count_mouse + 1;
    end
    
    % Save the results in SWR_DA_strength_all under the current structure name
    Ratio_strength_all.(curr_structure_name) = Ratio_strength_dF_all;
end


%%  Matrix to table 
tbl = table(matrix_valswr(:,1),matrix_valswr(:,2),matrix_valswr(:,3),matrix_valswr(:,4),matrix_valswr(:,5),'VariableNames',{'Mouse','Session','PrePost','SWRDA','sub'});

%% Add early late to the model 
list = zeros(height(tbl),1); 
% early sessions 1-4 = 1
I_early = find(tbl.Session < 5);  % Find indices where 'condition_row' is positive
list(I_early) = 1; 

% late sessions 5-6 = 2
I_late = find(tbl.Session > 4);  % Find indices where 'condition_row' is positive
list(I_late) = 2; 

% append list to table 
tbl.("EarlyLate") = list;

%% 
%tbl2 = tbl;
%tbl2.ratio = sqrt(tbl2.ratio);

%% Base model: 
lmebase = fitlme(tbl,'SWRDA ~ 1 + (1|Mouse) + (1|Session)'); % used to do 1|Session:Mouse but didn't make the model better...
disp(lmebase)

%% Full model: 
lmefull = fitlme(tbl,'SWRDA ~ PrePost + EarlyLate + sub + (1|Mouse) + (1|Session)'); % used to do 1|Session:Mouse but didn't make the model better...
disp(lmefull)

%% Full model - ratio: 
lmefullnosub = fitlme(tbl,'SWRDA ~ PrePost + EarlyLate + (1|Mouse) + (1|Session)'); % used to do 1|Session:Mouse but didn't make the model better...
disp(lmefullnosub)

%% ratio
sub = fitlme(tbl,'SWRDA ~ sub + (1|Mouse) + (1|Session)');
disp(sub)

% compared to the full model
[results,siminfo] = compare(lmefullnosub, lmefull,'nsim',1000)
[results,siminfo] = compare(lmebase, sub,'nsim',1000)


%% alternative models

% session as a factor, random intercepts for mouse 
lme1 = fitlme(tbl,'SWRDA ~ dFValue + Session + PrePost + (1|Mouse)');
disp(lme1)
% AIC: 220.89
% nothing is significant.

lme1_1 = fitlme(tbl,'SWRDA ~ Session + PrePost + (1|Mouse)');
disp(lme1_1)

compare(lme1,lme1_1)

lme1_2 = fitlme(tbl,'SWRDA ~ dFValue + Session + (1|Mouse)');
disp(lme1_2)
% AIC: 221.01

% random intercepts for mouse and session
lme2 = fitlme(tbl,'SWRDA ~ dFValue + PrePost + (1|Session) + (1|Mouse)');
disp(lme2)
% AIC: 221.57
% nothing is significant.

% random intercepts for mouse and session nested within mouse
lme3 = fitlme(tbl,'SWRDA ~ dFValue + PrePost + (1|Mouse) + (1|Session:Mouse)');
disp(lme3)
% AIC: 219.4
% nothing is significant. pre and post is almost
% makes the most sense.

% e.g. Horsepower|EngineType) session and mouse are correlated random effects
lme4 = fitlme(tbl,'SWRDA ~ dFValue + PrePost + (Session|Mouse)');
disp(lme4)
% AIC 212.28- best model!

lme5 = fitlme(tbl,'SWRDA ~ PrePost + (Session|Mouse)');
disp(lme5)
% AIC: 212.89

compare(lme4,lme5, 'nsim',1000) % You must use this test to test for both fixed- and random-effect terms. Note that both models are fit using the default fitting method, M

%
lme3_v2 = fitlme(tbl,'SWRDA ~ PrePost + (1|Mouse) + (1|Session:Mouse)');
disp(lme3_v2)

compare(lme3, lme3_v2, 'nsim',1000)


%% PRE TASK REST
structure_names = {'M433rat', 'M452rat','M453rat', 'M460rat','M533rat','M545rat','M547rat','M548rat'}; % Add all your structure names here
structure_names_swrda = {'M433sd', 'M452sd', 'M453sd','M460sd', 'M533sd','M545sd','M547sd','M548sd'}; % Add all your structure names here

% Initialize arrays to hold all pre values
SWR_DA_pre_all = [];
SWR_DA_post_all = [];
Ratio_pre_all = [];
Ratio_post_all = [];
group_labels = []; % To store group labels for coloring the plot

% Iterate through each structure
for s = 1:length(structure_names_swrda)
    curr_structure_name = structure_names_swrda{s}; % Get current structure name as a string
    curr_structure_name_rat = structure_names{s}; % Get current structure name as a string

    % Extract the SWR_DA_strength and RPE_strength for this structure
    SWR_DA_strength = SWR_DA_strength_all.(curr_structure_name) % Assuming you've already computed this
    Ratio_strength_0 = Ratio_strength_all.(curr_structure_name_rat); % Assuming this is computed similarly
    Ratio_strength.pre = (Ratio_strength_0.pre(Ratio_strength_0.pre ~= 0)); % Assuming this is computed similarly
    Ratio_strength.post = (Ratio_strength_0.post(Ratio_strength_0.post ~= 0)); % Assuming this is computed similarly
    Ratio_strength
    % Append the pre values to the arrays
    SWR_DA_pre_all = [SWR_DA_pre_all, SWR_DA_strength.pre]; % Append SWR_DA pre values
    SWR_DA_post_all = [SWR_DA_post_all, SWR_DA_strength.post]; % Append SWR_DA pre values

    Ratio_pre_all = [Ratio_pre_all, Ratio_strength.pre]; % Append SWR_DA pre values
    Ratio_post_all = [Ratio_post_all, Ratio_strength.post]; % Append SWR_DA pre values

    % Append the group labels for coloring
    group_labels = [group_labels, repmat(s, 1, length(SWR_DA_strength.pre))];
end


%% remove outliers: 
% Fit initial linear model
Ratio_col = Ratio_pre_all(:);
SWR_col = SWR_DA_pre_all(:);

mdl = fitlm(Ratio_col, SWR_col);

% Calculate Cook's distance
cookD = mdl.Diagnostics.CooksDistance;
cookD = cookD(:);

% Common cutoff
threshold = 4 / length(Ratio_col);

% Identify potentially influential observations
outliers = cookD > threshold;

% Display which observations are influential
disp('Potential influential observations:');
disp(find(outliers));

% Display their values
disp(table(Ratio_col(outliers), ...
           SWR_col(outliers), ...
           cookD(outliers), ...
           'VariableNames', {'Ratio','SWR_DA','CooksDistance'}));

% Remove influential observations
Ratio_clean = Ratio_col(~outliers);
SWR_DA_clean = SWR_col(~outliers);

% Refit model without influential observations
mdl_clean = fitlm(Ratio_clean, SWR_DA_clean);

% Display cleaned model
disp('Original model:');
disp(mdl);

disp('Cleaned model:');
disp(mdl_clean);

disp('Cleaned model coefficients:');
disp(mdl_clean.Coefficients);

disp('Cleaned model ANOVA:');
disp(anova(mdl_clean,'summary'));

disp('Cleaned adjusted R squared:');
disp(mdl_clean.Rsquared.Adjusted);


%% Pre-Task Rest Scatter Plot - Linear Model 
figure (1);
hold on;

% Optionally, fit a regression line to all data
mdl = fitlm(Ratio_pre_all,SWR_DA_pre_all);
h = plot(mdl); % Plot the fit line
disp(mdl.Coefficients);
disp(anova(mdl,'summary'));
disp(mdl)
for s = 1:length(structure_names)
    % Get the logical indices for the current group
    idx = group_labels == s;
    
    % Check if there are any points for this group
    if sum(idx) > 0
        % Plot each structure's values with different colors and make points semi-transparent
        scatter(Ratio_pre_all(idx),SWR_DA_pre_all(idx), ...
            100, colors(s,:), 'filled', 'DisplayName', structure_names{s}, ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7); % Transparency level of 0.6
    end
end


delete(h(1))
legend('hide')

adj_R_squared = mdl.Rsquared.Adjusted; % Extract adjusted R^2
text(max(SWR_DA_pre_all) * 0.5, max(Ratio_pre_all) * 0.8, ...
    ['Adjusted R^2 = ', num2str(adj_R_squared, '%.3f')], ...
    'FontSize', 16, 'Color', 'k');
disp('adjusted R squared:'); disp(adj_R_squared); 

set(gca,'fontsize', 16)
%set(gcf, 'color','none');
%set(gca,'color','none');
set(gcf, 'renderer','painters');
%fontname("AvenirNext LT Pro Regular");

%ylim([-1 5]);
xlim([0 50]);

ylabel('SWR-DA Strength');
xlabel('Ratio Strength');
title('Pre SWR-DA Strength vs Ratio Strength');
hold off;
% 
%cd 'C:\Users\mimia\Desktop\supp_figures2'
%exportgraphics(gcf,'ValuevsSWRDA_Pre_dF.eps','ContentType','vector'); 

%% Post-Task Rest Scatter Plot - Linear Model 
figure (2);
hold on;
mdl = fitlm( Ratio_post_all,SWR_DA_post_all);
h = plot(mdl); % Plot the fit line

adj_R_squared = mdl.Rsquared.Adjusted; % Extract adjusted R^2
text(max(SWR_DA_post_all) * 0.5, max(Ratio_post_all) * 0.8, ...
    ['Adjusted R^2 = ', num2str(adj_R_squared, '%.3f')], ...
    'FontSize', 16, 'Color', 'k');
legend('hide')

disp(mdl.Coefficients);
disp(anova(mdl,'summary'));
disp(mdl)

for s = 1:length(structure_names)
    % Get the logical indices for the current group
    idx = group_labels == s;
    
    % Check if there are any points for this group
    if sum(idx) > 0
        % Plot each structure's values with different colors and make points semi-transparent
        scatter(Ratio_post_all(idx),SWR_DA_post_all(idx), ...
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
%fontname("AvenirNext LT Pro Regular");

%xlim([0 20]);
%ylim([-1 5]);

% Title add labels
%legend('hide')
ylabel('SWR-DA Strength');
xlabel('Ratio Strength');
title('Post SWR-DA Strength vs Ratio Strength');
hold off;
