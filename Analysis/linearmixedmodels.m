%% Mixed Linear Effects Model: Value Strenght - SWR_DA 
% subject, session, pre/post, SWR-DA, dF-value RPE 

clear; clc;
cd 'D:\Mouse_avg'
load ('colors.mat')

matrix_valswr = zeros(86,5);

%%  Populate Matrix with Mouse Names 

% ~~~~~~~~~~~~~~ SWR-DA ~~~~~~~~~~~~~~
cd 'D:\M433\avg_data\avg_data'
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

cd 'D:\M453\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M453sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 2; % M453 is mouse 2 
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 2; % M453 is mouse 2 
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

cd 'D:\M460\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M460sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 3; % M453 is mouse 3
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 3; % M453 is mouse 3
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

cd 'D:\M533\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M533sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 4; % M453 is mouse 4
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 4; % M453 is mouse 4
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

cd 'D:\M534\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M534sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 5; % M453 is mouse 5
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 5; % M453 is mouse 5
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

cd 'D:\M545\avg_data\avg_data'
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

cd 'D:\M547\avg_data\avg_data'
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

cd 'D:\M548\avg_data\avg_data'
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
sessions = [1,1,2,2,3,3,5,5,6,6,7,7,8,8,2,2,3,3,4,4,5,5,8,8,1,1,5,5,6,6,7,7,1,1,2,2,3,3,4,4,5,5,6,6,7,7,1,1,3,3,4,4,4,4,5,5,6,6,7,7,1,1,2,2,3,3,4,4,5,5,6,6,1,1,2,2,3,3,4,4,5,5,6,6,7,7];

matrix_valswr(:,2) = sessions';

%% SAVE SWR-DA Strength 
SWR_DA_strength = {}; 
count_mouse = 1; 
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

        
        mu_pre = mean(curr_structure.(session).circ_avg_pre(first_half+1:end));
        mu_post = mean(curr_structure.(session).circ_avg_post(first_half+1:end));
        std_pre = mean(curr_structure.(session).circ_std_pre(first_half+1:end))/2; 
        std_post = mean(curr_structure.(session).circ_std_post(first_half+1:end))/2; 


        SWR_DA_strength.pre(1,i) = max((curr_structure.(session).avg_fiber_pre(first_half+1:end)-mu_pre)/std_pre);  % currently this is dividing by 2 sd. So I want to divide by one
        matrix_valswr(count_mouse,4) = max((curr_structure.(session).avg_fiber_pre(first_half+1:end)-mu_pre)/std_pre); % pre is 1
        count_mouse = count_mouse + 1;
        % did i z-score first for correlation plots?
        % can't compare the size of these two peaks because I'm z-scoring
        % based on two different means and stds. consider changing this in
        % the future. 

        % Post-condition SWR-DA strength
        SWR_DA_strength.post(1,i) = max((curr_structure.(session).avg_fiber_post(first_half+1:end)-mu_post)/std_post);
        matrix_valswr(count_mouse,4) =max((curr_structure.(session).avg_fiber_post(first_half+1:end)-mu_post)/std_post);
        count_mouse = count_mouse + 1;


    end
    
    % Save the results in SWR_DA_strength_all under the current structure name
    SWR_DA_strength_all.(curr_structure_name) = SWR_DA_strength;
end

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
count_mouse = 1; 

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
        RPE_strength_dF_all(1,i) = curr_structure.(session).dF_tstats.tstat;%dF_tstats.tstat); %p;
        matrix_valswr(count_mouse,5) = curr_structure.(session).dF_tstats.tstat; 
        count_mouse = count_mouse + 1;
        matrix_valswr(count_mouse,5) = curr_structure.(session).dF_tstats.tstat; 
        count_mouse = count_mouse + 1;
    end
    
    % Save the results in SWR_DA_strength_all under the current structure name
    RPE_strength_all.(curr_structure_name) = RPE_strength_dF_all;
end


%%  Matrix to table 

tbl = table(matrix_valswr(:,1),matrix_valswr(:,2),matrix_valswr(:,3),matrix_valswr(:,4),matrix_valswr(:,5),'VariableNames',{'Mouse','Session','PrePost','SWRDA','dFValue'});

% define Mouse, Session and PrePost a categorical variable 
%tbl.Mouse = nominal(tbl.Mouse);
%tbl.PrePost = nominal(tbl.PrePost);

lme_everything = fitlme(tbl,'dFValue ~ SWRDA + Session + PrePost + (1|Mouse)');
% AIC = 429.23
% BIC = 443.95

lme_everything2 = fitlme(tbl,'dFValue ~ SWRDA + PrePost + (1|Mouse)');

lme_everything3 = fitlme(tbl,'dFValue ~ SWRDA + (1|Mouse)');

lme_swrda = fitlme(tbl,'SWRDA ~ dFValue + PrePost + (1|Mouse)');
% AIC = 194.74
% BIC = 209.47

lme_swrda_value = fitlme(tbl,'SWRDA ~ dFValue + (1|Mouse)');

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



%%

F = fitted(lme_swrda);
R = residuals(lme_swrda);

plot(F,R,'bx')
xlabel('Fitted Values')
ylabel('Residuals')

figure()
gscatter(F,R,Session)
