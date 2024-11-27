%% code to find the percent sessions that [DA] exceeded a shuffle at any point 
% subject, session, pre/post, SWR-DA, motivation 
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

        % Pre-condition SWR-DA strength
        % subtract the mean and divide by sd and then find the max value. 
        SWR_DA_strength.pre(1,i) = max(curr_structure.(session).avg_fiber_pre(first_half+1:end))-max((curr_structure.(session).circ_std_pre(first_half+1:end))); 
        matrix_valswr(count_mouse,4) = max((curr_structure.(session).avg_fiber_pre(first_half+1:end))-max(curr_structure.(session).circ_std_pre(first_half+1:end))); % pre is 1
        count_mouse = count_mouse + 1;

        % Post-condition SWR-DA strength
        SWR_DA_strength.post(1,i) = max((curr_structure.(session).avg_fiber_post(first_half+1:end))- max(curr_structure.(session).circ_std_post(first_half+1:end)));
        matrix_valswr(count_mouse,4) = max((curr_structure.(session).avg_fiber_post(first_half+1:end))-max(curr_structure.(session).circ_std_post(first_half+1:end)));
        count_mouse = count_mouse + 1;

    end
    
    % Save the results in SWR_DA_strength_all under the current structure name
    SWR_DA_strength_all.(curr_structure_name) = SWR_DA_strength;
end



%%  Matrix to table 

tbl = table(matrix_valswr(:,1),matrix_valswr(:,2),matrix_valswr(:,3),matrix_valswr(:,4),'VariableNames',{'Mouse','Session','PrePost','SWRDA',});


%% Find percentage of sessions that are positive (meaning significant). 
ind_pos = find(tbl.SWRDA>0);
percent_sig_sess = (length(ind_pos) / length(tbl.SWRDA))*100; 


%% Pre sessions
pre = tbl(tbl.PrePost == 1,:); % 1 is pre; 2 is post. 
percent_sig_pre = (length(find(pre.SWRDA>0)) / length(pre.SWRDA))*100; 
disp(percent_sig_pre)
disp(length(find(pre.SWRDA>0)))

%% Post sessions
post = tbl(tbl.PrePost == 2,:); % 1 is pre; 2 is post. 
percent_sig_post = (length(find(post.SWRDA>0)) / length(post.SWRDA))*100; 
disp(percent_sig_post)
disp(length(find(post.SWRDA>0)))


%% SAME THING BUT FOR NEGATIVE SESSIONS: 
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

        % Pre-condition SWR-DA strength
        % subtract the mean and divide by sd and then find the max value. 
        SWR_DA_strength.pre(1,i) = min(curr_structure.(session).avg_fiber_pre(first_half+1:end))+max((curr_structure.(session).circ_std_pre(first_half+1:end))); 
        matrix_valswr(count_mouse,4) = min((curr_structure.(session).avg_fiber_pre(first_half+1:end))+max(curr_structure.(session).circ_std_pre(first_half+1:end))); % pre is 1
        count_mouse = count_mouse + 1;

        % Post-condition SWR-DA strength
        SWR_DA_strength.post(1,i) = min((curr_structure.(session).avg_fiber_post(first_half+1:end))+ max(curr_structure.(session).circ_std_post(first_half+1:end)));
        matrix_valswr(count_mouse,4) = min((curr_structure.(session).avg_fiber_post(first_half+1:end))+max(curr_structure.(session).circ_std_post(first_half+1:end)));
        count_mouse = count_mouse + 1;

    end
    
    % Save the results in SWR_DA_strength_all under the current structure name
    SWR_DA_strength_all.(curr_structure_name) = SWR_DA_strength;
end

%%
tbl = table(matrix_valswr(:,1),matrix_valswr(:,2),matrix_valswr(:,3),matrix_valswr(:,4),'VariableNames',{'Mouse','Session','PrePost','SWRDA',});

%% Find percentage of sessions that are negative (meaning significant). 
ind_neg = find(tbl.SWRDA<0);
disp('number of sig negative sessions')
disp(length(ind_neg))
percent_sig_sess = (length(ind_neg) / length(tbl.SWRDA))*100; 

disp('percent of sig neg sessions')
disp(percent_sig_sess)

%% Pre sessions
pre = tbl(tbl.PrePost == 1,:); % 1 is pre; 2 is post. 
percent_sig_pre = (length(find(pre.SWRDA<0)) / length(pre.SWRDA))*100; 
disp(percent_sig_pre)
disp(length(find(pre.SWRDA<0)))

%% Post sessions
post = tbl(tbl.PrePost == 2,:); % 1 is pre; 2 is post. 
percent_sig_post = (length(find(post.SWRDA<0)) / length(post.SWRDA))*100; 
disp(percent_sig_post)
disp(length(find(post.SWRDA<0)))
