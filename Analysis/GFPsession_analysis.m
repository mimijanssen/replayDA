%% Saving GFP Data for Session Analysis 
% subject, session, pre/post, SWR-DA
clear; clc;
cd 'D:\GFP_avg'
%load ('colors.mat')

matrix_valswr = zeros(48,4);

%%  Populate Matrix with Mouse Names  

% ~~~~~~~~~~~~~~ SWR-DA ~~~~~~~~~~~~~~
cd 'D:\M556\avg_data\avg_data'
Files=dir('*.*');
count_mouse = 1; 
for k=3:length(Files)
   FileNames=Files(k).name;
   M556sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 9; % 9 is M556 
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 9; % M556 is mouse 9
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

cd 'D:\M578\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M578sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 10; % M578 is mouse 10
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 10; % M578 is mouse 10
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

cd 'D:\M600\avg_data\avg_data'
Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   M600sd.(['sess',num2str(k-2)]) = load(FileNames);
   matrix_valswr(count_mouse,1) = 11; % M600 is mouse 11
   matrix_valswr(count_mouse,3) = 1; % pre is 1
   count_mouse = count_mouse + 1;
   matrix_valswr(count_mouse,1) = 11; % M600 is mouse 11
   matrix_valswr(count_mouse,3) = 2; % post is 2 
   count_mouse = count_mouse + 1;
end

% session needs to be hard coded. 
sessions = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8];

matrix_valswr(:,2) = sessions';

%% SAVE SWR-DA Strength 
SWR_DA_strength = {}; 
count_mouse = 1; 
% List of structure names
structure_names = {'M556sd',  'M578sd', 'M600sd'}; % Add all your structure names here

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
        
        % recalculate std and mean for whole session.
        mu_pre = mean(curr_structure.(session).circ_avg_pre(first_half+1:end));
        mu_post = mean(curr_structure.(session).circ_avg_post(first_half+1:end));
        std_pre = mean(curr_structure.(session).circ_std_pre(first_half+1:end))/2; 
        std_post = mean(curr_structure.(session).circ_std_post(first_half+1:end))/2; 

        % Pre-condition SWR-DA strength
        % subtract the mean and divide by sd and then find the max value.
        % -- actually I think we need to find the max absolute value....
        SWR_DA_strength.pre(1,i) = max((curr_structure.(session).avg_fiber_pre(first_half+1:end)-mu_pre)/std_pre);  % currently this is dividing by 2 sd. So I want to divide by one
        matrix_valswr(count_mouse,4) = max((curr_structure.(session).avg_fiber_pre(first_half+1:end)-mu_pre)/std_pre); % pre is 1
        count_mouse = count_mouse + 1;

        % Post-condition SWR-DA strength
        SWR_DA_strength.post(1,i) = max((curr_structure.(session).avg_fiber_post(first_half+1:end)-mu_post)/std_post);
        matrix_valswr(count_mouse,4) =max((curr_structure.(session).avg_fiber_post(first_half+1:end)-mu_post)/std_post);
        count_mouse = count_mouse + 1;
    end
    
    % Save the results in SWR_DA_strength_all under the current structure name
    SWR_DA_strength_all.(curr_structure_name) = SWR_DA_strength;
end


