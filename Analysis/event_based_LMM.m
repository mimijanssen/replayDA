%% load session data into a large matrix 

cd F:\SWR_DA_MegaMatrix_1s_withSWRinfo
%allTables = load('MegaMatrixALLDATA.mat');
%load('MegaMatrixALLDATA.mat');
%%
allTables = []; % Initialize an empty array for concatenation

Files=dir('*.*');
for k=3:length(Files)
   FileNames=Files(k).name;
   loadedData = load(FileNames); % Load the .mat file

    % Assuming your table is saved as 'matrix_sess' in each file
    if isfield(loadedData, 'matrix_sess')
        sessionTable = loadedData.matrix_sess;

        % Concatenate tables vertically
        if isempty(allTables)
            allTables = sessionTable; % Initialize with the first table
        else
            allTables = [allTables; sessionTable]; % Append subsequent tables
        end
    else
        fprintf('Warning: %s does not contain a table named matrix_sess.\n', FileName);
    end
end



%% Long format with before and after ...  for Before and After Peak (2 seconds) 
ProcPeakTbl = stack(allTables,{'TwosBeforePeak','TwosAfterPeak'},'NewDataVariableName','Peak','IndexVariableName','BeforeAfter');

%%
I_pre = find(allTables.PrePost == 1);  % Find indices where 'condition_row' is positive
test = allTables{I_pre, 'TwosAfterPeak'};
mean_pre = mean(allTables{I_pre, 'TwosAfterPeak'}, 'omitnan');  
mean_pre_time = mean(allTables{I_pre, 'TimeAfterPeak'}, 'omitnan'); 

std_pre = std(allTables{I_pre, 'TwosAfterPeak'}, 'omitnan');  
std_pre_time = std(allTables{I_pre, 'TimeAfterPeak'}, 'omitnan'); 

% Post rest session 
I_post = find(allTables.PrePost == 2);  % Find indices where 'condition_row' is positive
mean_post = mean(allTables{I_post, 'TwosAfterPeak'}, 'omitnan');  
mean_post_time = mean(allTables{I_post, 'TimeAfterPeak'}, 'omitnan'); 

std_post = std(allTables{I_post, 'TwosAfterPeak'}, 'omitnan');  
std_post_time = std(allTables{I_post, 'TimeAfterPeak'}, 'omitnan'); 