%% load session data into a large matrix 

cd D:\SWR_DA_MegaMatrix
%allTables = load('MegaMatrixALLDATA.mat')

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

%% Edit the table so that it only saves 1 second of data before swrs and after swrs 
% also edit the peak values 
Ones_table = allTables; 

% Edit the Ones_table.FoursPreRaw{1,1}.signal and
% Ones_table.FoursPreRaw{1,1}.time 
for i = 1:1:height(Ones_table)
    if Ones_table.PrePost{i} == 1
        Ones_table.FoursPreRaw{i,1}.signal(1:3000)= [];
        Ones_table.FoursPreRaw{i,1}.tvec(1:3000)= [];
        Ones_table.FoursPreProc{i,1}.signal(1:3000)= [];
        Ones_table.FoursPreProc{i,1}.tvec(1:3000)= [];

        % change the max values 
        Ones_table.TwosBeforePeak{i} = max(Ones_table.FoursPreProc{i,1}.signal);

end

% T = renamevars(T,["Var1","Var2"],["Nums1","Nums2"])



%% Long format with before and after ...  for Before and After Peak (2 seconds) 
ProcPeakTbl = stack(allTables,{'TwosBeforePeak','TwosAfterPeak'},'NewDataVariableName','Peak','IndexVariableName','BeforeAfter');
