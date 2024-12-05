%% load session data into a large matrix 

cd D:\SWR_DA_MegaMatrix

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

%% Linear Mixed Effects Model 
lme_swrevents = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');

%% preprocessed area under the curve
ProcAUCTbl = stack(allTables,{'TwosBeforeAUC','TwosAfterAUC'},'NewDataVariableName','AUC','IndexVariableName','BeforeAfter');
lme_swrevents_auc = fitlme(ProcAUCTbl,'AUC ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');

%% raw data
ProcAUCTbl_raw = stack(allTables,{'TwosBeforeAUCRAW','TwosAfterAUCRAW'},'NewDataVariableName','AUCRAW','IndexVariableName','BeforeAfter');
lme_swrevents_aucRAW = fitlme(ProcAUCTbl_raw,'AUCRAW ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');

%%
PeakTbl_raw = stack(allTables,{'TwosBeforePeakRAW','TwosAfterPeakRAW'},'NewDataVariableName','PeakRAW','IndexVariableName','BeforeAfter');
lme_swrevents_peakRAW = fitlme(PeakTbl_raw,'PeakRAW ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');


%% Post Peak Time... what is the average time and sd for the peak post event
% SAVED THIS WRONG!
% I save the index for  the Pre Peak time... would be interesting though
% to see if the standard deviation decreases

%% plotting things
figure(1)
histogram(allTables.('TwosAfterPeak'));
hold on
histogram(allTables.('TwosBeforePeak'));
legend('after','before')
hold off