cd D:\
load('MegaMatrix1s.mat')

%% mega table 1 s 
% average peak value ~
% Pre rest session
I_pre = find(allTables.PrePost == 1);  % Find indices where 'condition_row' is positive
test = allTables{I_pre, 'OnesAfterPeak'};
mean_pre = mean(allTables{I_pre, 'OnesAfterPeak'}, 'omitnan');  
mean_pre_time = mean(allTables{I_pre, 'TimeAfterPeak'}, 'omitnan'); 

std_pre = std(allTables{I_pre, 'OnesAfterPeak'}, 'omitnan');  
std_pre_time = std(allTables{I_pre, 'TimeAfterPeak'}, 'omitnan'); 

% Post rest session 
I_post = find(allTables.PrePost == 2);  % Find indices where 'condition_row' is positive
mean_post = mean(allTables{I_post, 'OnesAfterPeak'}, 'omitnan');  
mean_post_time = mean(allTables{I_post, 'TimeAfterPeak'}, 'omitnan'); 

std_post = std(allTables{I_post, 'OnesAfterPeak'}, 'omitnan');  
std_post_time = std(allTables{I_post, 'TimeAfterPeak'}, 'omitnan'); 

%%
ProcPeakTbl = stack(allTables,{'OnesBeforePeak','OnesAfterPeak'},'NewDataVariableName','Peak','IndexVariableName','BeforeAfter');

%% LMS-- BEFORE AND AFTER
lme_swrevt = fitlme(ProcPeakTbl,'Peak ~ PrePost + swrID+ (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt)
% AIC: 1.0814e+05 

lme_swrevt4 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + PrePost + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt4)
% AIC: 1.0802e+05

compare(lme_swrevt, lme_swrevt4, 'nsim',1000)

%% PRE AND POST 

lme_swrevt2 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + swrID+ (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt2)
% AIC: 1.0814e+05 

lme_swrevt5 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + PrePost + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt5)
% AIC: 1.0802e+05

compare(lme_swrevt2, lme_swrevt5,'nsim',1000)

%% Late vs. Early Sessions
list = zeros(height(ProcPeakTbl),1); 
% early sessions 1-4 = 1
I_early = find(ProcPeakTbl.sess < 5);  % Find indices where 'condition_row' is positive
list(I_early) = 1; 

% late sessions 5-6 = 2
I_late = find(ProcPeakTbl.sess > 4);  % Find indices where 'condition_row' is positive
list(I_late) = 2; 

% append list to table 
ProcPeakTbl.("EarlyLate") = list;

%% LMM
lme_swrevt7 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + PrePost + swrID+ (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt7)
% AIC: 1.0802e+05 

lme_swrevt8 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + PrePost + EarlyLate + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt8)
% AIC: 1.0802e+05

compare(lme_swrevt7, lme_swrevt8, 'nsim',1000)

%% OTHER GLMS TO COMPARE AIC 
% lme_swrevt1 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
% disp(lme_swrevt1)
% % NOTES: 
% % SWRID is significant.
% % AIC: 1.084e+05
% % random effects do not include 0
% 
% % 2) model 2: random intercepts for both mouse and session with session nested within
% % mouse 
% lme_swrevt2 =  fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess:mouseID)');
% disp(lme_swrevt2)
% % NOTES: 
% % AIC: 1.0802e+05
% % SWR ID is significant 
% % random effects do not include 0
% 
% % model 3: random intercepts and slopes for session nested within mouse
% lme_swrevt3 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1+ sess|mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
% disp(lme_swrevt3)
% AIC: 1.0834e+05

%% Plot all raw data 
list = zeros(height(ProcPeakTbl),1); 
% early sessions 1-4 = 1
I_pre = find(ProcPeakTbl.PrePost == 1);  % Find indices where 'condition_row' is positive
Pre_table = allTables(I_pre)



%% pre-track rest
% Initialize an empty array for storing valid signals

time = ProcPeakTbl.OnesPreProc{1,1}.tvec- ProcPeakTbl.OnesPreProc{1,1}.tvec(1); 

validSignals = [];

% Loop through the table entries
for i = 1:height(ProcPeakTbl)
    if ~isempty(ProcPeakTbl.OnesPreProc{i,1}) && isstruct(ProcPeakTbl.OnesPreProc{i,1}) ...
            && isfield(ProcPeakTbl.OnesPreProc{i,1}, 'signal')
        % Extract the signal
        signal = ProcPeakTbl.OnesPreProc{i,1}.signal;
        
        % Ensure signal is a column vector
        if isrow(signal)
            signal = signal'; % Convert to column vector if needed
        end
        
        % Store the signal in an array (assuming all signals have the same length)
        validSignals = [validSignals, signal]; % Concatenates column-wise
    end
end

avgSignal = mean(validSignals, 2); % Average across columns (each row is a time point)
stdSignal = std(validSignals, 0,2); 

    % Plot the average signal
figure (2);
%plot(avgSignal, 'LineWidth', 2);
%shadedErrorBar(time,avgSignal,stdSignal,'lineProps','-k','transparent',1) % subtract the circ mean here 
hold on
plot(time,avgSignal,'LineWidth',2,'Color','k') % subtract the circ mean here 

    %shadedErrorBar(time_extract_pre(1,:),circ_avg_fiber_pre,circ_std_fiber_pre,'lineProps','-k','transparent',1) % subtract the circ mean here 
xlabel('Time from SWR (s)');
ylabel('Fiber Signal (z-score)');
xticks([0 0.5 1])
xticklabels({'-0.5','0','0.5'})
title('Pre-task Rest');
hold off

%% Post task
validSignals = [];

% Loop through the table entries
for i = 1:height(ProcPeakTbl)
    if ~isempty(ProcPeakTbl.OnesPostProc{i,1}) && isstruct(ProcPeakTbl.OnesPostProc{i,1}) ...
            && isfield(ProcPeakTbl.OnesPostProc{i,1}, 'signal')
        % Extract the signal
        signal = ProcPeakTbl.OnesPostProc{i,1}.signal;
        
        % Ensure signal is a column vector
        if isrow(signal)
            signal = signal'; % Convert to column vector if needed
        end
        
        % Store the signal in an array (assuming all signals have the same length)
        validSignals = [validSignals, signal]; % Concatenates column-wise
    end
end

avgSignal = mean(validSignals, 2); % Average across columns (each row is a time point)
stdSignal = std(validSignals, 0,2); 

    % Plot the average signal
figure (3);
%plot(avgSignal, 'LineWidth', 2);
shadedErrorBar(time,avgSignal,stdSignal,'lineProps','-k','transparent',1) % subtract the circ mean here 
hold on
plot(time,avgSignal,'LineWidth',2,'Color','k') % subtract the circ mean here 

%shadedErrorBar(time_extract_pre(1,:),circ_avg_fiber_pre,circ_std_fiber_pre,'lineProps','-k','transparent',1) % subtract the circ mean here 
xlabel('Time from SWR (s)');
ylabel('Fiber Signal (z-score)');
xticks([0 0.5 1])
xticklabels({'-0.5','0','0.5'})
title('Post-task Rest');
hold off