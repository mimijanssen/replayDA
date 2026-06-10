cd F:\
load('MegaMatrix1s.mat')

%% mega table 1 s 
% average peak value ~
% Pre rest session
I_pre = find(allTables.PrePost == 1);  % Find indices where 'condition_row' is positive
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

% change I_pre to pre and I_post to post


%%
ProcPeakTbl = stack(allTables,{'OnesBeforePeak','OnesAfterPeak'},'NewDataVariableName','Peak','IndexVariableName','BeforeAfter');

%% convert post and pre to categorical 
ProcPeakTbl.PrePost = categorical(ProcPeakTbl.PrePost);
ProcPeakTbl.mouseID = categorical(ProcPeakTbl.mouseID);
ProcPeakTbl.BeforeAfter = categorical(ProcPeakTbl.BeforeAfter);

%% LMS-- BEFORE AND AFTER
lme_swrevt = fitlme(ProcPeakTbl,'Peak ~ 1 + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt)
% AIC: 1.0814e+05 

lme_swrevt4 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt4)
% AIC: 1.0802e+05

%%
ProcPeakTbl.PrePost = reordercats(ProcPeakTbl.PrePost, ...
    ['2'; setdiff(categories(ProcPeakTbl.PrePost), {'2'})]);
lme_swrevt5 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter*PrePost + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt5)

compare(lme_swrevt, lme_swrevt4, 'nsim',1000)


%% PRE AND POST 

lme_swrevt2 = fitlme(ProcPeakTbl,'Peak ~ 1 + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt2)
% AIC: 1.0814e+05 

lme_swrevt5 = fitlme(ProcPeakTbl,'Peak ~ PrePost + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt5)
% AIC: 1.0802e+05

compare(lme_swrevt2, lme_swrevt5,'nsim',1000)

% for these it is important to make sure the beforeafter and prepost are
% categorical variables. MATLAB defines the first as the reference category
% and the second as the comparative category. If there is a positive beta,
% that means the dependent variable is higher for this level compared to
% the reference. 
%% Late vs. Early Sessions
list = zeros(height(ProcPeakTbl),1); 
% early sessions 1-4 = 1
I_early = find(ProcPeakTbl.sess < 4);  % Find indices where 'condition_row' is positive
list(I_early) = 1; 

% late sessions 5-6 = 2
I_late = find(ProcPeakTbl.sess > 5);  % Find indices where 'condition_row' is positive
list(I_late) = 2; 

% append list to table 
ProcPeakTbl.("EarlyLate") = list;

%% LMM
lme_swrevt7 = fitlme(ProcPeakTbl,'Peak ~ 1+ (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt7)
% AIC: 1.0802e+05 

lme_swrevt8 = fitlme(ProcPeakTbl,'Peak ~  EarlyLate + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt8)
% AIC: 1.0802e+05

compare(lme_swrevt7, lme_swrevt8, 'nsim',1000)

%% Full model 
% no diff between 1|sess and 1|sess:mouse
% swrID doesn't add anything to the full model. 
full2 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + PrePost + EarlyLate + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(full2)

full3 = fitlme(ProcPeakTbl,'Peak ~ PrePost + EarlyLate + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(full3)

full4 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + EarlyLate + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(full4)

full5 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + PrePost + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(full5)

compare(full2, full3, 'nsim',1000)
compare(full2,full4,'NSim',1000) % no prepost
compare(full2,full5,'NSim',1000) % no early late


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

%% Plot the event based average: 
% Define time vector (assuming all should be the same length, e.g., 1001 points)
common_tvec = linspace(-1, 1, 1001); % Modify this to match the expected range

% Extract unique sessions and mice
unique_sess = unique(ProcPeakTbl.sess);
unique_mice = unique(ProcPeakTbl.mouseID);

% Initialize storage for Pre and Post

allPreProc = [];
allPostProc = [];

% Loop through each row in the table
for i = 1:height(ProcPeakTbl)
    if ProcPeakTbl.PrePost(i) == 1 % Pre session
        allPreProc = [allPreProc; ProcPeakTbl.OnesPreProc{i,1}.signal'];
    elseif ProcPeakTbl.PrePost(i) == 2 % Post session
        allPostProc = [allPostProc; ProcPeakTbl.OnesPostProc{i,1}.signal'];
    end
end

% Compute averages
meanPreProc = mean(allPreProc, 1, 'omitnan');
meanPostProc = mean(allPostProc, 1, 'omitnan');

% Plot overall averages
figure(1);
plot(common_tvec, meanPreProc, 'b', 'LineWidth', 1.5); hold on;
plot(common_tvec, meanPostProc, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Processed Signal'); title('Overall: Pre vs Post');
legend({'Pre', 'Post'});

% Plot per session
figure(2)
hold on;
for s = 1:length(unique_sess)
    sess_idx = (ProcPeakTbl.sess == unique_sess(s));
    preSess = []; postSess = [];
    
    for i = find(sess_idx)'
        if ProcPeakTbl.PrePost(i) == 1
            preSess = [preSess; ProcPeakTbl.OnesPreProc{i,1}.signal'];
        elseif ProcPeakTbl.PrePost(i) == 2
            postSess = [postSess; ProcPeakTbl.OnesPostProc{i,1}.signal'];
        end
    end
    
    plot(common_tvec, mean(preSess, 2, 'omitnan')', 'b');
    plot(common_tvec, mean(postSess, 2, 'omitnan')', 'r');
end
xlabel('Time (s)'); ylabel('Proc Signal'); title('Per Session Average');
legend({'Pre', 'Post'});

% Plot per mouse
subplot(2,2,4); hold on;
for m = 1:length(unique_mice)
    mouse_idx = (ProcPeakTbl.mouseID == unique_mice(m));
    preMouse = []; postMouse = [];
    
    for i = find(mouse_idx)'
        if ProcPeakTbl.PrePost(i) == 1
            preMouse = [preMouse, interp1(ProcPeakTbl.OnesPreProc{i,1}.tvec, ProcPeakTbl.OnesPreProc{i,1}.signal, common_tvec, 'linear', 'extrap')];
        elseif ProcPeakTbl.PrePost(i) == 2
            postMouse = [postMouse, interp1(ProcPeakTbl.OnesPosProc{i,1}.tvec, ProcPeakTbl.OnesPostProc{i,1}.signal, common_tvec, 'linear', 'extrap')];
        end
    end
    
    plot(common_tvec, mean(preMouse, 2, 'omitnan'), 'b');
    plot(common_tvec, mean(postMouse, 2, 'omitnan'), 'r');
end
xlabel('Time (s)'); ylabel('Proc Signal'); title('Per Mouse Average');
legend({'Pre', 'Post'});

sgtitle('Comparison of Pre vs Post Sessions Across Sessions and Mice');

