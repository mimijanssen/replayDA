%% load session data into a large matrix 

cd F:\SWR_DA_MegaMatrix_1s
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


%% Long format with before and after ...  for Before and After Peak (2 seconds) 
ProcPeakTbl = stack(allTables,{'OnesBeforePeak','OnesAfterPeak'},'NewDataVariableName','Peak','IndexVariableName','BeforeAfter');

%% Linear Mixed Effects Model 
% 1) model 1: random intercepts for both mouse and session 
lme_swrevt1 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt1)
% NOTES: 
% SWRID is not significant. - take it out
% AIC: 1.1287e+05 
% BIC: 1.1293e+05
% random effects do not include 0

% 2) model 2: random intercepts for both mouse and session with session nested within
% mouse 
lme_swrevt2 =  fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess:mouseID)');
disp(lme_swrevt2)
% NOTES: 
% AIC: 1.122e+05 
% BIC: 1.1226e+05 -- nested intercpet for sess is better 
% SWR ID is not significant 
% random effects do not include 0

% model 3: random intercepts and slopes for session nested within mouse
lme_swrevt3 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1+ sess|mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt3)
% NOTES: 
% AIC: 1.126e+05
% BIC: 1.1268e+05 -- nested slope for sess is worse
% SWRID is not significant 

% model 4: no swrID
lme_swrevt4 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + PrePost + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt4)
% AIC: 1.122e+05
% BIC: 1.1225e+05
% * best model* 

% no random effects for session 
lme_swrevt5 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter +  PrePost + (1|mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt5)
% AIC: 1.131e+05
% BIC: 1.1314e+05 -- taking out sess made the AIC worse. 
% random effects do not include 0

% no random effects
lme_swrevt6 = fitlme(ProcPeakTbl,'Peak ~ BeforeAfter + PrePost');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(lme_swrevt6)
% AIC:  1.1489e+05
% BIC: 1.1493e+05 -- even worse without mosue. 

compare(lme_swrevt2, lme_swrevt4)
% 4 has a lower BIC but just barely 

% GOING OVER MODELS FOR AUC AND COMPARING AIC BIC: 
lme_swrevt7 = fitlme(ProcAUCTbl,'AUC ~ BeforeAfter + PrePost + (1|mouseID) + (1|sess:mouseID)');
% AIC: 8.1615e+05
% BIC: 8.162e+05

compare(lme_swrevt4,lme_swrevt7)
% swrevt7 additions significantly makes the model better (AUC over peak
% values) 

% GOING OVER PEAK RAW DATA AND COMPARING THE AIC BIC: 
lme_swrevt8 = fitlme(PeakTbl_raw,'PeakRAW ~ BeforeAfter + PrePost + (1|mouseID) + (1|sess:mouseID)');
%AIC: -1.5811e+05; 
%BIC: -1.5806e+05;

lme_swrevt9 = fitlme(PeakTbl_raw,'PeakRAW ~ BeforeAfter + swrID + PrePost + swrID*PrePost + (1|mouseID) + (1|sess:mouseID)');
disp(lme_swrevt9)
% AIC: -1.5976e+05 % This has the lower AIC And BIC so it is better fit
% BIC: -1.5969e+05 

lme_swrevt10 = fitlme(PeakTbl_raw,'PeakRAW ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess:mouseID)');
disp(lme_swrevt10)
% AIC:  -1.5829e+05 -- worse fit...
% BIC: -1.5823e+05

compare(lme_swrevt4,lme_swrevt8)
% using raw values significantly mekes the model perform better
compare(lme_swrevt8,lme_swrevt9)
% adding swrID significantly makes the model perform better 

% GOING OVER AUC RAW DATA AND COMPARING 
lme_swrevt11 = fitlme(ProcAUCTbl_raw,'AUCRAW ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess:mouseID)');
disp(lme_swrevt11)
% AIC: 5.6949e+05 -- really high AIC 
% BIC: 5.6955e+05

% if confidence intervals of random effects include zero, then it is not
% significant. if it includes 0 then these terms are significant and needed
% in your model? Formally test random effects using the compare function. 

% AIC BIC for these models:

% idenitfy the significane of random effects by inspecting the random
% effects covariance parameters 

% perform ANOVA to assess fixed effects
anovaResutls = anova(lme_swrevt9); 
disp(anovaResutls)
    % 
    % Term                     FStat     DF1    DF2      pValue     
    % {'(Intercept)'  }        261.35    1      48515     1.2399e-58
    % {'swrID'        }        1644.2    1      48515              0
    % {'PrePost'      }         13250    1      48515              0
    % {'BeforeAfter'  }        40.315    1      48515     2.1804e-10
    % {'swrID:PrePost'}        1491.8    1      48515    1.8182e-321
    % 

%% Other LMEs... 
% %% preprocessed area under the curve
ProcAUCTbl = stack(allTables,{'TwosBeforeAUC','TwosAfterAUC'},'NewDataVariableName','AUC','IndexVariableName','BeforeAfter');
lme_swrevents_auc = fitlme(ProcAUCTbl,'AUC ~ BeforeAfter + swrID + PrePost + (sess|mouseID)');
% 
% %% raw data
ProcAUCTbl_raw = stack(allTables,{'TwosBeforeAUCRAW','TwosAfterAUCRAW'},'NewDataVariableName','AUCRAW','IndexVariableName','BeforeAfter');
% lme_swrevents_auc = fitlme(ProcAUCTbl_raw,'AUCRAW ~ BeforeAfter + swrID + PrePost + (sess|mouseID)');
% 
% %%
PeakTbl_raw = stack(allTables,{'TwosBeforePeakRAW','TwosAfterPeakRAW'},'NewDataVariableName','PeakRAW','IndexVariableName','BeforeAfter');
% lme_swrevents_peakRAW = fitlme(PeakTbl_raw,'PeakRAW ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');


%% ~ ASSUMPTIONS ~

% Testing Linearity : Plotting fitted resposne vs. obsereved response and
% residuals 
F = fitted(lme_swrevt9);
R = response(lme_swrevt9);
figure(1);
plot(R,F,'rx')
xlabel('Response')
ylabel('Fitted')
% Preprocessed: evt 4 -- This does not look like a linear fit. The fitted vs. observed resposne
% values do not form a 45 degree angle which would indicate good fit. 
% --> RAW peak values: evt 9 does look like a linear fit. 

% Testing Heteroscedasticity : Plotting residuals vs. fitted values. 
figure(3)
plotResiduals(lme_swrevt9,'fitted') 
% no immediate signs of heteroscedasticity 
% --> RAW peak values: do have some outliers

% --> grouped by mouse
%figure(3);
%gscatter(F,R,PeakTbl_raw.('mouseID')); % looks better for auc

%figure(3);
%gscatter(F,R,PeakTbl_raw.('sess')); % looks better for auc
% m548 has some outliers

% Testing if Residual Errors are normally distributed --need to do the same for Random Effects
figure(4) 
plotResiduals(lme_swrevt9,'probability') 
% --> Raw peak values: mostly normal, but tons of outliers 
high_outliers = find(residuals(lme_swrevt9)>0.15);
low_outliers = find(residuals(lme_swrevt9)<-0.2);

figure(6)
mouse_outliers_high = PeakTbl_raw.('mouseID')(high_outliers);
histogram(mouse_outliers_high); 
% some for each mouse, but dominated by Mouse 6 - M545

mouse_outliers_low = PeakTbl_raw.('mouseID')(low_outliers);
histogram(mouse_outliers_low);
% only mouse 8 - M548!

% any specific session? 
mouse_outliers_low2 = PeakTbl_raw.('sess')(low_outliers);
histogram(mouse_outliers_low2);
% ALL SESSION 6!!! Maybe consider removing this? 

% Testing to see if residuals are correlated
figure(5)
plotResiduals(lme_swrevt9,'lagged')
% most residuals do not seem to have a mild positive correlation
% there is one bunch that does seem to be. - M548

% r = residuals(lme_swrevents);
% pr = residuals(lme_swrevents,'ResidualType','Pearson');
% st = residuals(lme_swrevents,'ResidualType','Standardized');
% X = [r pr st];
% boxplot(X,'labels',{'Raw','Pearson','Standardized'})
% 
% % histogram of mouse outliers... 
% mouse_outliers = ProcPeakTbl.('mouseID')(outliers);
% histogram(mouse_outliers); % mouse 2 seems to have the mouse outliers. (M453)
% 
% % histogram of mouse non outliers
% ProcPeakTbl2 = ProcPeakTbl;
% ProcPeakTbl2(outliers,:)=[];
% mouse_outliers = ProcPeakTbl2.('mouseID');
% histogram(mouse_outliers); % mouse 2 seems to have the mouse outliers. (M453)
% % mouse 1, 2, and 8 make up a lot of the normal data

 figure(1)
 histogram(allTables.('TwosBeforeAUCRAW'));
 hold on
 histogram(allTables.('TwosAfterAUCRAW'));
 legend('before','after')
 ylabel('Frequency')
 xlabel('Peak')
 hold off
% % data seems to be normal
% 
% figure(2) 
% plotResiduals(lme_swrevents)

%% plot regular linear regression 
mdl = fitlm(PeakTbl_raw,'PeakRAW ~ BeforeAfter + swrID + PrePost + swrID*PrePost'); 

%% Fitting GLMS 
glme_raw = fitglme(PeakTbl_raw,'PeakRAW ~ BeforeAfter + swrID + PrePost + swrID*PrePost + (1|mouseID) + (1|sess:mouseID)'); 
glme_proc =fitglme(ProcPeakTbl,'Peak ~ BeforeAfter + PrePost + (1|mouseID) + (1|sess:mouseID)'); 

disp(glme_raw)
disp(glme_proc)
% raw is still better. 

%% OK NOW TRY TO Z-score the data and run this again???


%% Find average time 
% pre
pre_time_tbl = ProcPeakTbl(ProcPeakTbl.PrePost ==1,:);
avg_pre_time = mean(pre_time_tbl.TimeAfterPeak);
std_pre_time = std(pre_time_tbl.TimeAfterPeak);

% post
post_time_tbl = ProcPeakTbl(ProcPeakTbl.PrePost ==2,:);
avg_post_time = mean(post_time_tbl.TimeAfterPeak);
std_post_time = std(post_time_tbl.TimeAfterPeak);