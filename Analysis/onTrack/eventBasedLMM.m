%% Add the Early Late variable
allTables2 = Tw; 
list = zeros(height(allTables2),1); 
% early sessions 1-3 = 1
I_early = find(allTables2.sess < 4);  % Find indices where 'condition_row' is positive
list(I_early) = 1; 

% late sessions 6-8 = 2
I_late = find(allTables2.sess > 5);  % Find indices where 'condition_row' is positive
list(I_late) = 2; 

% append list to table 
allTables2.("EarlyLate") = list;

%% change things to categorical
allTables2.PrePost = categorical(allTables2.PrePost);
allTables2.mouseID = categorical(allTables2.mouseID);
allTables2.EarlyLate = categorical(allTables2.EarlyLate);
%allTables2.sleep = categorical(allTables2.sleep);

allTables2.PrePost = reordercats(allTables2.PrePost, ...
    ['2'; setdiff(categories(allTables2.PrePost), {'2'})]);
%
%% Modeling time 
% Full Model 
full = fitlme(allTables2,'PeakOneSecDiff ~ PrePost + EarlyLate + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(full)
% BIC: 54479 
% swrID is sig, PrePost is sig but that's it. 

% Base Model 
base = fitlme(allTables2,'PeakOneSecDiff ~ 1 + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(base)
% BIC 54442

% ~ PREPOST ~
prepost = fitlme(allTables2,'PeakOneSecDiff ~ PrePost + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(prepost)
noprepost = fitlme(allTables2,'PeakOneSecDiff ~ EarlyLate + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noprepost);
% base vs. prepost
compare(base, prepost,'nsim',1000)
% full vs. noprepost 
compare(noprepost,full,'nsim',1000)

% ~ EarlyLate ~
earlylate = fitlme(allTables2,'PeakOneSecDiff ~ EarlyLate + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(earlylate)
noearlylate = fitlme(allTables2,'PeakOneSecDiff ~ PrePost + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noearlylate);
% base vs. prepost
compare(base, earlylate,'nsim',1000)
% full vs. noprepost 
compare(noearlylate,full,'nsim',1000)


% ~ duration ~
dur = fitlme(allTables2,'PeakOneSecDiff ~ SWRdur + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(dur)
nodur = fitlme(allTables2,'PeakOneSecDiff ~ PrePost + EarlyLate + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nodur);
% base vs. prepost
compare(base, dur,'nsim',1000)
% full vs. noprepost 
compare(nodur,full,'nsim',1000)

% ~ swrpower ~
power = fitlme(allTables2,'PeakOneSecDiff ~ SWRpower + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(power)
nopower = fitlme(allTables2,'PeakOneSecDiff ~ PrePost + EarlyLate + SWRdur + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nopower);
% base vs. prepost
compare(base, power,'nsim',1000)
% full vs. noprepost 
compare(nopower,full,'nsim',1000)

% ~ swrID ~
id = fitlme(allTables2,'PeakOneSecDiff ~ swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(id)
noid = fitlme(allTables2,'PeakOneSecDiff ~ PrePost + EarlyLate + SWRdur + SWRpower + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noid);
% base vs. prepost
compare(base, id,'nsim',1000)
% full vs. noprepost 
compare(noid, full,'nsim',1000)

%%

% Full Model 
full = fitlme(allTables2,'MeanOneSecDiff ~ PrePost + EarlyLate + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(full)
% BIC: 54479 
% swrID is sig, PrePost is sig but that's it. 

% Base Model 
base = fitlme(allTables2,'MeanOneSecDiff ~ 1 + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(base)
% BIC 54442

% ~ PREPOST ~
prepost = fitlme(allTables2,'MeanOneSecDiff ~ PrePost + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(prepost)
noprepost = fitlme(allTables2,'MeanOneSecDiff ~ EarlyLate + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noprepost);
% base vs. prepost
compare(base, prepost,'nsim',1000)
% full vs. noprepost 
compare(noprepost,full,'nsim',1000)

% ~ EarlyLate ~
earlylate = fitlme(allTables2,'MeanOneSecDiff ~ EarlyLate + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(earlylate)
noearlylate = fitlme(allTables2,'MeanOneSecDiff ~ PrePost + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noearlylate);
% base vs. prepost
compare(base, earlylate,'nsim',1000)
% full vs. noprepost 
compare(noearlylate,full,'nsim',1000)

% ~ duration ~
dur = fitlme(allTables2,'MeanOneSecDiff ~ SWRdur + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(dur)
nodur = fitlme(allTables2,'MeanOneSecDiff ~ PrePost + EarlyLate + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nodur);
% base vs. prepost
compare(base, dur,'nsim',1000)
% full vs. noprepost 
compare(nodur,full,'nsim',1000)

% ~ swrpower ~
power = fitlme(allTables2,'MeanOneSecDiff ~ SWRpower + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(power)
nopower = fitlme(allTables2,'MeanOneSecDiff ~ PrePost + EarlyLate + SWRdur + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nopower);
% base vs. prepost
compare(base, power,'nsim',1000)
% full vs. noprepost 
compare(nopower,full,'nsim',1000)

% ~ swrID ~
id = fitlme(allTables2,'MeanOneSecDiff ~ swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(id)
noid = fitlme(allTables2,'MeanOneSecDiff ~ PrePost + EarlyLate + SWRdur + SWRpower + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noid);
% base vs. prepost
compare(base, id,'nsim',1000)
% full vs. noprepost 
compare(noid,full,'nsim',1000)