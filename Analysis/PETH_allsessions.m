%%
clc 
clear
rng(10) %NOTE FORGOT TO RUN THIS FOR EARLY 1-3 and LATE 6-8 sessions
%%
cd ('D:\SWR_DA_MegaMatrix_4s')

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

%% Table with NREM and Wake
% 24277 
% has a better SWR power as well so use that! 
% allTables has 24277 
allTables2 = allTables; 
%allTables2.sleep = sleep.allTables.sleep;
%allTables2.betterpower = sleep.allTables.swrp;

%% Find base_peak 
% AFTER SWR 
allTables2.PostSWRPeak = NaN(height(allTables2),1);
allTables2.PostSWRmean = NaN(height(allTables2),1);
for i = 1:height(allTables2)
    if allTables2.PrePost(i) == 1
        signal = allTables2.TwosPreProc{i}.signal;
    elseif allTables2.PrePost(i) == 2
        signal = allTables2.TwosPostProc{i}.signal;
    else
        continue
    end
    allTables2.PostSWRPeak(i) = max(signal(4000:5001)); % NEED TO CHANGE THIS
    allTables2.PostSWRmean(i) = mean(signal(4000:5001));
end

% Before SWR
allTables2.PreSWRPeak = NaN(height(allTables2),1);
allTables2.PreSWRmean = NaN(height(allTables2),1);
for i = 1:height(allTables2)
    if allTables2.PrePost(i) == 1
        signal = allTables2.TwosPreProc{i}.signal;
    elseif allTables2.PrePost(i) == 2
        signal = allTables2.TwosPostProc{i}.signal;
    else
        continue
    end
    allTables2.PreSWRPeak(i) = max(signal(3000:4001));
    allTables2.PreSWRmean(i) = mean(signal(3000:4001));
end

allTables2.Base_peak = allTables2.PostSWRPeak - allTables2.PreSWRPeak;
allTables2.Base_mean = allTables2.PostSWRmean - allTables2.PreSWRmean;

%% Add the Early Late variable 
list = zeros(height(allTables2),1); 
% early sessions 1-3 = 1
I_early = find(allTables2.sess < 5);  % Find indices where 'condition_row' is positive
list(I_early) = 1; 

% late sessions 6-8 = 2
I_late = find(allTables2.sess > 4);  % Find indices where 'condition_row' is positive
list(I_late) = 2; 

% append list to table 
allTables2.("EarlyLate") = list;

%% Loop through and make session averages? 

uni_mouse = unique(allTables2.mouseID);
uni_sess = unique(allTables2.sess);
uni_prepost = unique(allTables2.PrePost);  % Should be [1,2] for Pre/Post
uni_sleep = unique(allTables2.sleep);  
uni_earlylate = unique(allTables2.EarlyLate);  


% Initialize sess_avg_tbl with the same column names but no data
mouseID = 0;
sess = 0;
PrePost = 0;
sleep = 0;
earlylate = 0;
avg_signal = [];
temp_signal = [];

%sess_avg_tbl = table(mouseID, sess, PrePost, Peak_one_sec2); 
sess_avg_tbl = [];

% Loop over unique mice, sessions, sleep, earlylate, and Pre/Post categories
for i_mouse = 1:length(uni_mouse)
    for i_sess = 1:length(uni_sess)
        list_allsess = find(allTables2.sess == uni_sess(i_sess) & ...
        allTables2.mouseID == uni_mouse(i_mouse)); % create a list of all swrs within a session for a mouse.
        for i_prepost = 1:length(uni_prepost) % loops through both pre and post 
            for i_sleep = 1:length(uni_sleep) 
                for i_earlylate = 1:length(uni_earlylate)
                % Find indices where mouseID, session,PrePost, and beforeafter match
                    try
                        list = find(allTables2.sess == uni_sess(i_sess) & ...
                        allTables2.mouseID == uni_mouse(i_mouse) & ...
                        allTables2.PrePost == uni_prepost(i_prepost) &... 
                        allTables2.sleep == uni_sleep(i_sleep) &...
                        allTables2.EarlyLate == uni_earlylate(i_earlylate));
                        if ~isempty(list) % if that session exists. 
                            for i_list = 1:1:length(list) % iterate through all structures 
                                if i_prepost == 1
                                    signal_save = allTables2.TwosPreProc{list(i_list),1}.signal; % signals are saved in different columns but I'm saving everything 
                                    temp_signal = [temp_signal; signal_save'];
                                end
                                if i_prepost == 2
                                    signal_save = allTables2.TwosPostProc{list(i_list),1}.signal;
                                    temp_signal = [temp_signal; signal_save'];
                                end
                            end
                        end
                    % Create a new row to append
                    %new_row = {uni_mouse(i_mouse), uni_sess(i_sess), uni_prepost(i_prepost), avg_peak}; % for main analysis, extend this to also include before ... save avg signal here 
                    % Append the row to sess_avg_tbl
                    %sess_avg_tbl = [sess_avg_tbl; new_row];
                    catch
                        disp('no combo')
                    end
                    
                    avg_signal = mean(temp_signal); % average over all signal 

                    if ~isempty(list_allsess) % if that session exists 
                        avg_base_peak = mean(allTables2.Base_peak(list));
                        avg_base_mean = mean(allTables2.Base_mean(list));
                        avg_time = mean(allTables2.TimeAfterPeak(list));
                        new_row = {uni_mouse(i_mouse), uni_sess(i_sess), uni_prepost(i_prepost),uni_sleep(i_sleep),uni_earlylate(i_earlylate), avg_base_peak, avg_base_mean, avg_time, avg_signal}; 
                        sess_avg_tbl = [sess_avg_tbl; new_row];
                        temp_signal = [];
                    end
                end
            end
        end
    end
end

sess_avg_tbl = cell2table(sess_avg_tbl, 'variablenames',{'mouse','sess','prepost','sleep','earlylate','peak','mean','time','signal'});
% LOTS OF NO COMBOS
% might have to do a basic one with separating out each variable at a time.
sess_avg_tbl2 = sess_avg_tbl;
%% LMMS - make things categorical
sess_avg_tbl2.earlylate = categorical(sess_avg_tbl2.earlylate);
sess_avg_tbl2.sleep = categorical(sess_avg_tbl2.sleep);
sess_avg_tbl2.prepost = categorical(sess_avg_tbl2.prepost);
sess_avg_tbl2.mouse = categorical(sess_avg_tbl2.mouse);
sess_avg_tbl2.sess = categorical(sess_avg_tbl2.sess);

%% Remove nans: 
include=~isnan(sess_avg_tbl2.peak(:,1));
nonan_Tbl = sess_avg_tbl2(include,:);

%% Need to remove 0s for sleep 
nonan_Tbl(nonan_Tbl.sleep == '0', :) = [];
nonan_Tbl.sleep = removecats(nonan_Tbl.sleep);
nonan_Tbl.sleep = reordercats(nonan_Tbl.sleep, {'1','2'});

% need to remove 0 
nonan_Tbl(nonan_Tbl.earlylate == '0', :) = [];
nonan_Tbl.earlylate = removecats(nonan_Tbl.earlylate);
nonan_Tbl.earlylate = reordercats(nonan_Tbl.earlylate, {'1','2'});

%% PEAK analysis 
% intercept
base = fitlme(nonan_Tbl,'peak ~ 1 + (1|mouse) + (1|sess)'); 
disp(base)

% Full model
full = fitlme(nonan_Tbl,'peak ~ prepost + sleep + (1|mouse) + (1|sess)'); 
disp(full)

% ~ PREPOST ~
prepost = fitlme(nonan_Tbl,'peak ~ prepost + (1|mouse) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(prepost)
noprepost = fitlme(nonan_Tbl,'peak ~ sleep + earlylate + (1|mouse) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noprepost);
% base vs. prepost
compare(base, prepost,'nsim',1000)
% full vs. noprepost 
compare(noprepost,full,'nsim',1000)

% ~ sleep ~
sleep = fitlme(nonan_Tbl,'peak ~ sleep + (1|mouse) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(sleep)
nosleep = fitlme(nonan_Tbl,'peak ~ prepost + earlylate + (1|mouse) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nosleep);
% base vs. prepost
compare(base, sleep,'nsim',1000)
% full vs. noprepost 
compare(nosleep,full,'nsim',1000)

% ~ earlylate ~
earlylate = fitlme(nonan_Tbl,'peak ~ earlylate + (1|mouse) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(earlylate)
noearlylate = fitlme(nonan_Tbl,'peak ~ prepost + sleep + (1|mouse) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noearlylate);
% base vs. prepost
compare(base, earlylate,'nsim',1000)
% full vs. noprepost 
compare(noearlylate,full,'nsim',1000)


%% MEAN analysis 
% intercept
base = fitlme(nonan_Tbl,'mean ~ 1 + (1|mouse) + (1|sess)'); 
disp(base)

% Full model
full = fitlme(nonan_Tbl,'mean ~ prepost + sleep + earlylate + (1|mouse) + (1|sess)'); 
disp(full)

% ~ PREPOST ~
prepost = fitlme(nonan_Tbl,'mean ~ prepost + (1|mouse) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(prepost)
noprepost = fitlme(nonan_Tbl,'mean ~ sleep + earlylate + (1|mouse) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noprepost);
% base vs. prepost
compare(base, prepost,'nsim',1000)
% full vs. noprepost 
compare(noprepost,full,'nsim',1000)

% ~ sleep ~
sleep = fitlme(nonan_Tbl,'mean ~ sleep + (1|mouse) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(sleep)
nosleep = fitlme(nonan_Tbl,'mean ~ prepost + earlylate + (1|mouse) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nosleep);
% base vs. prepost
compare(base, sleep,'nsim',1000)
% full vs. noprepost 
compare(nosleep,full,'nsim',1000)

% ~ earlylate ~
earlylate = fitlme(nonan_Tbl,'mean ~ earlylate + (1|mouse) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(earlylate)
noearlylate = fitlme(nonan_Tbl,'mean ~ prepost + sleep + (1|mouse) + (1|sess)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noearlylate);
% base vs. prepost
compare(base, earlylate,'nsim',1000)
% full vs. noprepost 
compare(noearlylate,full,'nsim',1000)
