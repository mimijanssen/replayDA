%% Plotting lone SWRs. 
%cd F:\SWR_DA_MegaMatrix
%allTables = load('MegaMatrixALLDATA.mat');
cd F:\
allTables = load('MegaMatrix1s_sleep_basicswrs.mat');

% load color
cd C:\Users\mimia\Documents\GitHub\replayDA\Analysis
load('colors.mat')

%%
% allTables = []; % Initialize an empty array for concatenation
% 
% Files=dir('*.*');
% for k=3:length(Files)
%    FileNames=Files(k).name;
%    loadedData = load(FileNames); % Load the .mat file
% 
%     % Assuming your table is saved as 'matrix_sess' in each file
%     if isfield(loadedData, 'matrix_sess')
%         sessionTable = loadedData.matrix_sess;
% 
%         % Concatenate tables vertically
%         if isempty(allTables)
%             allTables = sessionTable; % Initialize with the first table
%         else
%             allTables = [allTables; sessionTable]; % Append subsequent tables
%         end
%     else
%         fprintf('Warning: %s does not contain a table named matrix_sess.\n', FileName);
%     end
% end

%%
ProcPeakTbl = stack(allTables,{'OnesBeforePeak','OnesAfterPeak'},'NewDataVariableName','Peak','IndexVariableName','BeforeAfter');

%% Find the max value within 1 second Before & After
% within 1 second of a SWR save: 
% 1) peak fiber value
% 2) time of peak fiber value
% 3) time of SWR

Peak_one_sec = [];
Time_one_sec = []; 


% If GFP Mice use this: 
%x1 = 1600;
%x2 = 3200; 

%x3 = 3201; 
%x4 = 4801; 

% IF regular MICE use this: 
%x1 = 1000; 
%x2 = 2000; 
%x3 = 2001; 
%x4 = 3001; 

% x1-x2 : 1-1000; 1001- 2001 mostly 
%1000:2000 % for GFP : 1600:3200 % finds the max signal one second before an swr 

for i = 1:1:height(ProcPeakTbl)
    if ProcPeakTbl.PrePost(i) == 1 % pre session 
        if ProcPeakTbl.BeforeAfter(i) == 'OnesBeforePeak' % before swr
            % find one second time range
            x2 = (length(ProcPeakTbl.TwosPreRaw{i,1}.signal) - 1)/2;
            x1 = 1;
            % parse out period of signal you are interested in 
            subset_signal = ProcPeakTbl.TwosPreProc{i,1}.signal(x1:x2);
            % find max signal during that time
            [max_before, I_before] = max(subset_signal); 
            % initialize time vector based on swr time! (x2)
            subset_time = ProcPeakTbl.TwosPreProc{i,1}.tvec(x1:x2)- ProcPeakTbl.TwosPreProc{i,1}.tvec(x2+1); % initialized
            % update matrix. 
            Peak_one_sec = [Peak_one_sec; max_before];
            Time_one_sec = [Time_one_sec; subset_time(I_before)];
        elseif ProcPeakTbl.BeforeAfter(i) == 'OnesAfterPeak' % after swr
            x3 = ((length(ProcPeakTbl.TwosPreRaw{i,1}.signal) - 1)/2);
            x4 = length(ProcPeakTbl.TwosPreRaw{i,1}.signal);
            subset_signal = ProcPeakTbl.TwosPreProc{i,1}.signal(x3:x4);
            subset_time = ProcPeakTbl.TwosPreProc{i,1}.tvec(x3:x4) - ProcPeakTbl.TwosPreProc{i,1}.tvec(x2+1); 
            [max_after, I_after] = max(subset_signal); %2001:3001 % for GFP 3201:4801 % finds the max signal one second before an swr 
            Peak_one_sec = [Peak_one_sec; max_after];
            Time_one_sec = [Time_one_sec; subset_time(I_after)];
        end
    elseif ProcPeakTbl.PrePost(i) == 2 % post session 
        if ProcPeakTbl.BeforeAfter(i) == 'OnesBeforePeak' % before swr
            x2 = (length(ProcPeakTbl.TwosPostRaw{i,1}.signal) - 1)/2;
            x1 = 1;
            subset_signal = ProcPeakTbl.TwosPostProc{i,1}.signal(x1:x2); 
            subset_time = ProcPeakTbl.TwosPostProc{i,1}.tvec(x1:x2) - ProcPeakTbl.TwosPostProc{i,1}.tvec(x2+1);
            [max_before, I_before] = max(subset_signal);
            Peak_one_sec = [Peak_one_sec; max_before];
            Time_one_sec = [Time_one_sec; subset_time(I_before)];
        elseif ProcPeakTbl.BeforeAfter(i) == 'OnesAfterPeak' % after swr
            x3 = ((length(ProcPeakTbl.TwosPostRaw{i,1}.signal) - 1)/2);
            x4 = length(ProcPeakTbl.TwosPostRaw{i,1}.signal);
            subset_signal = ProcPeakTbl.TwosPostProc{i,1}.signal(x3:x4);
            subset_time =  ProcPeakTbl.TwosPostProc{i,1}.tvec(x3:x4) - ProcPeakTbl.TwosPostProc{i,1}.tvec(x2+1); 
            [max_after, I_after] = max(subset_signal); % finds the max signal one second before an swr 
            Peak_one_sec = [Peak_one_sec; max_after];
            Time_one_sec = [Time_one_sec; subset_time(I_after)];
        end
    end
end

%% append table 

%Post_swr_peak_t=array2table((Time_swr)-(Time_one_sec));
%ProcPeakTbl = [ProcPeakTbl Post_swr_peak_t];

Peak_one_sec= array2table(Peak_one_sec);
Time_one_sec= array2table(Time_one_sec);
%Time_swr= array2table(Time_swr);
ProcPeakTbl = [ProcPeakTbl Peak_one_sec];
ProcPeakTbl = [ProcPeakTbl Time_one_sec];
%ProcPeakTbl = [ProcPeakTbl Time_swr];

Tbl = ProcPeakTbl;

% %% OK saving a new table because I'm doing some wrangling 
% Tbl = ProcPeakTbl; 
% 
% % delete before peaks 
% After_Tbl = ProcPeakTbl; 
% toDelete = After_Tbl.BeforeAfter == 'TwosBeforePeak';
% After_Tbl(toDelete,:) = [];
% 
% %% Find all lone SWRs and save in a new table 
% % 1) iterate through each row 
% % 2) check if the i and i+ are the same session (TO DO...) I might skip
% % this because the times will be so different. 
% % 3) if the same session, then if the swr are < 1 second away, delete them (0.001 s time step) 
% toDelete = [];
% for i_swr = 2:1:height(After_Tbl) % skips the first swr
%     if After_Tbl.Time_swr(i_swr) - After_Tbl.Time_swr(i_swr-1) < 2 % if the swr before is less than 1 second
%         toDelete = [toDelete; i_swr-1];
%     end
% end
% 
% % delete those
% After_Tbl(toDelete,:) = []; % save a list of indicies 
% 
% %% all SWRs
% % delete before peaks 
% After_Tbl2 = ProcPeakTbl; 
% toDelete = After_Tbl2.BeforeAfter == 'TwosBeforePeak';
% After_Tbl2(toDelete,:) = [];

%% AVERAGE OVER SESSION FOR PREPOST AND BEFORE AFTER. 

uni_mouse = unique(Tbl.mouseID);
uni_sess = unique(Tbl.sess);
uni_swr = unique(Tbl.BeforeAfter); 
uni_prepost = unique(Tbl.PrePost);  % Should be [1,2] for Pre/Post

% Initialize sess_avg_tbl with the same column names but no data
mouseID = 0;
sess = 0;
PrePost = 0;
Peak_one_sec = 0;
BeforeAfter = 0;
avg_signal = [];
temp_signal = [];

%sess_avg_tbl = table(mouseID, sess, PrePost, Peak_one_sec2); 
sess_avg_tbl = [];
% Loop over unique mice, sessions, and Pre/Post categories
for i_mouse = 1:length(uni_mouse)
    for i_sess = 1:length(uni_sess)
        list_allsess = find(Tbl.sess == uni_sess(i_sess) & ...
        Tbl.mouseID == uni_mouse(i_mouse)); % create a list of all swrs within a session for a mouse.
        for i_prepost = 1:length(uni_prepost) % loops through both pre and post 
            for i_beforeafter = 1:length(uni_swr) % loops through before first and then after...  
                % Find indices where mouseID, session,PrePost, and beforeafter match
                try
                list = find(Tbl.sess == uni_sess(i_sess) & ...
                       Tbl.mouseID == uni_mouse(i_mouse) & ...
                            Tbl.PrePost == uni_prepost(i_prepost) &... 
                            Tbl.BeforeAfter == uni_swr(i_beforeafter));

                if ~isempty(list) % if that session exists. 
                    for i_list = 1:1:length(list) % iterate through all structures 
                        if i_prepost == 1
                            signal_save = Tbl.TwosPreProc{list(i_list),1}.signal; % signals are saved in different columns but I'm saving everything 
                            temp_signal = [temp_signal; signal_save'];
                        end
                        if i_prepost == 2
                            signal_save = Tbl.TwosPostProc{list(i_list),1}.signal;
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
            avg_peak = mean(Tbl.Peak_one_sec(list));
            avg_time = mean(Tbl.Time_one_sec(list));
            new_row = {uni_mouse(i_mouse), uni_sess(i_sess), uni_swr(i_beforeafter), uni_prepost(i_prepost), avg_peak, avg_time, avg_signal}; 
            sess_avg_tbl = [sess_avg_tbl; new_row];
            temp_signal = [];
        end
        end
    end
    end
end

sess_avg_tbl2 = cell2table(sess_avg_tbl, 'variablenames',{'mouse','sess','beforeafter','prepost','peak','time','signal'});

%% Make things categorical 
sess_avg_tbl2.prepost = categorical(sess_avg_tbl2.prepost);
sess_avg_tbl2.mouse = categorical(sess_avg_tbl2.mouse);
sess_avg_tbl2.sess = categorical(sess_avg_tbl2.sess);
sess_avg_tbl2.beforeafter = categorical(sess_avg_tbl2.beforeafter);

%% Session Based LMM 
% Figuring out best base model! Do we want random intercept and slopes for
% mouse? or nested session? Let's see.

% MODEL 1: random intercept for mouse 
lme_base = fitlme(sess_avg_tbl2,'peak ~ 1 + (1|mouse)'); 
disp(lme_base)
% AIC: -158.9
% BIC: -149.2

% Model 2: intercept and session slope vary by mouse
lme_base2 = fitlme(sess_avg_tbl2,'peak ~ 1 + (sess|mouse)'); 
disp(lme_base2)
% AIC: -154.91
% BIC: -31.924 (maybe over fit)
% technically better model when using compare, but it uses a lot of
% parameters so the BIC is not the best. 

% Model 3: random intercept for session nested within mouse
lme_base3 = fitlme(sess_avg_tbl2,'peak ~ 1 + (1|sess:mouse)');
disp(lme_base3)
% AIC: -178.75
% BIC: -169.04
% ~~~~~ Lowest BIC! ~~~~~~

% Model 4: random intercept for mouse and sess 
lme_base4 = fitlme(sess_avg_tbl2,'peak ~ 1 + (1|mouse) + (1|sess)'); 
disp(lme_base4)
% AIC: -159.38
% BIC: -146.44

compare(lme_base, lme_base3, 'nsim',1000)
%lmebase3 is significanlty better base model! 

compare(lme_base3, lme_base4, 'nsim',1000) 
% lme_base4 is marginally better... by p = 0.025 --> p < 0.05. 

% Will use model 4 as the base model, with random intercepts for mouse and
% session :). 

%% Base model with added either prepost or before after 
% using lme_base2 
lme_beforeafter = fitlme(sess_avg_tbl2,'peak ~ beforeafter + (1|mouse) + (1|sess)'); 
disp(lme_beforeafter)

[results,siminfo] = compare(lme_base4, lme_beforeafter, 'nsim',1000) 


%%  New saved variables without before peaks...
sess_avg_tbl3 = sess_avg_tbl2; 
toDelete = sess_avg_tbl3.beforeafter == 'OnesBeforePeak';
sess_avg_tbl3(toDelete,:) = [];

%prepost base
lme_prepost_base = fitlme(sess_avg_tbl3,'peak ~ 1 + (1|mouse) + (1|sess)'); 
disp(lme_prepost_base)

lme_prepost = fitlme(sess_avg_tbl3,'peak ~ prepost + (1|mouse) + (1|sess)'); 
disp(lme_prepost)

[results,siminfo] = compare(lme_prepost_base, lme_prepost, 'nsim',1000) 


%% Early vs. Late

list = zeros(height(sess_avg_tbl3),1); 
% early sessions 1-4 = 1
I_early = find(double(sess_avg_tbl3.sess) < 5);  % Find indices where 'condition_row' is positive
list(I_early) = 1; 

% late sessions 5-6 = 2
I_late = find(double(sess_avg_tbl3.sess) > 4);  % Find indices where 'condition_row' is positive
list(I_late) = 2; 

% append list to table 
sess_avg_tbl3.("EarlyLate") = list;

%% Model for Early vs. Late

lme_earlylate = fitlme(sess_avg_tbl3,'peak ~ EarlyLate + (1|mouse) + (1|sess)'); 
disp(lme_earlylate)

[results,siminfo] = compare(lme_prepost_base, lme_earlylate, 'nsim',1000) 


%% Can't fit the right base model/idk how to fit a base model without a predictor. 
% 
% lme_beforeafter1 = fitlme(sess_avg_tbl2,'peak ~ beforeafter + (beforeafter|mouse) + (beforeafter|mouse:sess)'); % I used peak value here instead of the circ shift swr strength...
% disp(lme_beforeafter1)
% % this should be random slopes and intercepts for mouse and session 
% % AIC : -190.35
% 
% lme_beforeafter2 = fitlme(sess_avg_tbl2,'peak ~ beforeafter + (1|mouse) + (1|mouse:sess)'); % I used peak value here instead of the circ shift swr strength...
% disp(lme_beforeafter2)
% %random intercepts
% % AIC = -189.19 
% 
% compare(lme_beforeafter1, lme_beforeafter2, 'nsim',1000) % lme_base2 is a better model (p = 0.02697)
% % comparable 
% 
% lme_beforeafter3 = fitlme(sess_avg_tbl2,'peak ~ beforeafter + (1|mouse)');
% disp(lme_beforeafter3)
% % AIC = -153.02
% 
% compare(lme_beforeafter1, lme_beforeafter3, 'nsim',1000) % lme_base2 is a better model (p = 0.02697)
% % no difference... 
% 
% lme_beforeafter4 = fitlme(sess_avg_tbl2,'peak ~ beforeafter + (1|sess)');
% disp(lme_beforeafter4)
% %AIC = -105
% 
% compare(lme_beforeafter1, lme_beforeafter4, 'nsim',1000) % lme_base2 is a better model (p = 0.02697)
% % eh no difference
% 
% lme_beforeafter5 = fitlme(sess_avg_tbl2,'peak ~ beforeafter + (beforeafter|mouse)');
% disp(lme_beforeafter5)
% % AIC -152
% 
% compare(lme_beforeafter1, lme_beforeafter5, 'nsim',1000) % lme_base2 is a better model (p = 0.02697)
% % no difference 
% 
% lme_beforeafter6 = fitlme(sess_avg_tbl2,'peak ~ beforeafter + (beforeafter|mouse)');
% disp(lme_beforeafter6)
% % AIC -152
% 
% compare(lme_beforeafter1, lme_beforeafter6, 'nsim',1000) % lme_base2 is a better model (p = 0.02697)
% % no difference

%%
% with mouse ? 
% lme_mouse = fitlme(sess_avg_tbl2,'peak ~ mouse + (beforeafter|mouse) + (beforeafter|mouse:sess)'); % I used peak value here instead of the circ shift swr strength...
% disp(lme_mouse)
% % with before and after information 
% 
% lme_beforeafter1 = fitlme(sess_avg_tbl2,'peak ~ mouse + beforeafter + (beforeafter|mouse) + (beforeafter|mouse:sess)'); % I used peak value here instead of the circ shift swr strength...
% disp(lme_beforeafter1)
% 
% lme_beforeafter = fitlme(sess_avg_tbl2,'peak ~ beforeafter + (beforeafter|mouse) + (beforeafter|mouse:sess)');
% disp(lme_beforeafter)
% 
% compare(lme_mouse, lme_beforeafter1, 'nsim',1000) 
% % mouse vs. beforeafter1
% % duh beforeafter information improves the model. 
% 
% compare(lme_beforeafter1, lme_beforeafter, 'nsim',1000) % lme_base2 is a better model (p = 0.02697)
% %beforeafter1 (with mouse as dependent) vs. before after
% % no mouse dependent variable is better. 
% 
% lme_beforeafter_timeint = fitlme(sess_avg_tbl2,'peak ~ beforeafter + beforeafter*time + (beforeafter|mouse) + (beforeafter|mouse:sess)');
% disp(lme_beforeafter_timeint)
% % is a significant interaction 
% 
% compare(lme_beforeafter, lme_beforeafter_timeint, 'nsim',1000) % lme_base2 is a better model (p = 0.02697)
% TIMING TELLS YOU SOMETHING!!! 

% % PREPOST 
% %lme_beforeafter = fitlme(sess_avg_tbl2,'peak ~ beforeafter + mouse + (1|sess:mouse)'); % I used peak value here instead of the circ shift swr strength...
% %disp(lme_beforeafter)
% % AIC: -183.58 
% % beforeafter is significant ! VERY 
% 
% %compare(lme_base2, lme_beforeafter, 'nsim',1000) % lme_base2 is a better model (p = 0.02697)
% % LMEBefore and after is better :) 
% 
% % make prepost categorical 
% lme_prepost_mouse = fitlme(sess_avg_tbl2,'peak ~ prepost +  mouse + (beforeafter|mouse) + (beforeafter|mouse:sess)');
% disp(lme_prepost)
% 
% compare(lme_mouse, lme_prepost_mouse, 'nsim',1000) 
% % Just barely better (prepost)
% 
% lme_prepost = fitlme(sess_avg_tbl2,'peak ~ prepost +  (beforeafter|mouse) + (beforeafter|mouse:sess)');
% disp(lme_prepost)
% compare(lme_prepost, lme_prepost_mouse, 'nsim',1000) 
% 
% fullmodel = fitlme(sess_avg_tbl2,'peak~ beforeafter + prepost + (beforeafter|mouse) + (beforeafter|mouse:sess)');
% disp(fullmodel)
% 
% % might be also interesting to take out before data and see if time says
% % anything? 
% compare(fullmodel, lme_prepost, 'nsim',1000) % lme_base2 is a better model (p = 0.02697)
% compare(fullmodel, lme_beforeafter, 'nsim',1000) % lme_base2 is a better model (p = 0.02697)


%% All SWRs average for sessions
% 
% uni_mouse = unique(After_Tbl.mouseID);
% uni_sess = unique(After_Tbl.sess);
% uni_prepost = unique(After_Tbl.PrePost);  % Should be [1,2] for Pre/Post
% 
% % Initialize sess_avg_tbl with the same column names but no data
% mouseID = 0;
% sess = 0;
% PrePost = 0;
% Peak_one_sec = 0;
% avg_signal = [];
% temp_signal = [];
% 
% %sess_avg_tbl = table(mouseID, sess, PrePost, Peak_one_sec2); 
% sess_avg_tbl_all = [];
% % Loop over unique mice, sessions, and Pre/Post categories
% for i_mouse = 1:length(uni_mouse)
%     for i_sess = 1:length(uni_sess)
%         list_allsess = find(After_Tbl2.sess == uni_sess(i_sess) & ...
%         After_Tbl2.mouseID == uni_mouse(i_mouse)); % create a list of all swrs within a session for a mouse. 
%         for i_prepost = 1:length(uni_prepost) % loops through both pre and post  
%                 % Find indices where mouseID, session, and PrePost match
%                 try
%                 list = find(After_Tbl2.sess == uni_sess(i_sess) & ...
%                        After_Tbl2.mouseID == uni_mouse(i_mouse) & ...
%                             After_Tbl2.PrePost == uni_prepost(i_prepost));
% 
%                 if ~isempty(list) % if that session exists. 
%                     for i_list = 1:1:length(list) % iterate through all structures 
%                         if i_prepost == 1
%                             signal_save = After_Tbl2.FoursPreProc{list(i_list),1}.signal;
%                             temp_signal = [temp_signal; signal_save'];
%                         end
%                         if i_prepost == 2
%                             signal_save = After_Tbl2.FoursPostProc{list(i_list),1}.signal;
%                             temp_signal = [temp_signal; signal_save'];
%                         end
%                     end
%                 end
%                     % Create a new row to append
%                     %new_row = {uni_mouse(i_mouse), uni_sess(i_sess), uni_prepost(i_prepost), avg_peak}; % for main analysis, extend this to also include before ... save avg signal here 
%                     % Append the row to sess_avg_tbl
%                     %sess_avg_tbl = [sess_avg_tbl; new_row];
%                 catch
%                     disp('no combo')
%                 end
%                 avg_signal = mean(temp_signal);
%         end
% 
%         if ~isempty(list_allsess) % if that session exists 
%             avg_peak = mean(After_Tbl2.Peak_one_sec(list));
%             new_row = {uni_mouse(i_mouse), uni_sess(i_sess), avg_peak, avg_signal}; 
%             sess_avg_tbl_all = [sess_avg_tbl_all; new_row];
%             temp_signal = [];
%         end
%     end
% end
% 
% sess_avg_tbl_all2 = cell2table(sess_avg_tbl_all, 'variablenames',{'mouse','sess','peak','signal'});

%% Average for each mouse and then plot that!!! 
%new table
plottbl = sess_avg_tbl3;
grand_signal = []; 
grand_signal_m5 = [];
plottbl.mouse = double(plottbl.mouse);

std_grand = [];
std_grand_m5 = [];

% Loop through each mouse and collect their mean signals 'variablenames'
for m = 1:length(uni_mouse) % iterate over mice
    if m ~= 5 
        mouseI = ismember(plottbl.mouse, m); % find all those guys index 
        % Take only those indicies 
        temptable = sess_avg_tbl3(mouseI,:);
        signal_matrix = cell2mat(temptable.signal);  % results in an N x 2001 matrix
        mean_signal = mean(signal_matrix, 1);  % 1 x 2001
        grand_signal = [grand_signal; mean_signal]; % Collect pre means
        std_grand = [std_grand; std(signal_matrix,1)];
    else
        mouseI = ismember(plottbl.mouse, m); % find all those guys index 
        temptable = sess_avg_tbl3(mouseI,:);
        signal_matrix = cell2mat(temptable.signal);  % results in an N x 2001 matrix
        mean_signal = mean(signal_matrix, 1);  % 1 x 2001
        grand_signal_m5 = [grand_signal_m5; mean_signal]; % Collect pre means
        std_grand_m5 = [std_grand_m5 std(signal_matrix,1)];
    end
end

%% Decimate the m5 signal to be the same frequency as the rest of the data... 
% might need to go back and replot the other data.... (RPE and Main plot
% becuase the ddata is different size)

% t = 0:.00025:1;  % Time vector
% x = sin(2*pi*30*t) + sin(2*pi*60*t);
% y = decimate(x,4);
% subplot(1,2,1);
% stem(x(1:120)), axis([0 120 -2 2])   % Original signal
% title('Original Signal')
% subplot(1,2,2);
% stem(y(1:30))                        % Decimated signal
% title('Decimated Signal')

%dsf = /1000; 
%FP.data = decimate(FP.data,dsf);
%FP.tvec = downsample(FP.tvec,dsf);
%FP.cfg.hdr{1}.SamplingFrequency = FP.cfg.hdr{1}.SamplingFrequency/dsf;

% Use resample (requires Signal Processing Toolbox)
signal_2001 = resample(grand_signal_m5, 2001, 3001);  % 1x2001
signal_std = resample(std_grand_m5, 2001, 3001);
grand_signal = [grand_signal; signal_2001]; % good
std_signal = [std_grand; signal_std];


%% PLOT THAT WOOHOOO
common_tvec = linspace(-1, 1, 3001); % Modify based on your data
grand_mean_signal = mean(grand_signal);
grand_SEM_signal = std(grand_signal)/sqrt(7); 
figure(1)
shadedErrorBar(common_tvec, grand_mean_signal, grand_SEM_signal, 'lineprops', {'Color', colors(m, :), 'LineWidth', 1.5});
xlabel('time from SWR (s)')
ylabel('Mean [DA] (z-score)')

%% plot everything 
figure(2)
for m = 1:1:length(colors)
    shadedErrorBar(common_tvec, grand_signal(m,:),std_signal(m,:), 'lineprops', {'Color', colors(m,:), 'LineWidth', 1.5});
    hold on
end
xlabel('time from SWR (s)')
ylabel('Mean [DA] (z-score)')
legend()
%%
shadedErrorBar(common_tvec, grand_signal_m5, std(grand_signal_m5))
%% To compare, I want to plot the same way with all SWRS on the same graph... 
figure(2)
xl = xline(0,'--k',{'SWR'});
hold on
shadedErrorBar(common_tvec, grand_mean_signal, grand_SEM_signal, 'lineprops', {'Color', colors(1, :), 'LineWidth', 1.5});
hold on
shadedErrorBar(common_tvec, grand_mean_signal_all, grand_SEM_signal_all, 'lineprops', {'Color', colors(2, :), 'LineWidth', 1.5});
legend('','lone SWRs', 'all SWRs','Location','northwest')
legend 'boxoff'
xlabel('time from SWR (s)')
ylabel('Mean [DA] (z-score)')


%% Try separating pre and post? 
uni_mouse = unique(After_Tbl.mouseID);
uni_sess = unique(After_Tbl.sess);
uni_prepost = unique(After_Tbl.PrePost);  % Should be [1,2] for Pre/Post

% Initialize sess_avg_tbl with the same column names but no data
mouseID = 0;
sess = 0;
PrePost = 0;
Peak_one_sec = 0;
avg_signal = [];
temp_signal = [];

%sess_avg_tbl = table(mouseID, sess, PrePost, Peak_one_sec2); 
sess_avg_tbl_prepost = [];
% Loop over unique mice, sessions, and Pre/Post categories
for i_mouse = 1:length(uni_mouse)
    for i_sess = 1:length(uni_sess)
        for i_prepost = 1:length(uni_prepost) % loops through both pre and post  
                % Find indices where mouseID, session, and PrePost match
                list = find(After_Tbl2.sess == uni_sess(i_sess) & ...
                       After_Tbl2.mouseID == uni_mouse(i_mouse) & ...
                            After_Tbl2.PrePost == uni_prepost(i_prepost));

                if ~isempty(list) % if that session exists. 
                    for i_list = 1:1:length(list) % iterate through all structures 
                        if i_prepost == 1
                            signal_save = After_Tbl2.FoursPreProc{list(i_list),1}.signal;
                            temp_signal = [temp_signal; signal_save'];
                        end
                        if i_prepost == 2
                            signal_save = After_Tbl2.FoursPostProc{list(i_list),1}.signal;
                            temp_signal = [temp_signal; signal_save'];
                        end
                    end
                   
                    avg_signal = mean(temp_signal);
                    temp_signal = [];

                    avg_peak = mean(After_Tbl2.Peak(list));
                    new_row = {uni_mouse(i_mouse), uni_sess(i_sess), uni_prepost(i_prepost), avg_peak, avg_signal}; 
                    sess_avg_tbl_prepost = [sess_avg_tbl_prepost; new_row]; %there should be 48 here
                end
                clear list
    

        end
    end
end

sess_avg_tbl_prepost2 = cell2table(sess_avg_tbl_prepost, 'variablenames',{'mouse','sess','prepost','peak','signal'});

%% WHY NOT WORKING
uni_mouse = unique(After_Tbl2.mouseID);
uni_prepost = unique(After_Tbl2.PrePost);
sess_avg_tbl_prepost = {};

for i_mouse = 1:length(uni_mouse)
    mouse_sessions = unique(After_Tbl2.sess(After_Tbl2.mouseID == uni_mouse(i_mouse)));
    
    for i_sess = 1:length(mouse_sessions)
        for i_prepost = 1:length(uni_prepost)
            list = find(After_Tbl2.sess == mouse_sessions(i_sess) & ...
                        After_Tbl2.mouseID == uni_mouse(i_mouse) & ...
                        After_Tbl2.PrePost == uni_prepost(i_prepost));

            if ~isempty(list)
                temp_signal = [];
                for i_list = 1:length(list)
                    if i_prepost == 1
                        signal_save = After_Tbl2.FoursPreProc{i_list}.signal;
                    else
                        signal_save = After_Tbl2.FoursPostProc{i_list}.signal;
                    end
                    temp_signal = [temp_signal; signal_save'];
                end

                avg_signal = mean(temp_signal);
                avg_peak = mean(After_Tbl2.Peak(list));
                new_row = {uni_mouse(i_mouse), mouse_sessions(i_sess), uni_prepost(i_prepost), avg_peak, avg_signal}; 
                sess_avg_tbl_prepost = [sess_avg_tbl_prepost; new_row];
            end
        end
    end
end

sess_avg_tbl_prepost2 = cell2table(sess_avg_tbl_prepost, 'variablenames',{'mouse','sess','prepost','peak','signal'});


%% average the pre sessions for each mouse
grand_signal_pre = []; 
% Loop through each mouse and collect their mean signals'variablenames'
for m = 1:length(uni_mouse) % iterate over mice
    mouseI = ismember(sess_avg_tbl2.mouse, m); % find all those guys index 
    grand_signal = [grand_signal; mean(sess_avg_tbl2.signal(mouseI,:))]; % Collect pre means
end


