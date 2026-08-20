%%
clc 
clear
rng(10)

%%

cd ('D:\SWR_DA_MegaMatrix_4s_GFP')

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
    allTables2.PostSWRPeak(i) = max(signal(6400:8001)- mean(signal(1:4800))); % NEED TO CHANGE THIS
    allTables2.PostSWRmean(i) = mean(signal(6400:8001)- mean(signal(1:4800)));
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
    allTables2.PreSWRPeak(i) = max(signal(4800:6401)- mean(signal(1:4800)));
    allTables2.PreSWRmean(i) = mean(signal(4800:6401)- mean(signal(1:4800)));
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
uni_earlylate = unique(allTables2.EarlyLate);  


% Initialize sess_avg_tbl with the same column names but no data
mouseID = 0;
sess = 0;
PrePost = 0;
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
                for i_earlylate = 1:length(uni_earlylate)
                % Find indices where mouseID, session,PrePost, and beforeafter match
                    try
                        list = find(allTables2.sess == uni_sess(i_sess) & ...
                        allTables2.mouseID == uni_mouse(i_mouse) & ...
                        allTables2.PrePost == uni_prepost(i_prepost) &... 
                        allTables2.EarlyLate == uni_earlylate(i_earlylate));
                        if ~isempty(list) % if that session exists. 
                            for i_list = 1:1:length(list) % iterate through all structures 
                                if i_prepost == 1
                                    signal_save = allTables2.TwosPreProc{list(i_list),1}.signal - mean(allTables2.TwosPreProc{list(i_list),1}.signal(1:4800)); %allTables2.PreSWRmean(i_list); % signals are saved in different columns but I'm saving everything 
                                    temp_signal = [temp_signal; signal_save'];
                                end
                                if i_prepost == 2
                                    signal_save = allTables2.TwosPostProc{list(i_list),1}.signal- - mean(allTables2.TwosPostProc{list(i_list),1}.signal(1:4800)); %allTables2.PreSWRmean(i_list);
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
                        new_row = {uni_mouse(i_mouse), uni_sess(i_sess), uni_prepost(i_prepost),uni_earlylate(i_earlylate), avg_base_peak, avg_base_mean, avg_time, avg_signal}; 
                        sess_avg_tbl = [sess_avg_tbl; new_row];
                        temp_signal = [];
                    end
                end
        end
    end
end

sess_avg_tbl = cell2table(sess_avg_tbl, 'variablenames',{'mouse','sess','prepost','earlylate','peak','mean','time','signal'});
% LOTS OF NO COMBOS
% might have to do a basic one with separating out each variable at a time.
sess_avg_tbl2 = sess_avg_tbl;
%% LMMS - make things categorical
sess_avg_tbl2.earlylate = categorical(sess_avg_tbl2.earlylate);
sess_avg_tbl2.prepost = categorical(sess_avg_tbl2.prepost);
sess_avg_tbl2.mouse = categorical(sess_avg_tbl2.mouse);
sess_avg_tbl2.sess = categorical(sess_avg_tbl2.sess);

%% Remove nans: 
include=~isnan(sess_avg_tbl2.peak(:,1));
nonan_Tbl = sess_avg_tbl2(include,:);

%% Need to remove 0s for sleep 


% need to remove 0 
nonan_Tbl(nonan_Tbl.earlylate == '0', :) = [];
nonan_Tbl.earlylate = removecats(nonan_Tbl.earlylate);
nonan_Tbl.earlylate = reordercats(nonan_Tbl.earlylate, {'1','2'});
%% Plotting 
% Collapse to ONE row per mouse x session x prepost
% sleep and earlylate are averaged out
% prepost = 1 and 2 are kept separate
%% Remove the two problematic rows

%% Remove the two problematic rows

nonan_clean = nonan_Tbl;

%% Collapse to ONE row per mouse x session x prepost

mouseIDs = unique(nonan_clean.mouse);
prepostIDs = categories(nonan_clean.prepost);

nSamples = length(nonan_clean.signal{1});

collapsed_signal = {};
collapsed_mouse = [];
collapsed_sess = [];
collapsed_prepost = categorical();

row = 1;

for m = 1:length(mouseIDs)

    theseSessions = unique(nonan_clean.sess( ...
        nonan_clean.mouse == mouseIDs(m)));

    for s = 1:length(theseSessions)

        for p = 1:length(prepostIDs)

            idx = nonan_clean.mouse == mouseIDs(m) & ...
                  nonan_clean.sess == theseSessions(s) & ...
                  nonan_clean.prepost == prepostIDs{p};

            if ~any(idx)
                continue
            end

            signals = nonan_clean.signal(idx);

            % Make all signals row vectors
            signals = cellfun(@(x) x(:)', signals, ...
                              'UniformOutput', false);

            % Combine signals
            signals = vertcat(signals{:});

            % Average sleep/earlylate
            avg_signal = mean(signals, 1, 'omitnan');

            % Store
            collapsed_mouse(row,1) = mouseIDs(m);
            collapsed_sess(row,1) = theseSessions(s);
            collapsed_prepost(row,1) = categorical(prepostIDs(p));
            collapsed_signal{row,1} = avg_signal;

            row = row + 1;

        end
    end
end

collapsed_Tbl = table( ...
    collapsed_mouse, ...
    collapsed_sess, ...
    collapsed_prepost, ...
    collapsed_signal, ...
    'VariableNames', {'mouse','sess','prepost','signal'});

%%
%% GRAND PETH FROM collapsed_Tbl
%
% collapsed_Tbl contains:
%   mouse
%   sess
%   prepost
%   signal
%
% Pre and post are averaged together.
% Each session gets equal weight.
% Each mouse gets equal weight.
% SEM is calculated across mice.

% ---------------------------------------------------------
% SETTINGS
% ---------------------------------------------------------

nSamples = length(collapsed_Tbl.signal{1});

mouseIDs = unique(collapsed_Tbl.mouse);
nMice = length(mouseIDs);

% Store one trace per mouse
mouse_traces = nan(nMice, nSamples);
% ---------------------------------------------------------
% AVERAGE PRE + POST WITHIN EACH SESSION,
% THEN AVERAGE SESSIONS WITHIN EACH MOUSE
% ---------------------------------------------------------

for m = 1:nMice

    thisMouse = mouseIDs(m);

    % Get this mouse's data
    mouse_idx = collapsed_Tbl.mouse == thisMouse;
    T_mouse = collapsed_Tbl(mouse_idx,:);

    % Find sessions for this mouse
    sessionIDs = unique(T_mouse.sess);
    nSessions = length(sessionIDs);

    % Store one trace per session
    session_traces = nan(nSessions, nSamples);

    for s = 1:nSessions

        thisSession = sessionIDs(s);

        % Get Pre + Post for this session
        sess_idx = T_mouse.sess == thisSession;
        T_session = T_mouse(sess_idx,:);

        % Get signals
        signals = T_session.signal;

        % Make sure all are row vectors
        signals = cellfun(@(x) x(:)', signals, ...
                          'UniformOutput', false);

        % Combine Pre + Post
        signals = vertcat(signals{:});

        % Average Pre + Post
        session_traces(s,:) = mean(signals, 1, 'omitnan');

    end

    % Average all sessions for this mouse
    mouse_traces(m,:) = mean(session_traces, 1, 'omitnan');

end

% ---------------------------------------------------------
% GRAND AVERAGE ACROSS MICE
% ---------------------------------------------------------

grand_trace = mean(mouse_traces, 1, 'omitnan');

% ---------------------------------------------------------
% SEM ACROSS MICE
% ---------------------------------------------------------

sem_trace = std(mouse_traces, 0, 1, 'omitnan') ./ ...
            sqrt(sum(~isnan(mouse_traces), 1));

% ---------------------------------------------------------
% TIME VECTOR
% ---------------------------------------------------------

% 8001 points from -4 to +4 seconds
time = linspace(-4, 4, nSamples);

% ---------------------------------------------------------
% PLOT GRAND PETH
% ---------------------------------------------------------

figure;
hold on;

% SEM shading
fill([time fliplr(time)], ...
     [grand_trace + sem_trace, ...
      fliplr(grand_trace - sem_trace)], ...
     [0.3 0.7 0.3], ...
     'FaceAlpha', 0.25, ...
     'EdgeColor', 'none');

% Grand average trace
plot(time, grand_trace, ...
     'Color', [0 0.5 0], ...
     'LineWidth', 3);

% Event at time 0
xline(0, 'k-', 'LineWidth', 1.5);

% Formatting

xlim([-4 4]);
xticks([-4 -2 0 2 4]);
yticks([-0.04 0 0.04 0.08 0.12])
ylim([-0.04 0.16]);

xlabel('Time from event (s)');
ylabel('Signal');

title('Grand Average PETH');

box off;
set(gca, 'FontSize', 12);

hold off;

%% Magnitude & Timing + SEM 
% Time vector

nSamples = size(mouse_traces,2);
time = linspace(-4,4,nSamples);

nMice = size(mouse_traces,1);

% Preallocate

peak_magnitude = nan(nMice,1);
peak_time = nan(nMice,1);

% Find peak for each mouse

for m = 1:nMice

    this_trace = mouse_traces(m,:);

    % Find maximum value
    [peak_magnitude(m), peak_idx] = max(this_trace);

    % Corresponding time
    peak_time(m) = time(peak_idx);

end

% Grand mean peak magnitude

mean_peak_magnitude = mean(peak_magnitude,'omitnan');

% SEM peak magnitude

sem_peak_magnitude = std(peak_magnitude,0,'omitnan') / ...
                     sqrt(sum(~isnan(peak_magnitude)));
% Grand mean peak timing

mean_peak_time = mean(peak_time,'omitnan');

% SEM peak timing

sem_peak_time = std(peak_time,0,'omitnan') / ...
                sqrt(sum(~isnan(peak_time)));

% Display results

fprintf('\n============================\n');
fprintf('GRAND PETH PEAK RESULTS\n');
fprintf('============================\n');

fprintf('Peak magnitude: %.4f ± %.4f SEM\n', ...
        mean_peak_magnitude, sem_peak_magnitude);

fprintf('Peak timing: %.4f ± %.4f s SEM\n', ...
        mean_peak_time, sem_peak_time);

fprintf('Number of mice: %d\n', nMice);

%% Just basic analysis 
base_peak = fitlme(nonan_Tbl,'peak ~ 1 + (1|mouse)'); 
disp(base_peak)

base_mean = fitlme(nonan_Tbl,'mean ~ 1 + (1|mouse)'); 
disp(base_mean)

%% PEAK analysis 
% intercept
base = fitlme(nonan_Tbl,'peak ~ 1 + (1|mouse)'); 
disp(base)

% Full model
full = fitlme(nonan_Tbl,'peak ~ prepost + sleep + earlylate + (1|mouse)'); % removed 1|session to be consistent with event based
disp(full)

% ~ PREPOST ~
prepost = fitlme(nonan_Tbl,'peak ~ prepost + (1|mouse)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(prepost)
noprepost = fitlme(nonan_Tbl,'peak ~ sleep + earlylate + (1|mouse)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noprepost);
% base vs. prepost
compare(base, prepost,'nsim',1000)
% full vs. noprepost 
compare(noprepost,full,'nsim',1000)

% ~ sleep ~
sleep = fitlme(nonan_Tbl,'peak ~ sleep + (1|mouse)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(sleep)
nosleep = fitlme(nonan_Tbl,'peak ~ prepost + earlylate + (1|mouse)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nosleep);
% base vs. prepost
compare(base, sleep,'nsim',1000)
% full vs. noprepost 
compare(nosleep,full,'nsim',1000)

% ~ earlylate ~
earlylate = fitlme(nonan_Tbl,'peak ~ earlylate + (1|mouse)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(earlylate)
noearlylate = fitlme(nonan_Tbl,'peak ~ prepost + sleep + (1|mouse)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noearlylate);
% base vs. prepost
compare(base, earlylate,'nsim',1000)
% full vs. noprepost 
compare(noearlylate,full,'nsim',1000)


%% MEAN analysis 
% intercept
base = fitlme(nonan_Tbl,'mean ~ 1 + (1|mouse)'); 
disp(base)

% Full model
full = fitlme(nonan_Tbl,'mean ~ prepost + sleep + earlylate + (1|mouse)'); 
disp(full)

% ~ PREPOST ~
prepost = fitlme(nonan_Tbl,'mean ~ prepost + (1|mouse)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(prepost)
noprepost = fitlme(nonan_Tbl,'mean ~ sleep + earlylate + (1|mouse)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noprepost);
% base vs. prepost
compare(base, prepost,'nsim',1000)
% full vs. noprepost 
compare(noprepost,full,'nsim',1000)

% ~ sleep ~
sleep = fitlme(nonan_Tbl,'mean ~ sleep + (1|mouse)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(sleep)
nosleep = fitlme(nonan_Tbl,'mean ~ prepost + earlylate + (1|mouse)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nosleep);
% base vs. prepost
compare(base, sleep,'nsim',1000)
% full vs. noprepost 
compare(nosleep,full,'nsim',1000)

% ~ earlylate ~
earlylate = fitlme(nonan_Tbl,'mean ~ earlylate + (1|mouse)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(earlylate)
noearlylate = fitlme(nonan_Tbl,'mean ~ prepost + sleep + (1|mouse)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noearlylate);
% base vs. prepost
compare(base, earlylate,'nsim',1000)
% full vs. noprepost 
compare(noearlylate,full,'nsim',1000)
