%%
rng(10)
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

% ALSO LOAD THE OTHER DATA TABLE . 
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
for i = 1:height(allTables2)
    if allTables2.PrePost(i) == 1
        signal = allTables2.TwosPreProc{i}.signal;
    elseif allTables2.PrePost(i) == 2
        signal = allTables2.TwosPostProc{i}.signal;
    else
        continue
    end
    allTables2.PostSWRPeak(i) = max(signal(4001:end)); % NEED TO CHANGE THIS
end

% Before SWR
allTables2.PreSWRPeak = NaN(height(allTables2),1);
for i = 1:height(allTables2)
    if allTables2.PrePost(i) == 1
        signal = allTables2.TwosPreProc{i}.signal;
    elseif allTables2.PrePost(i) == 2
        signal = allTables2.TwosPostProc{i}.signal;
    else
        continue
    end
    allTables2.PreSWRPeak(i) = max(signal(3000:4000));
end

allTables2.Base_peak = allTables2.PostSWRPeak - allTables2.PreSWRPeak;

%% Add the Early Late variable 
list = zeros(height(allTables2),1); 
% early sessions 1-3 = 1
I_early = find(allTables2.sess < 4);  % Find indices where 'condition_row' is positive
list(I_early) = 1; 

% late sessions 6-8 = 2
I_late = find(allTables2.sess > 5);  % Find indices where 'condition_row' is positive
list(I_late) = 2; 

% append list to table 
allTables2.("EarlyLate") = list;

%% SOME PLOTS 
% Setup
Fs        = 1000;  % adjust if needed
n_samples = 8001;
tvec      = linspace(-4, 4, n_samples);


% Define groups
groups = {
    'All Events',     (allTables2.PrePost == 1 | allTables2.PrePost == 2);
    'Pre-Track',      (allTables2.PrePost == 1);
    'Post-Track',     (allTables2.PrePost == 2);
    'NREM',           (allTables2.sleep == 1);
    'Wake',           (allTables2.sleep == 2);
    'Early Sessions', (allTables2.EarlyLate == 1);
    'Late Sessions',  (allTables2.EarlyLate == 2);
};

% Plot all groups
figure('Position', [100 100 1800 900]);
for g = 1:size(groups, 1)
    label = groups{g, 1};
    mask  = groups{g, 2};
    [grand_mean, grand_sd, n_events] = compute_peth(allTables2, mask, n_samples);
    subplot(2, 4, g);
    hold on;

    % Shaded SD
    %fill([tvec, fliplr(tvec)], ...
    %     [grand_mean + grand_sd, fliplr(grand_mean - grand_sd)], ...
    %     [0.6 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    % Mean trace
    plot(tvec, grand_mean, 'b-', 'LineWidth', 2);

    % Reference lines
    xline(0, 'k--', 'LineWidth', 1.2);
    yline(0, 'k:',  'LineWidth', 0.8);

    xlabel('Time relative to SWR (s)');
    ylabel('DA signal (z-score)');
    title(sprintf('%s (n = %d events)', label, n_events));
    xlim([-4 4]);
    ylim([-0.06 0.09])
    box off;
    set(gca, 'FontSize', 12);
end

sgtitle('SWR-triggered DA signal', 'FontSize', 16);
%set(gcf, 'color', 'none');
set(gcf, 'renderer', 'painters');
%fontname("AvenirNext LT Pro Regular");

%% change things to categorical
allTables2.PrePost = categorical(allTables2.PrePost);
allTables2.mouseID = categorical(allTables2.mouseID);
allTables2.EarlyLate = categorical(allTables2.EarlyLate);
allTables2.sleep = categorical(allTables2.sleep);

%%
% need to remove 0 
allTables2(allTables2.sleep == '0', :) = [];
allTables2.sleep = removecats(allTables2.sleep);

allTables2.sleep = reordercats(allTables2.sleep, {'1','2'});
%% Modeling time 
% Full Model 
full = fitlme(allTables2,'Base_peak ~ PrePost + EarlyLate + sleep + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(full)
% BIC: 54479 
% swrID is sig, PrePost is sig but that's it. 

% Base Model 
base = fitlme(allTables2,'Base_peak ~ 1 + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(base)
% BIC 54442

% ~ PREPOST ~
prepost = fitlme(allTables2,'Base_peak ~ PrePost + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(prepost)
noprepost = fitlme(allTables2,'Base_peak ~ EarlyLate + sleep + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noprepost);
% base vs. prepost
compare(base, prepost,'nsim',1000)
% full vs. noprepost 
compare(noprepost,full,'nsim',1000)

% ~ EarlyLate ~
earlylate = fitlme(allTables2,'Base_peak ~ EarlyLate + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(earlylate)
noearlylate = fitlme(allTables2,'Base_peak ~ PrePost + sleep + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noearlylate);
% base vs. prepost
compare(base, earlylate,'nsim',1000)
% full vs. noprepost 
compare(noearlylate,full,'nsim',1000)

% ~ WakeNrem ~
sleep = fitlme(allTables2,'Base_peak ~ sleep + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(sleep)
nosleep = fitlme(allTables2,'Base_peak ~ PrePost + EarlyLate + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nosleep);
% base vs. prepost
compare(base, sleep,'nsim',1000)
% full vs. noprepost 
compare(nosleep,full,'nsim',1000)

% ~ duration ~
dur = fitlme(allTables2,'Base_peak ~ SWRdur + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(dur)
nodur = fitlme(allTables2,'Base_peak ~ PrePost + EarlyLate + sleep + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nodur);
% base vs. prepost
compare(base, dur,'nsim',1000)
% full vs. noprepost 
compare(nodur,full,'nsim',1000)

% ~ swrpower ~
power = fitlme(allTables2,'Base_peak ~ SWRpower + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(power)
nopower = fitlme(allTables2,'Base_peak ~ PrePost + EarlyLate + sleep + SWRdur + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nopower);
% base vs. prepost
compare(base, power,'nsim',1000)
% full vs. noprepost 
compare(nopower,full,'nsim',1000)

% ~ swrID ~
id = fitlme(allTables2,'Base_peak ~ swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(id)
noid = fitlme(allTables2,'Base_peak ~ PrePost + EarlyLate + sleep + SWRdur + SWRpower + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noid);
% base vs. prepost
compare(base, id,'nsim',1000)
% full vs. noprepost 
compare(noid, full,'nsim',1000)

%% Find base_AUC
% AFTER SWR
allTables2.PostSWRAUC = NaN(height(allTables2),1);
for i = 1:height(allTables2)
    if allTables2.PrePost(i) == '1'
        signal = allTables2.TwosPreProc{i}.signal;
    elseif allTables2.PrePost(i) == '2'
        signal = allTables2.TwosPostProc{i}.signal;
    else
        continue
    end
    allTables2.PostSWRAUC(i) = trapz(signal(1002:end));
end

% Before SWR
allTables2.PreSWRAUC = NaN(height(allTables2),1);
for i = 1:height(allTables2)
    if allTables2.PrePost(i) == '1'
        signal = allTables2.TwosPreProc{i}.signal;
    elseif allTables2.PrePost(i) == '2'
        signal = allTables2.TwosPostProc{i}.signal;
    else
        continue
    end
    allTables2.PreSWRAUC(i) = trapz(signal(1:1001));
end

allTables2.Base_AUC = allTables2.PostSWRAUC - allTables2.PreSWRAUC;


%% AUC Modeling time 
% Full Model 
full = fitlme(allTables2,'Base_AUC ~ PrePost + EarlyLate + sleep + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(full)
% BIC: 54479 
% swrID is sig, PrePost is sig but that's it. 

% Base Model 
base = fitlme(allTables2,'Base_AUC ~ 1 + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(base)
% BIC 54442

% ~ PREPOST ~
prepost = fitlme(allTables2,'Base_AUC ~ PrePost + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(prepost)
noprepost = fitlme(allTables2,'Base_AUC ~ EarlyLate + sleep + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noprepost);
% base vs. prepost
compare(base, prepost,'nsim',1000)
% full vs. noprepost 
compare(noprepost,full,'nsim',1000)

% ~ EarlyLate ~
earlylate = fitlme(allTables2,'Base_AUC ~ EarlyLate + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(earlylate)
noearlylate = fitlme(allTables2,'Base_AUC ~ PrePost + sleep + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noearlylate);
% base vs. prepost
compare(base, earlylate,'nsim',1000)
% full vs. noprepost 
compare(noearlylate,full,'nsim',1000)

% ~ WakeNrem ~
sleep = fitlme(allTables2,'Base_AUC ~ sleep + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(sleep)
nosleep = fitlme(allTables2,'Base_AUC ~ PrePost + EarlyLate + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nosleep);
% base vs. prepost
compare(base, sleep,'nsim',1000)
% full vs. noprepost 
compare(nosleep,full,'nsim',1000)

% ~ duration ~
dur = fitlme(allTables2,'Base_AUC ~ SWRdur + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(dur)
nodur = fitlme(allTables2,'Base_AUC ~ PrePost + EarlyLate + sleep + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nodur);
% base vs. prepost
compare(base, dur,'nsim',1000)
% full vs. noprepost 
compare(nodur,full,'nsim',1000)

% ~ swrpower ~
power = fitlme(allTables2,'Base_AUC ~ SWRpower + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(power)
nopower = fitlme(allTables2,'Base_AUC ~ PrePost + EarlyLate + sleep + SWRdur + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nopower);
% base vs. prepost
compare(base, power,'nsim',1000)
% full vs. noprepost 
compare(nopower,full,'nsim',1000)

% ~ swrID ~
id = fitlme(allTables2,'Base_AUC ~ swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(id)
noid = fitlme(allTables2,'Base_AUC ~ PrePost + EarlyLate + sleep + SWRdur + SWRpower + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noid);
% base vs. prepost
compare(base, id,'nsim',1000)
% full vs. noprepost 
compare(noid,full,'nsim',1000)


%% MEAN 
% AFTER SWR
allTables2.PostSWRmean = NaN(height(allTables2),1);
for i = 1:height(allTables2)
    if allTables2.PrePost(i) == '1'
        signal = allTables2.TwosPreProc{i}.signal;
    elseif allTables2.PrePost(i) == '2'
        signal = allTables2.TwosPostProc{i}.signal;
    else
        continue
    end
    allTables2.PostSWRmean(i) = mean(signal(1002:end));
end

% Before SWR
allTables2.PreSWRmean = NaN(height(allTables2),1);
for i = 1:height(allTables2)
    if allTables2.PrePost(i) == '1'
        signal = allTables2.TwosPreProc{i}.signal;
    elseif allTables2.PrePost(i) == '2'
        signal = allTables2.TwosPostProc{i}.signal;
    else
        continue
    end
    allTables2.PreSWRmean(i) = mean(signal(1:1001));
end

allTables2.Base_mean = allTables2.PostSWRmean - allTables2.PreSWRmean;

%%

%%
% Full Model 
full = fitlme(allTables2,'Base_mean ~ PrePost + EarlyLate + sleep + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(full)
% BIC: 54479 
% swrID is sig, PrePost is sig but that's it. 

% Base Model 
base = fitlme(allTables2,'Base_mean ~ 1 + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(base)
% BIC 54442

% ~ PREPOST ~
prepost = fitlme(allTables2,'Base_mean ~ PrePost + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(prepost)
noprepost = fitlme(allTables2,'Base_mean ~ EarlyLate + sleep + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noprepost);
% base vs. prepost
compare(base, prepost,'nsim',1000)
% full vs. noprepost 
compare(noprepost,full,'nsim',1000)

% ~ EarlyLate ~
earlylate = fitlme(allTables2,'Base_mean ~ EarlyLate + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(earlylate)
noearlylate = fitlme(allTables2,'Base_mean ~ PrePost + sleep + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noearlylate);
% base vs. prepost
compare(base, earlylate,'nsim',1000)
% full vs. noprepost 
compare(noearlylate,full,'nsim',1000)



% ~ WakeNrem ~
sleep = fitlme(allTables2,'Base_mean ~ sleep + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(sleep)
nosleep = fitlme(allTables2,'Base_mean ~ PrePost + EarlyLate + SWRdur + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nosleep);
% base vs. prepost
compare(base, sleep,'nsim',1000)
% full vs. noprepost 
compare(nosleep,full,'nsim',1000)

% ~ duration ~
dur = fitlme(allTables2,'Base_mean ~ SWRdur + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(dur)
nodur = fitlme(allTables2,'Base_mean ~ PrePost + EarlyLate + sleep + SWRpower + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nodur);
% base vs. prepost
compare(base, dur,'nsim',1000)
% full vs. noprepost 
compare(nodur,full,'nsim',1000)

% ~ swrpower ~
power = fitlme(allTables2,'Base_mean ~ SWRpower + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(power)
nopower = fitlme(allTables2,'Base_mean ~ PrePost + EarlyLate + sleep + SWRdur + swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(nopower);
% base vs. prepost
compare(base, power,'nsim',1000)
% full vs. noprepost 
compare(nopower,full,'nsim',1000)

% ~ swrID ~
id = fitlme(allTables2,'Base_mean ~ swrID + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(id)
noid = fitlme(allTables2,'Base_mean ~ PrePost + EarlyLate + sleep + SWRdur + SWRpower + (1|mouseID) + (1|sess:mouseID)');%(ProcPeakTbl,'Peak ~ BeforeAfter + swrID + PrePost + (1|mouseID) + (1|sess)');
disp(noid);
% base vs. prepost
compare(base, id,'nsim',1000)
% full vs. noprepost 
compare(noid,full,'nsim',1000)

%% Accumulator function (on the fly, no big matrix)
function [grand_mean, grand_sd, n_events] = compute_peth(allTables2, mask, n_samples)
    n_acc    = zeros(1, n_samples);
    sum_acc  = zeros(1, n_samples);
    sum2_acc = zeros(1, n_samples);

    idx = find(mask);
    for i = 1:length(idx)
        row = idx(i);
        if allTables2.PrePost(row) == 1
            sig = single(allTables2.TwosPreProc{row}.signal(:)');   % single precision
        elseif allTables2.PrePost(row) == 2
            sig = single(allTables2.TwosPostProc{row}.signal(:)');  % single precision
        else
            continue
        end
        if length(sig) == n_samples
            valid         = ~isnan(sig);
            n_acc         = n_acc    + valid;
            sum_acc       = sum_acc  + sig;
            sum2_acc      = sum2_acc + sig.^2;
        end
    end

    grand_mean = sum_acc ./ n_acc;
    grand_sd   = sqrt((sum2_acc - (sum_acc.^2) ./ n_acc) ./ (n_acc - 1));
    n_events   = max(n_acc);
end