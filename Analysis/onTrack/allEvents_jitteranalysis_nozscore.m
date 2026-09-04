%% ============================================================
% PETH + FAST JITTER + PEAK-CORRECTED EMPIRICAL P-VALUE
%
% IMPORTANT:
% The input DA signals are ALREADY z-scored.
% Analysis:
%
% 1. Calculate the REAL mean PETH for each condition.
% 2. Find the REAL peak between 0 and <1 second.
% 3. Record the time and amplitude of that peak.
% 4. Generate 1000 jittered mean PETHs.
% 5. For every jitter:
%       find the maximum between 0 and <1 second.
% 6. Compare the REAL peak to the distribution of
%    1000 jitter peaks.
% 7. Calculate an empirical, one-sided p-value.
% ============================================================

% PARAMETERS
Fs = 1600;
n_samples = 12801;
tvec = linspace(-4,4,n_samples);
n_jitter = 1000;
jitter_sec = 1;       % +/- 1 second


% FIGURE
figure('Position',[100 100 700 500]);
hold on;

% GROUP DEFINITIONS
group_defs = {
    'Pre-Track',  (allTables2.PrePost == categorical(1)), [0.4 0.6 0.8];
    'Track',      (allTables2.PrePost == categorical(2)), [0.9 0.4 0.3];
    'Post-Track', (allTables2.PrePost == categorical(3)), [0.2 0.7 0.3];
};


% STORAGE FOR RESULTS
results = struct();
% LOOP THROUGH GROUPS

for g = 1:size(group_defs,1)

    label = group_defs{g,1};
    mask = group_defs{g,2};
    col = group_defs{g,3};

    % REAL PETH
    [grand_mean, grand_sd, n_events] = compute_peth(allTables2,mask,n_samples);


    [jitter_mean, jitter_sem, ...
     jitter_null_mean, jitter_null_sd, ...
     jitter_peths] = ...
        fast_jitter_peth(allTables2,mask,n_samples,...
                         Fs,n_jitter,jitter_sec);

    fill([tvec fliplr(tvec)], ...
         [jitter_mean + jitter_sem ...
          fliplr(jitter_mean - jitter_sem)], ...
         col, ...
         'FaceAlpha',0.10, ...
         'EdgeColor','none', ...
         'HandleVisibility','off');


    plot(tvec,jitter_mean,'--', ...
         'Color',col, ...
         'LineWidth',1.5, ...
         'DisplayName',sprintf('%s jitter',label));

    plot(tvec,grand_mean,'-', ...
         'Color',col, ...
         'LineWidth',2, ...
         'DisplayName',sprintf('%s (n=%d)',label,n_events));


    % FIND REAL PEAK
    %
    % Search only from:
    %
    %       0 <= time < 1 second
    %
    % The peak is found from the REAL mean PETH.
    %
    % IMPORTANT:
    % This is a peak in the already-z-scored DA signal.
    % No additional normalization is performed.

    post_idx = tvec >= 0 & tvec < 1;
    post_times = tvec(post_idx);
    real_post = grand_mean(post_idx);

    % Find real peak

    [real_peak,local_idx] = max(real_post);


    % Peak time relative to SWR

    peak_time = post_times(local_idx);


    % Convert local index to full time-vector index

    full_idx = find(post_idx);

    peak_idx = full_idx(local_idx);


    %% ========================================================
    % VALUE OF EVERY JITTER PETH AT THE REAL PEAK TIME
    %
    % This is useful as a secondary analysis.
    %
    % We take the value from every jitter PETH at the EXACT
    % time where the real PETH peaked.
    %
    % Result:
    %
    %       1000 jitter values
    %
    % =========================================================

    null_at_real_peak = jitter_peths(:,peak_idx);
    null_at_real_peak = null_at_real_peak(~isnan(null_at_real_peak));


    %% ========================================================
    % FIXED-TIME EMPIRICAL P-VALUE
    %
    % Question:
    %
    % "At the time where the real PETH peaked, is the real
    % value larger than expected from jitter?"
    %
    % One-sided test.
    % =========================================================

    n_null_fixed = length(null_at_real_peak);
    if n_null_fixed > 0
        p_fixed = (sum(null_at_real_peak >= real_peak) + 1) / (n_null_fixed + 1);
    else
        p_fixed = NaN;
    end

    % FIND PEAK OF EVERY JITTER PETH
    null_peak_values = NaN(n_jitter,1);
    null_peak_times = NaN(n_jitter,1);

    for j = 1:n_jitter

        this_jitter = jitter_peths(j,post_idx);


        % Skip if entirely NaN

        if all(isnan(this_jitter))
            continue
        end


        % Find maximum for this jitter

        [this_peak,this_idx] = max(this_jitter);


        % Store peak amplitude

        null_peak_values(j) = this_peak;


        % Store peak time

        null_peak_times(j) = post_times(this_idx);

    end


    %% ========================================================
    % REMOVE INVALID NULL PEAKS
    % =========================================================

    valid_null_peaks = ~isnan(null_peak_values);

    null_peak_values = ...
        null_peak_values(valid_null_peaks);

    null_peak_times = ...
        null_peak_times(valid_null_peaks);


    %% ========================================================
    % PEAK-CORRECTED EMPIRICAL P-VALUE
    %
    % PRIMARY STATISTICAL TEST
    %
    % Compare:
    %
    %       REAL peak
    %
    % against:
    %
    %       maximum peak from each jitter realization
    %
    % This accounts for searching across the entire 0-1 sec
    % window for the peak.
    % =========================================================

    n_null_peaks = length(null_peak_values);

    if n_null_peaks > 0

        p_peak = (sum(null_peak_values >= real_peak) + 1) /(n_null_peaks + 1);
    else
        p_peak = NaN;
    end


    % DESCRIPTIVE JITTER VALUES AT REAL PEAK TIME

    null_mean_at_real_peak = mean(null_at_real_peak,'omitnan');
    null_sd_at_real_peak = std(null_at_real_peak,'omitnan');

    % STORE RESULTS

    results(g).label = label;

    results(g).n_events = n_events;

    results(g).real_peth = grand_mean;

    results(g).real_sd = grand_sd;

    results(g).peak_time = peak_time;

    results(g).peak_idx = peak_idx;

    results(g).real_peak = real_peak;

    results(g).jitter_peths = jitter_peths;

    results(g).jitter_mean = jitter_mean;

    results(g).jitter_sem = jitter_sem;

    results(g).jitter_null_mean = jitter_null_mean;

    results(g).jitter_null_sd = jitter_null_sd;

    results(g).null_at_real_peak = null_at_real_peak;

    results(g).null_peak_values = null_peak_values;

    results(g).null_peak_times = null_peak_times;

    results(g).null_mean_at_real_peak = ...
        null_mean_at_real_peak;

    results(g).null_sd_at_real_peak = ...
        null_sd_at_real_peak;

    results(g).p_fixed = p_fixed;

    results(g).p_peak = p_peak;


    %% ========================================================
    % DISPLAY RESULTS
    % =========================================================

    fprintf('\n%s\n',label);

    fprintf('--------------------------------------\n');

    fprintf('N events:             %d\n',n_events);

    fprintf('Real peak time:       %.4f s\n',peak_time);

    fprintf('Real peak:            %.4f z\n',real_peak);

    fprintf('Jitter mean at time:  %.4f z\n', ...
        null_mean_at_real_peak);

    fprintf('Jitter SD at time:    %.4f z\n', ...
        null_sd_at_real_peak);

    fprintf('Fixed-time p:%.5f\n',p_fixed);

    fprintf('Peak-corrected p: %.5f\n',p_peak);


    if ~isnan(p_peak)

        if p_peak < 0.001

            fprintf('RESULT: Significant (p < 0.001)\n');

        elseif p_peak < 0.01

            fprintf('RESULT: Significant (p < 0.01)\n');

        elseif p_peak < 0.05

            fprintf('RESULT: Significant (p < 0.05)\n');

        else

            fprintf('RESULT: Not significant (p >= 0.05)\n');

        end

    end

end


%% ============================================================
% FORMATTING
% ============================================================

xline(0,'k--','LineWidth',1.2);

yline(0,'k:','LineWidth',0.8);

xlabel('Time relative to SWR (s)');

ylabel('DA signal (z-score)');

title('Pre-Track vs Track vs Post-Track');

xlim([-4 4]);

legend('Location','northwest','FontSize',9);

box off;

set(gca,'FontSize',12);

set(gcf,'renderer','painters');


%% ============================================================
% SUMMARY TABLE
% ============================================================

Summary = table( ...
    {results.label}', ...
    [results.n_events]', ...
    [results.peak_time]', ...
    [results.real_peak]', ...
    [results.null_mean_at_real_peak]', ...
    [results.null_sd_at_real_peak]', ...
    [results.p_fixed]', ...
    [results.p_peak]', ...
    'VariableNames', { ...
    'Condition', ...
    'N_Events', ...
    'Peak_Time_s', ...
    'Real_Peak_z', ...
    'Jitter_Mean_at_Real_Peak', ...
    'Jitter_SD_at_Real_Peak', ...
    'Fixed_Time_P', ...
    'Peak_Corrected_P'});


disp(Summary);


%% FUNCTION: COMPUTE REAL PETH
function [grand_mean,grand_sd,n_events] = compute_peth(allTables2,mask,n_samples)

    n_acc = zeros(1,n_samples);

    sum_acc = zeros(1,n_samples);

    sum2_acc = zeros(1,n_samples);

    idx = find(mask);

    n_used = 0;

    for i = 1:length(idx)

        row = idx(i);

        if allTables2.PrePost(row) == '1'
            sig = single(allTables2.TwosPreProc{row}.signal(:)');
        elseif allTables2.PrePost(row) == '2'
            sig = single(allTables2.TwosTrackProc{row}.signal(:)');
        elseif allTables2.PrePost(row) == '3'
            sig = single(allTables2.TwosPostProc{row}.signal(:)');

        else
            continue
        end


        if length(sig) ~= n_samples
            continue
        end

        baseline = mean(sig(1:4800),'omitnan');
        sig = sig - baseline;

        valid = ~isnan(sig);

        sig(~valid) = 0;

        n_acc = n_acc + valid;
        sum_acc = sum_acc + sig;
        sum2_acc = sum2_acc + sig.^2;


        n_used = n_used + 1;

    end

    grand_mean = NaN(1,n_samples);
    valid_mean = n_acc > 0;
    grand_mean(valid_mean) = sum_acc(valid_mean) ./ n_acc(valid_mean);
    grand_sd = NaN(1,n_samples);
    valid_n = n_acc > 1;
    grand_sd(valid_n) = sqrt((sum2_acc(valid_n) - (sum_acc(valid_n).^2 ./ n_acc(valid_n))) ./ (n_acc(valid_n) - 1));
    n_events = n_used;

end


%% ============================================================
% FUNCTION: FAST JITTER PETH

function [jitter_mean,jitter_sem, ...
          jitter_null_mean,jitter_null_sd, ...
          jitter_peths] = ...
          fast_jitter_peth(allTables2,mask,n_samples,...
                           Fs,n_jitter,jitter_sec)

    idx = find(mask);
    n_events = length(idx);
    event_signals = NaN(n_events,n_samples);
    n_used = 0;

    for i = 1:n_events
        row = idx(i);
        if allTables2.PrePost(row) == '1'
            sig = single( allTables2.TwosPreProc{row}.signal(:)');

        elseif allTables2.PrePost(row) == '2'
            sig = single( allTables2.TwosTrackProc{row}.signal(:)');

        elseif allTables2.PrePost(row) == '3'
            sig = single(allTables2.TwosPostProc{row}.signal(:)');
        else
            continue
        end

        if length(sig) ~= n_samples
            continue
        end

        n_used = n_used + 1;
        event_signals(n_used,:) = sig;

    end

    event_signals = event_signals(1:n_used,:);

    jitter_peths = NaN(n_jitter,n_samples);

    for j = 1:n_jitter
        shifts_sec = (2*rand(n_used,1)-1) * jitter_sec;
        shifts_samples = round(shifts_sec * Fs);
        jittered = NaN(n_used,n_samples);
        for e = 1:n_used
            sig = event_signals(e,:);
            sig = circshift( sig, shifts_samples(e));
            baseline = mean(sig(1:4800),'omitnan');
            jittered(e,:) = sig-baseline;
        end
        jitter_peths(j,:) = mean(jittered,1,'omitnan');
    end

    jitter_null_mean = mean(jitter_peths,1,'omitnan');
    jitter_null_sd = std(jitter_peths,0,1,'omitnan');
    n_valid = sum(~isnan(jitter_peths),1);

    jitter_sem = jitter_null_sd ./ sqrt(n_valid);
    jitter_mean = jitter_null_mean;

end



%% What is the peak at 380ms and waht is the pvalue of the permutation test between them. 
% Pre 
real_peak_pre = results(1).real_peth(6780)
real_peak_track = results(2).real_peth(6780)
real_peak_post = results(3).real_peth(6780)

null_peak_values_pre = results(1).jitter_peths(:,6780);
null_peak_values_track = results(2).jitter_peths(:,6780);
null_peak_values_post = results(3).jitter_peths(:,6780);

mean(null_peak_values_pre)
mean(null_peak_values_track)
mean(null_peak_values_post)

p_peak_pre = (sum(null_peak_values_pre >= real_peak_pre) + 1) /(1000 + 1)
p_peak_track = (sum(null_peak_values_track >= real_peak_track) + 1) /(1000 + 1)
p_peak_post = (sum(null_peak_values_post >= real_peak_post) + 1) /(1000 + 1)



