%% ============================================================
%  PETH + FAST JITTER + PEAK-BASED EMPIRICAL P-VALUE
%
%  Real PETH:
%       - Mean across all events in each condition
%       - Find peak from 0 to +1 sec
%
%  Jitter null:
%       - 1000 jittered mean PETHs
%       - Same peak-search procedure applied to every jitter
%
%  Outputs:
%       1. Real peak time
%       2. Real peak value
%       3. Jitter distribution at real peak time
%       4. Fixed-time empirical p-value
%       5. Jitter peak distribution
%       6. Peak-corrected empirical p-value
%       7. Z-score at real peak
%
% =============================================================

Fs = 1600;
n_samples = 12801;
tvec = linspace(-4,4,n_samples);

n_jitter = 1000;
jitter_sec = 1;       % +/- 1 second

figure('Position',[100 100 700 500]);
hold on;

group_defs = {
    'Pre-Track',  (allTables2.PrePost == categorical(1)), [0.4 0.6 0.8];
    'Track',      (allTables2.PrePost == categorical(2)), [0.9 0.4 0.3];
    'Post-Track', (allTables2.PrePost == categorical(3)), [0.2 0.7 0.3];
};


results = struct();

for g = 1:size(group_defs,1)

    label = group_defs{g,1};
    mask  = group_defs{g,2};
    col   = group_defs{g,3};

    [grand_mean, grand_sd, n_events] = ...
        compute_peth(allTables2,mask,n_samples);

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

    
    z_trace = (grand_mean - jitter_null_mean) ./ jitter_null_sd;


    % FIND REAL PEAK
    %
    % 0 <= time < 1 second
    %
    % IMPORTANT:
    % We find the maximum of the REAL PETH.

    post_idx = tvec >= 0 & tvec < 1;

    post_times = tvec(post_idx);
    real_post = grand_mean(post_idx);

    [real_peak, local_idx] = max(real_post);

    peak_time = post_times(local_idx);

    % Convert local index to index in full time vector
    full_idx = find(post_idx);
    peak_idx = full_idx(local_idx);


    % REAL Z-SCORE AT PEAK

    peak_z = z_trace(peak_idx);


    % NULL DISTRIBUTION AT THE REAL PEAK TIME
    %
    % Every jitter PETH is evaluated at the SAME time
    % as the observed real peak.
    %
    % Result:
    %
    %       1000 x 1 distribution
    %

    null_at_peak = jitter_peths(:,peak_idx);

    null_at_peak = null_at_peak(~isnan(null_at_peak));


    % FIXED-TIME EMPIRICAL P-VALUE
    %
    % Tests:
    %
    %   Is the real value at the observed peak time
    %   larger than expected from jitter?
    %
    % One-sided test.

    n_null = length(null_at_peak);

    p_fixed = ...
        (sum(null_at_peak >= real_peak) + 1) / ...
        (n_null + 1);


    % NULL PEAK DISTRIBUTION
    %
    % For EACH jitter:
    %
    %   find its maximum between 0 and 1 sec
    %
    % This accounts for the fact that the real peak was
    % selected by searching across the entire 0-1 sec window.

    null_peak_values = NaN(n_jitter,1);
    null_peak_times  = NaN(n_jitter,1);

    for j = 1:n_jitter

        this_jitter = jitter_peths(j,post_idx);

        if all(isnan(this_jitter))
            continue
        end

        [this_peak, this_idx] = max(this_jitter);

        null_peak_values(j) = this_peak;
        null_peak_times(j) = post_times(this_idx);

    end


    % Remove invalid jitter peaks
    valid_null_peaks = ~isnan(null_peak_values);

    null_peak_values = null_peak_values(valid_null_peaks);
    null_peak_times  = null_peak_times(valid_null_peaks);


    % PEAK-CORRECTED EMPIRICAL P-VALUE
    %
    % Compare:
    %
    %       REAL peak
    %
    % against
    %
    %       maximum peak from each jitter realization
    %

    n_null_peaks = length(null_peak_values);

    p_peak = ...
        (sum(null_peak_values >= real_peak) + 1) / ...
        (n_null_peaks + 1);


    % NULL STATISTICS AT REAL PEAK TIME

    null_mean_at_peak = mean(null_at_peak,'omitnan');
    null_sd_at_peak   = std(null_at_peak,'omitnan');


    % =========================================================
    % Z-SCORE CALCULATED DIRECTLY FROM PEAK DISTRIBUTION
    % =========================================================

    peak_z_direct = ...
        (real_peak - null_mean_at_peak) / null_sd_at_peak;


    % STORE RESULTS

    results(g).label = label;

    results(g).real_peth = grand_mean;

    results(g).real_sd = grand_sd;

    results(g).n_events = n_events;

    results(g).peak_time = peak_time;

    results(g).peak_idx = peak_idx;

    results(g).real_peak = real_peak;

    results(g).peak_z = peak_z;

    results(g).peak_z_direct = peak_z_direct;

    results(g).jitter_peths = jitter_peths;

    results(g).null_at_peak = null_at_peak;

    results(g).null_peak_values = null_peak_values;

    results(g).null_peak_times = null_peak_times;

    results(g).p_fixed = p_fixed;

    results(g).p_peak = p_peak;

    results(g).null_mean_at_peak = null_mean_at_peak;

    results(g).null_sd_at_peak = null_sd_at_peak;


    % =========================================================
    % DISPLAY RESULTS
    % =========================================================

    fprintf('\n%s\n',label);
    fprintf('--------------------------------------\n');

    fprintf('N events:         %d\n',n_events);

    fprintf('Peak time:        %.4f s\n',peak_time);

    fprintf('Real peak:        %.4f\n',real_peak);

    fprintf('Jitter mean:      %.4f\n',null_mean_at_peak);

    fprintf('Jitter SD:        %.4f\n',null_sd_at_peak);

    fprintf('Peak Z-score:     %.3f\n',peak_z_direct);

    fprintf('Fixed-time p:     %.5f\n',p_fixed);

    fprintf('Peak-corrected p: %.5f\n',p_peak);

    fprintf('\n');


    if p_peak < 0.05
        fprintf('RESULT:            Significant peak (p < 0.05)\n');
    else
        fprintf('RESULT:            Not significant peak (p >= 0.05)\n');
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
% OPTIONAL: SUMMARY TABLE
% ============================================================

Summary = table( ...
    {results.label}', ...
    [results.n_events]', ...
    [results.peak_time]', ...
    [results.real_peak]', ...
    [results.null_mean_at_peak]', ...
    [results.null_sd_at_peak]', ...
    [results.peak_z_direct]', ...
    [results.p_fixed]', ...
    [results.p_peak]', ...
    'VariableNames', { ...
    'Condition', ...
    'N_Events', ...
    'Peak_Time_s', ...
    'Real_Peak', ...
    'Jitter_Mean', ...
    'Jitter_SD', ...
    'Peak_Z', ...
    'Fixed_Time_P', ...
    'Peak_Corrected_P'});


disp(Summary);


%% ============================================================
% FUNCTION: COMPUTE REAL PETH
% ============================================================

function [grand_mean, grand_sd, n_events] = ...
    compute_peth(allTables2, mask, n_samples)

    n_acc    = zeros(1,n_samples);
    sum_acc  = zeros(1,n_samples);
    sum2_acc = zeros(1,n_samples);

    idx = find(mask);

    fs = 1600;

    n_used = 0;

    for i = 1:length(idx)

        row = idx(i);


        % -----------------------------------------------------
        % Get signal
        % -----------------------------------------------------

        if allTables2.PrePost(row) == '1'

            sig = single( ...
                allTables2.TwosPreProc{row}.signal(:)');

        elseif allTables2.PrePost(row) == '2'

            sig = single( ...
                allTables2.TwosTrackProc{row}.signal(:)');

        elseif allTables2.PrePost(row) == '3'

            sig = single( ...
                allTables2.TwosPostProc{row}.signal(:)');

        else
            continue
        end


        % -----------------------------------------------------
        % Check signal length
        % -----------------------------------------------------

        if length(sig) ~= n_samples
            continue
        end


        % -----------------------------------------------------
        % Baseline
        %
        % IMPORTANT:
        % No jitter is applied to the real PETH.
        % -----------------------------------------------------

        baseline = mean(sig(1:4800),'omitnan');

        sig = sig - baseline;


        % -----------------------------------------------------
        % Handle NaNs
        % -----------------------------------------------------

        valid = ~isnan(sig);

        sig(~valid) = 0;

        n_acc    = n_acc + valid;

        sum_acc  = sum_acc + sig;

        sum2_acc = sum2_acc + sig.^2;

        n_used = n_used + 1;

    end


    % ---------------------------------------------------------
    % Grand mean
    % ---------------------------------------------------------

    grand_mean = NaN(1,n_samples);

    valid_mean = n_acc > 0;

    grand_mean(valid_mean) = ...
        sum_acc(valid_mean) ./ n_acc(valid_mean);


    % ---------------------------------------------------------
    % Grand SD
    % ---------------------------------------------------------

    grand_sd = NaN(1,n_samples);

    valid_n = n_acc > 1;

    grand_sd(valid_n) = sqrt( ...
        (sum2_acc(valid_n) - ...
        (sum_acc(valid_n).^2 ./ n_acc(valid_n))) ./ ...
        (n_acc(valid_n) - 1));


    % ---------------------------------------------------------
    % Number of events
    % ---------------------------------------------------------

    n_events = n_used;

end


%% ============================================================
% FUNCTION: FAST JITTER PETH
%
% Returns:
%
% jitter_mean
% jitter_sem
% jitter_null_mean
% jitter_null_sd
% jitter_peths
%
% jitter_peths is:
%
%       n_jitter x n_samples
%
% This is what allows us to calculate the empirical p-value.
% ============================================================

function [jitter_mean, jitter_sem, ...
          jitter_null_mean, jitter_null_sd, ...
          jitter_peths] = ...
          fast_jitter_peth(allTables2,mask,n_samples,...
                           Fs,n_jitter,jitter_sec)


    % ---------------------------------------------------------
    % Find events
    % ---------------------------------------------------------

    idx = find(mask);

    n_events = length(idx);


    % ---------------------------------------------------------
    % Preallocate event signals
    % ---------------------------------------------------------

    event_signals = NaN(n_events,n_samples);


    % ---------------------------------------------------------
    % Load all signals once
    % ---------------------------------------------------------

    n_used = 0;

    for i = 1:n_events

        row = idx(i);


        if allTables2.PrePost(row) == '1'

            sig = single( ...
                allTables2.TwosPreProc{row}.signal(:)');

        elseif allTables2.PrePost(row) == '2'

            sig = single( ...
                allTables2.TwosTrackProc{row}.signal(:)');

        elseif allTables2.PrePost(row) == '3'

            sig = single( ...
                allTables2.TwosPostProc{row}.signal(:)');

        else
            continue
        end


        % Check length

        if length(sig) ~= n_samples
            continue
        end


        % Store

        n_used = n_used + 1;

        event_signals(n_used,:) = sig;

    end


    % Remove unused rows

    event_signals = event_signals(1:n_used,:);


    % ---------------------------------------------------------
    % Preallocate jitter PETH matrix
    %
    % Rows = jitter iterations
    % Columns = time
    % ---------------------------------------------------------

    jitter_peths = NaN(n_jitter,n_samples);


    % ---------------------------------------------------------
    % Generate jittered PETHs
    % ---------------------------------------------------------

    for j = 1:n_jitter


        % -----------------------------------------------------
        % Random shift for EACH EVENT
        %
        % Uniformly distributed between:
        %
        %       -jitter_sec and +jitter_sec
        %
        % -----------------------------------------------------

        shifts_sec = ...
            (2*rand(n_used,1)-1) * jitter_sec;

        shifts_samples = ...
            round(shifts_sec * Fs);


        % -----------------------------------------------------
        % Jitter every event
        % -----------------------------------------------------

        jittered = NaN(n_used,n_samples);


        for e = 1:n_used

            sig = event_signals(e,:);

            sig = circshift(sig,shifts_samples(e));


            % -------------------------------------------------
            % Baseline AFTER jitter
            % -------------------------------------------------

            baseline = mean(sig(1:4800),'omitnan');

            sig = sig - baseline;

            jittered(e,:) = sig;

        end


        % -----------------------------------------------------
        % Mean across events
        % -----------------------------------------------------

        jitter_peths(j,:) = ...
            mean(jittered,1,'omitnan');

    end


    % ---------------------------------------------------------
    % Mean of the 1000 jitter PETHs
    % ---------------------------------------------------------

    jitter_null_mean = ...
        mean(jitter_peths,1,'omitnan');


    % ---------------------------------------------------------
    % SD of the 1000 jitter PETHs
    % ---------------------------------------------------------

    jitter_null_sd = ...
        std(jitter_peths,0,1,'omitnan');


    % ---------------------------------------------------------
    % SEM across jitter iterations
    % ---------------------------------------------------------

    n_valid = sum(~isnan(jitter_peths),1);

    jitter_sem = ...
        jitter_null_sd ./ sqrt(n_valid);


    % ---------------------------------------------------------
    % Same mean returned as jitter_mean for compatibility
    % with your original plotting code.
    % ---------------------------------------------------------

    jitter_mean = jitter_null_mean;

end
