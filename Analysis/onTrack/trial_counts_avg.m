%% Trial count analysis across mice
% For each mouse, loops through session folders containing "recording"
% in the folder name, loads events, and computes trial_count per session.
% Then computes: mean trial_count per mouse (average across sessions),
% and mean/std/sem across mice (average of mouse-level averages).

mice = {'M646', 'M648', 'M650', 'M654'};
data_root = 'D:\';

mouse_means = nan(1, numel(mice));
all_trial_counts = cell(1, numel(mice)); % keep session-level data per mouse, in case you want it

for m = 1:numel(mice)
    mouse_id = mice{m};
    mouse_dir = fullfile(data_root, mouse_id);

    d = dir(mouse_dir);
    d = d([d.isdir]);                       % keep only folders
    session_names = {d.name};
    is_recording = contains(session_names, 'recording');
    session_names = session_names(is_recording);

    if isempty(session_names)
        warning('No "recording" session folders found for %s', mouse_id);
    end

    trial_counts = nan(1, numel(session_names));

    for s = 1:numel(session_names)
        session_path = fullfile(mouse_dir, session_names{s});
        cd(session_path);

        clear ExpKeys evt3
        LoadExpKeys();

        cfg_evt = [];
        cfg_evt.eventList = ExpKeys.eventList;
        cfg_evt.eventLabel = ExpKeys.eventLabel;
        evt3 = LoadEvents(cfg_evt);

        evt_ordered = sort([evt3.t{1}, evt3.t{2}]);

        if length(evt_ordered) < 60
            warning('%s has only %d events (<60) — skipping trim', ...
                session_names{s}, length(evt_ordered));
        else
            evt_ordered = evt_ordered(1:60); % only want 60 of them
        end

        trial_counts(s) = length(evt_ordered);

        fprintf('%s | %s | trial_count = %d\n', ...
            mouse_id, session_names{s}, trial_counts(s));
    end

    all_trial_counts{m} = trial_counts;
    mouse_means(m) = mean(trial_counts, 'omitnan');
end

%% Step 2: average over mice (using each mouse's session-average)
grand_mean = mean(mouse_means, 'omitnan');
grand_std  = std(mouse_means, 'omitnan');
grand_sem  = grand_std / sqrt(sum(~isnan(mouse_means)));

fprintf('\n--- Per-mouse averages ---\n');
for m = 1:numel(mice)
    fprintf('%s: mean trial_count = %.2f (n sessions = %d)\n', ...
        mice{m}, mouse_means(m), numel(all_trial_counts{m}));
end

fprintf('\n--- Across mice ---\n');
fprintf('Mean: %.2f\n', grand_mean);
fprintf('Std: %.2f\n', grand_std);
fprintf('SEM: %.2f\n', grand_sem);