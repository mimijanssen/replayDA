%% Step 1: Average SWRpower per session per mouse per PrePost
% (session-level averages first)
[G, mouse_ids, sess_ids, prepost_ids] = findgroups(...
    allTables.mouseID, allTables.sess, allTables.PrePost);

sess_avg = splitapply(@mean, allTables.SWRpower, G);

% Build session-level table
sess_tbl = table(mouse_ids, sess_ids, prepost_ids, sess_avg, ...
    'VariableNames', {'mouseID','sess','PrePost','meanSWRpower'});

%% Step 2: Average across sessions per mouse per PrePost
% (mouse-level averages)
[G2, mouse_ids2, prepost_ids2] = findgroups(...
    sess_tbl.mouseID, sess_tbl.PrePost);

mouse_avg = splitapply(@mean, sess_tbl.meanSWRpower, G2);

mouse_tbl = table(mouse_ids2, prepost_ids2, mouse_avg, ...
    'VariableNames', {'mouseID','PrePost','meanSWRpower'});

%% Step 3: Grand average and SEM per PrePost
pre_vals  = mouse_tbl.meanSWRpower(mouse_tbl.PrePost == 1);
post_vals = mouse_tbl.meanSWRpower(mouse_tbl.PrePost == 2);

grand_mean = [mean(pre_vals),  mean(post_vals)];
grand_sem  = [std(pre_vals)/sqrt(length(pre_vals)), ...
              std(post_vals)/sqrt(length(post_vals))];

%% Step 4: Plot
figure; clf;

bar_colors = [0.6 0.6 0.6; 0.2 0.2 0.2]; % light gray = pre, dark gray = post
x_pos = [1 2];
jitter_width = 0.08;

hold on;

% Draw bars
for i = 1:2
    b = bar(x_pos(i), grand_mean(i), 0.5);
    b.FaceColor = bar_colors(i,:);
    b.EdgeColor = 'none';
end

% Error bars (SEM)
errorbar(x_pos, grand_mean, grand_sem, grand_sem, ...
    'k.', 'LineWidth', 1.5, 'CapSize', 8);

% Jittered individual mouse points
vals = {pre_vals, post_vals};
for i = 1:2
    jitter = (rand(length(vals{i}),1) - 0.5) * jitter_width;
    scatter(x_pos(i) + jitter, vals{i}, ...
        50, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
end

% Connect paired points across pre and post
for i = 1:length(pre_vals)
    jitter_pre  = (rand - 0.5) * jitter_width;
    jitter_post = (rand - 0.5) * jitter_width;
    plot([1 + jitter_pre, 2 + jitter_post], ...
         [pre_vals(i), post_vals(i)], ...
         'k-', 'LineWidth', 0.8, 'Color', [0 0 0 0.3]);
end

% Labels
ylabel('Mean SWR Power (z-score)')
xticks([1 2])
xticklabels({'Pre', 'Post'})
xlim([0.5 2.5])
title('SWR Power: Pre vs Post')
box off

%% Step 5: Stats (Wilcoxon signed-rank, paired)
[p, ~, stats] = signrank(pre_vals, post_vals);
fprintf('Wilcoxon signed-rank (paired): p = %.4f, stat = %.3f\n', ...
    p, stats.signedrank);

% Add p-value to plot
y_max = max([pre_vals; post_vals]) * 1.1;
plot([1 2], [y_max y_max], 'k-', 'LineWidth', 1)
text(1.5, y_max * 1.02, sprintf('p = %.3f', p), ...
    'HorizontalAlignment', 'center', 'FontSize', 10)