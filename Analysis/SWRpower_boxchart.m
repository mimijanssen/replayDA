%% load session data into a large matrix 

cd  F:\SWR_DA_MegaMatrix_1s_withSWRinfo
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



%% Step 1: Average SWRpower per session per mouse per PrePost
[G, mouse_ids, sess_ids, prepost_ids] = findgroups(...
    allTables.mouseID, allTables.sess, allTables.PrePost);
sess_avg = splitapply(@mean, allTables.SWRpower, G);
sess_tbl = table(mouse_ids, sess_ids, prepost_ids, sess_avg, ...
    'VariableNames', {'mouseID','sess','PrePost','meanSWRpower'});

%% Step 2: Mouse-level averages
[G2, mouse_ids2, prepost_ids2] = findgroups(...
    sess_tbl.mouseID, sess_tbl.PrePost);
mouse_avg = splitapply(@mean, sess_tbl.meanSWRpower, G2);
mouse_tbl = table(mouse_ids2, prepost_ids2, mouse_avg, ...
    'VariableNames', {'mouseID','PrePost','meanSWRpower'});

%% Step 3: Separate pre and post
pre_vals  = mouse_tbl.meanSWRpower(mouse_tbl.PrePost == 1);
post_vals = mouse_tbl.meanSWRpower(mouse_tbl.PrePost == 2);

%%
[h, p, ci, stats] = ttest2(pre_vals,post_vals)
mean(pre_vals)
std(pre_vals)
mean(post_vals)
std(post_vals)
%% Step 4: Plot
% define your colors to match your existing figures
color2  = [0.5 0.5 0.5];   % box color
low_c   = [0.3 0.3 0.8];   % scatter point color -- change to match your palette

figure;
h = boxchart([pre_vals, post_vals]);
hold on;
h.BoxFaceColor = color2;
h.WhiskerLineColor = 'k';
h.MarkerColor = 'k'; % outlier marker color

% Jittered scatter points
xDataNumeric = [1 2];
scatter(repmat(xDataNumeric(1), 1, numel(pre_vals))  + randn(1, numel(pre_vals))  * 0.05, pre_vals', ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');
scatter(repmat(xDataNumeric(2), 1, numel(post_vals)) + randn(1, numel(post_vals)) * 0.05, post_vals', ...
    60, 'MarkerFaceColor', low_c, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o');

% Connect paired points
% for i = 1:length(pre_vals)
%     plot([1 + randn*0.05, 2 + randn*0.05], ...
%          [pre_vals(i), post_vals(i)], ...
%          'k-', 'LineWidth', 0.8, 'Color', [0 0 0 0.3]);
% end

% Stats
% [p, ~, stats] = signrank(pre_vals, post_vals);
% fprintf('Wilcoxon signed-rank (paired): p = %.4f, stat = %.3f\n', p, stats.signedrank);

% Significance bracket
% y_max = max([pre_vals; post_vals]) * 1.15;
% plot([1 2], [y_max y_max], 'k-', 'LineWidth', 1)
% text(1.5, y_max * 1.02, sprintf('p = %.3f', p), ...
%     'HorizontalAlignment', 'center', 'FontSize', 14)

% Labels and formatting
xticklabels(["Pre-Track" "Post-Track"])
ylabel('Mean SWR Power (z-score)')
xlabel('Rest Session')
title('SWR Power: Pre vs Post')
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
set(gca, 'fontsize', 18)
set(gcf, 'renderer', 'painters');
fontname("AvenirNext LT Pro Regular");

% Export
%cd('C:\Users\mimia\Desktop\descriptive stats swrs')
%exportgraphics(gcf, 'SWRpower_prepost.png', 'ContentType', 'vector');
hold off