%% Make RPE barchart 
% based on dF or AUC for high and low trials 
%% Load Data  
folders = {'F:\M433\avg_data\RPE_ttest',
    'F:\M452\avg_data\RPE_ttest',
    'F:\M453\avg_data\RPE_Ttest',
    'F:\M460\avg_data\RPE_Ttest',
    'F:\M533\avg_data\RPE_Ttest',
    'F:\M545\avg_data\RPE_Ttest',
    'F:\M547\avg_data\RPE_Ttest',
    'F:\M548\avg_data\RPE_Ttest'};

%% Loop through each folder and load all .mat files
allData = struct();  % store everything here

for i = 1:length(folders)
    folderPath = folders{i};
    % Get subject ID from folder path (e.g., 'M433')
    parts = strsplit(folderPath, '\');
    subjectID = parts{2};  % grabs the MXX part
    % Find all .mat files in the folder
    files = dir(fullfile(folderPath, '*.mat'));
    if isempty(files)
        warning('No .mat files found in %s', folderPath);
        continue;
    end
      % Load each file
    for j = 1:length(files)
        filePath = fullfile(files(j).folder, files(j).name);
        fprintf('Loading: %s\n', filePath);
        loaded = load(filePath);
        % Store under subject ID + file name (without extension)
        [~, fname, ~] = fileparts(files(j).name);
        allData.(subjectID).(fname) = loaded;
    end
end

disp('All data loaded successfully.')

%% Take average for each session and then average for each mouse AUC
mouseIDs = {'M433', 'M452', 'M453', 'M460', 'M533','M545', 'M547', 'M548'};

% Use cell arrays to store session-level data per mouse
sess_high_AUC = cell(1, length(mouseIDs));
sess_low_AUC  = cell(1, length(mouseIDs));
sess_dF_high  = cell(1, length(mouseIDs));
sess_dF_low   = cell(1, length(mouseIDs));
sess_mean_high  = cell(1, length(mouseIDs));
sess_mean_low   = cell(1, length(mouseIDs));
sess_mean_med = cell(1, length(mouseIDs));

% Keep the mouse-level averages too
avg_high_AUC = zeros(1, length(mouseIDs));
avg_low_AUC  = zeros(1, length(mouseIDs));
avg_dF_high  = zeros(1, length(mouseIDs));
avg_dF_low   = zeros(1, length(mouseIDs));
avg_mean_high  = zeros(1, length(mouseIDs));
avg_mean_low   = zeros(1, length(mouseIDs));
avg_mean_med   = zeros(1, length(mouseIDs));

for i = 1:length(mouseIDs)
    mouseID = mouseIDs{i};
    sessions = fieldnames(allData.(mouseID));

    tmp_high_AUC = zeros(1, length(sessions));
    tmp_low_AUC  = zeros(1, length(sessions));
    tmp_dF_high  = zeros(1, length(sessions));
    tmp_dF_low   = zeros(1, length(sessions));
    tmp_mean_high  = zeros(1, length(sessions));
    tmp_mean_low   = zeros(1, length(sessions));
    tmp_mean_med   = zeros(1, length(sessions));


    for j = 1:length(sessions)
        sess = sessions{j};
        tmp_high_AUC(j) = mean(allData.(mouseID).(sess).AUC_high);
        tmp_low_AUC(j)  = mean(allData.(mouseID).(sess).AUC_low);
        tmp_dF_high(j)  = mean(allData.(mouseID).(sess).dF_high);
        tmp_dF_low(j)   = mean(allData.(mouseID).(sess).dF_low);
        tmp_mean_high(j)  = mean(allData.(mouseID).(sess).mean_high);
        tmp_mean_low(j)   = mean(allData.(mouseID).(sess).mean_low);
        tmp_mean_med(j)   = mean(allData.(mouseID).(sess).mean_med);

    end

    % Store session-level arrays in cell arrays
    sess_high_AUC{i} = tmp_high_AUC;
    sess_low_AUC{i}  = tmp_low_AUC;
    sess_dF_high{i}  = tmp_dF_high;
    sess_dF_low{i}   = tmp_dF_low;
    sess_mean_high{i}  = tmp_mean_high;
    sess_mean_low{i}   = tmp_mean_low;
    sess_mean_med{i}   = tmp_mean_med;

    % Also store mouse-level averages
    avg_high_AUC(i) = mean(tmp_high_AUC);
    avg_low_AUC(i)  = mean(tmp_low_AUC);
    avg_dF_high(i)  = mean(tmp_dF_high);
    avg_dF_low(i)   = mean(tmp_dF_low);
    avg_mean_high(i)  = mean(tmp_mean_high);
    avg_mean_low(i)   = mean(tmp_mean_low);
    avg_mean_med(i)   = mean(tmp_mean_med);

end

disp('Done.')

%% MEANS 

grand_mean_low = mean(avg_mean_low,'omitnan');
grand_mean_med = mean(avg_mean_med,'omitnan');
grand_mean_high = mean(avg_mean_high,'omitnan');


fprintf('Omission Trials Mean = %.4f\n',grand_mean_low);
fprintf('Medium Trials Mean  = %.4f\n',grand_mean_med);
fprintf('High-Vol Trials Mean  = %.4f\n',grand_mean_high);


sem_low = std(avg_mean_low,'omitnan') / ...
            sqrt(sum(~isnan(avg_mean_low)))
sem_med = std(avg_mean_med,'omitnan') / ...
            sqrt(sum(~isnan(avg_mean_med)))
sem_high = std(avg_mean_high,'omitnan') / ...
            sqrt(sum(~isnan(avg_mean_high)))

%% FIGURE TIME : Mean
figure(2);
hold on;

low_c = [157,130,187]./255;
high_c = [255, 165,0]./255;

high_positions = zeros(1, length(mouseIDs));
low_positions  = zeros(1, length(mouseIDs));
high_medians   = zeros(1, length(mouseIDs));
low_medians    = zeros(1, length(mouseIDs));

count = 0;
for i = 1:length(mouseIDs)
    % --- HIGH AUC ---
    count = count + 1;
    high_positions(i) = count;
    data_high = sess_mean_high{i}';
    b_high = boxchart(count*ones(size(data_high)), data_high, ...
        'BoxFaceColor', high_c, 'MarkerStyle', 'none');

    % Extend whiskers to min/max
    minVal = min(data_high);
    maxVal = max(data_high);
    line([count-0.2, count+0.2], [minVal, minVal], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [maxVal, maxVal], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [minVal, prctile(data_high, 25)], 'Color', 'k', 'LineStyle', '-');
    line([count, count], [maxVal, prctile(data_high, 75)], 'Color', 'k', 'LineStyle', '-');
    
    high_medians(i) = median(data_high);

    % --- LOW AUC ---
    count = count + 1;
    low_positions(i) = count;

    data_low = sess_mean_low{i}';
    b_low = boxchart(count*ones(size(data_low)), data_low, ...
        'BoxFaceColor', low_c, 'MarkerStyle', 'none');
    
    % Extend whiskers to min/max
    minVal = min(data_low);
    maxVal = max(data_low);
    line([count-0.2, count+0.2], [minVal, minVal], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [maxVal, maxVal], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [minVal, prctile(data_low, 25)], 'Color', 'k', 'LineStyle', '--');
    line([count, count], [maxVal, prctile(data_low, 75)], 'Color', 'k', 'LineStyle', '--');

    low_medians(i) = median(data_low);
end 

for i = 1:length(mouseIDs)

    % Connecting line
    plot([high_positions(i), low_positions(i)], ...
         [high_medians(i),   low_medians(i)], ...
         'k-', 'LineWidth', 1.5);

    % Per-mouse significance
    [p_mouse, ~] = signrank(sess_mean_high{i}', sess_mean_low{i}')

    if p_mouse < 0.001
        sig_str = '***';
    elseif p_mouse < 0.01
        sig_str = '**';
    elseif p_mouse < 0.05
        sig_str = '*';
    else
        sig_str = '';
    end

    if ~isempty(sig_str)
        x_mid = mean([high_positions(i), low_positions(i)]);
        y_sig = max([high_medians(i), low_medians(i)]) + ...
                range([high_medians, low_medians]) * 0.08;
        text(x_mid, y_sig, sig_str, ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 16, 'FontWeight', 'bold');
    end
end

% X ticks centered between each low/high pair
xticks([1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5, 15.5]);
xticklabels(mouseIDs);
ylabel('Mean DA (z-score)');

%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
set(gca, 'fontsize', 18);
%cd ('C:\Users\mimia\Desktop\UpdatedFigs')
%exportgraphics(gcf, 'Amp_mice_recolored2.eps', 'ContentType','vector');  % Export as PDF
hold off;

%% t-tests? 
for i = 1:length(mouseIDs)
    disp(mouseIDs(i));
    low  = sess_low_AUC{1,i};
    high = sess_high_AUC{1,i};
    [h, p, ci, stats] = ttest2(low, high);
    fprintf('Mouse %d: t(%d) = %.3f, p = %.4f\n', i, stats.df, stats.tstat, p);
end

avg_high_AUC
avg_low_AUC

%% FIGURE TIME: DF
figure(2);
hold on;
count = 0;
for i = 1:length(mouseIDs)

    % --- dF LOW ---
    count = count + 1;
    data_low = sess_dF_low{i}';
    boxchart(count*ones(size(data_low)), data_low, ...
        'BoxFaceColor', low_c, 'MarkerStyle', 'none');

    minVal = min(data_low);
    maxVal = max(data_low);
    line([count-0.2, count+0.2], [minVal, minVal], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [maxVal, maxVal], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [minVal, prctile(data_low, 25)], 'Color', 'k', 'LineStyle', '-');
    line([count, count], [maxVal, prctile(data_low, 75)], 'Color', 'k', 'LineStyle', '-');

    % --- dF HIGH ---
    count = count + 1;
    data_high = sess_dF_high{i}';
    boxchart(count*ones(size(data_high)), data_high, ...
        'BoxFaceColor', high_c, 'MarkerStyle', 'none');

    minVal = min(data_high);
    maxVal = max(data_high);
    line([count-0.2, count+0.2], [minVal, minVal], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [maxVal, maxVal], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [minVal, prctile(data_high, 25)], 'Color', 'k', 'LineStyle', '--');
    line([count, count], [maxVal, prctile(data_high, 75)], 'Color', 'k', 'LineStyle', '--');

end

xticks([1.5,3.5,5.5,7.5,9.5,11.5,13.5,15.5]);
xticklabels({"M433","M452","M453","M460","M533","M545","M547","M548"});
ylabel("Mean dF (z-score)");
set(gca, 'fontsize', 18);
set(gcf, 'renderer', 'painters');