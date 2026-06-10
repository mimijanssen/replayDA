%% Make RPE barchart 
% based on dF or AUC for high and low trials 
%% Load Data  
folders = {'F:\M433\avg_data\RPE',
    'F:\M452\avg_data\RPE',
    'F:\M453\avg_data\RPE',
    'F:\M460\avg_data\RPE',
    'F:\M533\avg_data\RPE',
    'F:\M545\avg_data\RPE',
    'F:\M547\avg_data\RPE',
    'F:\M548\avg_data\RPE'};

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

%% take points 8000- 
x_values = 8001:1:12000; % 4 seconds after photobeam break


%% Take average for each session and then average for each mouse AUC
mouseIDs = {'M433', 'M452', 'M453', 'M460', 'M533','M545', 'M547', 'M548'};

% Use cell arrays to store session-level data per mouse
sess_high_avg = cell(1, length(mouseIDs));
sess_low_avg  = cell(1, length(mouseIDs));
sess_med_avg  = cell(1, length(mouseIDs));
sess_high_con   = cell(1, length(mouseIDs));
sess_med_con   = cell(1, length(mouseIDs));
sess_low_con   = cell(1, length(mouseIDs));

% Keep the mouse-level averages too
avg_high_avg = zeros(1, length(mouseIDs));
avg_low_avg  = zeros(1, length(mouseIDs));
avg_med_avg  = zeros(1, length(mouseIDs));
avg_high_con  = zeros(1, length(mouseIDs));
avg_med_con   = zeros(1, length(mouseIDs));
avg_low_con   = zeros(1, length(mouseIDs));

for i = 1:length(mouseIDs)
    mouseID = mouseIDs{i};
    sessions = fieldnames(allData.(mouseID));

    tmp_high_avg = zeros(1, length(sessions));
    tmp_low_avg  = zeros(1, length(sessions));
    tmp_med_avg  = zeros(1, length(sessions));
    tmp_high_con = zeros(1, length(sessions));
    tmp_low_con  = zeros(1, length(sessions));
    tmp_med_con  = zeros(1, length(sessions));


    for j = 1:length(sessions)
        sess = sessions{j};
        tmp_high_avg(j) = max(allData.(mouseID).(sess).avg_high(x_values));
        tmp_low_avg(j)  = min(allData.(mouseID).(sess).avg_low(x_values));
        tmp_med_avg(j)  = max(allData.(mouseID).(sess).avg_med(x_values));
        tmp_low_con(j)  = max(allData.(mouseID).(sess).avg_low(x_values));
        tmp_med_con(j)  = min(allData.(mouseID).(sess).avg_med(x_values));
        tmp_high_con(j) = min(allData.(mouseID).(sess).avg_high(x_values));
    end

    % Store session-level arrays in cell arrays

    sess_high_avg{i} = tmp_high_avg;
    sess_low_avg{i}  = tmp_low_avg; %cell(1, length(mouseIDs));
    sess_med_avg{i}  = tmp_med_avg; %cell(1, length(mouseIDs));
    sess_high_con{i}   = tmp_high_con;
    sess_med_con{i}   = tmp_med_con;
    sess_low_con{i}   = tmp_low_con;

    % Also store mouse-level averages
    avg_high_avg(i) = mean(tmp_high_avg);
    avg_low_avg(i)  = mean(tmp_low_avg);
    avg_med_avg(i)  = mean(tmp_med_avg);
    avg_high_con(i)   = mean(tmp_high_con);
    avg_low_con(i)   = mean(tmp_low_con);
    avg_med_con(i)   = mean(tmp_med_con);

end

disp('Done.')

%% New fig  with control boxes??

figure(3);
hold on;

% Colors - paired, cons are grey
high_c     = [0.85 0.15 0.15];  % red
med_c      = [0.15 0.15 0.85];  % blue
low_c      = [0.15 0.75 0.15];  % green
high_con_c = [0.6 0.6 0.6];     % dark grey
med_con_c  = [0.6 0.6 0.6];     % dark grey
low_con_c  = [0.6 0.6 0.6];     % dark grey

count = 0;
for i = 1:length(mouseIDs)

    % --- HIGH AVG + HIGH CON (paired) ---
    count = count + 1;
    data = sess_high_avg{i}';
    boxchart(count*ones(size(data)), data, 'BoxFaceColor', high_c, 'MarkerStyle', 'none');
    line([count-0.2, count+0.2], [min(data), min(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [max(data), max(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [min(data), prctile(data,25)], 'Color', 'k', 'LineStyle', '--');
    line([count, count], [max(data), prctile(data,75)], 'Color', 'k', 'LineStyle', '--');

    count = count + 1;
    data = sess_high_con{i}';
    boxchart(count*ones(size(data)), data, 'BoxFaceColor', high_con_c, 'MarkerStyle', 'none');
    line([count-0.2, count+0.2], [min(data), min(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [max(data), max(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [min(data), prctile(data,25)], 'Color', 'k', 'LineStyle', '--');
    line([count, count], [max(data), prctile(data,75)], 'Color', 'k', 'LineStyle', '--');

    % --- MED AVG + MED CON (paired) ---
    count = count + 1;
    data = sess_med_avg{i}';
    boxchart(count*ones(size(data)), data, 'BoxFaceColor', med_c, 'MarkerStyle', 'none');
    line([count-0.2, count+0.2], [min(data), min(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [max(data), max(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [min(data), prctile(data,25)], 'Color', 'k', 'LineStyle', '--');
    line([count, count], [max(data), prctile(data,75)], 'Color', 'k', 'LineStyle', '--');

    count = count + 1;
    data = sess_med_con{i}';
    boxchart(count*ones(size(data)), data, 'BoxFaceColor', med_con_c, 'MarkerStyle', 'none');
    line([count-0.2, count+0.2], [min(data), min(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [max(data), max(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [min(data), prctile(data,25)], 'Color', 'k', 'LineStyle', '--');
    line([count, count], [max(data), prctile(data,75)], 'Color', 'k', 'LineStyle', '--');

    % --- LOW AVG + LOW CON (paired) ---
    count = count + 1;
    data = sess_low_avg{i}';
    boxchart(count*ones(size(data)), data, 'BoxFaceColor', low_c, 'MarkerStyle', 'none');
    line([count-0.2, count+0.2], [min(data), min(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [max(data), max(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [min(data), prctile(data,25)], 'Color', 'k', 'LineStyle', '--');
    line([count, count], [max(data), prctile(data,75)], 'Color', 'k', 'LineStyle', '--');

    count = count + 1;
    data = sess_low_con{i}';
    boxchart(count*ones(size(data)), data, 'BoxFaceColor', low_con_c, 'MarkerStyle', 'none');
    line([count-0.2, count+0.2], [min(data), min(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [max(data), max(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [min(data), prctile(data,25)], 'Color', 'k', 'LineStyle', '--');
    line([count, count], [max(data), prctile(data,75)], 'Color', 'k', 'LineStyle', '--');

end

% X ticks centered in the middle of each group of 6
xticks([3.5, 9.5, 15.5, 21.5, 27.5, 33.5, 39.5, 45.5]);
xticklabels({"M433","M452","M453","M460","M533","M545","M547","M548"});
ylabel('Mean Avg');
set(gca, 'fontsize', 18);
set(gcf, 'renderer', 'painters');

% Legend - grab one representative boxchart handle per group
h1 = boxchart(nan, nan, 'BoxFaceColor', high_c);
h2 = boxchart(nan, nan, 'BoxFaceColor', med_c);
h3 = boxchart(nan, nan, 'BoxFaceColor', low_c);
h4 = boxchart(nan, nan, 'BoxFaceColor', high_con_c);
legend([h1 h2 h3 h4], {'High','Med','Low','Control'}, 'Location', 'best');

%%
figure(4);
hold on;

high_c = [255, 165,0]./255;  % orange
med_c  = [104,187,225]./255;  % green
low_c  = [0,104,87]./255;  % purple

count = 0;
for i = 1:length(mouseIDs)

    % --- HIGH AVG ---
    count = count + 1;
    data = sess_high_avg{i}';
    boxchart(count*ones(size(data)), data, 'BoxFaceColor', high_c, 'MarkerStyle', 'none');
    line([count-0.2, count+0.2], [min(data), min(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [max(data), max(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [min(data), prctile(data,25)], 'Color', 'k', 'LineStyle', '-');
    line([count, count], [max(data), prctile(data,75)], 'Color', 'k', 'LineStyle', '-');

    % --- MED AVG ---
    count = count + 1;
    data = sess_med_avg{i}';
    boxchart(count*ones(size(data)), data, 'BoxFaceColor', med_c, 'MarkerStyle', 'none');
    line([count-0.2, count+0.2], [min(data), min(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [max(data), max(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [min(data), prctile(data,25)], 'Color', 'k', 'LineStyle', '-');
    line([count, count], [max(data), prctile(data,75)], 'Color', 'k', 'LineStyle', '-');

    % --- LOW AVG ---
    count = count + 1;
    data = sess_low_avg{i}';
    boxchart(count*ones(size(data)), data, 'BoxFaceColor', low_c, 'MarkerStyle', 'none');
    line([count-0.2, count+0.2], [min(data), min(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [max(data), max(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [min(data), prctile(data,25)], 'Color', 'k', 'LineStyle', '-');
    line([count, count], [max(data), prctile(data,75)], 'Color', 'k', 'LineStyle', '-');

end

% X ticks centered in the middle of each group of 3
xticks([2, 5, 8, 11, 14, 17, 20, 23]);
xticklabels({"M433","M452","M453","M460","M533","M545","M547","M548"});
ylabel('Peak or Trough DA (z-score)');
set(gca, 'fontsize', 18);
set(gcf, 'renderer', 'painters');

% Legend
h1 = boxchart(nan, nan, 'BoxFaceColor', high_c);
h2 = boxchart(nan, nan, 'BoxFaceColor', med_c);
h3 = boxchart(nan, nan, 'BoxFaceColor', low_c);
%legend([h1 h2 h3], {'High','Med','Low'}, 'Location', 'best');

%%
%% Print t-test results for High vs Low per mouse
fprintf('\n========================================\n');
fprintf('  T-test Results: High vs Low per Mouse\n');
fprintf('========================================\n');
fprintf('%-8s  %-10s  %-10s  %-10s\n', 'Mouse', 't-score', 'p-value', 'Sig');
fprintf('----------------------------------------\n');

for i = 1:length(mouseIDs)
    dH = sess_high_avg{i}';
    dL = sess_low_avg{i}';

    [~, p, ~, stats] = ttest2(dH, dL)
    t = stats.tstat;

    % Significance label
    if p < 0.001
        sig = '***';
    elseif p < 0.01
        sig = '**';
    elseif p < 0.05
        sig = '*';
    else
        sig = 'ns';
    end

    fprintf('%-8s  %-10.4f  %-10.4f  %-10s\n', mouseIDs{i}, t, p, sig);
end

fprintf('========================================\n');
fprintf('Significance: * p<0.05  ** p<0.01  *** p<0.001\n\n');
%% Significance added 

figure(4);
hold on;

high_c = [255, 165,0]./255;  % orange
med_c  = [104,187,225]./255;  % green
low_c  = [0,104,87]./255;  % purple


count = 0;
x_positions = zeros(length(mouseIDs), 3); % store x pos for each box [high, med, low]

for i = 1:length(mouseIDs)

    % --- HIGH AVG ---
    count = count + 1;
    x_positions(i,1) = count;
    data = sess_high_avg{i}';
    boxchart(count*ones(size(data)), data, 'BoxFaceColor', high_c, 'MarkerStyle', 'none');
    line([count-0.2, count+0.2], [min(data), min(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [max(data), max(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [min(data), prctile(data,25)], 'Color', 'k', 'LineStyle', '--');
    line([count, count], [max(data), prctile(data,75)], 'Color', 'k', 'LineStyle', '--');

    % --- MED AVG ---
    count = count + 1;
    x_positions(i,2) = count;
    data = sess_med_avg{i}';
    boxchart(count*ones(size(data)), data, 'BoxFaceColor', med_c, 'MarkerStyle', 'none');
    line([count-0.2, count+0.2], [min(data), min(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [max(data), max(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [min(data), prctile(data,25)], 'Color', 'k', 'LineStyle', '--');
    line([count, count], [max(data), prctile(data,75)], 'Color', 'k', 'LineStyle', '--');

    % --- LOW AVG ---
    count = count + 1;
    x_positions(i,3) = count;
    data = sess_low_avg{i}';
    boxchart(count*ones(size(data)), data, 'BoxFaceColor', low_c, 'MarkerStyle', 'none');
    line([count-0.2, count+0.2], [min(data), min(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [max(data), max(data)], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [min(data), prctile(data,25)], 'Color', 'k', 'LineStyle', '--');
    line([count, count], [max(data), prctile(data,75)], 'Color', 'k', 'LineStyle', '--');

end

% --- Significance stars ---
% Helper function to convert p to stars
% --- Significance stars ---
for i = 1:length(mouseIDs)
    xH = x_positions(i,1);
    xM = x_positions(i,2);
    xL = x_positions(i,3);

    dH = sess_high_avg{i}';
    dM = sess_med_avg{i}';
    dL = sess_low_avg{i}';

    % Run t-tests
    [~, p_HM] = ttest2(dH, dM);
    [~, p_HL] = ttest2(dH, dL);
    [~, p_ML] = ttest2(dM, dL);

    % Stagger bracket heights
    y1 = y_max + bracket_gap;
    y2 = y_max + bracket_gap*2.5;
    y3 = y_max + bracket_gap*4;

    % High vs Med
    line([xH, xH, xM, xM], [y1, y1+bracket_gap*0.4, y1+bracket_gap*0.4, y1], 'Color', 'k', 'LineWidth', 1);
    text((xH+xM)/2, y1+bracket_gap*0.5, getStars(p_HM), 'HorizontalAlignment', 'center', 'FontSize', 14);

    % Med vs Low
    line([xM, xM, xL, xL], [y2, y2+bracket_gap*0.4, y2+bracket_gap*0.4, y2], 'Color', 'k', 'LineWidth', 1);
    text((xM+xL)/2, y2+bracket_gap*0.5, getStars(p_ML), 'HorizontalAlignment', 'center', 'FontSize', 14);

    % High vs Low
    line([xH, xH, xL, xL], [y3, y3+bracket_gap*0.4, y3+bracket_gap*0.4, y3], 'Color', 'k', 'LineWidth', 1);
    text((xH+xL)/2, y3+bracket_gap*0.5, getStars(p_HL), 'HorizontalAlignment', 'center', 'FontSize', 14);
end

% % Get y axis limit to stack brackets above data
% y_max = max(ylim);
% bracket_gap = (max(ylim) - min(ylim)) * 0.05; % spacing between brackets
% 
% for i = 1:length(mouseIDs)
%     xH = x_positions(i,1);
%     xM = x_positions(i,2);
%     xL = x_positions(i,3);
% 
%     dH = sess_high_avg{i}';
%     dM = sess_med_avg{i}';
%     dL = sess_low_avg{i}';
% 
%     % Run t-tests
%     [~, p_HM] = ttest2(dH, dM);
%     [~, p_HL] = ttest2(dH, dL);
%     [~, p_ML] = ttest2(dM, dL);
% 
%     % Stagger bracket heights
%     y1 = y_max + bracket_gap;
%     y2 = y_max + bracket_gap*2.5;
%     y3 = y_max + bracket_gap*4;
% 
%     % High vs Med
%     line([xH, xH, xM, xM], [y1, y1+bracket_gap*0.4, y1+bracket_gap*0.4, y1], 'Color', 'k', 'LineWidth', 1);
%     text((xH+xM)/2, y1+bracket_gap*0.5, pToStars(p_HM), 'HorizontalAlignment', 'center', 'FontSize', 14);
% 
%     % Med vs Low
%     line([xM, xM, xL, xL], [y2, y2+bracket_gap*0.4, y2+bracket_gap*0.4, y2], 'Color', 'k', 'LineWidth', 1);
%     text((xM+xL)/2, y2+bracket_gap*0.5, pToStars(p_ML), 'HorizontalAlignment', 'center', 'FontSize', 14);
% 
%     % High vs Low (highest bracket)
%     line([xH, xH, xL, xL], [y3, y3+bracket_gap*0.4, y3+bracket_gap*0.4, y3], 'Color', 'k', 'LineWidth', 1);
%     text((xH+xL)/2, y3+bracket_gap*0.5, pToStars(p_HL), 'HorizontalAlignment', 'center', 'FontSize', 14);
% end

xticks([2, 5, 8, 11, 14, 17, 20, 23]);
xticklabels({"M433","M452","M453","M460","M533","M545","M547","M548"});
ylabel('Mean Avg');
set(gca, 'fontsize', 18)


%% FIGURE TIME : AUC 
figure(2);
hold on;

low_c  = [0.3 0.6 1.0];
high_c = [1.0 0.3 0.3];

count = 0;
for i = 1:length(mouseIDs)
    
    % --- LOW AUC ---
    count = count + 1;
    data_low = sess_low_AUC{i}';
    b_low = boxchart(count*ones(size(data_low)), data_low, ...
        'BoxFaceColor', low_c, 'MarkerStyle', 'none');
    
    % Extend whiskers to min/max
    minVal = min(data_low);
    maxVal = max(data_low);
    line([count-0.2, count+0.2], [minVal, minVal], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [maxVal, maxVal], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [minVal, prctile(data_low, 25)], 'Color', 'k', 'LineStyle', '--');
    line([count, count], [maxVal, prctile(data_low, 75)], 'Color', 'k', 'LineStyle', '--');

    % --- HIGH AUC ---
    count = count + 1;
    data_high = sess_high_AUC{i}';
    b_high = boxchart(count*ones(size(data_high)), data_high, ...
        'BoxFaceColor', high_c, 'MarkerStyle', 'none');

    % Extend whiskers to min/max
    minVal = min(data_high);
    maxVal = max(data_high);
    line([count-0.2, count+0.2], [minVal, minVal], 'Color', 'k', 'LineWidth', 1.5);
    line([count-0.2, count+0.2], [maxVal, maxVal], 'Color', 'k', 'LineWidth', 1.5);
    line([count, count], [minVal, prctile(data_high, 25)], 'Color', 'k', 'LineStyle', '-');
    line([count, count], [maxVal, prctile(data_high, 75)], 'Color', 'k', 'LineStyle', '-');

end

% X ticks centered between each low/high pair
xticks([1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5, 15.5]);
xticklabels(mouseIDs);
ylabel('Mean AUC');

%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
set(gca, 'fontsize', 18);
%cd ('C:\Users\mimia\OneDrive\Desktop\figures')
%exportgraphics(gcf, 'Amp_mice_recolored2.png', 'ContentType','vector');  % Export as PDF
%hold off;

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
%%
function stars = getStars(p)
    if p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = 'ns';
    end
end