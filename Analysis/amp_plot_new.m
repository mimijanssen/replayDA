%% Initialize
clear; clc;
addpath('C:\Users\mimia\OneDrive\Documents\Toolboxes\raacampbell-shadedErrorBar'); % path to shaded error bar 

% Define folder paths
salineFolder = 'D:\Saline';
amphetamineFolder = 'D:\Amphetamine';

% Get list of mouse folders
salineMice = dir(fullfile(salineFolder, 'M*'));
amphetamineMice = dir(fullfile(amphetamineFolder, 'M*'));

%% 1) save all session data (-4 min to 11 min from injection time) 
% --> and find longest injection time 

% initialize variables
% 25 min = 1500 s ; 1500 s * 1000 samples/s = 1500000 samples
% double check what I expect ^
salineMice = dir(fullfile(salineFolder, 'M*'));
ampMice = dir(fullfile(amphetamineFolder, 'M*'));
referenceMouse = 'M433'; % Update with your choice
referenceFile = 'M433_2023_10_18processed'; % Update with your choice
refData = load(fullfile(salineFolder, referenceMouse, referenceFile));
minutes = 15; 

r = 100; % decimate by this amount. 


samples = minutes*60*refData.cfg.hdr{1,1}.SamplingFrequency/r;
samples_before = 4*60*refData.cfg.hdr{1,1}.SamplingFrequency/r;
samples_after = 11*60*refData.cfg.hdr{1,1}.SamplingFrequency/r;

Saline_sess = {}; 
Amp_sess = {}; 
Saline_time = {};
Amp_time = {};
Sal_inj_time = {}; 
Amp_inj_time = {}; 

% find injection time : avg between drug admin and admin done 
% --- Process saline group ---
for i = 1:length(salineMice) % iterate through the mice
    mouseFolder = fullfile(salineFolder, salineMice(i).name);
    dataFiles = dir(fullfile(mouseFolder, '*processed*')); % list of sessions for that mouse
    mouseName = salineMice(i).name;
    numSessions = length(dataFiles);
    sessionData = zeros(numSessions, samples); % assuming all sessions have 1500000 samples
    sess_salineTime = zeros(numSessions, samples);
    sess_injTime = zeros(numSessions,1);

    for j = 1:numSessions % now open each of those sessions 

        data = load(fullfile(mouseFolder, dataFiles(j).name)); % data.tvec; data.F_zscored_base for detrended and z_scored data
        keys = (fullfile(mouseFolder, strrep(dataFiles(j).name, 'processed.mat', '')));
        cd (keys);
        LoadExpKeys();
        cfg_evt = [];
        cfg_evt.eventList = ExpKeys.eventList;
        cfg_evt.eventLabel = ExpKeys.eventLabel;
        evt = LoadEvents(cfg_evt); % 7 is injection start ; 8 is injection end

        % ok, I did not save the original tvec. so now I need to extract
        % these injection times and then Load the original data and
        % subtract the original tvec. Then I will take the average. 
        inj_init = evt.t{1,7};
        inj_end = evt.t{1,8}; 

        raw_data_path = fullfile(strrep(fileparts(mouseFolder), 'Saline', salineMice(i).name), strrep(dataFiles(j).name, 'processed.mat', ''));
        cd(raw_data_path)
        cfg.fc = {'CSC30.ncs'};
        csc_photo = LoadCSC(cfg);

        inj_init = inj_init - csc_photo.tvec(1); 
        inj_end = inj_end - csc_photo.tvec(1); 
        
        inj_time = (inj_init + inj_end)/2; 

        deci_time = decimate(data.tvec,r); 
        index_time = nearest_idx3(inj_time, deci_time); % values, look-up values

        init_bound = index_time - samples_before; 
        end_bound = index_time + samples_after; 
        %check = end_bound - init_bound % should be number of samples 
        
        deci_data = decimate(data.dF_base, r); 

        % take -10 to 15 min of that and put it in a structure 
        sessionData(j, :) = deci_data(init_bound:end_bound-1); %data.F_zscored_base
        sess_salineTime(j,:) = deci_time(init_bound:end_bound-1); % data.tvec

        % save injection time
        sess_injTime(j,:) = inj_time; 
    end
    
    Saline_sess.(mouseName) = sessionData;
    Saline_time.(mouseName) = sess_salineTime;
    Sal_inj_time.(mouseName) = sess_injTime; 

end

% --- Process amphetamine group ---
for i = 1:length(ampMice) % iterate through the mice
    mouseFolder = fullfile(amphetamineFolder, ampMice(i).name);
    dataFiles = dir(fullfile(mouseFolder, '*processed*')); % list of sessions for that mouse
    mouseName = ampMice(i).name;
    numSessions = length(dataFiles);
    sessionData = zeros(numSessions, samples); % assuming all sessions have 1500000 samples
    sess_ampTime = zeros(numSessions, samples);
    sess_injTime = zeros(numSessions,1);

    for j = 1:numSessions % now open each of those sessions 

        data = load(fullfile(mouseFolder, dataFiles(j).name)); % data.tvec; data.F_zscored_base for detrended and z_scored data
        keys = (fullfile(mouseFolder, strrep(dataFiles(j).name, 'processed.mat', '')));
        cd (keys);
        LoadExpKeys();
        cfg_evt = [];
        cfg_evt.eventList = ExpKeys.eventList;
        cfg_evt.eventLabel = ExpKeys.eventLabel;
        evt = LoadEvents(cfg_evt); % 7 is injection start ; 8 is injection end

        % ok, I did not save the original tvec. so now I need to extract
        % these injection times and then Load the original data and
        % subtract the original tvec. Then I will take the average. 
        inj_init = evt.t{1,7}(1);
        inj_end = evt.t{1,8}(1); 

        raw_data_path = fullfile(strrep(fileparts(mouseFolder), 'Amphetamine', ampMice(i).name), strrep(dataFiles(j).name, 'processed.mat', ''));
        cd(raw_data_path)
        cfg.fc = {'CSC30.ncs'};
        csc_photo = LoadCSC(cfg);

        inj_init = inj_init - csc_photo.tvec(1); 
        inj_end = inj_end - csc_photo.tvec(1); 
        
        inj_time = (inj_init + inj_end)/2; 

        deci_time = decimate(data.tvec,r); 
        index_time = nearest_idx3(inj_time, deci_time); % values, look-up values

        init_bound = index_time - samples_before; 
        end_bound = index_time + samples_after; 
        %check = end_bound - init_bound % should be number of samples 

        deci_data = decimate(data.dF_base, r); 
        
        % take -10 to 15 min of that and put it in a structure 
        sessionData(j, :) = deci_data(init_bound:end_bound-1); % dF_base 
        sess_ampTime(j,:) = deci_time(init_bound:end_bound-1); 

        % save injection time
        sess_injTime(j,:) = inj_time; 
    end
    
    Amp_sess.(mouseName) = sessionData;
    Amp_time.(mouseName) = sess_ampTime;
    Amp_inj_time.(mouseName) = sess_injTime; 

end

% find max injection time for both saline and amphetamine group 

% find max mouse for saline sessions: 
%[M,I] = max([max(Sal_inj_time.M433),max(Sal_inj_time.M453),max(Sal_inj_time.M460),max(Sal_inj_time.M533),max(Sal_inj_time.M534),max(Sal_inj_time.M545),max(Sal_inj_time.M547),max(Sal_inj_time.M548)]); 
% 808.3070 ; mouse 6 
%[M2,I2] = max([max(Amp_inj_time.M433),max(Amp_inj_time.M453),max(Amp_inj_time.M460),max(Amp_inj_time.M533),max(Amp_inj_time.M534),max(Amp_inj_time.M545),max(Amp_inj_time.M547),max(Amp_inj_time.M548)]); 
% 786 ; mouse 6 
% Saline session of M545 is the longest injection time 

% find session: 
%[M3,I3] = max(Sal_inj_time.M545); 
% session 1. 

%% 2) save all average mouse data (-10 to 14 min from max injection time)
% average each field and save in an overall matrix 

% Get the names of the fields in the structure
fields_saline = fieldnames(Saline_sess);
fields_amp = fieldnames(Amp_sess); 

% Initialize an empty matrix to store the averages
avg_sal = [];
avg_amp = []; 

% Loop through each field
for iter = 1:length(fields_saline) % iterate through each mouse 
    % Get the data for the current field
    data_saline_sess = Saline_sess.(fields_saline{iter});
    % time_saline_sess = Saline_time.(fields_saline{iter});  
    % inj_saline_sess = Sal_inj_time.(fields_saline{iter});
    % Calculate the average of the current field
    field_average = mean(data_saline_sess,1);
    % Append the average to the matrix as a new row
    avg_sal = [avg_sal; field_average];
end

for iter = 1:length(fields_amp)
    data_amp_sess = Amp_sess.(fields_amp{iter});
    field_average = mean(data_amp_sess,1);
    avg_amp = [avg_amp; field_average];
end


%% 3) find mean and SEM for plotting 
% average and SEM 
meanSal = mean(avg_sal,1); 
meanAmp = mean(avg_amp,1); 
SEMSal = std(avg_sal)/sqrt(size(avg_sal, 1));
SEMAmp = std(avg_amp)/sqrt(size(avg_amp, 1)); 

%% 4) plotting 
time_for_plotting = Saline_time.M433(1,:) - Saline_time.M433(1,1); % intialized. 
% This means that the injection variable should be + samples before + 1
% this is for 1,000 Hz. How would we plot minutes? can relabel the x axis. 
time_for_plotting = time_for_plotting/(60); % in minutes now? 
%injection_time_index = 600001;

%% colors
high_c = [238,119,51]./255 ; %[255, 165,0]./255; % rgb(204, 204, 255) periwinkle
med_c = [0, 153, 136]./255; %[78,178,101]./255;%[104,187,225]./255; % rgb(167, 199, 231) rgb(255, 165, 0)  blue: rgb(104,187,227)


%%
figure (1);

% Saline
shadedErrorBar(time_for_plotting, meanSal, SEMSal, 'lineprops',{'-','color',med_c,'MarkerFaceColor',med_c});
hold on;

% Amphetamine
shadedErrorBar(time_for_plotting, meanAmp, SEMAmp, 'lineprops',{'-','color',high_c,'MarkerFaceColor',high_c});
xline(time_for_plotting(samples_before+1), '--k', 'LineWidth', 1);

xlim([0 15])
xticks([0 2 4 6 8 10 12 14]) % 4 is about the time where injection is. 
xticklabels({'-4','','0','','4','','','10'})

% Labels and legend
xlabel('Time from injection (min)');
ylabel('Mean [DA] (dF/F)');
legend({'Saline', 'Amphetamine'}, 'Location', 'best');
%title('Average Fiber Photometry Traces');
set(gca,'fontsize', 18)
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');

set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\OneDrive\Desktop\figures')
exportgraphics(gcf, 'Amp_recolored.png', 'ContentType','vector');  % Export as PDF

hold off;


%% statistics dF 
% change in dF : for each mouse, 1) find the max value in the post injection period , 2)
% find the min value in the post injection period, 3) take the difference,
% 4) enter that into a t-test 

% 1) find the max value in the post injection period : all samples after
% samples_before 

max_saline = max(avg_sal(:,samples_before+1:end),[],2); 
max_amp =  max(avg_amp(:,samples_before+1:end),[],2); 

% 2) find the min value in the post injection period: all samples after
% samples before
min_saline = min(avg_sal(:,samples_before+1:end),[],2); 
min_amp =  min(avg_amp(:,samples_before+1:end),[],2); 

% 3) take the difference 
dF_saline = max_saline - min_saline;
dF_amp = max_amp - min_amp; 

% 4) t-test
[dF_h, dF_p,~, dF_tstats] = ttest(dF_saline, dF_amp); 

%% separate out for each mouse 
% Amp_sess.MXXX (each row is a new session ex. 2x90000)

% 1) combine data into a column vector for boxchart input 
% first group of rows is each session for one mouse for saline 
% second group of rows is each session for one mouse for amp
% repeat for all mice 

% along with this, you need a categories variable that is the same size, 
% but a different number for each box. 

combinedData = [];
categories = [];
count = 0; 

% Loop through each field
for iter = 1:length(fields_saline) % iterate through each mouse -- same for amphetamine...
   
    % get session data 
    count = count + 1; 
    data_sal_sess = Saline_sess.(fields_saline{iter}); % 2x9000
    % find the min and max values 
    sal_sess_dF = max(data_sal_sess(:,samples_before+1:end),[],2) - min(data_sal_sess(:,samples_before+1:end),[],2); 
    cat_sal_sess = count*ones(1,numel(sal_sess_dF));

    count = count + 1;
    data_amp_sess = Amp_sess.(fields_saline{iter}); 
    amp_sess_dF = max(data_amp_sess(:,samples_before+1:end),[],2) - min(data_amp_sess(:,samples_before+1:end),[],2); 
    cat_amp_sess = count*ones(1,numel(amp_sess_dF)); 

    combinedData = [combinedData; sal_sess_dF; amp_sess_dF];
    categories = [categories, cat_sal_sess, cat_amp_sess]; 

end

% as one vector? that's crazy
combinedData2 = reshape(combinedData,1,[]);  % convert matrix to row vector

%%
% 2) FIGURE TIME 
figure(2);
%h = boxchart(categories, combinedData2);  % Use categories for x grouping
hold on;

% loop for boxchart color 
boxchart_size = [];
count = 0; 
% Loop through each field
for iter = 1:length(fields_saline) % iterate through each mouse -- same for amphetamine...
    count = count + 1; 
    data_sal_sess = Saline_sess.(fields_saline{iter}); % 2x9000
    sal_sess_dF = max(data_sal_sess(:,samples_before+1:end),[],2) - min(data_sal_sess(:,samples_before+1:end),[],2); 
    boxchart(count*ones(1,numel(sal_sess_dF)), reshape(sal_sess_dF,1,[]), 'BoxFaceColor', med_c);

    count = count + 1;
    data_amp_sess = Amp_sess.(fields_saline{iter}); 
    amp_sess_dF = max(data_amp_sess(:,samples_before+1:end),[],2) - min(data_amp_sess(:,samples_before+1:end),[],2); 
    boxchart(count*ones(1,numel(amp_sess_dF)), reshape(amp_sess_dF,1,[]), 'BoxFaceColor', high_c);
end

% Set x-tick labels and axis properties
 xticks([1.5,3.5,5.5,7.5,9.5,11.5,13.5,15.5]);
xticklabels({"M433", "M453", "M460", "M533","M534","M545","M547","M548"});
ylabel("Mean \Delta [DA] (dF/F)");
%  Set figure properties
 set(gcf, 'color', 'none');
 set(gca, 'color', 'none');
 set(gca, 'fontsize', 18);

set(gcf, 'renderer', 'painters');
cd ('C:\Users\mimia\OneDrive\Desktop\figures')
exportgraphics(gcf, 'Amp_mice_recolored.png', 'ContentType','vector');  % Export as PDF
hold off;

%% Stats for each mice! 
% because we had a significant group mean difference between amphetamine and salien groups (shown by the t-test)
% . I want to break it out by mice to see if all mice have a significant
% difference in session. 

% I stored the data in a really odd way for the boxchart function to work.
% So let's re-run that in a more intuitive way for the t-test, where each
% structure is a mouse and each row is a session. 

amp_dF_mice = {};
sal_dF_mice = {}; 
% Loop through each field
for iter = 1:length(fields_saline) % iterate through each mouse -- same for amphetamine...
   
    % get session data 
    data_sal_sess = Saline_sess.(fields_saline{iter}); % 2x9000
    name = fields_saline{iter}; 
    % find the min and max values 
    sal_sess_dF = max(data_sal_sess(:,samples_before+1:end),[],2) - min(data_sal_sess(:,samples_before+1:end),[],2); % 2x1
    sal_dF_mice.(name) = sal_sess_dF; 
    
    data_amp_sess = Amp_sess.(fields_saline{iter}); 
    amp_sess_dF = max(data_amp_sess(:,samples_before+1:end),[],2) - min(data_amp_sess(:,samples_before+1:end),[],2); 
    amp_dF_mice.(name) = amp_sess_dF;

end

% Loop and display t-tests in command window 
for iter = 1:length(fields_saline) % iterate through each mouse -- same for amphetamine...
     name = fields_saline{iter}; 
     disp(name);
     [dF_m_h, dF_m_p,~, dF_m_tstats] = ttest2(sal_dF_mice.(name), amp_dF_mice.(name)); % paired needs to be the same size 
     disp(dF_m_h);
     disp(dF_m_p);
     disp(dF_m_tstats);
end


%% Stats but AUC... 
amp_auc_mice = {};
sal_auc_mice = {}; 
% Loop through each field
for iter = 1:length(fields_saline) % iterate through each mouse -- same for amphetamine...
   
    % get session data 
    data_sal_sess = Saline_sess.(fields_saline{iter}); % 2x9000
    name = fields_saline{iter}; 
    % find the min and max values 
    sal_sess_auc = trapz(samples_before+1:90000,data_sal_sess(:,samples_before+1:90000),2); % 2x1
    % trapz(x1, mouse_fiber_post(i_post,x1));
    sal_auc_mice.(name) = sal_sess_auc; 
    
    data_amp_sess = Amp_sess.(fields_saline{iter}); 
    amp_sess_auc = trapz(samples_before+1:90000,data_amp_sess(:,samples_before+1:90000),2); % 2x1
    amp_auc_mice.(name) = amp_sess_auc;

end

% Loop and display t-tests in command window 
for iter = 1:length(fields_saline) % iterate through each mouse -- same for amphetamine...
     name = fields_saline{iter}; 
     disp(name);
     [dF_m_h, dF_m_p,~, dF_m_tstats] = ttest2(sal_auc_mice.(name), amp_auc_mice.(name)); % paired needs to be the same size 
     disp(dF_m_h);
     disp(dF_m_p);
     disp(dF_m_tstats);
end

%% TRASH CODE 
% Injection window shading...
%injectStart = refData.tvec(refAlignIdx);
%injectEnd = refData.tvec(refAlignIdx + 500); % Example: adjust based on injection duration
%patch([injectStart, injectEnd, injectEnd, injectStart], ...
    %  [min(meanSaline), min(meanSaline), max(meanAmphetamine), max(meanAmphetamine)], ...
    %  [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Dashed line for injection onset
% xline(injection_time_index, '--k', 'LineWidth', 1);
% 
% 
% x_ticks = get(gca, 'XTick');  % Get current x-tick positions
% x_ticks_in_minutes = (x_ticks - injection_time_index) / 60*1000;  % Convert to minutes
% 
% % Set the new x-tick positions and labels
% xticks(x_ticks);  % Set x-tick positions
% xticklabels(arrayfun(@(x) sprintf('%d', x), x_ticks_in_minutes, 'UniformOutput', false));  % Set x-tick labels


% additional code for if I want to attempt adding individual session
% points. 

% Set x-axis positions for each group with consistent jitter
% jitterAmount = 0.05;

% x_jitter_early_before = 1 + randn(1, numel(dF1_pre_early)) * jitterAmount;
% x_jitter_early_after = 2 + randn(1, numel(dF2_pre_early)) * jitterAmount;
% x_jitter_late_before = 3 + randn(1, numel(dF1_pre_late)) * jitterAmount;
% x_jitter_late_after = 4 + randn(1, numel(dF2_pre_late)) * jitterAmount;
% 
% % Draw connecting lines for each subject between before and after within early and late sessions
% for i = 1:numel(dF1_pre_early)
%     plot([x_jitter_early_before(i), x_jitter_early_after(i)], [dF1_pre_early(i), dF2_pre_early(i)], ...
%         'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey lines for early
% end
% 
% for i = 1:numel(dF1_pre_late)
%     plot([x_jitter_late_before(i), x_jitter_late_after(i)], [dF1_pre_late(i), dF2_pre_late(i)], ...
%         'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % Light grey lines for late
% end
% 
% % Scatter points for each group
% scatter(x_jitter_early_before, dF1_pre_early, 60, 'MarkerFaceColor', lightp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);
% scatter(x_jitter_early_after, dF2_pre_early, 60, 'MarkerFaceColor', lightp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);
% scatter(x_jitter_late_before, dF1_pre_late, 60, 'MarkerFaceColor', darkerp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);
% scatter(x_jitter_late_after, dF2_pre_late, 60, 'MarkerFaceColor', darkerp, 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Marker', 'o','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8);
% 
