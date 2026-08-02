%% Mouse permutation means 
% load all the mouse files
clear; clc;
cd ('F:\M548\avg_data\RPE_mean')
RPEFolder = 'F:\M548\avg_data\RPE_mean';
datafiles= dir(fullfile(RPEFolder, 'M*'));
filename = 'M548';

%%
nSessions = length(datafiles);
nPerms = 1000;
perm_all = zeros(nPerms,nSessions);
real_all = zeros(nSessions,1);
for s = 1:nSessions
    D = load(fullfile(RPEFolder,datafiles(s).name));
    perm_all(:,s) = D.perm_stats;
    real_all(s) = D.real_stat;
end

%% average
real_mouse = mean(real_all);
perm_mouse = mean(perm_all,2);

%% P value
p_one = (sum(perm_mouse >= real_mouse)+1)/(nPerms+1);

fprintf('Observed statistic = %.3f\n',real_mouse)
fprintf('Permutation mean   = %.3f\n',mean(perm_mouse))
fprintf('p = %.4f\n',p_one)

%% Plot

figure
histogram(perm_mouse,50,...
    'FaceColor',[0.7 .7 .7],...
    'EdgeColor','none')
hold on
xline(real_mouse,'r','LineWidth',2)
xline(mean(perm_mouse),'k--','LineWidth',2)
xlabel('High - Low response')
ylabel('Count')

title(sprintf('%s  p = %.4f',RPEFolder,p_one))

% %% Take all sessions and mea\n it
% % initalize data 
% n_perms = 1000; 
% n = size(datafiles,1); 
% 
% perm_stats_high = zeros(1000,n); % each session is a column 
% perm_stats_low = zeros(1000,n); 
% 
% real_stats_high = zeros(1,n);
% real_stats_low = zeros(1,n);
% 
% % load all data
% for session = 1:1:length(datafiles)
%     sessName = load(datafiles(session).name);
%     % store data in a large table n x 1000
%     perm_stats_high(:,session) = sessName.perm_stats_high; 
%     perm_stats_low(:,session) = sessName.perm_stats_low; 
%     real_stats_high(1,session) = sessName.sess_mean_high;
%     real_stats_low(1,session) = sessName.sess_mean_low;
% end
% 
% perm_stats_high_avg = mean(perm_stats_high,2); 
% perm_stats_low_avg = mean(perm_stats_low,2); 
% real_stats_high_avg= mean(real_stats_high);
% real_stats_low_avg= mean(real_stats_low);
% 
% %%
% % subtract high and low 
% real_stat = real_stats_high_avg - real_stats_low_avg; 
% perm_stats = perm_stats_high_avg - perm_stats_low_avg;
% 
% % Find perm stat? 
% p_one = (sum(perm_stats>=real_stat)+1)/(n_perms+1); % upper tailed test. 
% 
% fprintf('Observed statistic = %.4f\n',real_stat);
% fprintf('Permutation mean   = %.4f\n',mean(perm_stats));
% fprintf('Upper-Tail p       = %.4f\n',p_one);
% 
% % plot 
% 
% figure(1)
% histogram(perm_stats,200,...
%         'FaceColor',[0.6 0.6 0.6],...
%         'EdgeColor','none');
% hold on 
% 
% xline(real_stat,'r','LineWidth',2);
% 
% xline(mean(perm_stats),'k--','LineWidth',2);
% 
% xlabel('Mean high trial DA - Mean low trial DA')
% ylabel('Count')
% title(sprintf(filename))
% 
% set(gca,'FontSize',12)
% box off

%% Store everything 
%cd 'F:\RPE_means';
%filename = append(filename, "RPE_m_avg.mat");
%save(filename, 'real_stats_high_avg','real_stats_low_avg','perm_stats_high_avg','perm_stats_low_avg');





