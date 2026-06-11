% Initialize storage
allPre2  = [];
allPost2 = [];

% Loop through each row and collect signals
for i = 1:height(ProcAUCTbl)
    if ProcAUCTbl.PrePost(i) == 1  % Pre
        allPre2  = [allPre2;  ProcAUCTbl.AUC(i)];   % ensure row vector before concatenating
        
    elseif ProcAUCTbl.PrePost(i) == 2  % Post
        allPost2 = [allPost2; ProcAUCTbl.AUC(i)];
    end
end

% % Compute mean and SEM
% mean_pre  = mean(allPre2,  1);
% mean_post = mean(allPost2, 1);
% sem_pre   = std(allPre,  0, 1) / sqrt(size(allPre,  1));
% sem_post  = std(allPost, 0, 1) / sqrt(size(allPost, 1));

% Time vector -- adjust to match your actual signal
% tvec = linspace(-1, 1, length(mean_pre));
%% Baseline the AUC 
allTables4 = allTables;

allTables4.PostMinusPreAUC = NaN(height(allTables4),1);
allTables4.PostMinusPreMean = NaN(height(allTables4),1);

for i = 1:height(allTables4)

    if allTables4.PrePost(i) == 1
        signal = allTables4.TwosPreProc{i}.signal;
    elseif allTables4.PrePost(i) == 2
        signal = allTables4.TwosPostProc{i}.signal;
    else
        continue
    end

    preSignal  = signal(1:1001);
    postSignal = signal(1002:end);

    % AUC difference (post - pre)
    allTables4.PostMinusPreAUC(i) = ...
        trapz(postSignal) - trapz(preSignal);

    % Mean difference (post - pre)
    allTables4.PostMinusPreMean(i) = ...
        mean(postSignal) - mean(preSignal);

end


%%

% Initialize storage
allPre2AUC  = [];
allPost2AUC = [];

% Loop through each row and collect signals
for i = 1:height(allTables4)
    if allTables4.PrePost(i) == 1  % Pre
        allPre2AUC  = [allPre2AUC;  allTables4.PostMinusPreAUC(i)];   % ensure row vector before concatenating
    elseif allTables4.PrePost(i) == 2  % Post
        allPost2AUC = [allPost2AUC; allTables4.PostMinusPreAUC(i)];
    end
end

%% baselined mean

allPre2Mean  = [];
allPost2Mean = [];

% Loop through each row and collect signals
for i = 1:height(allTables4)
    if allTables4.PrePost(i) == 1  % Pre
        allPre2Mean  = [allPre2Mean;  allTables4.PostMinusPreMean(i)];   % ensure row vector before concatenating
    elseif allTables4.PrePost(i) == 2  % Post
        allPost2Mean = [allPost2Mean; allTables4.PostMinusPreMean(i)];
    end
end

%%

figure(2); 
boxplot(allPre2,allPost2)

data = [allPre2; allPost2];
group = [repmat({'X'},length(allPre2),1);
         repmat({'Y'},length(allPost2),1)];

boxplot(data, group)
ylabel('Value')
set(gca,'XTickLabel',{'Pre','Post'})
ylabel('SWR-DA (z-score)')
title('Peak comparison')

%%
figure

subplot(1,2,1)
histogram(allPre2AUC,50,'FaceColor',[0.2 0.4 0.8],...
    'FaceAlpha',0.5)
hold on
histogram(allPost2AUC,50,'FaceColor',[0.8 0.2 0.2],...
    'FaceAlpha',0.5)
xlabel('SWR-DA AUC (z-score)')
ylabel('SWR Count')
legend('Pre-track rest','Post-track rest')
title('Baselined AUC distributions')

subplot(1,2,2)
boxplot([allPre2AUC; allPost2AUC], ...
    [ones(size(allPre2AUC)); 2*ones(size(allPost2AUC))])
set(gca,'XTickLabel',{'Pre','Post'})
ylabel('SWR-DA AUC (z-score)')
title('Baselined AUC comparison')

median(allPre2AUC)
mean(allPre2AUC)
std(allPre2AUC)
disp('------')
median(allPost2AUC)
mean(allPost2AUC)
std(allPost2AUC)

%% figure for mean

figure

subplot(1,2,1)
histogram(allPre2Mean,50,'FaceColor',[0.2 0.4 0.8],...
    'FaceAlpha',0.5)
hold on
histogram(allPost2Mean,50,'FaceColor',[0.8 0.2 0.2],...
    'FaceAlpha',0.5)
xlabel('SWR-DA Mean (z-score)')
ylabel('SWR Count')
legend('Pre-track rest','Post-track rest')
title('Baselined Mean distributions')

subplot(1,2,2)
boxplot([allPre2Mean; allPost2Mean], ...
    [ones(size(allPre2Mean)); 2*ones(size(allPost2Mean))])
set(gca,'XTickLabel',{'Pre','Post'})
ylabel('SWR-DA Mean (z-score)')
title('Baselined Mean comparison')

median(allPre2Mean)
mean(allPre2Mean)
std(allPre2Mean)
disp('------')
median(allPost2Mean)
mean(allPost2Mean)
std(allPost2Mean)
