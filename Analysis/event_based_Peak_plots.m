%% Peak in 1 second after SWR from stored traces
allTables2 = allTables;

allTables2.PostSWRPeak = NaN(height(allTables2),1);

for i = 1:height(allTables2)

    if allTables2.PrePost(i) == 1
        signal = allTables2.TwosPreProc{i}.signal;

    elseif allTables2.PrePost(i) == 2
        signal = allTables2.TwosPostProc{i}.signal;

    else
        continue
    end

    % 2001-point trace:
    % 1:1000      = 1 s before SWR
    % 1001        = SWR center
    % 1002:2001   = 1 s after SWR

    allTables2.PostSWRPeak(i) = max(signal(1002:end));

end


%% BEFORE SWR PEAK VALUES 

%% Peak in 1 second after SWR from stored traces
allTables2.PreSWRPeak = NaN(height(allTables2),1);

for i = 1:height(allTables2)

    if allTables2.PrePost(i) == 1
        signal = allTables2.TwosPreProc{i}.signal;

    elseif allTables2.PrePost(i) == 2
        signal = allTables2.TwosPostProc{i}.signal;

    else
        continue
    end

    % 2001-point trace:
    % 1:1000      = 1 s before SWR
    % 1001        = SWR center
    % 1002:2001   = 1 s after SWR

    allTables2.PreSWRPeak(i) = max(signal(1:1001));

end

%%

% Separate peak values by sleep epoch
pre_peaks  = allTables2.PostSWRPeak(allTables.PrePost == 1)-allTables2.PreSWRPeak(allTables.PrePost == 1);
post_peaks = allTables2.PostSWRPeak(allTables.PrePost == 2)-allTables2.PreSWRPeak(allTables.PrePost == 2);

figure

subplot(1,2,1)
histogram(pre_peaks,50,'FaceColor',[0.2 0.4 0.8],...
    'FaceAlpha',0.5)
hold on
histogram(post_peaks,50,'FaceColor',[0.8 0.2 0.2],...
    'FaceAlpha',0.5)
xlabel('SWR-DA (z-score)')
ylabel('SWR Count')
legend('Pre-track rest','Post-track rest')
title('Peak distributions')

subplot(1,2,2)
boxplot([pre_peaks; post_peaks], ...
    [ones(size(pre_peaks)); 2*ones(size(post_peaks))])
set(gca,'XTickLabel',{'Pre','Post'})
ylabel('SWR-DA (z-score)')
title('Peak comparison')
%%
median(pre_peaks)
mean(pre_peaks)
std(pre_peaks)
disp('------')
median(post_peaks)
mean(post_peaks)
std(post_peaks)
%%
figure

scatter(allTables.OnesAfterPeak,...
        allTables2.PostSWRPeak,...
        10,'filled')

hold on

lims = [
    min([allTables.OnesAfterPeak; allTables2.PostSWRPeak]), ...
    max([allTables.OnesAfterPeak; allTables2.PostSWRPeak])
];

plot(lims,lims,'r--','LineWidth',2)

xlabel('Stored OnesAfterPeak')
ylabel('Recomputed PostSWRPeak')
title('Peak verification')
axis square


%%

% Separate peak values by sleep epoch
pre_peaks  = allTables2.PostSWRPeak(allTables.PrePost == 1)-allTables2.PreSWRPeak(allTables.PrePost == 1);
post_peaks = allTables2.PostSWRPeak(allTables.PrePost == 2)-allTables2.PreSWRPeak(allTables.PrePost == 2);

figure

subplot(1,2,1)
histogram(pre_peaks,50,'FaceColor',[0.2 0.4 0.8],...
    'FaceAlpha',0.5)
hold on
histogram(post_peaks,50,'FaceColor',[0.8 0.2 0.2],...
    'FaceAlpha',0.5)
xlabel('SWR-DA (z-score)')
ylabel('SWR Count')
legend('Pre-track rest','Post-track rest')
title('Peak Baseline distributions')

subplot(1,2,2)
boxplot([pre_peaks; post_peaks], ...
    [ones(size(pre_peaks)); 2*ones(size(post_peaks))])
set(gca,'XTickLabel',{'Pre','Post'})
ylabel('SWR-DA (z-score)')
title('Peak Baseline comparison')

%%
median(pre_peaks)
mean(pre_peaks)
std(pre_peaks)
disp('------')
median(post_peaks)
mean(post_peaks)
std(post_peaks)