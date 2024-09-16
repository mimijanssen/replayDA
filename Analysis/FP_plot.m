%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT PHOTOMETRY SIGNAL NON-ISOSBESTIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% From NLX (CW)
Fs = FP_data.final.Fs;
time = FP.tvec;
FP = FP.data;

%% From NLX (FM)
rawFs = FP_data.acq.Fs;
excDemod = excDemod.data;
isoDemod = isoDemod.data;
excRef = csc_ref470.data';
isoRef = csc_ref405.data';

%% From WS (CW)
Fs = data.gen.Fs;
time = data.final.time;
FP = data.final.FP{1};

%% From WS (FM)
rawFs = data.gen.acqFs;
rawTime = data.acq.time;

Fs = data.gen.Fs;
time = data.final.time;
FP = data.final.FP{1};

excDemod = data.final.nbFP{1};
isoDemod = data.final.iso{1};
excRef = data.acq.refSig{1};
isoRef = data.acq.refSig{2};

%%
sessionTitle = 'CW_';
time_ranges = [110, 200, 1100]; % 500 previously 500s (too long for this session) [10, 150,500]
to_remove = 100; %added
% try time ranges 110, 200, 1100 and subtract 100
% previously 10, 100,1000
for t_i = 1:length(time_ranges)
    t_range = 1:Fs*time_ranges(t_i);
    subplot(3, 1, t_i);
    %nothing here, but im creating x_hund
    x_hund = (1:Fs*to_remove); %added
    %FP_rem = FP(x_hund);
    %FP_rem = [];
    plot(time(length(x_hund):length(t_range)), FP(length(x_hund):length(t_range)), 'Color', [0 0.5 0])
    %plot(time(t_range), FP(t_range), 'Color', [0 0.5 0])
    title([sessionTitle, num2str(time_ranges(t_i))], 'Interpreter','none')
    ylabel('Fluorescence (dF/F)'); xlabel('Time (s)');
end

%%
sessionTitle = 'FM_';
time_ranges = [10, 30, 60];

for t_i = 1:length(time_ranges)
    t_range = 1:Fs*time_ranges(t_i);
    FP_demods = {FP, excDemod, isoDemod};
    FP_titles = {'signal', 'excitation', 'isosbestic'};
    figure;
    for m_i = 1:length(FP_demods)
        subplot(3, 1, m_i);
        plot(time(t_range), FP_demods{m_i}(t_range), 'Color', [0 0.5 0])
        title([sessionTitle, FP_titles{m_i}], 'Interpreter','none')
        ylabel('Fluorescence (dF/F)'); xlabel('Time (s)');
    end
end

%%
sessionTitle = 'FM_';
time_ranges = [0.1];

for t_i = 1:length(time_ranges)
    t_range = 1:rawFs*time_ranges(t_i);
    FP_demods = {excRef, isoRef};
    FP_titles = {'Excitation', 'Isosbestic'};
    figure;
    for m_i = 1:length(FP_demods)
        subplot(2, 1, m_i);
        plot(rawTime(t_range), FP_demods{m_i}(t_range), 'Color', [0 0.5 0])
        title([sessionTitle, FP_titles{m_i}], 'Interpreter','none')
        ylabel('V'); xlabel('Time (0.1 s)');
    end
end

%% PETH
plot(left_end_peth.tvec, left_end_peth.data); 
hold on;
plot(right_end_peth.tvec, right_end_peth.data); % some sessions will not plot because data is NaN 
legend('Left end', 'Right end')