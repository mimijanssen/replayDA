cd 'C:\Users\mimia\Documents'
data = readcell('DG_Spike.xlsx');  % or readtable
col = data(:,1);  % assuming it's all in column 1

splitData = cellfun(@(x) strsplit(x, ','), col, 'UniformOutput', false);

timestamps = cellfun(@(x) str2double(x{1}), splitData);
values     = cellfun(@(x) str2double(x{2}), splitData);
labels     = cellfun(@(x) string(x{3}), splitData);

idx.prerecord = labels == "prerecord";
idx.task      = labels == "task";
idx.postrecord= labels == "postrecord";

T = table(timestamps, values, labels);

SWR_evt = load('2023-09-22_M433_recording3_5detectedSWRs.mat');

%%

LoadExpKeys;
cfg = []; cfg.fc = {'CSC11.ncs'};
lfp = LoadCSC(cfg);

%% DG iv struct 

DG = iv; 
DG.tstart = T.timestamps;
DG.tend = T.timestamps;

fg = [];
cfg.display = 'iv';
cfg.mode = 'center';
cfg.fgcol = 'k';
 
PlotTSDfromIV(cfg,DG,lfp);

%%
cfg_temp = []; cfg_temp.getRatings = 0; cfg_temp.load_questionable_cells = cfg.load_questionable_cells; cfg_temp.verbose = cfg.verbose;
%S = LoadSpikes(cfg_temp);
cfg_temp.fc = ExpKeys.goodSpikes;
S = LoadNST(cfg_temp); 
%%
cfg_plot.lfp(1) = lfp; 
%spd_zscored = spd;
%spd_zscored.data = zscore(spd_zbscored.data);
%cfg_plot.lfp(5) = spd_zscored;
cfg_plot.evt = DG;
MultiRaster(cfg_plot,[])

