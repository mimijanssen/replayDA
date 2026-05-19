load('M648_2026_02_20_recording2-manualIV.mat')
ncfs = SWRfreak([], evt, CSC);
figure; plot(ncfs.freqs1);
LoadMetadata();
metadata.SWRfreqs = ncfs;
save(FindFile('*Metadata.mat'), 'metadata')