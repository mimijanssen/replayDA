load('M453-2024-02-07_ontrack4-manualIV.mat')
ncfs = SWRfreak([], evt, CSC);
figure; plot(ncfs.freqs1);
LoadMetadata();
metadata.SWRfreqs = ncfs;
save(FindFile('*Metadata.mat'), 'metadata')