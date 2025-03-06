load('M578_2025_01_14_recording1-manualIV.mat')
ncfs = SWRfreak([], evt, CSC);
figure; plot(ncfs.freqs1);
LoadMetadata();
metadata.SWRfreqs = ncfs;
save(FindFile('*Metadata.mat'), 'metadata')