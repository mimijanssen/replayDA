# replayDA
Fiber Photometry Analysis Code 
Inspired by Thomas Akam's code (https://github.com/ThomasAkam/photometry_preprocessing/blob/master/Photometry%20data%20preprocessing.ipynb)

Step one: Run Preprocessing Code (preprocessing.m)
input: 1) constant wavelength neuralynx fiber photometry data & 2) ExpKeys
output: preprocessing figures and data saved as your "file_nameprocessed"
notes: 
-  rename the desired path and directory in the code (folder that has the fiber data (mine is labeled CSC30.ncs)
-  rename the recording files name in the code
-  this code decimates to 1000 Hz.

Step two: Run Plotting Linear Track Code (plotting_prob_lineartrack.m)
input: 1) constant wavelength neuralynx fiber photometry data, 2) ExpKeys preprocessed data, 3) preprocessed data, 4) pseudo probability data (the volumes of each reward delivered)
outputs: plots individual trials over each other, plots averaged signal for each reward volume, and saves this averaged data for additional plots averaging session data ("file_name_4avgsess.mat") 
notes: 
- rename paths and directories
- rename filename

Step three: Run Plotting Average Session Data (plotting_avg_session_data.m)
input: 1) file_name_4avgsess.mat aquired from step two
output: plots from averaging sessions 
notes: 
- requires shaded error bar: https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar

  

