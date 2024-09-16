%% Load 
posdata = LoadPos([]); 

%% Plot x against y : quick check
plot(posdata.data(1,:),posdata.data(2,:),'.');

%% More general plot x against y if you don't know which dimension is which 
plot(getd(posdata,'x'),getd(posdata,'y'),'.');

%% Video data
% load video data (make sure the VT1.zip file is unzipped first and now present in MATLAB's working folder!)
[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, [] );
% The 1s say load everything 
% help Nlx2MatVT for the details --> MEX file 
whos 
plot(X,Y);
% The lines mean that neuralynx is missing data. (samples on which no
% position data could be acquired). 

%% Plot data again without missing data 
% Define a variable keep_idx that contains the indices of those samples
% which you want to keep (that are not (0,0)) 
% find(~A) finds the index of a 0 if find(~X) == find(~Y) then set to 0 
keep_idx = [find(X)',find(Y)'];
%coordinates for variables to keep. find(X) returns linear indicies of nonzero elements 
% This just deletes all 0s ... so it is wrong

%% plot video data -- use a new cell so that you can rerun this without also reloading the data
fh = figure; %set(fh,'Color',[0 0 0]);
plot(X(keep_idx),Y(keep_idx),'.','Color',[0 0 0],'MarkerSize',1);
xlabel('x position (pixels)','Color','k')
ylabel('y position (pixels)','Color','k')
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
exportgraphics(gca,'plot.pdf','BackgroundColor','none')


%% Timestamps by plotting X data as a function of time 
plot(Timestamps(keep_idx),X(keep_idx),'.k','MarkerSize',3)
xlabel('t (us)','Color','k')
ylabel('x position (pixels)','Color','k')
box off;
set(gca,'FontSize',24);

%% Plotting X as a function of s not ms us * 1s/1000000 us

%find the index of keepidx and only keep those timesteps 
Timestamps_sec = (Timestamps/1000000);
Timestamps_init = [Timestamps_sec - Timestamps(1)/1000000 ]; 
plot(Timestamps_init(keep_idx),X(keep_idx),'.b','MarkerSize',3)
xlabel('t (s)','Color','k')
ylabel('x position (pixels)','Color','k')
box off;
set(gca,'FontSize',24);
exportgraphics(gca,'plot.pdf','BackgroundColor','none')


%% Video Tracking sampling rate
% delete gaps in data? 
X_sec_diff = diff(X/1000000);
datapoints = nonzeros(X_sec_diff); %43851 samples per 131898 seconds 0.33 samples per second? It is actually 30 Hz so figure out what went wrong here 
plot(Timestamps(keep_idx),X_sec(keep_idx),'.r','MarkerSize',3)
box off;
set(gca,'FontSize',24);