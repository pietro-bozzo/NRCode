% Visualize ultra-slow avalanches in Nucleus Reuniens
% requires FMAToolbox, Regions

% load session
session = '/mnt/hubel-data-131/perceval/Rat003_20231222/Rat003_20231222.xml';
[filebase,basename] = fileparts(session);
R = regions(session,regions=70);

% see spikes
start = 2000; % s
stop = 3000; % s
R.plotSpikeRaster(start,stop);

% load ultra-slow avalanches
us_intervals = readmatrix(fullfile(filebase,'UltraSlowRythm',basename+".slownr"),FileType='text');
us_avals = readmatrix(fullfile(filebase,'UltraSlowRythm',basename+".slowavalnr"),FileType='text');
us_avals = us_avals(:,1:2);

% see population activity
R.plotFiringRates(start,stop,step=5,smooth=45);
PlotIntervals(us_intervals,'legend','US')
PlotIntervals(us_avals,'color',[0.8,0.2,0.2],'legend','avalanches')

% TO DO
% 1. use GetWidebandData (from FMAT) to load raw data from 1 or few channels
% 2. filter it with correct high-pass to detect spikes (ask Pascale about the frequency and how to count spikes)
% 3. see if resulting signal is higher during US avalanches