session = '/mnt/hubel-data-139/perceval/Rat003_20231226/Rat003_20231226.xml';
[filebase,basename] = fileparts(session);
R = regions(session);
us_intervals = readmatrix(fullfile(filebase,'UltraSlowRythm',basename+".slownr"),FileType='text');

% see US rhythm
R.plotFiringRates(regions=70,smooth=100)
event_int = R.eventIntervals('sleepm');
[unit_firing,time] = R.unitFiringRates('all',70,window=0.15,smooth=25,zscore=true); % CAN TRY WITH zscore AS WELL
unit_firing = Restrict([time,unit_firing],event_int);
time = unit_firing(:,1);
unit_firing = unit_firing(:,2:end);

%% SHOULD RESTRICT TO SWS AND ONLY KEEP NEURONS WHICH CONTRIBUTE TO ISR!!

% see raw data
figure
plot(time,unit_firing(:,1))

% spikes = R.spikes('all',70);
% start = 500; stop = 5500; step = 0.15;
% raster = spikeRaster(spikes,start,stop,'step',step,'mode','int');
% time = start : step : stop-step/3;
% 
% % see raw data
% figure
% plot_data = smoothdata(raster,2,'gaussian',20);
% plot(time,plot_data(1,:))

% on smoothed firing rates
[reduction, umap, clusterIdentifiers, extras]=run_umap(unit_firing,'min_dist',0.6,'n_neighbors',50,'n_components',3); % 'metric','cosine',

figure
hold on
[~,ind] = Restrict(time,us_intervals);
us_ind = false(size(reduction,1),1);
us_ind(ind) = true;
subs_factor = 3;
plot_data = reduction(1:subs_factor:end,:);
plot_ind = us_ind(1:subs_factor:end);
scatter3(plot_data(plot_ind,1),plot_data(plot_ind,2),plot_data(plot_ind,3),[],myColors(1));
scatter3(plot_data(~plot_ind,1),plot_data(~plot_ind,2),plot_data(~plot_ind,3),[],myColors(2));
clear subs_factor us_ind plot_data plot_ind

figure
plot_data = reduction(1:10:end,:);
scatter3(plot_data(:,1),plot_data(:,2),plot_data(:,3),[],time(1:10:end)); % linspace(start,stop,size(plot_data,1)))
colorbar