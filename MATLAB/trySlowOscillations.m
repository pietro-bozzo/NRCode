%session = '/mnt/hubel-data-143/peter/Rat005_20240529/Rat005_20240529.xml';
%session = '/mnt/hubel-data-145/peter/Rat005_20240603/Rat005_20240603.xml';
%session = '/mnt/hubel-data-145/peter/Rat005_20240610/Rat005_20240610.xml';
%session = '/mnt/hubel-data-149/Rat012/Rat012_2025-12-09/Rat012_2025-12-09.xml';
session = '/mnt/hubel-data-131/perceval/Rat003_20231222/Rat003_20231222.xml';
%session = '/mnt/hubel-data-139/perceval/Rat003_20231228/Rat003_20231228.xml';
%session = '/mnt/hubel-data-140/karadoc/Rat004_20240313/Rat004_20240313.xml';
%session = '/mnt/hubel-data-140/karadoc/Rat004_20240313/Rat004_20240313.xml';

R = regions(session,'states',["sws","rem"],'events',["InfraSlowRhythm/slownr","InfraSlowRhythm/slowavalnr"]);

[filebase,basename] = fileparts(session);
load(fullfile(filebase,basename+".deltaWaves.events.mat"))
load(fullfile(filebase,basename+".spikeDeltaWaves.events.mat"))

R.plotFiringRates(0,10000,'region',["pfc","nr"],'smooth',100,'scale',400); ylim([-1,350])

PlotIntervals(R.eventIntervals('slownr'),'legend','isr')
PlotIntervals(R.eventIntervals('slowavalnr'),'legend','on isr','alpha',0.7)

RasterPlot([deltaWaves.peaks,mean(yLim)*ones(size(deltaWaves.timestamps,1),1)],diff(yLim),'Color','#0fbd5d','DisplayName','delta peak');

R.plotSpikeRaster(2600,2800,'region','nr')

start = 2100; stop = 2150;
[fig,axs] = makeFigure('test','',[3,1]);
set(axs(1:2),'XTick',[])
% PFC
R.plotFiringRates(start,stop,0.01,'region',"pfc",'smooth',10,'ax',axs(1)); ylim(axs(1),[0,300])
RasterPlot([deltaWaves.peaks,mean(yLim)*ones(size(deltaWaves.peaks,1),1)],diff(yLim),'Color','#0fbd5d','DisplayName','Δ peak',ax=axs(1));
RasterPlot([spikeDeltaWaves.peaks,mean(yLim)*ones(size(spikeDeltaWaves.peaks,1),1)],diff(yLim),'Color','#074f00','DisplayName','spike-Δ peak',ax=axs(1));
% NR
R.plotFiringRates(start,stop,'region',"nr",'smooth',10,'ax',axs(2));
PlotIntervals(R.eventIntervals('slowavalnr'),'legend','on isr','alpha',0.7,'ax',axs(2))
% NR raster
R.plotSpikeRaster(start,stop,'region','nr','ax',axs(3))
xlabel(axs(1:2),'')

saveFig(fig,'/mnt/hubel-data-103/Pietro/ReuSlowRythm/Results/FinalFigures/SupplDeltas/nested_deltas',["png","svg"])