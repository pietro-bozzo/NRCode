%session = '/mnt/hubel-data-125/Rat386-20180919/Rat386-20180919.xml';
%session = '/mnt/hubel-data-139/perceval/Rat003_20231226/Rat003_20231226.xml';
%session = '/mnt/hubel-data-139/perceval/Rat003_20231229/Rat003_20231229.xml';
session = '/mnt/hubel-data-139/perceval/Rat003_20231228/Rat003_20231228.xml';
%session = '/mnt/hubel-data-139/karadoc/Rat004_20240229/Rat004_20240229.xml';

events = ["sleepm";"tachem";"sleepn"];
%states = ["trials","successtrials","sws","rem"];
states = ["sws","rem"];
regs = [12;32;71];
[names,stamps] = loadEvents(session);
task = [stamps{2}(1,1),stamps{2}(1,2)];

%%
R = regions(session,events="all",states=states);
R = R.loadSpikes();

%%
fig = R.plotSpikeRaster(0,3820);

%%
saveas(fig,append(PietroPath(),'/AsmbRegions/Results/spike_raster_task1'),'svg')
print(fig,append(PietroPath(),'/AsmbRegions/Results/spike_raster_reuTRY'),'-dsvg','-vector')

%%
fig = R.plotFiringRates(0,0,smooth=10,states=["all",states]);

%%
R = regions(session,phases=names(1:3),states=states,regions=regs);
R = R.loadAssemblies('ICA',0.03,'all','sws');

fig = R.plotAsmbRaster(0,0,states=["all",states]);
task = [stamps{2}(1,1),stamps{2}(1,2)];
xline(task,'r',HandleVisibility='off',LineWidth=2)

%%
saveas(fig,append(PietroPath(),'/AsmbRegions/Results/asmb_raster_ICA_sws'),'svg')


%%
R1 = regions(session,phases=names(1:3),states=states,regions=regs);
R1 = R1.loadAssemblies('ISAC',0.03,'all','all');

fig1 = R1.plotAsmbRaster(0,0,states=["all",states]);
reference = (stamps{2}(1,1) + stamps{2}(1,2)) / 2;
task = [stamps{2}(1,1),stamps{2}(1,2)];
xline(task,'r',HandleVisibility='off',LineWidth=2)

saveas(fig1,append(PietroPath(),'/AsmbRegions/Results/asmb_raster_ISAC'),'svg')

%% load activations, get pre / post task, make plot
stamps_s1 = stamps{1};
stamps_s2 = stamps{3};
act = R1.asmbActivations('all');
assemblies = unique(act(:,2));
act_s1 = Restrict(act,stamps_s1);
act_s2 = Restrict(act,stamps_s2);
% compute activation rates
[n_act_s1,ind1] = groupcounts(act_s1(:,2));
[n_act_s2,ind2] = groupcounts(act_s2(:,2));
% handle assemblies which never activate
activ_rate = zeros(numel(assemblies),2);
activ_rate(ind1,1) = n_act_s1 / sum(stamps_s1(:,2)-stamps_s1(:,1)); % normalize by time duration
activ_rate(ind2,2) = n_act_s2 / sum(stamps_s2(:,2)-stamps_s2(:,1));
% keep most active assemblies
%activ_rate = activ_rate(any(activ_rate > 0.2,2),:);
diff = activ_rate(:,2) - activ_rate(:,1);
[h,p,c] = ttest(diff);
% pair plot
figure(Position=get(0,'ScreenSize'),Name='pair',NumberTitle='off'); hold on
groups = repmat([1;2;NaN],size(activ_rate,1),1);
values = [activ_rate,nan(size(activ_rate,1),1)].';
plot(groups,values(:),color=[0.5,0.5,0.5]);
plot([1,2],mean(activ_rate),color=[0,0,0]);
xlim([0.5,2.5]); ylim([-0.25,2]); set(gca,TickDir='out',Box='off')

figure(Position=get(0,'ScreenSize'),Name='box',NumberTitle='off'); hold on
groups = [ones(size(activ_rate,1),1);repmat(2,size(activ_rate,1),1)];
boxplot([activ_rate(:,1);activ_rate(:,2)],groups);
title('Whole-brain-assemblies'' activation rate',FontSize=17,FontWeight='normal')
ylabel('activation rate (Hz)',FontSize=14)
set(gca,TickDir='out',Box='off')
text(0.5,1.5,append('ttest p: ',string(p),' n: ',string(numel(diff))))
xticklabels(["pre-task sleep","post-task sleep"])