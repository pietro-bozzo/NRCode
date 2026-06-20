% Figure 1f, g: Nucleus Reuniens infra-slow rhythm concentrates in nREM sleep

file_root = fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename)));

% set up batch
batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_states.batch';
states = ["sws","rem","other"];
args = {states,'event_restrict','sleepm'};
saveFlag = false;

[isr_f_perstate,state_f,ref_f,~,animal] = runBatch(batch_file,@slowIntervalStates_,args);
isr_f_perstate = reverseCellStruct(isr_f_perstate);
state_f = reverseCellStruct(state_f);
ref_f = reverseCellStruct(ref_f);
animal = vertcat(animal{:});
[~,~,animal] = unique(animal,'stable');
clear args


%% 1f. fraction of ISR time spent in each state
makeFigure('state_ref','','size',[4.5,3],'format','poster');
plot([-5,105],[-5,105],'Color',[0.8,0.8,0.8],'LineStyle','--');
data_2d = []; groups = []; centroids = zeros(numel(states),2);
for j = 1 : numel(states)
  a = isr_f_perstate.(states(j));
  b = ref_f.(states(j));
  groups = [groups;j*ones(size(a))];
  data_2d = [data_2d;a,b];
  centroids(j,:) = mean(data_2d(groups==j,:)) * 100;
  h(j) = scatter((b(animal==1))*100,a(animal==1)*100,35,'MarkerFaceColor',paperColors(j+4),'MarkerFaceAlpha',0.7,'MarkerEdgeColor',mean([paperColors(j+4);.7*ones(1,3)]));
  scatter((b(animal==2))*100,a(animal==2)*100,35,'diamond','MarkerFaceColor',paperColors(j+4),'MarkerFaceAlpha',0.7,'MarkerEdgeColor',mean([paperColors(j+4);.7*ones(1,3)]));
end
xlabel('time spent in state (%)'), ylabel('time spent in ISR (%)'), xlim([-5,105]), ylim([-5,105]);
legend(h,["nREM","REM","wake"],'Location','west')

% bootsrap centroid distance
distance = 0;
for j = 1 : numel(states)
  distance = distance + norm(centroids(j,:)-centroids(mod(j,numel(states))+1,:));
end
N = 1000; centroids(:) = 0; distance_sh = zeros(N,1);
for i = 1 : N
  groups_sh = groups(randperm(numel(groups)));
  for j = 1 : numel(states)
    centroids(j,:) = mean(data_2d(groups_sh==j,:)) * 100;
  end
  for j = 1 : numel(states)
    distance_sh(i) = distance_sh(i) + norm(centroids(j,:)-centroids(mod(j,numel(states))+1,:));
  end
end
p = percentRank(distance_sh,distance,"down");

makeInset(0.75,0.4,0.3,0.35,'Color','none','XLim',[0,250],'XTick',[],'YTick',[]);
plotDistr(distance_sh,'color',[0.7,0.7,0.7],'name','\langledistance\rangle','unit','a.u.')
xline(distance,'r','LineWidth',1.2)
%xlabel('\langledistance\rangle (a.u.)','FontSize',7)

saveFlag && saveFig(gcf,fullfile(file_root,'Results/PosterFigures/states'),"svg");
clear a b j h scale p N i pos centroids distance distance_sh groups groups_sh data_2d box_data tit


%% 1g. fraction of each state spent in slow rhythm
box_data = [];
for j = 1 : numel(states)
  box_data = [box_data,state_f.(states(j))];
end

tit = "medians: "+num2str(median(box_data(:,1))*100)+" %, "+num2str(median(box_data(:,2))*100)+" %, "+num2str(median(box_data(:,3))*100)+" %";
makeFigure('state_inv',tit,size=[4,3]);
v = violinplot(ones(numel(box_data),1),box_data(:)*100,'GroupByColor',repelem((1:3).',size(box_data,1),1));
for i = 1 : numel(v)
  v(i).FaceColor = paperColors(4+i);
end

set(gca,'YLim',[-1,95],'YTick',[0,25,50,75],'XTick',[0.78,1,1.22],'XTicklabel',["nREM","REM","wake"])
ylabel('state time (%)')
[~,h] = ANOVATests(box_data,[],'parametric',false,'paired',false,'alpha',[0,0.05]);
pBar(h.h1,xticks,-1,'dy',2)

saveFlag && saveFig(gcf,fullfile(file_root,'Figures/Figure1/Fig1g_state_inv'),["png","svg"]);
clear j box_data h tit v