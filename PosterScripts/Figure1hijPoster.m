% Figure 1h, i, j: Nucleus Reuniens infra-slow rhythm dynamics during nREM sleep

file_root = fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename)));

% set up batch
batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_states.batch';
states = ["sws","rem","other"];
args = {states,'sleepm'};
saveFlag = false;

[IS_f,animal] = runBatch(batch_file,@ISStateDuration_,args);
f = reverseCellStruct(IS_f);
animal = vertcat(animal{:});
[~,~,animal] = unique(animal,'stable');
for a = 1 : 2
  animal_f{a} = reverseCellStruct(IS_f(animal==a));
end
clear args a


%% 1h. + S1g. ISR % as a f of state duration
[fig,ax] = makeFigure('dur','','size',[3,3],'format','poster');
drawnow
[fig_suppl,axs] = makeFigure('dur','',[1,3],'size',[10,3]);

bin_size = [10,2.5,20];
x_lim = [400,60,1000];
for i = 1 : numel(states)

  bins = 0 : bin_size(i) : max(f.(states(i))(:,2))+bin_size(i)*2;
  bin_ind = discretize(f.(states(i))(:,2),bins);
  time = (bins(1:end-1)+bins(2:end)) / 2;

  y = [];
  for j = 1 : numel(bins)-1
    y = inhomogeneousHorzcat(y,f.(states(i))(bin_ind==j,1),'pad',NaN);
  end

  % repeat per animal
  for a = 1 : 2
    animal_bin_ind{a} = discretize(animal_f{a}.(states(i))(:,2),bins);
  end
  [animal_y{[1,2]}] = deal([]);
  for j = 1 : numel(bins)-1
    for a = 1 : 2
      animal_y{a} = inhomogeneousHorzcat(animal_y{a},animal_f{a}.(states(i))(animal_bin_ind{a}==j,1),'pad',NaN);
    end
  end

  % main figure
  if i == 1
    for a = 1 : 2
      plot(ax,time,mean(animal_y{a},'omitmissing')*100,'Color',mean([paperColors(4+i);1,1,1]))
    end
    semplot(time,y*100,paperColors(4+i),'patchProp',{'FaceAlpha',0.3},'lineProp',{'Marker','o','MarkerFaceColor','w','MarkerSize',1},'ax',ax)
    xlabel(ax,'interval duration (s)'), xlim(ax,[0,x_lim(i)]), ylim(ax,[0,100]), yticks(ax,[25,50,75,100])
    % test for increase
    test_time = time(time<200);
    reference = y(:,time>=200);
    p_vals = zeros(size(test_time));
    for t = 1 : numel(test_time)
      p_vals(t) = ranksum(y(:,t),reference(:),'tail','left');
    end
    p_vals = holmBonferroni(p_vals);
    p_vals(p_vals==0) = NaN;
    p_vals = [p_vals,p_vals].';
    test_time = [test_time-bin_size(i)/2;test_time+bin_size(i)/2];
    plot(ax,test_time(:),p_vals(:)*90,'k')
  end

  % supplementary figure
  for a = 1 : 2
    plot(axs(i),time,mean(animal_y{a},'omitmissing')*100,'Color',mean([paperColors(4+i);1,1,1]))
  end
  semplot(time,y*100,paperColors(4+i),'patchProp',{'FaceAlpha',0.3},'lineProp',{'Marker','o','MarkerFaceColor','w','MarkerSize',1},'ax',axs(i))
  xlabel(axs(i),'interval duration (s)'), xlim(axs(i),[0,x_lim(i)]), ylim(axs(i),[0,100]), yticks(axs(i),[25,50,75,100])

end
ylabel([ax,axs(1)],'ISR percentage (%)')

saveFlag && saveFig(fig,fullfile(file_root,'Results/PosterFigures/nrem_dur'),"svg");
saveFlag && saveFig(fig_suppl,fullfile(file_root,'Figures/FigureS1/FigS1g_state_dur'),["png","svg"]);
clearvars bin_size x_lim bins bin_ind time y animal_bin_ind animal_y a axs i j t test_time reference ax p_vals fig fig_suppl



%% ISR occurrence inside state intervals

clearvars f IS_f animal_f time

window = [50,25,50];
for i = 1 : numel(states)
  args = {states,states(i),'sleepm','min_len',window(i)};
  [f,t] = runBatch(batch_file,@ISRStartStop_,args);
  isr_f.(states(i)) = reverseCellStruct(f,@(x) mean(x,1));
  for a = 1 : 2
    animal_f{a}.(states(i)) = reverseCellStruct(f(animal==a),@(x) mean(x,1));
  end
  time.(states(i)) = t{1};
end
clear args i t f a

%% 1i. + S1h. ISR % inside intervals start and stop
[fig,ax] = makeFigure('dur','','size',[3,3],'format','poster');
drawnow
[fig_suppl,axs] = makeFigure('dur','',[1,3],'size',[10,3]);
colors = [5,6,7];
for i = 1 : numel(states)

  idx = find(isnan(time.(states(i))));
  x_ticks = [0;time.(states(i))([idx-1,idx+1,numel(time.(states(i)))])];
  x_label = [0;window(i)/2;-window(i)/2;0];

  % main figure
  if i == 1
    for a = 1 : 2
      plot(ax,time.(states(i)),mean(animal_f{a}.(states(i)).sleepm,'omitmissing')*100,'Color',sqrt(paperColors(colors(i))))
    end
    semplot(time.(states(i)),isr_f.(states(i)).sleepm*100,paperColors(colors(i)),'patchProp',{'FaceAlpha',0.3},'ax',ax)
    xlabel(ax,'time (s)'), set(ax,'XTick',x_ticks,'XTickLabel',x_label,'YLim',[0,100],'YTick',[25,50,75,100])
  end

  % supplementary figure
  for a = 1 : 2
    plot(axs(i),time.(states(i)),mean(animal_f{a}.(states(i)).sleepm,'omitmissing')*100,'Color',sqrt(paperColors(colors(i))))
  end
  semplot(time.(states(i)),isr_f.(states(i)).sleepm*100,paperColors(colors(i)),'patchProp',{'FaceAlpha',0.3},'ax',axs(i))
  xlabel(axs(i),'time (s)'), set(axs(i),'XTick',x_ticks,'XTickLabel',x_label,'YLim',[0,100],'YTick',[25,50,75,100])

end
ylabel(axs(1),'ISR percentage (%)')
saveFlag && saveFig(fig,fullfile(file_root,'Results/PosterFigures/start_stop'),"svg");
saveFlag && saveFig(fig_suppl,fullfile(file_root,'Figures/FigureS1/FigS1h_start_stop'),["png","svg"]);
clearvars a axs colors i ax fig fig_suppl colors idx x_ticks x_label



%% ISR occurrence inside intervals rescaled in time

clearvars f IS_f animal_f time

for i = 1 : numel(states)
  args = {states,states(i),"sleepm"};
  [IS_f_t.(states(i)),time.(states(i))] = runBatch(batch_file,@ISInInterval_,args);
  f_t.(states(i)) = reverseCellStruct(IS_f_t.(states(i)),@(x) mean(x,1));
  for a = 1 : 2
    animal_f_t{a}.(states(i)) = reverseCellStruct(IS_f_t.(states(i))(animal==a),@(x) mean(x,1));
  end
  time.(states(i)) = time.(states(i)){1};
end
clear args i

%% S1i. ISR % inside intervals rescaled in time
[fig,ax] = makeFigure('dur','','size',[3,3],'format','poster');
drawnow
[fig_suppl,axs] = makeFigure('dur','',[1,3],'size',[10,3]);
colors = [5,6,7];
for i = 1 : numel(states)

  % main figure
  if i == 1
    for a = 1 : 2
      plot(ax,time.(states(i)),mean(animal_f_t{a}.(states(i)).sleepm,'omitmissing')*100,'Color',sqrt(paperColors(colors(i))))
    end
    semplot(time.(states(i)),f_t.(states(i)).sleepm*100,paperColors(colors(i)),'patchProp',{'FaceAlpha',0.3},'ax',ax)
    xlabel(ax,'rescaled time'), xlim(ax,[0,1]), ylim(ax,[0,100]), yticks(ax,[25,50,75,100])
  end

  % supplementary figure
  for a = 1 : 2
    plot(axs(i),time.(states(i)),mean(animal_f_t{a}.(states(i)).sleepm,'omitmissing')*100,'Color',sqrt(paperColors(colors(i))))
  end
  semplot(time.(states(i)),f_t.(states(i)).sleepm*100,paperColors(colors(i)),'patchProp',{'FaceAlpha',0.3},'ax',axs(i))
  xlabel(axs(i),'rescaled time'), xlim(axs(i),[0,1]), ylim(axs(i),[0,100]), yticks(axs(i),[25,50,75,100])

end
ylabel(axs(1),'ISR percentage (%)')
%saveFlag && saveFig(fig,fullfile(file_root,'Figures/Figure1/Fig1i_rescaled'),["png","svg"]);
saveFlag && saveFig(fig_suppl,fullfile(file_root,'Figures/FigureS1/FigS1i_rescaled'),["png","svg"]);
clearvars a axs colors i ax fig fig_suppl colors



%% 2. ISR occurrence over shifted state time

clearvars animal_f_t f_t IS_f_t time

window = [60,25,60];
for i = 1 : numel(states)
  args = {states,states(i),"sleepm",'window',window(i)};
  [IS_f_t.(states(i)),time.(states(i))] = runBatch(batch_file,@ISShiftedTime_,args);
  f_t.(states(i)) = reverseCellStruct(IS_f_t.(states(i)),@(x) x.');
  for a = 1 : 2
    animal_f_t{a}.(states(i)) = reverseCellStruct(IS_f_t.(states(i))(animal==a),@(x) x.');
  end
  time.(states(i)) = time.(states(i)){1};
end
clear args i


%% 1j. + S1j. ISR % over shifted state time
[fig,ax] = makeFigure('dur','','size',[3,3],'format','poster');
drawnow
[fig_suppl,axs] = makeFigure('dur','',[1,3],'size',[10,3]);
x_lim = [4700,1000,5500];
colors = [5,6,7];
for i = 1 : numel(states)

  % main figure
  if i == 1
    for a = 1 : 2
      plot(ax,time.(states(i)),mean(animal_f_t{a}.(states(i)).sleepm,'omitmissing')*100,'Color',sqrt(paperColors(colors(i))))
    end
    semplot(time.(states(i)),f_t.(states(i)).sleepm*100,paperColors(colors(i)),'patchProp',{'FaceAlpha',0.3},'lineProp',{'Marker','o','MarkerFaceColor','w','MarkerSize',1},'ax',ax)
    xlabel(ax,'shifted time (s)'), xlim(ax,[0,x_lim(i)]), ylim(ax,[0,100]), yticks(ax,[25,50,75,100])
    % test for decrease
    is_test = time.(states(i)) < 0 | time.(states(i)) > 1500;
    test_time = time.(states(i))(is_test);
    reference = f_t.(states(i)).sleepm(:,~is_test);
    p_vals = zeros(size(test_time));
    m = 1;
    for t = find(is_test).'
      p_vals(m) = ranksum(f_t.(states(i)).sleepm(:,t),reference(:),'tail','left');
      m = m + 1;
    end
    p_vals = holmBonferroni(p_vals);
    p_vals(p_vals==0) = NaN;
    p_vals = [p_vals,p_vals].';
    test_time = [test_time-window(i)/2;test_time+window(i)/2];
    plot(ax,test_time(:),p_vals(:)*90,'k')
  end

  % supplementary figure
  for a = 1 : 2
    plot(axs(i),time.(states(i)),mean(animal_f_t{a}.(states(i)).sleepm,'omitmissing')*100,'Color',sqrt(paperColors(colors(i))))
  end
  semplot(time.(states(i)),f_t.(states(i)).sleepm*100,paperColors(colors(i)),'patchProp',{'FaceAlpha',0.3},'lineProp',{'Marker','o','MarkerFaceColor','w','MarkerSize',1},'ax',axs(i))
  xlabel(axs(i),'shifted time (s)'), xlim(axs(i),[0,x_lim(i)]), ylim(axs(i),[0,100]), yticks(axs(i),[25,50,75,100])

end
ylabel(axs(1),'ISR percentage (%)')
saveFlag && saveFig(fig,fullfile(file_root,'Results/PosterFigures/shifted'),"svg");
saveFlag && saveFig(fig_suppl,fullfile(file_root,'Figures/FigureS1/FigS1j_shifted'),["png","svg"]);
clearvars a axs colors i x_lim fig fig_suppl ax is_test m p_vals test_time t reference