%% Analyse performance of animals in PacMaze task ~

% set up batch
batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/PacMazePerformance.batch';
[perf,data] = runBatch(batch_file,@PacPerformance_,{},'verbose',true);
data = reverseCellStruct(data);
data.animal = repelem([1;2],21,1);

%% 1 performance across the session
for field = fieldnames(perf{1})'
  performance.(field{1}) = [];
  for i = 1 : numel(perf)
    performance.(field{1}) = inhomogeneousVertcat(performance.(field{1}),perf{i}.(field{1}).');
  end
end
clearvars field i

%% a) over all sessions
makeFigure('perf','Performance','size',[1200,600]);
semplot(1:size(performance.tachem,2),performance.tachem,myColors(1),'legend','task1')
semplot(1:size(performance.tachea,2),performance.tachea,myColors(2),'legend','task2')
xlabel('trials')
ylabel('successes (%)')
ylim([0,1.01])

%% b) divided by weeks
[~,axs] = makeFigure('perf','Performance',[1,3],'size',[1500,500]);
titles = "week " + ["1 (","2 (no ","3 ("] + "sleep)";

for w = 1 : 3

  is_week = data.week == w;

  for a = 1 : 2
    plot(axs(w),mean(performance.tachem(is_week & data.animal==a,:),'omitmissing'),'Color',sqrt(myColors(1)))
    plot(axs(w),mean(performance.tachea(is_week & data.animal==a,:),'omitmissing'),'Color',sqrt(myColors(2)))
  end

  h(1) = plot(axs(w),mean(performance.tachem(is_week,:),'omitmissing'),'Color',myColors(1),'DisplayName','task1');
  h(2) = plot(axs(w),mean(performance.tachea(is_week,:),'omitmissing'),'Color',myColors(2),'DisplayName','task2');
  legend(axs(w),h,'Location','southeast')
  title(axs(w),titles(w))

end

xlabel(axs,'trials')
ylabel(axs(1),'successes (%)')
ylim(axs,[0,1.01]), xlim(axs,[0,25])

clearvars a w h is_week

%% c) divided by days
[~,axs] = makeFigure('perf','Performance',[2,4],'size',[1500,700]);

for d = 1 : 7

  is_day = data.day == d;

  for a = 1 : 2
    plot(axs(d),mean(performance.tachem(is_day & data.animal==a,:),'omitmissing'),'Color',sqrt(myColors(1)))
    plot(axs(d),mean(performance.tachea(is_day & data.animal==a,:),'omitmissing'),'Color',sqrt(myColors(2)))
  end

  h(1) = plot(axs(d),mean(performance.tachem(is_day,:),'omitmissing'),'Color',myColors(1),'DisplayName','task1');
  h(2) = plot(axs(d),mean(performance.tachea(is_day,:),'omitmissing'),'Color',myColors(2),'DisplayName','task2');
  legend(axs(d),h,'Location','southeast')

end

xlabel(axs,'trials')
ylabel(axs(1),'successes (%)')
ylim(axs,[0,1.01]), xlim(axs,[0,25])

clearvars a d h is_day

%% d) performance as a coefficient

perf_coeff = zeros(2,3);

for a = 1 : 2

  for w = 1 : 3
    
    this_perf.tachem = mean(performance.tachem(data.week == w & data.animal == a,:),1,'omitmissing');
    this_perf.tachea = mean(performance.tachea(data.week == w & data.animal == a,:),1,'omitmissing');
    threshold = prctile(this_perf.tachem,70,2);

    trials.tachem = find(this_perf.tachem >= threshold);
    trials.tachea = find(this_perf.tachea >= threshold);
    disp(string(threshold)+" "+trials.tachem(find(diff(trials.tachem)==1,1))+" "+trials.tachea(find(diff(trials.tachea)==1,1)));

    perf_coeff(a,w) = trials.tachea(find(diff(trials.tachea)==1,1)); % trials.tachem(find(diff(trials.tachem)==1,1));

  end

end


%% e) observing day by day, every week performance is worse than previous one

perf_coeff = zeros(14,3);

for w = 1 : 3

  this_perf.tachem = performance.tachem(data.week == w,:);
  this_perf.tachem(isnan(this_perf.tachem)) = 0;
  this_perf.tachea = performance.tachea(data.week == w,:);
  this_perf.tachea(isnan(this_perf.tachea)) = 0;

  for i = 1 : size(this_perf.tachea,1)
    perf_coeff(i,w) = find(this_perf.tachea(i,1:end-1) & this_perf.tachea(i,2:end),1) / find(this_perf.tachem(i,1:end-1) & this_perf.tachem(i,2:end),1);
    disp(string(find(this_perf.tachem(i,1:end-1) & this_perf.tachem(i,2:end),1))+" "+find(this_perf.tachea(i,1:end-1) & this_perf.tachea(i,2:end),1))
  end
 

end

figure
distPlot(perf_coeff,'withinlines',1,'group2',data.animal(data.week==1)), ylim([0,1.3])



% potrei usare la velocita` dell'animale nel trial per avere un data point per trial?