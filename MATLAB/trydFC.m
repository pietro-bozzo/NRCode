batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_units.batch';
regs = ["hpc","nr","pfc","v1"];
saveFlag = false;

% sleepm
[dFC_p,sh,time,animal] = runBatchParallel(batch_file,@ISdFC_,{regs,'shuffle',1000});
time = time{1};
dFC_peth = reverseCellStruct(dFC_p);
shuffled = reverseCellStruct(sh,@(x) x,@(x) permute(cat(3,x{:}),[3,2,1]));

% sleepn
%[dFC_p,sh] = runBatchParallel(batch_file,@ISdFC_,{regs,'sleepn'});
%dFC_peth.sleepn = reverseCellStruct(dFC_p);
%shuffled.sleepn = reverseCellStruct(sh);
%has_sleepn = ~cellfun(@isempty,dFC_p);

%% figure with significance bar
[~,axs] = makeFigure('corr','',[1,3],'size',[12,4]);
fields = string(fieldnames(dFC_peth));
j = [1,2,4]; % fields to actually plot
for i = 1 : 3
  
  semplot(time,squeeze(mean(shuffled.(fields(j(i))),3,'omitmissing')),.7*[1,1,1],'alpha',0.3,'legend','shuffled','ax',axs(i),'lineProp',{'Color','none'}); 
  % CANNOT WORK WITH 3d SHUFFLED, INSTEAD PLOT: semplot( mean_sessions(shuffle) ) without mean line maybe + 2 ylines with thresholds from maxstat (converted back to data)
  semplot(time,dFC_peth.(fields(j(i))),myColors(i),'alpha',0.3,'legend',strrep(fields(j(i)),'_','-'),'ax',axs(i));
  xline(axs(i),0,'--')
  title(axs(i),upper(strrep(fields(j(i)),'_','-')))

  % test and plot significance
  [~,h] = maxStatisticTest(dFC_peth.(fields(j(i))),shuffled.(fields(j(i))));
  h = double(h);
  h(h==0) = NaN;
  h = [h;h];
  bin_size = time(2) - time(1);
  test_time = [time-bin_size/2;time+bin_size/2];
  yLim = ylim(axs(i));
  dy = diff(yLim);
  plot(axs(i),test_time(:),h(:)*yLim(2)-dy/10,'Color',myColors(i),'HandleVisibility','off','LineWidth',1.5)

end
xlabel(axs(2),'time from OFF-ON transition (s)'), ylabel(axs(1),'\langle|\rho_{ij}|\rangle')
clearvars hand i a fields axs h bin_size test_time yLim dy j

saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISUnits/3g_corr',["png","svg"]);

%% MAYBE MAKE COLORPLOTS INSTEAD? ONE PER TRANSITION
[~,axs] = makeFigure('corr','Average pairwise unit correlation',[1,3],'size',[10,3]);
fields = string(fieldnames(dFC_peth.sleepm));
j = [1,2,4]; % fields to actually plot
for i = 1 : 3
  
  semplot(time,shuffled.sleepm.(fields(j(i)))(has_sleepn,:),.7*[1,1,1],'alpha',0.3,'legend','shuffled','ax',axs(i));
  for a = 1 : 3
    plot(axs(i),time,mean(dFC_peth.sleepm.(fields(j(i)))(animal==a&has_sleepn,:),1,'omitmissing'),'Color',sqrt(paperColors(i+3)))
  end
  semplot(time,dFC_peth.sleepm.(fields(j(i)))(has_sleepn,:),paperColors(i+3),'alpha',0.3,'legend',strrep(fields(j(i)),'_','-'),'ax',axs(i));
  xline(axs(i),0,'--')
  title(axs(i),strrep(fields(j(i)),'_','-')+", sleepm")

  % semplot(time,shuffled.sleepn.(fields(j(i))),.7*[1,1,1],'alpha',0.3,'legend','shuffled','ax',axs(i+3));
  % for a = 1 : 3
  %   plot(axs(i+3),time,mean(dFC_peth.sleepn.(fields(j(i)))(animal(has_sleepn)==a,:),1,'omitmissing'),'Color',sqrt(paperColors(i+3)))
  % end
  % semplot(time,dFC_peth.sleepn.(fields(j(i))),paperColors(i+3),'alpha',0.3,'legend',strrep(fields(j(i)),'_','-'),'ax',axs(i+3));
  % xline(axs(i+3),0,'--')
  % title(axs(i+3),'sleepn')

end
xlabel(axs,'time from OFF-ON transition (s)'), ylabel(axs([1,4]),'\langleunit correlation\rangle'), ylim(axs([1,4]),[-0.005,0.011]), ylim(axs([2,5]),[0,0.03]), ylim(axs([3,6]),[0,0.11])
clearvars hand i a fields axs

saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISUnits/3g_corr',["png","svg"]);