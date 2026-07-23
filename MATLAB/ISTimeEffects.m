%% 1. E(ISR% | interval duration) says to analyze nREM intervals starting from 40s

batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/US_intervals_states.batch';
states = ["sws","rem","other"];
args = {states,'sleepm'};
IS_f = runBatch(batch_file,@ISStateDuration_,args,'verbose',true);
f = reverseCellStruct(IS_f);
animal = repelem((1:2).',[14,17],1);
for a = 1 : 2
  animal_f{a} = reverseCellStruct(IS_f(animal==a));
end
clear args a

%%
[~,axs] = makeFigure('dur','E(ISR% | interval duration)',[1,3],'size',[1500,450]);

bin_size = [10,2.5,20];
x_lim = [400,60,1000];
for i = 1 : 3

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

  for a = 1 : 2
    plot(axs(i),time,mean(animal_y{a},'omitmissing')*100,'Color',mean([paperColors(4+i);1,1,1]))
  end

  axes(axs(i))
  semplot(time,y*100,paperColors(4+i),'patchProp',{'FaceAlpha',0.3},'lineProp',{'Marker','o','MarkerFaceColor','w','MarkerSize',1})
  xlabel(axs(i),'interval duration (s)'), xlim(axs(i),[0,x_lim(i)]), ylim(axs(i),[0,100]), yticks(axs(i),[25,50,75,100])%, xticks(axs(i),[40,100,160,220])

end
ylabel(axs(1),'ISR percentage (%)')

clearvars bin_size x_lim bins bin_ind time y animal_bin_ind animal_y a axs i j
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISTimeEffect/2a_dur',["png","svg"]);




%% 2. percentage of slow rhythm over state time
states = ["sws","rem","other"];
window = [60,25,60];
for i = 1 : numel(states)
  args = {states,states(i),"sleepm",'window',window(i)};
  [IS_f_t.(states(i)),time.(states(i))] = runBatch(batch_file,@USTimeEffect_,args,verbose=true);
  f_t.(states(i)) = reverseCellStruct(IS_f_t.(states(i)),@(x) x.');
  for a = 1 : 2
    animal_f_t{a}.(states(i)) = reverseCellStruct(IS_f_t.(states(i))(animal==a),@(x) x.');
  end
  time.(states(i)) = time.(states(i)){1};
end
clear args i

%%
[~,axs] = makeFigure('dur','State time spent in ISR',[1,3],'size',[1500,450]);
x_lim = [4700,1000,5500];
colors = [5,6,7];
for i = 1 : 3

  for a = 1 : 2
    plot(axs(i),time.(states(i)),mean(animal_f_t{a}.(states(i)).sleepm,'omitmissing')*100,'Color',mean([paperColors(colors(i));1,1,1]))
  end

  axes(axs(i))
  semplot(time.(states(i)),f_t.(states(i)).sleepm*100,paperColors(colors(i)),'patchProp',{'FaceAlpha',0.3},'lineProp',{'Marker','o','MarkerFaceColor','w','MarkerSize',1})
  xlabel(axs(i),'shifted time (s)'), xlim(axs(i),[0,x_lim(i)]), ylim(axs(i),[0,100]), yticks(axs(i),[25,50,75,100])

end
ylabel(axs(1),'ISR percentage (%)')
clearvars a animal_f_t axs colors i x_lim
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISTimeEffect/2b_time',["png","svg"]);

% %%
% [r,p] = corrcoef(mean(ref_dur_pool.sws,'omitmissing'),mean(US_fraction_t.sws,'omitmissing'));
% r = r(2); p = p(2);
% 
% makeFigure('time_f',"Fraction of nREM sleep spent in slow rhythm, r: "+num2str(r,4)+", p: "+num2str(p,4),size=[1200,600]);
% semplot(1:size(ref_dur_pool.sws,2),ref_dur_pool.sws-60,[.7,.7,.7],'legend','\langlenREM\rangle','patchProp',{'EdgeColor',[.9,.9,.9],'EdgeAlpha',1,'LineWidth',1.3,'FaceColor','none'}, ...
%   'lineProp',{'LineStyle','--','Marker','o','MarkerFaceColor','w','MarkerSize',8});
% semplot(1:size(US_fraction_t.sws,2),US_fraction_t.sws*100,paperColors(5),'legend','US / nREM','patchProp',{'FaceAlpha',0.3},'lineProp',{'Marker','o','MarkerFaceColor','w','MarkerSize',8});
% xticks(1:10), xlabel('time windows (10 min each)'), xlim([0.5,10.5]), ylabel('\color{black} percentage (%)'), yticks([0,30,60,90]), ylim([-1,100]), legend
% yyaxis right
% set(gca,'YColor',[.7,.7,.7]), ylabel('\color{black} \langleduration\rangle (s)'), ylim([-1,100]), yticks([0,30,60,90]), yticklabels(["60","90","120","150"])
% set(gca, 'YTickLabel', cellfun(@(x) ['\color{black} ' x], get(gca,'YTickLabel'),'UniformOutput',false)); % make labels black
% yyaxis left
% set(gca,'YColor',paperColors(5),'YTickLabel',cellfun(@(x) ['\color{black} ' x], get(gca,'YTickLabel'),'UniformOutput',false)); % make labels black
% 
% makeInset(0.05,0.7,NaN,0.25,'XColor',[.7,.7,.7],'YColor',paperColors(5),'XTick',[60,90,120],'YTick',[0,30,60],'XTickLabel',["60","90","120"],'YTickLabel',["0","30","60"]);
% set(gca, 'XTickLabel', cellfun(@(x) ['\color{black} ' x], get(gca,'XTickLabel'),'UniformOutput',false)); % make labels black
% set(gca, 'YTickLabel', cellfun(@(x) ['\color{black} ' x], get(gca,'YTickLabel'),'UniformOutput',false)); % make labels black
% scatter(mean(ref_dur_pool.sws,'omitmissing'),mean(US_fraction_t.sws,'omitmissing')*100,20,'ok','MarkerFaceColor','k')
% xlim([60,120]), ylim([0,60])
% 
% %%
% [r,p] = corrcoef(mean(ref_dur_pool.other,'omitmissing'),mean(US_fraction_t.other,'omitmissing'));
% r = r(2); p = p(2);
% 
% makeFigure('time_f',"Fraction of awake sleep spent in slow rhythm, r: "+num2str(r,4)+", p: "+num2str(p,4),size=[1200,600]);
% offset = 120; scale = 0.05;
% semplot(1:size(ref_dur_pool.other,2),(ref_dur_pool.other-offset)*scale,[.7,.7,.7],'legend','\langleawake\rangle','patchProp',{'EdgeColor',[.9,.9,.9],'EdgeAlpha',1,'LineWidth',1.3,'FaceColor','none'}, ...
%   'lineProp',{'LineStyle','--','Marker','o','MarkerFaceColor','w','MarkerSize',8});
% semplot(1:size(US_fraction_t.other,2),US_fraction_t.other*100,paperColors(7),'legend','US / awake','patchProp',{'FaceAlpha',0.3},'lineProp',{'Marker','o','MarkerFaceColor','w','MarkerSize',8});
% xticks(1:10), xlabel('time windows (10 min each)'), xlim([0.5,10.5]), ylabel('\color{black} percentage (%)'), yticks([0,10,20]), ylim([-1,25]), legend
% yyaxis right
% set(gca,'YColor',[.7,.7,.7]), ylabel('\color{black} \langleduration\rangle (s)'), ylim([-1,25]), yticks([0,10,20]), yticklabels(string([0,10,20]/scale+offset))
% set(gca, 'YTickLabel', cellfun(@(x) ['\color{black} ' x], get(gca,'YTickLabel'),'UniformOutput',false)); % make labels black
% yyaxis left
% set(gca,'YColor',paperColors(7),'YTickLabel',cellfun(@(x) ['\color{black} ' x], get(gca,'YTickLabel'),'UniformOutput',false)); % make labels black
% 
% makeInset(0.15,0.7,NaN,0.25,'XColor',[.7,.7,.7],'YColor',paperColors(7),'XTick',[150,500],'YTick',[0,5,10],'XTickLabel',["150","500"],'YTickLabel',["0","5","10"]);
% set(gca, 'XTickLabel', cellfun(@(x) ['\color{black} ' x], get(gca,'XTickLabel'),'UniformOutput',false)); % make labels black
% set(gca, 'YTickLabel', cellfun(@(x) ['\color{black} ' x], get(gca,'YTickLabel'),'UniformOutput',false)); % make labels black
% scatter(mean(ref_dur_pool.other,'omitmissing'),mean(US_fraction_t.other,'omitmissing')*100,20,'ok','MarkerFaceColor','k')
% xlim([120,520]), ylim([0,10])




%% fGLS BAD MODEL, LINDEAR INETRACTION ARE NOT ENOUGH
X = [zeros(1,size(US_fraction_t,2)-1);diag(ones(size(US_fraction_t,2)-1,1))]; % categorical variables to define time point
X = repmat(X,size(US_fraction_t,1),1);
Y = US_fraction_t.';
[beta,se] = fgls(X,Y(:),ARLags=4,NumIter=100);

Y2 = Y;
Y2(isnan(Y2)) = 0;
[beta2,se2] = fgls(X,Y2(:),ARLags=4,NumIter=100);


errorbar(beta*100,se*100,'LineWidth',1.3,'Color',myColors(1))
errorbar(beta2*100,se2*100,'LineWidth',1.3,'Color',myColors(2))
xticks(1:10), xlabel('time windows (10 min each)'), ylabel('percentage (%)'), yticks([0,30,60,90]), xlim([0.5,10.5]), ylim([-1,100])



makeFigure('time_f','Fraction of nREM sleep spent in slow rhythm',size=[1200,600]);
distPlot(US_fraction_t*100,[],scatter=1,scale=0.5,salpha=0.7,ssize=35,colors=paperColors(5));
errorbar(beta*100,se*100)
xticks(1:10), xlabel('time windows (10 min each)'), ylabel('percentage (%)'), yticks([20,40,60]), xlim([0.5,10.5])

[beta,p,h,CI] = timeEffectGLS(Y(:),size(US_fraction_t,2));

clear i

% makeFigure('time_f','Fraction of nREM sleep spent in slow rhythm');
% semplot(1:10,US_fraction_t*100,myColors(1),'legend','all')
% xticks(1:10), xlabel('time windows (10 min each)'), ylabel('percentage (%)'), yticks([20,40,60]), xlim([0.5,10.5])



%% 3. same for every US intervals

states = ["sws","rem","other"];
for i = 1 : numel(states)
  args = {states,states(i),"sleepm"};
  [IS_f_t.(states(i)),time.(states(i))] = runBatch(batch_file,@ISInInterval_,args,'verbose',true);
  f_t.(states(i)) = reverseCellStruct(IS_f_t.(states(i)));
  for a = 1 : 2
    animal_f_t{a}.(states(i)) = reverseCellStruct(IS_f_t.(states(i))(animal==a));
  end
  time.(states(i)) = time.(states(i)){1};
end
clear args i

%%
[~,axs] = makeFigure('dur','Interval time spent in ISR',[1,3],'size',[1500,450]);
colors = [5,6,7];
for i = 1 : 3

  for a = 1 : 2
    plot(axs(i),time.(states(i)),mean(animal_f_t{a}.(states(i)).sleepm,'omitmissing')*100,'Color',mean([paperColors(colors(i));1,1,1]))
  end

  axes(axs(i))
  semplot(time.(states(i)),f_t.(states(i)).sleepm*100,paperColors(colors(i)),'patchProp',{'FaceAlpha',0.3})
  xlabel(axs(i),'rescaled time'), xlim(axs(i),[0,1]), ylim(axs(i),[0,100]), yticks(axs(i),[25,50,75,100])

end
ylabel(axs(1),'ISR percentage (%)')
clearvars a axs colors i x_lim
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISTimeEffect/2c_time_int',["png","svg"]);



%% 4. start and stop of ISR inside nREM bout



% MAKE F THAT FINDS WHEN ANIMAL IS AWAKE, THEN ELAPSED TIME BEFORE rhythm DISTR

% LOOk at REM: always inside? its length influences something