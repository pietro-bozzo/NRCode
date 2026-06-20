% Figure 1c, d, e: Quantification of Nucleus Reuniens infra-slow rhythm

file_root = fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename)));

% set up batch
batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_intervals.batch';
regs = ["pfc","hpc","nr","th","v1"];
fig_labels = ["mPFC","iHPC","NR","TH","V1"];
color_ind = [1,2,3,4,4,4,4,4];
saveFlag = false;

% compute or load intervals
args = {regs,'sleep#1','save',saveFlag,'load',true,'verbose',true};
[isr_intervals,isr_on,isr_fraction,animal] = runBatch(batch_file,@slowIntervals_,args);
isr_intervals = reverseCellStruct(isr_intervals,@(x) x,@(x) x);
isr_on = reverseCellStruct(isr_on,@(x) x(:,1:end-1),@(x) x);
isr_fraction = reverseCellStruct(isr_fraction);

animal = vertcat(animal{:});
[~,~,animal] = unique(animal,'stable');
animal = struct('pfc',animal);
for r = regs(2:end)
  animal.(r) = animal.pfc;
end

% add amygdala from fear conditioning dataset
batch_file_FC = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_intervals_FC.batch';
args = {["ldmPFC","rdmPFC","lvmPFC","rvmPFC","lAMY","rAMY"],'sleep1','legend',dataPath+"/bilateral.anat",'save',saveFlag,'load',true,'verbose',true};
[isr_intervals_FC,~,isr_fraction_FC,animal_FC] = runBatch(batch_file_FC,@slowIntervals_,args);
animal_FC = vertcat(animal_FC{:});
[~,~,animal_FC] = unique(animal_FC,'stable');
animal_FC = animal_FC + numel(unique(animal.pfc));
% concatenate data
isr_intervals_FC = reverseCellStruct(isr_intervals_FC,@(x) x,@(x) x);
isr_intervals.amy = [isr_intervals_FC.lamy;isr_intervals_FC.ramy];
isr_intervals.pfc = [isr_intervals.pfc;isr_intervals_FC.ldmpfc;isr_intervals_FC.rdmpfc;isr_intervals_FC.lvmpfc;isr_intervals_FC.rvmpfc];
isr_fraction_FC = reverseCellStruct(isr_fraction_FC);
isr_fraction.amy = [isr_fraction_FC.lamy;isr_fraction_FC.ramy];
isr_fraction.pfc = [isr_fraction.pfc;isr_fraction_FC.ldmpfc;isr_fraction_FC.rdmpfc;isr_fraction_FC.lvmpfc;isr_fraction_FC.rvmpfc];
regs = [regs,"amy"];
fig_labels = [fig_labels,"AMY"];
animal.amy = [animal_FC;animal_FC];
animal.pfc = [animal.pfc;repmat(animal_FC,4,1)];

% add dorsal and ventral HPC from fear conditioning dataset
batch_file_FC = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_intervals_dvHPC.batch';
args = {["lvHPC","rvHPC","ldHPC","rdHPC"],'sleep1','legend',dataPath+"/bilateral.anat",'save',saveFlag,'load',true,'verbose',true};
[isr_intervals_FC,~,isr_fraction_FC,animal_FC] = runBatch(batch_file_FC,@slowIntervals_,args);
animal_FC = vertcat(animal_FC{:});
[~,~,animal_FC] = unique(animal_FC,'stable');
animal_FC = animal_FC + numel(unique(animal.nr));
% concatenate data
isr_intervals_FC = reverseCellStruct(isr_intervals_FC,@(x) x,@(x) x);
isr_intervals.dhpc = [isr_intervals_FC.ldhpc;isr_intervals_FC.rdhpc];
isr_intervals.vhpc = [isr_intervals_FC.lvhpc;isr_intervals_FC.rvhpc];
isr_fraction_FC = reverseCellStruct(isr_fraction_FC);
isr_fraction.dhpc = [isr_fraction_FC.ldhpc;isr_fraction_FC.rdhpc];
isr_fraction.vhpc = [isr_fraction_FC.lvhpc;isr_fraction_FC.rvhpc];
regs = [regs,"dhpc","vhpc"];
fig_labels = [fig_labels,"dHPC","vHPC"];
[animal.dhpc,animal.vhpc] = deal([animal_FC;animal_FC]);

% pool results
isr_cycles.nr = cellfun(@(x) [x(1:end-1,1),x(2:end,1)],isr_on.nr,'UniformOutput',false);
for i = 1 : numel(isr_intervals.nr)
  [~,ind] = Restrict(isr_cycles.nr{i}(:,2)-1e-7,isr_intervals.nr{i},'verbose','off');
  isr_cycles.nr{i} = isr_cycles.nr{i}(ind,:);
  isr_off.nr{i,1} = SubtractIntervals(isr_intervals.nr{i},isr_on.nr{i});
end

clearvars i ind args batch_file_FC isr_intervals_FC isr_fraction_FC animal_FC r


%% 1c. fraction of time spent in slow rhythm
M = mean(isr_fraction.nr*100,'omitnan');
SE = nansem(isr_fraction.nr*100);
data = structCat(isr_fraction,"vert");
group = repelem((1:numel(regs)).',structfun(@numel,isr_fraction),1);
makeFigure('frac','','size',[3.6,3.5],'format','poster');
b = boxchart(ones(size(data)),data*100,'GroupByColor',group,'MarkerStyle','.','MarkerColor','k','BoxWidth',.7);
for i = 1 : numel(b)
  b(i).LineWidth = 1.5;
  b(i).BoxFaceColor = paperColors(color_ind(i));
end
x_ticks = linspace(.65,1.35,numel(regs));
set(gca,'XTick',x_ticks,'XTickLabel',fig_labels,'XLim',[.6,1.4],'YLim',[-5,45],'Ytick',[0,20,40])
ylabel('recording time (%)')

% test significance of difference
[p,h] = ANOVATests(data,group,'parametric',false,'paired',false,'precedence',2,'alpha',[0,0.05],'test',@ranksum);
pBar(h.h1,x_ticks,-1,'dy',2)

clearvars i data group p M SE b h p x_ticks
saveFlag && saveFig(gcf,fullfile(file_root,'Results/PosterFigures/isr_perc'),"svg");


%% 1d. ON+OFF cycles duration
dur = cellfun(@(x) diff(x,1,2),isr_cycles.nr,'UniformOutput',false);
dur_pool = vertcat(dur{:});
M = mean(dur_pool,'omitnan');
SE = nansem(dur_pool);
tit = num2str(M,4)+" ± "+num2str(SE,4)+" s (mean ± SEM), "+num2str(median(cellfun(@mean,dur),'omitnan'),4)+" (median of avrg)";
makeFigure('cycle','','size',[3.5,1.7],'format','poster');
plotDistr(dur_pool,'color',paperColors(3),'name','duration','unit','s')
set(gca,'YTick',[])
xlabel('cicle duration (s)')
clearvars dur dur_pool M SE tit
saveFlag && saveFig(gcf,fullfile(file_root,'Results/PosterFigures/cycles'),"svg");


%% 1e. ISR epochs duration
dur = cellfun(@(x) diff(x,1,2),isr_intervals.nr,'UniformOutput',false);
dur_pool = vertcat(dur{:});
M = mean(dur_pool,'omitnan');
SE = nansem(dur_pool);
tit = num2str(M,4)+" ± "+num2str(SE,4)+" s (mean ± SEM), "+num2str(median(cellfun(@mean,dur),'omitnan'),4)+" (median of avrg)";
makeFigure('interval','','size',[3.5,1.7],'format','poster');
plotDistr(dur_pool,'nbins',40,'color',paperColors(3),'name','duration','unit','s')
xlabel('ISR epoch duration (s)')
set(gca,'XLim',[0,300],'YTick',[])
saveFlag && saveFig(gcf,fullfile(file_root,'Results/PosterFigures/isr_epochs'),"svg");
clearvars dur dur_pool M SE tit


%%  --- Supplementary figures ---


% compute infra-slow intervals on spikes shuffled per unit
args = {regs,'sleepm','save',saveFlag,'load',false,'shuffle',1000,'verbose',false};
[intervals_sh,slow_avals_sh,isr_fraction_sh,animal] = runBatch(batch_file,@slowIntervals_,args);

% CONCATENATE DATA FROM OTHER DATASETS

% increase from shuffle
for i = 1 : numel(regs)
  isr_fraction_ifs.(regs(i)) = isr_fraction.(regs(i)) - cellfun(@(x) mean(cellfun(@(y) y.(regs(i)), x)), isr_fraction_sh);
  % p-value
  for j = 1 : numel(isr_fraction_sh)
    isr_fraction_p.(regs(i))(j,1) = percentRank(cellfun(@(x) x.(regs(i)), isr_fraction_sh{j}), isr_fraction.(regs(i))(j),'center');
  end
  isr_fraction_p.(regs(i)) = holmBonferroni(isr_fraction_p.(regs(i)));
end
clear i j args

%% S1a. increase from shuffle
M = mean(isr_fraction_ifs.nr*100,'omitnan');
SE = nansem(isr_fraction_ifs.nr*100);
data = structCat(isr_fraction_ifs,"vert");
group = repelem((1:numel(regs)).',structfun(@numel,isr_fraction_ifs),1);
makeFigure('frac',num2str(M,4)+" ± "+num2str(SE,4)+" % (mean "+char(177)+" SEM of NR)",size=[7,3.5]);
b = boxchart(ones(size(data)),data*100,'GroupByColor',group,'MarkerStyle','.','MarkerColor','k');
for i = 1 : numel(b)
  b(i).BoxFaceColor = paperColors(color_ind(i));
end
x_ticks = linspace(.65,1.35,numel(regs));
set(gca,'XTick',x_ticks,'XTickLabel',fig_labels,'XLim',[.6,1.4],'YLim',[-5,45])
ylabel('recording time (%)')

% test significance of difference of each distribution from zero
[p,h] = ANOVATests(data,group,'parametric',false,'paired',false,'precedence',2,'alpha',[0.05,0.05],'test',@ranksum);
pBar(h.h1,x_ticks,-1,'dy',4) % USELESS? JUST CHECK
for i = 1 : numel(h.h0)
  if h.h0(i)
    scatter(x_ticks(i),90,10,'k','*')
  end
end

clearvars i data group M SE b h p x_ticks
saveFlag && saveFig(gcf,fullfile(file_root,'Figures/FigureS1/FigS1a_perc'),["png","svg"]);






%% S1b. inter - ON center intervals duration
for i = 1 : numel(isr_cycles.nr)
  d = diff(mean(isr_on.nr{i},2));
  [~,ok_ind] = Restrict(isr_on.nr{i}(2:end,1)-1e-7,isr_intervals.nr{i},'verbose','off');
  dist{i} = d(ok_ind);
end
dist_pool = vertcat(dist{:});
M = mean(dist_pool,'omitnan');
SE = nansem(dist_pool);
tit = num2str(M,4)+" ± "+num2str(SE,4)+" s (mean ± SEM), "+num2str(median(cellfun(@mean,dist),'omitnan'),4)+" (median of avrg)";
makeFigure('interval',tit,size=[4,3]);
plotDistr(dist_pool,'color',paperColors(3),'name','distance','unit','s')
clearvars i d ok_ind dist dist_pool M SE tit
saveFlag && saveFig(gcf,fullfile(file_root,'Figures/FigureS1/FigS1b_on_center'),["png","svg"]);


%% S1c. ON intervals duration
dur = cellfun(@(x) diff(x(:,:),1,2),isr_on.nr,'UniformOutput',false);
dur_pool = vertcat(dur{:});
M = mean(dur_pool,'omitnan');
SE = nansem(dur_pool);
tit = num2str(M,4)+" ± "+num2str(SE,4)+" s (mean ± SEM), "+num2str(median(cellfun(@mean,dur),'omitnan'),4)+" (median of avrg)";
makeFigure('on',tit,size=[4,3]);
plotDistr(dur_pool,'color',paperColors(3),'name','duration','unit','s')
clearvars dur dur_pool M SE tit
saveFlag && saveFig(gcf,fullfile(file_root,'Figures/FigureS1/FigS1c_on'),["png","svg"]);


%% S1d. OFF intervals duration
dur = cellfun(@(x) diff(x,1,2),isr_off.nr,'UniformOutput',false);
dur_pool = vertcat(dur{:});
M = mean(dur_pool,'omitnan');
SE = nansem(dur_pool);
tit = num2str(M,4)+" ± "+num2str(SE,4)+" s (mean ± SEM), "+num2str(median(cellfun(@mean,dur),'omitnan'),4)+" (median of avrg)";
makeFigure('off',tit,size=[4,3]);
plotDistr(dur_pool,'color',paperColors(3),'name','duration','unit','s')
clearvars dur dur_pool M SE tit
saveFlag && saveFig(gcf,fullfile(file_root,'Figures/FigureS1/FigS1d_off'),["png","svg"]);


%% S1e. n of ISR epochs
n = cellfun(@(x) size(x(all(~isnan(x),2),:),1), isr_intervals.nr);
is_ok = cellfun(@(x) isempty(x) || ~all(isnan(x),'all'), isr_intervals.nr); % sessions with NR
n = n(is_ok);
makeFigure('n',num2str(median(n,'omitnan'),4)+" (median over sessions)",size=[4,3]); set(gca,'XColor','none')
distPlot(n,'colors',paperColors(3),'group2',animal.nr(is_ok));
ylabel('n')
set(gca,'YLim',[-1,43],'YTick',[0,10,20,30,40])
saveFlag && saveFig(gcf,fullfile(file_root,'Figures/FigureS1/FigS1e_n'),["png","svg"]);
clearvars n is_ok


%% S1f. n of ON intervals per ISR epoch
n = []; n_avrg = [];
for j = 1 : numel(isr_on.nr)
  avals = isr_on.nr{j}(all(~isnan(isr_on.nr{j}),2),1:2);
  if ~isempty(avals)
    % n avals
    [~,~,ind] = Restrict(avals, isr_intervals.nr{j});
    n_this = accumarray(ind,1);
    n = [n;n_this];
    n_avrg = [n_avrg;mean(n)];
  end
end
M = mean(n,'omitnan');
SE = nansem(n);
tit = num2str(M,4)+" ± "+num2str(SE,4)+" s (mean ± SEM), "+num2str(median(n_avrg,'omitnan'),4)+" (median of average)";
makeFigure('n_on',tit,size=[4,3]);
plotDistr(n,'color',paperColors(3),'name','n','nbins',35)
xlim([0,35])
clear n avals ind n_avrg M SE tit
saveFlag && saveFig(gcf,fullfile(file_root,'Figures/FigureS1/FigS1f_n_on_per_isr'),["png","svg"]);