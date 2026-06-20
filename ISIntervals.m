%% Find and analyse infra-slow rhythm intervals in Nucleus Reuniens population firing rate

% set up batch
batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_intervals.batch'; % SESSION FROM Rat004_20240228 TO Rat004_20240303 HAVE TH FOR NOW

% choose parameter values
regs = ["pfc","hpc","nr","th","v1"];
fig_labels = ["mPFC","iHPC","NR","TH","V1"];
saveFlag = false;

%% 1. describe slow-rhythm intervals
% load intervals
args = {regs,'sleep#1','save',saveFlag,'load',true,'verbose',true};
[is_intervals,is_on,is_fraction] = runBatch(batch_file,@slowIntervals_,args,'verbose',true);
is_intervals = reverseCellStruct(is_intervals,@(x) x,@(x) x);
is_on = reverseCellStruct(is_on,@(x) x(:,1:end-1),@(x) x);
is_fraction = reverseCellStruct(is_fraction);

% add amygdala from fear conditioning dataset
batch_file_FC = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_intervals_FC.batch';
args = {["ldmPFC","rdmPFC","lvmPFC","rvmPFC","lAMY","rAMY"],'sleep1','legend',dataPath+"/bilateral.anat",'save',saveFlag,'load',true,'verbose',true};
[is_intervals_FC,~,is_fraction_FC] = runBatch(batch_file_FC,@slowIntervals_,args,verbose=true);
% concatenate data
is_intervals_FC = reverseCellStruct(is_intervals_FC,@(x) x,@(x) x);
is_intervals.amy = [is_intervals_FC.lamy;is_intervals_FC.ramy];
is_intervals.pfc = [is_intervals.pfc;is_intervals_FC.ldmpfc;is_intervals_FC.rdmpfc;is_intervals_FC.lvmpfc;is_intervals_FC.rvmpfc];
is_fraction_FC = reverseCellStruct(is_fraction_FC);
is_fraction.amy = [is_fraction_FC.lamy;is_fraction_FC.ramy];
is_fraction.pfc = [is_fraction.pfc;is_fraction_FC.ldmpfc;is_fraction_FC.rdmpfc;is_fraction_FC.lvmpfc;is_fraction_FC.rvmpfc];
regs = [regs,"amy"];
fig_labels = [fig_labels,"AMY"];

% add dorsal and ventral HPC from fear conditioning dataset
batch_file_FC = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_intervals_dvHPC.batch';
args = {["lvHPC","rvHPC","ldHPC","rdHPC"],'sleep1','legend',dataPath+"/bilateral.anat",'save',saveFlag,'load',true,'verbose',true};
[is_intervals_FC,~,is_fraction_FC] = runBatch(batch_file_FC,@slowIntervals_,args,verbose=true);
% concatenate data
is_intervals_FC = reverseCellStruct(is_intervals_FC,@(x) x,@(x) x);
is_intervals.dhpc = [is_intervals_FC.ldhpc;is_intervals_FC.rdhpc];
is_intervals.vhpc = [is_intervals_FC.lvhpc;is_intervals_FC.rvhpc];
is_fraction_FC = reverseCellStruct(is_fraction_FC);
is_fraction.dhpc = [is_fraction_FC.ldhpc;is_fraction_FC.rdhpc];
is_fraction.vhpc = [is_fraction_FC.lvhpc;is_fraction_FC.rvhpc];
regs = [regs,"dhpc","vhpc"];
fig_labels = [fig_labels,"dHPC","vHPC"];

% pool results
is_cycles.nr = cellfun(@(x) [x(1:end-1,1),x(2:end,1)],is_on.nr,'UniformOutput',false);
for i = 1 : numel(is_intervals.nr)
  [~,ind] = Restrict(is_cycles.nr{i}(:,2)-1e-7,is_intervals.nr{i},'verbose','off');
  is_cycles.nr{i} = is_cycles.nr{i}(ind,:);
  is_off.nr{i,1} = SubtractIntervals(is_intervals.nr{i},is_on.nr{i});
end

% animal id
animal.pfc = repelem([1:8,5:8,5:8,5:8].',[14,17,12,10,4*ones(1,16)],1);
[animal.hpc,animal.nr,animal.th,animal.v1] = deal(repelem((1:4).',[14;17;12;10],1));
animal.amy = repelem([5:8,5:8].',4,1);
[animal.dhpc,animal.vhpc] = deal(repelem([5:8,5:8].',[3,4,4,3,3,4,4,3],1));

clear i j inter on_intervals ind args batch_file_FC is_intervals_FC is_fraction_FC
saveFlag = false;

%% a. fraction of time spent in slow rhythm
M = mean(is_fraction.nr*100,'omitnan');
SE = nansem(is_fraction.nr*100);
data = structCat(is_fraction,"vert");
group = repelem((1:numel(regs)).',structfun(@numel,is_fraction),1);
group2 = structCat(animal,"vert");
makeFigure('frac',"Percentage of recording time in ISR",size=[800,630]);
subtitle(num2str(M,4)+" "+char(177)+" "+num2str(SE,4)+" % (mean "+char(177)+" SEM of NR)")
set(gca,'TitleFontSizeMultiplier',1)
distPlot(data*100,group,'group2',group2,'salpha',0.7,'colors',paperColors([1:4,4,4,4,4]),'labels',fig_labels);
set(gca,YLim=[-1,45],YTick=[0,20,40])
ylabel('fraction of recording time (%)')
% test significance of difference MAYBE DO MANY RANKSUM TESTS
[p,h] = ANOVATests(data,group,'parametric',false,'paired',false,'precedence',2,'alpha',[0,0.05],'test',@ranksum);
% IMPLEMENT h!!
pBar(p.p1(h.h1(:,3)~=0,:),xticks)

clear i data group group2 p M SE
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISIntervals/1c_perc_recording',["png","svg"]);

%% b. cycle duration (ON+OFF)
dur = cellfun(@(x) diff(x,1,2),is_cycles.nr,'UniformOutput',false);
dur_pool = vertcat(dur{:});
M = mean(dur_pool,'omitnan');
SE = nansem(dur_pool);
makeFigure('cycle',"ISR cycles duration",size=[800,630]);
subtitle(num2str(M,4)+" ± "+num2str(SE,4)+" s (mean ± SEM), "+num2str(median(cellfun(@mean,dur),'omitnan'),4)+" (median of average)");
set(gca,'TitleFontSizeMultiplier',1)
plotDistr(dur_pool,'colors',paperColors(3),'name','duration','unit','s')
scale = 0.6;
pos = get(gca,'Position'); pos = pos + pos([3,4,3,4]).*scale.*[1,1,-1.1,-1.1];
axes('Position',pos), adjustAxes(gca,'Color','none','XColor','none')
distPlot(cellfun(@mean, dur),'group2',animal.nr,'colors',paperColors(3)); ylabel('\langleduration\rangle (s)')
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISIntervals/1d_cycles',["png","svg"]);

% c. inter-on-center distance
for i = 1 : numel(is_cycles.nr)
  d = diff(mean(is_on.nr{i},2));
  [~,ok_ind] = Restrict(is_on.nr{i}(2:end,1)-1e-7,is_intervals.nr{i},'verbose','off');
  dist{i} = d(ok_ind);
end
dist_pool = vertcat(dist{:});
M = mean(dist_pool,'omitnan');
SE = nansem(dist_pool);
makeFigure('interval',"ISR inter-on distance",size=[800,630]);
subtitle(num2str(M,4)+" ± "+num2str(SE,4)+" s (mean ± SEM), "+num2str(median(cellfun(@mean,dist),'omitnan'),4)+" (median of average)");
set(gca,'TitleFontSizeMultiplier',1)
plotDistr(dist_pool,'colors',paperColors(3),'name','distance','unit','s')
scale = 0.6;
pos = get(gca,'Position'); pos = pos + pos([3,4,3,4]).*scale.*[1,1,-1.1,-1.1];
axes('Position',pos), adjustAxes(gca,'Color','none','XColor','none')
distPlot(cellfun(@mean, dist).','group2',animal.nr,'colors',paperColors(3)); ylabel('\langledistance\rangle (s)')
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISIntervals/1e_ici',["png","svg"]);

% c. US rhythm intervals duration
dur = cellfun(@(x) diff(x,1,2),is_intervals.nr,'UniformOutput',false);
dur_pool = vertcat(dur{:});
M = mean(dur_pool,'omitnan');
SE = nansem(dur_pool);
makeFigure('interval',"ISR intervals duration",size=[800,630]);
subtitle(num2str(M,4)+" ± "+num2str(SE,4)+" s (mean ± SEM), "+num2str(median(cellfun(@mean,dur),'omitnan'),4)+" (median of average)");
set(gca,'TitleFontSizeMultiplier',1)
plotDistr(dur_pool,'colors',paperColors(3),'name','duration','unit','s')
scale = 0.6;
pos = get(gca,'Position'); pos = pos + pos([3,4,3,4]).*scale.*[1,1,-1.1,-1.1];
axes('Position',pos), adjustAxes(gca,'Color','none','XColor','none')
distPlot(cellfun(@mean, dur),'group2',animal.nr,'colors',paperColors(3)); ylabel('\langleduration\rangle (s)')
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISIntervals/1f_intervals',["png","svg"]);
clear dur dur_pool scale pos M SE dist dist_pool i d ok_ind

%% Supplementaries
% a. US rhythm ON intervals duration NEED SUPPLEMENTARY WHERE I CONTROL HOW MAY AVALANCHES THERE WOULD BE WITHOUT THE LONG AVAL CUTOFF
dur = cellfun(@(x) diff(x(:,:),1,2),is_on.nr,'UniformOutput',false);
dur_pool = vertcat(dur{:});
M = mean(dur_pool,'omitnan');
SE = nansem(dur_pool);
makeFigure('on','US rhythm ON intervals duration',size=[800,630]);
subtitle(num2str(M,4)+" ± "+num2str(SE,4)+" s (mean ± SEM), "+num2str(median(cellfun(@mean,dur),'omitnan'),4)+" (median of average)");
set(gca,'TitleFontSizeMultiplier',1)
plotDistr(dur_pool,'colors',paperColors(3),'name','duration','unit','s')
scale = 0.6;
pos = get(gca,'Position'); pos = pos + pos([3,4,3,4]).*scale.*[1,1,-1.1,-1.1];
axes('Position',pos), adjustAxes(gca,'Color','none','XColor','none')
distPlot(cellfun(@mean, dur),'group2',animal.nr,'colors',paperColors(3)); ylabel('\langleduration\rangle (s)')
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISIntervals/Supplementary/S1a_on',["png","svg"]);

% b. US rhythm OFF intervals duration  NEED SUPPLEMENTARY WHERE I CONTROL HOW MAY AVALANCHES THERE WOULD BE WITHOUT THE LONG SILENCE CUTOFF
dur = cellfun(@(x) diff(x,1,2),is_off.nr,'UniformOutput',false);
dur_pool = vertcat(dur{:});
M = mean(dur_pool,'omitnan');
SE = nansem(dur_pool);
makeFigure('off','US rhythm OFF intervals duration',size=[800,630]);
subtitle(num2str(M,4)+" ± "+num2str(SE,4)+" s (mean ± SEM), "+num2str(median(cellfun(@mean,dur),'omitnan'),4)+" (median of average)");
set(gca,'TitleFontSizeMultiplier',1)
plotDistr(dur_pool,'colors',paperColors(3),'name','duration','unit','s')
scale = 0.6;
pos = get(gca,'Position'); pos = pos + pos([3,4,3,4]).*scale.*[1,1,-1.1,-1.1];
axes('Position',pos), adjustAxes(gca,'Color','none','XColor','none')
distPlot(cellfun(@mean, dur),'group2',animal.nr,'colors',paperColors(3)); ylabel('\langleduration\rangle (s)')
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISIntervals/Supplementary/S1b_off',["png","svg"]);

% c. number of US inetrvals
n = cellfun(@(x) size(x(all(~isnan(x),2),:),1), is_intervals.nr);
n = n(cellfun(@(x) isempty(x) || ~all(isnan(x),'all'), is_intervals.nr)); % remove sessions without NR
makeFigure('n','number of US rhythm intervals',size=[800,630]); set(gca,'XColor','none')
subtitle(num2str(median(n,'omitnan'),4)+" (median over sessions)")
distPlot(n,'colors',paperColors(3)); ylabel('n')
set(gca,'TitleFontSizeMultiplier',1,'YLim',[-1,43],'YTick',[0,10,20,30,40])
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISIntervals/Supplementary/S1c_n_intervals',["png","svg"]);

% d. number of ON intervals per US interval
n = []; n_avrg = [];
for j = 1 : numel(is_on.nr)
  avals = is_on.nr{j}(all(~isnan(is_on.nr{j}),2),1:2);
  if ~isempty(avals)
    % n avals
    [~,~,ind] = Restrict(avals, is_intervals.nr{j});
    n_this = accumarray(ind,1);
    n = [n;n_this];
    n_avrg = [n_avrg;mean(n)];
  end
end
M = mean(n,'omitnan');
SE = nansem(n);
makeFigure('n_on','number of ON intervals per US rhythm interval',size=[800,630]);
subtitle(num2str(M,4)+" ± "+num2str(SE,4)+" s (mean ± SEM), "+num2str(median(n_avrg,'omitnan'),4)+" (median of average)");
set(gca,'TitleFontSizeMultiplier',1)
plotDistr(n,'colors',paperColors(3),'name','n','nbins',35)
scale = 0.6;
pos = get(gca,'Position'); pos = pos + pos([3,4,3,4]).*scale.*[1,1,-1.1,-1.1];
axes('Position',pos), adjustAxes(gca,'Color','none','XColor','none')
distPlot(n_avrg,'colors',paperColors(3)); ylabel('\langlen\rangle')
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISIntervals/Supplementary/S1d_n_on_per_interval',["png","svg"]);

clear dur dur_pool scale pos n n_avrg ind avals j n_this M SE




%% 2. average cycle shape
% set up batch
batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_intervals.batch';
regs = ["pfc","hpc","nr","v1"];
args = {regs};
[firing_on,firing_off,peth] = runBatch(batch_file,@ISFiringRate_,args,'ignore_args',true,'verbose',true);
f_on = reverseCellStruct(firing_on);
f_off = reverseCellStruct(firing_off);
p = reverseCellStruct(peth,@(x) mean(x,1,'omitmissing'));

animal = repelem((1:3).',[14,17,12],1);

%%
[~,axs] = makeFigure('shape','ISR average firing rate PETH',[1,4],'size',[1400,400]);
for i = 1 : 4
  title(axs(i),regs(i))
  for j = 1 : 3
    animal_peth = reverseCellStruct(peth(animal==j),@(x) mean(x,1,'omitmissing'));
    plot(axs(i),animal_peth.t(1,:),mean(animal_peth.(regs(i)),'omitmissing'),'Color',sqrt(paperColors(i)));
  end
  axes(axs(i))
  semplot(p.t(1,:),p.(regs(i)),paperColors(i));
end
ylabel(axs(1),'population firing rate (Hz)'), xlabel(axs([1,3]),'time from ON / OFF transition'), set(axs,'YTick',[1,2,3],'YLim',[0.5,3.5])
clear i j animal_peth axs

saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISIntervals/1f_peth',["png","svg"]);

%%
[~,axs] = makeFigure('shape','ISR average firing rate',[1,4],'size',[1400,400]);
OnOffAxes([0,3.5],axs)
phi = linspace(0,2*pi,1000);
for i = 1 : 4
  title(axs(i),regs(i))
  for j = 1 : 3
    animal_firing_on = reverseCellStruct(firing_on(animal==j));
    animal_firing_off = reverseCellStruct(firing_off(animal==j));
    plot(axs(i),phi,[mean(animal_firing_on.(regs(i)),'omitmissing'),mean(animal_firing_off.(regs(i)),'omitmissing')],'Color',mean([paperColors(i);1,1,1]));
  end
  axes(axs(i))
  semplot(phi(1:500),f_on.(regs(i)),paperColors(i));
  semplot(phi(501:1000),f_off.(regs(i)),paperColors(i));
  
  axs(i).XAxis.Label.Color = [0,0,0];
end
ylabel(axs(1),'population firing rate (Hz/unit)'), xlabel(axs,'normalized time'), set(axs,'XLim',[0,2*pi])
clear i j animal_firing_on animal_firing_off axs

%% FOR NOW I MANUALLY REMOVE EXTRA LABELS
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISIntervals/1g_traces',["png","svg"]);







%% 3. autocorrelation of population firing rate
batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_intervals.batch';
regs = ["pfc","hpc","nr","v1"];
args = {regs};
[firing_on,firing_off,peth] = runBatch(batch_file,@ISAutocorrelation_,args,'ignore_args',true,'verbose',true);







%% 4. state repartition of intervals
batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/US_intervals_states.batch';
states = ["sws","rem","other"];
args = {states,'event_restric','sleepm'};
[us_f_perstate,state_f,ref_f] = runBatch(batch_file,@slowIntervalStates_,args,verbose=true);
us_f_perstate = reverseCellStruct(us_f_perstate);
state_f = reverseCellStruct(state_f);
ref_f = reverseCellStruct(ref_f);
animal = repelem((1:2).',[14,17],1);
clear args

%% fraction of slow rhythm spent in each state
makeFigure('state_ref','Repartition of NR ISR between states',size=[800,630]);
plot([-5,105],[-5,105],'Color',[0.8,0.8,0.8],'LineStyle','--');
data_2d = []; groups = []; centroids = zeros(numel(states),2);
for j = 1 : numel(states)
  a = us_f_perstate.(states(j));
  b = ref_f.(states(j));
  groups = [groups;j*ones(size(a))];
  data_2d = [data_2d;a,b];
  centroids(j,:) = mean(data_2d(groups==j,:)) * 100;
  h(j) = scatter((b(animal==1))*100,a(animal==1)*100,35,'MarkerFaceColor',paperColors(j+4),'MarkerFaceAlpha',0.7,'MarkerEdgeColor',mean([paperColors(j+4);.5*ones(1,3)]));
  scatter((b(animal==2))*100,a(animal==2)*100,35,'diamond','MarkerFaceColor',paperColors(j+4),'MarkerFaceAlpha',0.7,'MarkerEdgeColor',mean([paperColors(j+4);.5*ones(1,3)]));
end
xlabel('time spent in state (%)'), ylabel('time spent in ISR (%)'), xlim([-5,105]), ylim([-5,105]);
legend(h,["nREM","REM","wake"],'Location','west')
subtitle("nREM center of mass: ["+num2str(centroids(1,2),4)+", "+num2str(centroids(1,1),4)+"] %"), set(gca,'TitleFontSizeMultiplier',1)

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
title("Repartition of NR ISR between states, p: "+num2str(p,4))

scale = 0.7;
pos = get(gca,'Position'); pos = pos + pos([3,4,3,4]).*scale.*[1,0.9,-1.02,-1.1];
axes('Position',pos), adjustAxes(gca,'XTick',[],'YTick',[],'XLim',[0,distance*1.05],'LabelFontSizeMultiplier',1)
plotDistr(distance_sh,'colors',[0.7,0.7,0.7],'name','\langledistance\rangle','unit','a.u.')
xline(distance,'r','LineWidth',1.2)

saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISIntervals/1h_states',["png","svg"]);
clear a b j h scale p N i pos centroids distance distance_sh groups groups_sh data_2d box_data

%% fraction of each state spent in slow rhythm
box_data = [];
for j = 1 : numel(states)
  box_data = [box_data,state_f.(states(j))];
end

makeFigure('state_inv','Percentage of each state spent in NR ISR',size=[800,630]);
subtitle("medians: "+num2str(median(box_data(:,1))*100)+" %, "+num2str(median(box_data(:,2))*100)+" %, "+num2str(median(box_data(:,3))*100)+" %")
distPlot(box_data*100,'group2',animal,'labels',["nREM","REM","awake"],colors=paperColors(5:7));
set(gca,YLim=[-1,90],YTick=[0,25,50,75,100])
ylabel('fraction of time (%)')
p = ANOVATests(box_data,[],'parametric',false,'paired',false,'alpha',[0,0.05]);
pBar(p.p1,xticks)

saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISIntervals/1i_state_inverse',["png","svg"]);
clear j box_data p











%% 5. repeat 1. with shuffled data
args = {regs,labels,'sleepm','save',saveFlag,'load',true,'shuffle',500,'verbose',false};
[intervals_sh,slow_avals_sh,slow_f_sh] = runBatch(batch_file,@slowIntervals_,args,verbose=true);
saveFlag = false; clear args

%%
for j = 1 : numel(labels)
  % slow_fraction
  % increase from shuffle
  us_f_ifs.(labels(j)) = us_fraction.(labels(j)) - cellfun(@(x) mean(cellfun(@(y) y.(labels(j)), x)), slow_f_sh);
  for i = 1 : numel(slow_f_sh)
    % p-value
    us_f_p.(labels(j))(i,1) = percentRank(cellfun(@(y) y.(labels(j)), slow_f_sh{i}), us_fraction.(labels(j))(i),'center');
  end

  % slow-rhythm interval average length
  % dur = cellfun(@(x) mean(cellfun(@(y) mean(diff(y.(labels(j)),1,2),'omitmissing'), x),'omitmissing'), intervals_sh);
  % dur(isnan(dur)) = 0;
  % avrg_dur = us_dur_avrg.(labels(j));
  % avrg_dur(isnan(avrg_dur)) = 0;
  % us_dur_ifs.(labels(j)) = avrg_dur - dur;
  % for i = 1 : numel(intervals_sh)
  %   % p-value
  %   us_dur_p.(labels(j))(i,1) = percentRank(cellfun(@(y) mean(diff(y.(labels(j)),1,2)), intervals_sh{i}), us_dur_avrg.(labels(j))(i),'center');
  % end
end
clear i j dur avrg_dur

%% IFS of US fraction
box_data = []; groups = [];
for i = 1 : 5
  box_data = [box_data;us_f_ifs.(labels(i))];
  groups = [groups;i*ones(size(us_f_ifs.(labels(i))))];
end

makeFigure('dur','Increase from shuffle of ISR time fraction',size=[800,600]);
distPlot(box_data*100,groups,'xtlabels',labels,scatter=1,scale=0.5,salpha=0.7,colors=paperColors([1:4,4]));
set(gca,YLim=[-3,45],YTick=[0,20,40]), ylabel('Δ percentage (%)')

[p,h] = ANOVATests(box_data,groups,'parametric',false,'paired',false,'precedence',2);
pBar(p.p1,1:5), text(1,20,'NOTE: only NR\newlinesignificantly ~= 0')

clear i box_data goups











%% 3. NR population firing rate concentrates power below 0.1 Hz DEPRECATED

args = {regs,labels,'f_min',0.005};
[psd,power_ratio,duration] = runBatch(batch_file,@slowPSD_,args,verbose=true);

% get frequency range
[~,i] = max(cellfun(@(x) numel(x.f),psd));
f_range = psd{i}.f;

% old code used when f differs across sessions
% f_min = min(cellfun(@(x) x.f(1),psd));
% f_max = max(cellfun(@(x) x.f(end),psd));
% f_N = round(mean(cellfun(@(x) numel(x.f),psd)));
% f_range = logspace(log10(f_min),log10(f_max),f_N);
% for i = 1 : numel(psd)
%   for field = labels
%     ind = discretize(psd{i}.f,f_range);
%     averr = accumarray(ind,psd{i}.(field),[],@mean);
%     sem = accumarray(ind,psd{i}.(field),[],@(x) std(x)/sqrt(length(x))); % CHECK
%     psd_sem.(field) = psd_sem.(field) + averr;
%   end
% end

% pool psd
for field = labels
  psd_sem.(field) = nan(numel(psd),numel(f_range));
  power_box_struct.(field) = [];
end
for i = 1 : numel(psd)
  fields = string(fieldnames(psd{i})).';
  for field = fields(fields~='f')
    % get psd
    ind = ismember(f_range,psd{i}.f);
    psd_sem.(field)(i,ind) = psd{i}.(field);
    % pool power ratio
    power_box_struct.(field) = [power_box_struct.(field);power_ratio{i}.(field)];
  end
end
power_box = [];
power_order = [3,1,2,4];
for i = 1 : numel(power_order)
  field = labels(power_order(i));
  values = power_box_struct.(field);
  power_box = [power_box;values,i*ones(size(values))];
end
clear i field fields power_box_struct values

% TRY AGAIN COMPUTING PSD WITH FIXED LOW F, SHOULD ALSO REMOVE SESSIONS
% WHERE DURATION IS NOT BIG ENOUGH TO ANALIZE TARGET F!!
%% plot psd
fig = makeFigure('psd','Wavelet power spectral density of population firing rate during slow rhythm');
semplot(f_range,psd_sem.nr,myColors(1),legend='nr')
semplot(f_range,psd_sem.pfc,myColors(2),legend='pfc')
semplot(f_range,psd_sem.hpc,myColors(3),legend='hpc')
semplot(f_range,psd_sem.th,myColors(4),legend='th')
xline(0.1,Color=myColors(5,'IBMcb'),LineWidth=1.8,DisplayName='target f')
legend()
adjustAxes(gca,'Xscale','log','YScale','log','XLim',[f_range(1),f_range(end)],'YLim',[500,10^6])
xlabel('time (s, log)');
ylabel('psd (W * s, log)');
saveFlag && saveFig(fig,PietroPath+"/ReuSlowRythm/Results/PSD/psd",'svg',pause=1);
saveFlag = false; clear fig