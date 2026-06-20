%% Analyse entrainment of units by US rhythm

file_root = fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename)));

% set up batch
batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_units.batch';
regs = ["nr","hpc","pfc","v1"];
color_ind = [3,2,1,4];
saveFlag = false;

%% 0. validate method
% shuffle sessions per animal
sessions = readBatchFile(batch_file);
animal = string(cellfun(@(y) y{4},cellfun(@(x) strsplit(x,'/'),sessions,'UniformOutput',false),'UniformOutput',false));
for a = unique(animal)'
  s = sessions(animal==a);
  sessions(animal==a) = s(randperm(numel(s)));
end
args = [cellstr(sessions),repmat({regs,labels,'sleepm','shuffle',500,'load',true,'save',false,'verbose',false},numel(sessions),1)];
stats = runBatch(batch_file,@validateCircEntr_,args,verbose=true);
clear a s args
% ERRORS IN Rat004_20240317 Rat004_20240306 Rat003_20231223
% Unable to perform assignment because the size of the left side is 4136-by-1 and the size of the right side is 4137-by-1.
% Error in circEntrainment (line 62)

%% test false positives 
h = reverseCellStruct(stats,@(x) x.h);
ok_ind = structfun(@(x) ~isnan(x),h,'UniformOutput',false);
R = reverseCellStruct(stats,@(x) x.R);

for j = 1 : numel(labels)
  h.(labels(j)) = logical(h.(labels(j))(ok_ind.(labels(j))));
  R.(labels(j)) = R.(labels(j))(ok_ind.(labels(j)));
  frac.(labels(j)) = mean(h.(labels(j)));
end

% acceptable 5% of false positives
for j = 1 : 3
  disp(labels(j)+": "+strjoin(string(binomialCredibleInt(sum(h.(labels(j))),numel(h.(labels(j))))),', '))
end
clear j ok_ind

% all R distributions are different...
% BUT Bayesian fctor proves that an overall model fits the data better in all combinations!
i=3;j=2; binomialBFTest(sum(h.(labels(i))),numel(h.(labels(i))),sum(h.(labels(j))),numel(h.(labels(j))))





%% 1. modulation of spikes by slow-rhythm intervals
args = {regs,'sleepm','shuffle',500,'load',true,'save',false,'verbose',false};
stat = runBatch(batch_file,@USUnitModulation_,args,'verbose',true);
clearvars args

is_ok = reverseCellStruct(stat,@(x) ~isnan(x.h),'fields',regs); % indeces of non-NaN units
h = reverseCellStruct(stat,@(x) x.h,'fields',regs);
h = structFun(@(x,y) logical(x(y)),h,is_ok);
R = reverseCellStruct(stat,@(x) x.R,'fields',regs);
R = structFun(@(x,y) x(y),R,is_ok);
phi = reverseCellStruct(stat,@(x) x.phi,'fields',regs);
phi = structFun(@(x,y) x(y),phi,is_ok);
frac = structfun(@mean,h,'UniformOutput',false);

%% get distributions and shuffled data
n_bins = 50;
args = {regs,'sleepm','n_bins',n_bins,'shuffle',1,'load',false,'verbose',false};
[~,distr,shuffled] = runBatchParallel(batch_file,@USUnitModulation_,args,'verbose',true);
clearvars args

for r = [regs,'bins']
  place_holder.(r) = nan(1,n_bins);
end
[distr{cellfun(@isempty,distr)}] = deal(place_holder);
distr = reverseCellStruct(distr);
phase_bins = distr.bins(1,:);
distr = rmfield(distr,'bins');
distr = structFun(@(x,y) x(y,:),distr,is_ok);
for r = regs
  place_holder.(r) = struct('phi',NaN);
end
[shuffled{cellfun(@isempty,shuffled)}] = deal(place_holder);
phi_sh = reverseCellStruct(shuffled,@(x) x.phi);
phi_sh = structFun(@(x,y) x(y),phi_sh,is_ok);
clearvars is_ok place_holder

%% a) colormap of units
[~,axs] = makeFigure('colormap','',[1,4],'size',[13,3],'format','poster');
OnOffAxes(2600,axs,'repeat',false)
phase_shift = 3/4;
phase_ind = [ceil(n_bins*phase_shift):n_bins,1:floor(n_bins*phase_shift)];
if phase_ind(1) == phase_ind(end)
  phase_ind(1) = [];
end
bins = unwrap(phase_bins(phase_ind));
labels = ["off","off","off","p(s) (z-score)"];
for i = 1 : 4

  % sort non-significantly modulated units
  ns_ind = find(~h.(regs(i)));
  [~,sort_ns] = sort(modulo(phi.(regs(i))(ns_ind),min(bins),max(bins)));
  ns_ind = ns_ind(sort_ns);

  % sort significantly modulated units
  sig_ind = find(h.(regs(i)));
  [~,sort_sig] = sort(modulo(phi.(regs(i))(sig_ind),min(bins),max(bins)));
  sig_ind = sig_ind(sort_sig);

  % z-score distributions
  sort_ind = [ns_ind;sig_ind];
  sorted_distr = zscore(circularSmooth(distr.(regs(i))(sort_ind,:),2),0,2);

  PlotColorMap(sorted_distr(:,phase_ind),'cutoffs',[0,1.6],'x',bins-bins(1),'bar',labels(i),'map','parula','ax',axs(i))
  x_lim = xlim;
  plot(axs(i),x_lim(1)*[1,1],numel(ns_ind)+[0,numel(sig_ind)],'r','LineWidth',1.6)

  title(axs(i),upper(regs(i)))

end
ylabel(axs(1),'units')
clearvars phase_shift phase_ind bins labels ns_ind sort_ns sig_ind sort_sig sort_ind sorted_distr axs i x_lim
saveFlag && saveFig(gcf,fullfile(file_root,'Results/PosterFigures/colormap'),"svg");

%% S3a) example of entrainment assessment
i = 1;
sig_ind = find(h.(regs(i)));
[~,axs] = makeFigure('example',"Example "+upper(regs(i))+" unit",[1,2],'size',[1000,400],'polar',[false,true]);

% linear
x = modulo(phase_bins,-pi/2,3*pi/2) + pi/2;
[~,sort_ind] = sort(x);
x = [x(sort_ind),max(x)+x(sort_ind)];
y = circularSmooth(distr.(regs(i))(sig_ind(17),sort_ind),'gaussian',5);
OnOffAxes(0.21,axs(1))
plot(axs(1),x,[y,y],'Color',paperColors(color_ind(i)))
ylabel(axs(1),'p(s)')

% polar
OnOffAxes(0.3,axs(2),"ticks",0.2)
polarplot(axs(2),phase_bins,circularSmooth(distr.(regs(i))(sig_ind(17),:),'gaussian',5),'Color',paperColors(color_ind(i)))
polarscatter(axs(2),phi.(regs(i))(sig_ind(17)),R.(regs(i))(sig_ind(17)),'filled','MarkerFaceColor',paperColors(color_ind(i+2)))

saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISUnits/S3a_entrain',["png","svg"]);
clearvars sig_ind axs i x y sort_ind

%% Supplementary: run on control sessions

%% b) percentage
makeFigure('signif','Percentage of significantly entrained units per animal','size',[800,500]);
for a = 1 : 3
  frac_animal(:,a) = structFun(@(x,y) mean(x(y==a)),h,animal);
end
plotPercent(frac_animal*100,'colors',repelem(paperColors(color_ind),3,1) .* [repmat([.8;1;1.2]*[1,1,1],3,1);1.05*ones(3,3)],'labels',upper(regs))
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISUnits/3b_perc',["png","svg"]);
clearvars a frac_animal

%% S3c) classification of subpopulations ADD COLORED LINES ON Y AXIS TO SHOW TWO SUBPOP
[~,axs] = makeFigure('subpop','Subpopulations of entrained units',[1,4],'size',[9,3],'format','poster');
OnOffAxes(2600,axs)
phase_shift = 3/4;
phase_ind = [ceil(n_bins*phase_shift):n_bins,1:floor(n_bins*phase_shift)];
if phase_ind(1) == phase_ind(end)
  phase_ind(1) = [];
end
bins = unwrap(phase_bins(phase_ind));
phi_thresh = pi/5;

for i = 1 : 4

  % sort non-significantly modulated units
  ns_ind = find(~h.(regs(i)));
  [~,sort_ns] = sort(modulo(phi.(regs(i))(ns_ind),min(bins),max(bins)));
  ns_ind = ns_ind(sort_ns);

  % sort significantly modulated units
  sig_ind = find(h.(regs(i)));
  [~,sort_sig] = sort(modulo(phi.(regs(i))(sig_ind),min(bins),max(bins)));
  sig_ind = sig_ind(sort_sig);

  x = modulo(phi.(regs(i))(ns_ind),min(bins),max(bins)) - bins(1);
  x = [x;x+x(end)];
  plot(axs(i),x,(1:2*numel(ns_ind))+0.5,'Color',paperColors(color_ind(i)));

  x = modulo(phi.(regs(i))(sig_ind),min(bins),max(bins)) - bins(1);
  x = [x;x+x(end)];
  plot(axs(i),x,(1:2*numel(sig_ind))+2*numel(ns_ind)+0.5,'Color',paperColors(color_ind(i)));

  xline(axs(i),[pi/2-phi_thresh,pi/2+phi_thresh,3*pi/2-phi_thresh,3*pi/2+phi_thresh],'--')

  title(axs(i),upper(regs(i))), ylim(axs(i),[0,2*(numel(ns_ind)+numel(sig_ind))])

end
ylabel(axs(1),'units')
clearvars phase_shift phase_ind bins ns_ind sort_ns sig_ind sort_sig axs i x
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISUnits/S3c_peaks',["png","svg"]);

%% c) status
[~,axs] = makeFigure('psi','',[1,4],'size',[10.5,2.7],'format','poster');
OnOffAxes(3,axs,'ticks',0.5)
for i = 1 : 4

  % distribution per animal
  for a = 1 : 3
    signif_ind = h.(regs(i)) & animal.(regs(i)) == a;
    plotDistr(phi.(regs(i))(signif_ind)+pi/2,'nbins',35,'polar',true,'color',sqrt(paperColors(color_ind(i))),'ax',axs(i))
  end

  % average distribution
  signif_ind = h.(regs(i));
  plotDistr(phi.(regs(i))(signif_ind)+pi/2,'nbins',35,'polar',true,'color',paperColors(color_ind(i)),'name','s','ax',axs(i))
  set(axs(i),'YLim',[0,2.6],'YTick',[0,1,2])

  % add inset with subpopulation percentages
  ax = makeInset(0.5,0.5,0.48,0.48,'ax',axs(i));
  b = plotPercent([mean(h.(regs(i)));mean(~h.(regs(i)))],'ax',ax,'colors',[sqrt(paperColors(color_ind(i)));0.7,0.7,0.7]);
  x = b.XData(1) + [-0.5,0.5]*b.BarWidth;

  prefer_on = phi.(regs(i)) > -phi_thresh & phi.(regs(i)) < phi_thresh;
  prefer_off = phi.(regs(i)) > pi-phi_thresh | phi.(regs(i)) < -pi+phi_thresh;
  perc_on = mean(h.(regs(i)) & prefer_on);
  perc_off = mean(h.(regs(i)) & prefer_off);

  fill(ax,[x,fliplr(x)],[0,0,perc_on,perc_on],paperColors(color_ind(i)))
  fill(ax,[x,fliplr(x)],perc_on+[0,0,perc_off,perc_off],log(1+paperColors(color_ind(i))))

  set(ax,'XColor','none','YLim',[0,1],'FontSize',10)
  if i == 1
    ylabel(ax,'fraction')
  else
    ylabel(ax,'')
  end

end
ylabel(axs(2:end),'')
saveFlag && saveFig(gcf,fullfile(file_root,'Results/PosterFigures/status'),"svg");
clear i a axs signif_ind ns_ind ax b prefer_on prefer_off perc_on perc_off b x ax

%% d) R OLD CODE
% [~,axs] = makeFigure('R','Unit entrainment strength',[1,4],'size',[1400,400]);
% for i = 1 : 4
% 
%   for a = 1 : 3
%     signif_ind = h.(regs(i)) & animal.(regs(i)) == a;
%     plotDistr(R.(regs(i))(signif_ind),'color',sqrt(paperColors(color_ind(i))),'ax',axs(i),'nbins',17)
%   end
% 
%   signif_ind = h.(regs(i));
%   plotDistr(R.(regs(i))(~signif_ind),R.(regs(i))(signif_ind),'color',[0.7,0.7,0.7;paperColors(color_ind(i))],'name','R','label',["n.s.","significant"],'ax',axs(i),'nbins',17)
%   disp([median(R.(regs(i))(~signif_ind)),median(R.(regs(i))(signif_ind))])
%   set(axs(i),'XLim',[0,0.75],'YLim',[0,16],'YTick',[0,6,12])
% 
% end
% ylabel(axs(2:end),'')

%% d) R
makeFigure('R','Unit entrainment strength','size',[1400,400]);
ind = structfun(@numel,R);
ind = repelem((1:numel(ind)).',ind,1);

for a = 1 : 3
  signif_ind = structCat(h,'vert') & structCat(animal,'vert') == a;
  signif_R = structCat(R,'vert');
  signif_R = signif_R(signif_ind);
  s_ind = ind(signif_ind);
  to_remove = ~ismember(1:numel(regs),unique(s_ind));
  for i = 1 : 4
    if to_remove(i)
      s_ind = [s_ind;i];
      signif_R = [signif_R;NaN];
    end
  end
  v = violinplot(s_ind,signif_R,'DensityDirection','negative','GroupByColor',s_ind,'ColorGroupWidth',1, ...
    'DensityWidth',4,'DensityScale','area','LineWidth',1);
  for i = 1 : 4
    v(i).FaceColor = 'none';
    v(i).EdgeColor = sqrt(paperColors(color_ind(i)));
    if to_remove(i)
      delete(v(i))
    end
  end
end

v = violinplot(ind(structCat(h,'vert')),structCat(structFun(@(x,y) x(y),R,h),'vert'),'DensityDirection','negative','GroupByColor',ind(structCat(h,'vert')),'ColorGroupWidth',1, ...
  'DensityWidth',4,'DensityScale','area','LineWidth',1);
for i = 1 : 4
  v(i).FaceColor = paperColors(color_ind(i));
end
v = violinplot(ind(~structCat(h,'vert')),structCat(structFun(@(x,y) x(~y),R,h),'vert'),'DensityDirection','positive','GroupByColor',ind(~structCat(h,'vert')),'ColorGroupWidth',1, ...
  'DensityWidth',4,'DensityScale','area','LineWidth',1);
for i = 1 : 4
  v(i).FaceColor = [.7,.7,.7];
end
x = (1:4)+(-1.5:1.5)*0.25;
set(gca,'XTick',x,'XTickLabel',upper(regs),'YLim',[-0.1,0.6])
ylabel('R')

ind = 2 * ind;
ind(structCat(h,'vert')) = ind(structCat(h,'vert')) - 1;
[p,H] = ANOVATests(structCat(R,'vert'),ind,'paired',false,'parametric',true);
x = x + [-0.2;0.2];
pBar(p.p1,x(:),'draw',[0,0,0,1],'dy',2)

saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISUnits/3d_R',["png","svg"],'pause',1);
clearvars ind a signif_ind signif_R s_ind to_remove i v x p H

%% e) Z OLD CODE
% [~,axs] = makeFigure('psi','Unit entrainment recap (Z = R exp(i*s))',[1,4],'size',[1400,400],'polar',true);
% symbols = ["o","diamond","^"];
% OnOffAxes(1,axs,'ticks',0.5)
% for i = 1 : 4
% 
%   for a = 1 : 3
%     ns_ind = ~h.(regs(i)) & animal.(regs(i)) == a;
%     hand(a) = polarscatter(axs(i),phi.(regs(i))(ns_ind),R.(regs(i))(ns_ind),30,symbols(a),'filled','MarkerFaceAlpha',0.3,'MarkerFaceColor',[0.7,0.7,0.7],'DisplayName','n.s.');
%   end
%   for a = 1 : 3
%     signif_ind = h.(regs(i)) & animal.(regs(i)) == a;
%     hand(a+3) = polarscatter(axs(i),phi.(regs(i))(signif_ind),R.(regs(i))(signif_ind),30,symbols(a),'filled','MarkerFaceAlpha',0.3,'MarkerFaceColor',paperColors(color_ind(i)), ...
%       'DisplayName','significant');
%   end
% 
%   legend(axs(i),hand([1,4]),'location','northwest')
% 
% end
% 
% saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISUnits/3d_Z',["png","svg"]);
% clearvars axs i a hand symbols ns_ind signif_ind


%% e) Z
[~,axs] = makeFigure('Z','Unit entrainment recap (Z = R exp(i*s))',[2,4],'size',[1400,600],'polar',true);
OnOffAxes(1,axs,'grid',false)
bin_size = 0.02;
for i = 1 : 4

  % significantly entrained units
  [x,y] = pol2cart(phi.(regs(i))(h.(regs(i))),R.(regs(i))(h.(regs(i))));
  [density,x_ax,y_ax] = histcounts2(x,y,-1:bin_size:1,-1:bin_size:1,'Normalization','pdf');
  x_ax = (x_ax(1:end-1)+x_ax(2:end)) / 2;
  y_ax = (y_ax(1:end-1)+y_ax(2:end)) / 2;
  [x_ax,y_ax] = meshgrid(x_ax,y_ax);
  [theta,rho] = cart2pol(x_ax,y_ax);

  surf(axs(i),theta,rho,smoothdata2(density,'gaussian',5),'EdgeColor','none')

  % n.s.
  [x,y] = pol2cart(phi.(regs(i))(~h.(regs(i))),R.(regs(i))(~h.(regs(i))));
  [density,x_ax,y_ax] = histcounts2(x,y,-1:bin_size:1,-1:bin_size:1,'Normalization','pdf');
  x_ax = (x_ax(1:end-1)+x_ax(2:end)) / 2;
  y_ax = (y_ax(1:end-1)+y_ax(2:end)) / 2;
  [x_ax,y_ax] = meshgrid(x_ax,y_ax);
  [theta,rho] = cart2pol(x_ax,y_ax);

  surf(axs(i+numel(regs)),theta,rho,smoothdata2(density,'gaussian',5),'EdgeColor','none')

end

saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISUnits/3e_Z',["png","svg"]);
clearvars axs bin_size i x y density x_ax y_ax theta rho

%% e) dividing by ON / OFF preference doesn't show tuning to specific transitions (distributions always center on ON, OFF) OLD CODE
% [~,axs] = makeFigure('dist','Average unit spiking distribution',[2,4],'size',[1400,800]);
% color_ind = [3,2,1,4];
% OnOffAxes([0.01,0.35],axs)
% bin_width = 2*pi / 250;
% phase = bin_width/2 : bin_width : 4*pi;
% phase_ind = [ceil(500*7/8):500,1:floor(500*7/8)];
% on_ind = structfun(@(x) (x>-pi/2)&(x<pi/2),phi,'UniformOutput',false);
% for i = 1 : 4
% 
%   signif_ind = h.(regs(i));
% 
%   for a = 1 : 3
%     a_ind = animal.(regs(i)) == a;
%     distr_avrg = repmat(circularSmooth(mean(distr.(regs(i))(~signif_ind & on_ind.(regs(i)) & a_ind,:))),1,2);
%     plot(axs(i),phase,distr_avrg(phase_ind),'Color',sqrt(.7)*[1,1,1],'LineWidth',1.3);
%     distr_avrg = repmat(circularSmooth(mean(distr.(regs(i))(signif_ind & on_ind.(regs(i)) & a_ind,:))),1,2);
%     plot(axs(i),phase,distr_avrg(phase_ind),'Color',sqrt(paperColors(color_ind(i))),'LineWidth',1.3);
% 
%     distr_avrg = repmat(circularSmooth(mean(distr.(regs(i))(~signif_ind & ~on_ind.(regs(i)) & a_ind,:))),1,2);
%     plot(axs(i+4),phase,distr_avrg(phase_ind),'Color',sqrt(.7)*[1,1,1],'LineWidth',1.3);
%     distr_avrg = repmat(circularSmooth(mean(distr.(regs(i))(signif_ind & ~on_ind.(regs(i)) & a_ind,:))),1,2);
%     plot(axs(i+4),phase,distr_avrg(phase_ind),'Color',sqrt(paperColors(color_ind(i))),'LineWidth',1.3);
%   end
% 
%   distr_avrg = repmat(circularSmooth(mean(distr.(regs(i))(~signif_ind & on_ind.(regs(i)),:))),1,2);
%   plot(axs(i),phase,distr_avrg(phase_ind),'Color',.7*[1,1,1],'LineWidth',1.3,'DisplayName','distribution');
%   distr_avrg = repmat(circularSmooth(mean(distr.(regs(i))(signif_ind & on_ind.(regs(i)),:))),1,2);
%   plot(axs(i),phase,distr_avrg(phase_ind),'Color',paperColors(color_ind(i)),'LineWidth',1.3,'DisplayName','distribution');
% 
%   distr_avrg = repmat(circularSmooth(mean(distr.(regs(i))(~signif_ind & ~on_ind.(regs(i)),:))),1,2);
%   plot(axs(i+4),phase,distr_avrg(phase_ind),'Color',.7*[1,1,1],'LineWidth',1.3,'DisplayName','distribution');
%   distr_avrg = repmat(circularSmooth(mean(distr.(regs(i))(signif_ind & ~on_ind.(regs(i)),:))),1,2);
%   plot(axs(i+4),phase,distr_avrg(phase_ind),'Color',paperColors(color_ind(i)),'LineWidth',1.3,'DisplayName','distribution');
% 
% end
% ylabel(axs(1),'\langlespike distr\rangle, prefer ON')
% ylabel(axs(5),'\langlespike distr\rangle, prefer OFF')
% saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISUnits/3e_dist',["png","svg"]);
% clearvars axs bin_width distr_avrg i on_ind phase signif_ind


%% load unit PETHs
args = {regs,'sleepm'};
[peth_on_off,peth_off_on,t] = runBatchParallel(batch_file,@ISPeth_,args,'verbose',true);
clearvars args

is_ok = reverseCellStruct(stat,@(x) ~isnan(x.h),'fields',regs); % indeces of non-NaN units
peth_on_off = reverseCellStruct(peth_on_off);
peth_on_off = structFun(@(x,y) x(y,:),peth_on_off,is_ok);
peth_off_on = reverseCellStruct(peth_off_on);
peth_off_on = structFun(@(x,y) x(y,:),peth_off_on,is_ok);
t = t{1};
clearvars is_ok

%% f) PETHs at transitions
% NOTE: PFC units are tuned to transition in 2 ways: INHIB, EXCIT, show that PREFER ON are actually PREFER OFF WITH INHIB?
% CAN MAYBE USE DIFFERENCE WITH SHUFFLE T OPROVE OR MAYBE SHOW THAT PEAKS ARE DIAGONAL< THROUGHS ARE VERTICAL
% BUT FIRST RE COMPUTE WITH ONLY ON OFF INTERVALS LONGER THAN 5 s
% Science Drieu 2018 fa sottrazione di matrici con bella palette di colori, probabilmente ci vuole smooth 2d

[~,axs] = makeFigure('peth','Unit transition PETH',[2,4],'size',[1200,800],'axProp',{'FontSize',10,'TitleFontSizeMultiplier',1});
%on_ind = structfun(@(x) (x>-phi_thresh)&(x<phi_thresh),phi,'UniformOutput',false);
%off_ind = structfun(@(x) (x<-pi+phi_thresh)|(x>pi-phi_thresh),phi,'UniformOutput',false);
labels = ["off","off","off","event rate (z-score)"];
smooth = 15;
limits = [-1.7,1.7];
% replicate sorting of a)
phase_shift = 3/4;
phase_ind = [ceil(n_bins*phase_shift):n_bins,1:floor(n_bins*phase_shift)];
if phase_ind(1) == phase_ind(end)
  phase_ind(1) = [];
end
bins = unwrap(phase_bins(phase_ind));

for i = 1 : 4

  signif_ind = find(h.(regs(i)));

  [~,sort_sig] = sort(modulo(phi.(regs(i))(signif_ind),min(bins),max(bins)));
  signif_ind = signif_ind(sort_sig);

  data = zscore(smoothdata(peth_on_off.(regs(i))(signif_ind,:),2,'gaussian',smooth),0,2);
  PlotColorMap(data,'cutoffs',limits,'x',t,'bar',labels(i),'map','parula','ax',axs(i),'barProp',{'FontSize',9});
  title(axs(i),upper(regs(i))+', ON ➞ OFF')

  data = zscore(smoothdata(peth_off_on.(regs(i))(signif_ind,:),2,'gaussian',smooth),0,2);
  PlotColorMap(data,'cutoffs',limits,'x',t,'bar',labels(i),'map','parula','ax',axs(i+numel(regs)),'barProp',{'FontSize',9});
  title(axs(i+numel(regs)),'OFF ➞ ON')

  % OLD
  % data = zscore(smoothdata(peth_on_off.(regs(i))(signif_ind & on_ind.(regs(i)),:),2,'gaussian',smooth),0,2);
  % [~,sort_ind] = sort(max(data,[],2));
  % PlotColorMap(data(sort_ind,:),'cutoffs',limits,'x',t,'bar',labels(i),'map','parula','ax',axs(i),'barProp',{'FontSize',9});
  % title(axs(i),upper(regs(i))+', ON ➞ OFF, pref ON')
  % 
  % data = zscore(smoothdata(peth_on_off.(regs(i))(signif_ind & off_ind.(regs(i)),:),2,'gaussian',smooth),0,2);
  % [~,sort_ind] = sort(max(data,[],2));
  % PlotColorMap(data(sort_ind,:),'cutoffs',limits,'x',t,'bar',labels(i),'map','parula','ax',axs(i+numel(regs)),'barProp',{'FontSize',9});
  % title(axs(i+numel(regs)),'ON ➞ OFF, pref OFF')
  % 
  % data = zscore(smoothdata(peth_off_on.(regs(i))(signif_ind & on_ind.(regs(i)),:),2,'gaussian',smooth),0,2);
  % [~,sort_ind] = sort(max(data,[],2));
  % PlotColorMap(data(sort_ind,:),'cutoffs',limits,'x',t,'bar',labels(i),'map','parula','ax',axs(i+2*numel(regs)),'barProp',{'FontSize',9});
  % title(axs(i+2*numel(regs)),'OFF ➞ ON, pref ON')
  % 
  % data = zscore(smoothdata(peth_off_on.(regs(i))(signif_ind & off_ind.(regs(i)),:),2,'gaussian',smooth),0,2);
  % [~,sort_ind] = sort(max(data,[],2));
  % PlotColorMap(data(sort_ind,:),'cutoffs',limits,'x',t,'bar',labels(i),'map','parula','ax',axs(i+3*numel(regs)),'barProp',{'FontSize',9});
  % title(axs(i+3*numel(regs)),'OFF ➞ ON, pref OFF')
  
end

xlabel(axs(end-4:end),'time from transition (s)'), ylabel(axs(1:4:end),'units')
clearvars signif_ind axs i labels smooth
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISUnits/3f_peth',["png","svg"]);



% Next up:
% a. units fire more than once per aval?
% b. units phase precess?
% c. all units escape spikes in down? p of escaping changes in different US periods?