function [stats,distr] = validateCircEntr_(session,other_session,regs,labels,events,opt)

arguments
  session (1,:) char
  other_session (1,:) char
  regs (:,1) {mustBeNumeric,mustBeInteger,mustBePositive}
  labels (:,1) string
  events (:,1) string = []
  opt.shuffle (1,1) {mustBeNumeric,mustBeInteger,mustBeNonnegative} = 0
  opt.save (1,1) {mustBeLogical} = false
  opt.load (1,1) {mustBeLogical} = true
  opt.verbose (1,1) {mustBeLogical} = false
end

if numel(regs) ~= numel(labels)
  error('validateCircEntr_:labelsSize','Arguments ''regs'' and ''labels'' must have the same number of elements')
end

% try loading data
[filebase,basename] = fileparts(session);
fname = fullfile(filebase,"UltraSlowRythm",basename+"_validateUnitModulation.mat");
if opt.load
  try
    load(fname,'stats');
    extra_regions = setdiff(labels,fieldnames(stats));
    for r = extra_regions'
      [stats.(r).R0,stats.(r).phi0,stats.(r).R,stats.(r).phi,stats.(r).p,stats.(r).h] = deal(NaN);
    end
    for label = labels'
      distr.(label) = [];
    end
    return
  catch
  end
end

% if loading failed, compute unit modulation statistics
R = regions(session,verbose=false);
event_intervals = R.eventIntervals(events);

% keep only regions found in data
[~,ind] = ismember(R.ids,regs);
ind = ind(ind~=0);
labels = labels(ind);
regs = regs(ind);

% slow-rhythm intervals from different session
[filebase,basename] = fileparts(other_session);
us_intervals = readmatrix(fullfile(filebase,'UltraSlowRythm',basename+".slownr"),FileType='text');
us_avals = readmatrix(fullfile(filebase,'UltraSlowRythm',basename+".slowavalnr"),FileType='text');
% restrict to event
us_intervals = IntersectIntervals(us_intervals,event_intervals);
us_avals = Restrict(us_avals,event_intervals);
% if no rhythm is detected
if all(isnan(us_intervals),'all') || isempty(us_intervals)
  for r = labels'
    [stats.(r).R0,stats.(r).phi0,stats.(r).R,stats.(r).phi,stats.(r).p,stats.(r).h] = deal(NaN);
  end
  if opt.save, save(fname,'stats'); end
  return
end

% phase of slow rhythm
phase = interpolatePhase(0.005,us_avals(:,1),3*pi/2,us_avals(:,2),pi/2);
% set up phase values at extremes and mid-points of avals OLD
% phase = [us_avals(:,1),mean(us_avals(:,1:2),2),us_avals(:,2),[mean([us_avals(1:end-1,2),us_avals(2:end,1)],2);NaN]].';
% phase = [phase(1:2,:);phase(2,:)+0.001;phase(3:4,:)];
% phase = [phase(:),repmat([3*pi/2;2*pi;0;pi/2;pi],size(us_avals,1),1)];
% phase = phase(1:end-1,:);
% % phase(t)
% phase = Interpolate(phase,0:0.005:phase(end,1));
phase = Restrict(phase,us_intervals);
% distribution of US phase in slow-rhythm intervals
n_bins = 250; % high number to estimate pdf correction
[dist_phase,bins_phase] = CircularDistribution(phase(:,2),'nBins',n_bins,'normalize','pdf'); % phase was linearly produced: bi-uniform distr

for j = 1 : numel(labels)

  % region data
  spikes = R.spikes('all',regs(j));
  spikes = Restrict(spikes,us_intervals);
  reg_units = R.units(regs(j));

  % declare variables
  distr.(labels(j)) = nan(numel(reg_units),n_bins);
  distr_orig.(labels(j)) = nan(numel(reg_units),n_bins);
  phase_unit.(labels(j)) = cell(numel(reg_units),1);
  for i = ["R0","phi0","R","phi","p","h"]
    stats.(labels(j)).(i) = nan(numel(reg_units),1);
  end

  % analyse
  for i = 1 : numel(reg_units)
    spikes_unit = spikes(spikes(:,2)==reg_units(i),1);
    if numel(spikes_unit) < 100
      continue
    end

    [unit_stats,distr.(labels(j))(i,:),distr_orig.(labels(j))(i,:),shuffled] = circEntrainment(spikes_unit,phase,'intervals',us_intervals, ...
      'shuffle',opt.shuffle);
    for k = ["R0","phi0","R","phi","p","h"]
      stats.(labels(j)).(k)(i) = unit_stats.(k);
    end

  end

  % overall statistics
  stats.(labels(j)).h = holmBonferroni(stats.(labels(j)).p);
  to_log = stats.(labels(j)).h(~isnan(stats.(labels(j)).h));
  opt.verbose && fprintf(1,labels(j)+": "+sum(to_log)+" out of "+numel(to_log)+"\n");

end

if opt.save
  save(fname,'stats')
end

end

% --- Extra code to plot examples in debug mode ---

function plotInDebug()
  % this function is not meant to be called, rather its code can be executed in debug mode to produce plots

  phase = phase; % to avoid MATLAB compilation error due to unknown variable

  % plot phase
  makeFigure('phase',"US phase, "+R.printBasename());
  PlotXY(phase)
  PlotIntervals(us_intervals);
  PlotIntervals(us_avals(:,1:2));

  % % phase PETH wrt spike times
  % [A,B,C] = PETH(phase,spikes_unit,'durations',[-5,5]);
  % figure
  % % SHOULD SORT before
  % PlotColorMap(A,'x',B,'bar','phase')
  % figure
  % plot(C)

  %% linear plot phase distribution for a unit SOME NR UNITS HAVE PREFERRED PHASE!!
  makeFigure('unit',"Phase distribution, unit "+string(i)+", "+labels(j));
  plot(bins_phase,dist_phase,'Color',myColors(1),'DisplayName','US phase')
  plot(bins_unit,distr_orig.(labels(j))(i,:),'Color',myColors(2),'DisplayName','unit phase')
  plot(bins_unit,distr.(labels(j))(i,:),'Color',myColors(4),'DisplayName','corrected unit phase')
  xlim([0,2*pi]), xticks([0,pi/2,pi,3*pi/2,2*pi]), xticklabels(["0","π/2",'π','3π/2','2π']), xlabel('phase (rad)'), ylabel('p(phase)'), legend

  %% plot phase cdf for a unit
  makeFigure('cdf');
  plot(bins_unit,cumtrapz(bins_unit,corrected_dist_unit),'Color',myColors(1),'DisplayName','empirical cdf')
  plot([0,2*pi],[0,1],'--','Color',[0.7,0.7,0.7],'DisplayName','uniform distribution')
  xlim([0,2*pi]), xticks([0,pi/2,pi,3*pi/2,2*pi]), xticklabels(["0","π/2",'π','3π/2','2π']), xlabel('phase (rad)'), ylabel('cdf(phase)'), ylim([0,1]), legend

  %% polar plot of original unit data for significance test
  makeFigure('signif',"Significance test, unit "+string(i)+", "+labels(j)+", signif: "+string(stats.(labels(j)).h(i)),polar=true);
  l = 1;
  for k = 1 : 5 : opt.shuffle
    hs(l) = polarhistogram(shuffled.phase(:,k),50,'Normalization','pdf','EdgeColor',myColors(2),'EdgeAlpha',0.5,'DisplayStyle','stairs','LineWidth',1.1,'DisplayName','shuffled'); l = l + 1; end
  for k = 1 : opt.shuffle
    polarscatter(shuffled.phi0(k),shuffled.R0(k),75,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor',myColors(2)); end
  h(1) = polarhistogram(phase(:,2),100,'Normalization','pdf','EdgeColor',[0.7,0.7,0.7],'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','rhythm');
  h(2) = polarhistogram(phase_unit.(labels(j)){i},50,'Normalization','pdf','EdgeColor',myColors(1),'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','uncorrected');
  polarscatter(stats.(labels(j)).phi0(i),stats.(labels(j)).R0(i),75,'filled','MarkerFaceAlpha',1,'MarkerFaceColor',myColors(1))
  legend([h(2),hs(1),h(1)]); clear k l h hs

  %% polar plot of corrected unit data for significance test
  makeFigure('signif',"Significance test, unit "+string(i)+", "+labels(j)+", signif: "+string(stats.(labels(j)).h(i)),polar=true);
  half_bin_width = (bins_phase(2)-bins_phase(1)) / 2;
  bin_edges = [bins_phase.'-half_bin_width,bins_phase(end)+half_bin_width];
  l = 1;
  for k = 1 : 5 : opt.shuffle
    hs(l) = polarhistogram('BinEdges',bin_edges,'BinCounts',movmean(shuffled.distr(:,k),5),'Normalization','pdf','EdgeColor',myColors(2),'EdgeAlpha',0.5,'DisplayStyle','stairs','LineWidth',1.1,'DisplayName','shuffled'); l = l + 1; end
  for k = 1 : opt.shuffle
    polarscatter(shuffled.phi(k),shuffled.R(k),75,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor',myColors(2)); end
  h(1) = polarhistogram(phase(:,2),100,'Normalization','pdf','EdgeColor',[0.7,0.7,0.7],'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','rhythm');
  h(2) = polarhistogram('BinEdges',bin_edges,'BinCounts',movmean(distr.(labels(j))(i,:),5),'Normalization','pdf','EdgeColor',myColors(1),'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','corrected');
  polarscatter(stats.(labels(j)).phi(i),stats.(labels(j)).R(i),75,'filled','MarkerFaceAlpha',1,'MarkerFaceColor',myColors(1))
  legend([h(2),hs(1),h(1)]); clear k l h hs half_bin_width bin_edges

  %% plot unit pdf (corrected and not)
  makeFigure('stats',"Phase pdf, unit "+string(i)+", "+labels(j)+", signif: "+string(stats.(labels(j)).h(i)),polar=true);
  half_bin_width = (bins_phase(2)-bins_phase(1)) / 2;
  bin_edges = [bins_phase.'-half_bin_width,bins_phase(end)+half_bin_width];
  h(1) = polarhistogram(phase(:,2),100,'Normalization','pdf','EdgeColor',[0.7,0.7,0.7],'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','rhythm');
  h(2) = polarhistogram(phase_unit.(labels(j)){i},50,'Normalization','pdf','EdgeColor',myColors(1),'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','uncorrected');
  polarscatter(stats.(labels(j)).phi0(i),stats.(labels(j)).R0(i),75,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor',myColors(1))
  h(3) = polarhistogram('BinEdges',bin_edges,'BinCounts',movmean(distr.(labels(j))(i,:),5),'Normalization','pdf','EdgeColor',myColors(2),'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','corrected');
  polarscatter(stats.(labels(j)).phi(i),stats.(labels(j)).R(i),75,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor',myColors(2))
  rticks([0.2,0.4]), rlim([0,0.5]), legend(h(3:-1:1))
  clear half_bin_width bin_edges

  %% color map of phase distribution for a region
  figure
  [~,ind] = sort(stats.(labels(j)).p);
  PlotColorMap(distr.(labels(j))(ind,:),'cutoffs',[NaN,0.8],'bar','pdf')

  %% scatter plot of original data for region
  makeFigure('signif',"Average entrainment per unit, "+labels(j),polar=true);
  polarscatter(stats.(labels(j)).phi0(stats.(labels(j)).h),stats.(labels(j)).R0(stats.(labels(j)).h),75,'filled','MarkerFaceColor',paperColors(j),'MarkerFaceAlpha',0.5,'DisplayName','significant')
  polarscatter(stats.(labels(j)).phi0(~stats.(labels(j)).h),stats.(labels(j)).R0(~stats.(labels(j)).h),75,'filled','MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5,'DisplayName','n. s.')
  legend

  %% scatter plot of corrected Z for region
  makeFigure('signif',"Corrected average entrainment per unit, "+labels(j),polar=true);
  polarscatter(stats.(labels(j)).phi(stats.(labels(j)).h),stats.(labels(j)).R(stats.(labels(j)).h),75,'filled','MarkerFaceColor',paperColors(j),'MarkerFaceAlpha',0.5,'DisplayName','significant')
  polarscatter(stats.(labels(j)).phi(~stats.(labels(j)).h),stats.(labels(j)).R(~stats.(labels(j)).h),75,'filled','MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5,'DisplayName','n. s.')
  rticks([0.4,0.8]), legend

  %% average phase distribution per region
  makeFigure('average','Average distribution of US phase')
  for k = 1 : numel(labels)
    plot(bins_phase,smoothdata(mean(distr.(labels(k))),'gaussian',10),'DisplayName',labels(j))
  end
  xlim([0,2*pi]), xticks([0,pi/2,pi,3*pi/2,2*pi]), xticklabels(["0","π/2",'π','3π/2','2π']), xlabel('phase (rad)'), ylabel('<p(phase)>'), legend

end