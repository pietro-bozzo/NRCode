function [stats,distr,distr_orig] = USRipples_(session,opt)
% USRipples_ Analize ripple entrainment w.r.t. slow-rythm oscillations

arguments
  session (1,:) char
  opt.type (1,:) string = "ripples" % type of events to analize
  opt.shuffle (1,1) {mustBeNumeric,mustBeInteger,mustBeNonnegative} = 0
end

[filebase,basename] = fileparts(session);

% slow-rhythm intervals
us_intervals = readmatrix(fullfile(filebase,'UltraSlowRythm',basename+".slownr"),FileType='text');
us_avals = readmatrix(fullfile(filebase,'UltraSlowRythm',basename+".slowavalnr"),FileType='text');

% phase of slow rhythm
% set up phase values at extremes and mid-points of avals
phase = [us_avals(:,1),mean(us_avals(:,1:2),2),us_avals(:,2),[mean([us_avals(1:end-1,2),us_avals(2:end,1)],2);NaN]].';
phase = [phase(1:2,:);phase(2,:)+0.001;phase(3:4,:)];
phase = [phase(:),repmat([3*pi/2;2*pi;0;pi/2;pi],size(us_avals,1),1)];
phase = phase(1:end-1,:);
% phase(t)
phase = Interpolate(phase,0:0.005:phase(end,1));
phase = Restrict(phase,us_intervals);
% distribution of US phase in slow-rhythm intervals
n_bins = 250; % high number to estimate pdf correction
[dist_phase,bins_phase] = CircularDistribution(phase(:,2),'nBins',n_bins,'normalize','pdf'); % phase was linearly produced: bi-uniform distr

for type = opt.type
  % default values
  [stats.(type).R0,stats.(type).phi0,stats.(type).R,stats.(type).phi,stats.(type).p,stats.(type).h] = deal(NaN);
  [distr.(type),distr_orig.(type)] = deal(nan(n_bins,1));

  % try loading events
  try
    event = readmatrix(fullfile(filebase,basename+"."+type),FileType='text');
  catch
    continue
  end
  peak_t = event(:,2);
  % restrict in US intervals
  peak_t = Restrict(peak_t,us_intervals);
  % phase of each ripple
  phase_event = Interpolate(phase,peak_t);

  % distribution of ripple phase values, corrected by prevalence of every phase bin
  [distr_orig.(type),stats.(type).R0,stats.(type).phi0,distr.(type),stats.(type).R,stats.(type).phi] = correctDistr(phase_event(:,2),dist_phase);

  % repeat on shuffled ripples
  shuffled.time = zeros(numel(peak_t),opt.shuffle);
  shuffled.phase = zeros(numel(peak_t),opt.shuffle);
  shuffled.distr = zeros(n_bins,opt.shuffle);
  shuffled.R0 = zeros(opt.shuffle,1);
  shuffled.phi0 = zeros(opt.shuffle,1);
  shuffled.R = zeros(opt.shuffle,1);
  shuffled.phi = zeros(opt.shuffle,1);
  for k = 1 : opt.shuffle
    time_sh = Restrict(peak_t,us_intervals,'shift','on');
    shuffled.time(:,k) = Unshift(shuffleSpikes(time_sh),us_intervals);
    phase_sh = Interpolate(phase,shuffled.time(:,k),'trim','off');
    shuffled.phase(:,k) = phase_sh(:,2);
    [~,shuffled.R0(k),shuffled.phi0(k),shuffled.distr(:,k),shuffled.R(k),shuffled.phi(k)] = correctDistr(shuffled.phase(:,k),dist_phase);
  end
  % H0: shuffled spikes can produce R as high as observed
  stats.(type).p = 1 - percentRank(shuffled.R,stats.(type).R);
  stats.(type).h = stats.(type).p < 0.05;
end

end

% --- helper functions ---

function [distr0,R0,phi0,distr,R,phi] = correctDistr(phase,reference_distr)
  % distribution of unit-spike phase values
  [distr0,bins_unit,statistics] = CircularDistribution(phase,'nBins',numel(reference_distr),'normalize','pdf');
  R0 = statistics.r;
  phi0 = statistics.m;

  % correct distribution by prevalence of every phase bin
  distr = distr0 ./ reference_distr;
  distr = distr / trapz(bins_unit,distr); % normalize to get pdf

  % estimate Z from pdf
  Z = trapz(bins_unit,exp(1i*bins_unit).*distr);
  R = abs(Z);
  phi = angle(Z);
end

% --- Extra code to plot examples in debug mode ---

function plotInDebug()
  % this function is not meant to be called, rather its code can be executed in debug mode to produce plots

  phase = phase; % to avoid MATLAB compilation error due to unknown variable

  %%
  R = regions(session,load_spikes=false,verbose=false);

  %% plot phase
  makeFigure('phase',"US phase, "+R.printBasename());
  PlotXY(phase)
  PlotIntervals(us_intervals);
  PlotIntervals(us_avals(:,1:2));

  %% plot phase distribution for events
  makeFigure('event',R.printBasename+", "+type+" phase distribution");
  plot(bins_phase,dist_phase,'Color',[0.7,0.7,0.7],'DisplayName','US phase')
  plot(bins_phase,distr_orig.(type),'Color',myColors(1),'DisplayName','ripple phase')
  plot(bins_phase,distr.(type),'Color',myColors(2),'DisplayName','corrected ripple phase')
  xlim([0,2*pi]), xticks([0,pi/2,pi,3*pi/2,2*pi]), xticklabels(["0","π/2",'π','3π/2','2π']), xlabel('phase (rad)'), ylabel('p(phase)'), legend

  %% polar plot of original ripple data for significance test
  makeFigure('signif',"Significance test, p: "+string(stats.(type).p),polar=true);
  l = 1;
  for k = 1 : 5 : opt.shuffle
    hs(l) = polarhistogram(shuffled.phase(:,k),50,'Normalization','pdf','EdgeColor',myColors(2),'EdgeAlpha',0.5,'DisplayStyle','stairs','LineWidth',1.1,'DisplayName','shuffled'); l = l + 1; end
  for k = 1 : opt.shuffle
    polarscatter(shuffled.phi0(k),shuffled.R0(k),75,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor',myColors(2)); end
  h(1) = polarhistogram(phase(:,2),100,'Normalization','pdf','EdgeColor',[0.7,0.7,0.7],'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','rhythm');
  h(2) = polarhistogram(phase_event,50,'Normalization','pdf','EdgeColor',myColors(1),'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','uncorrected');
  polarscatter(stats.(type).phi0,stats.(type).R0,75,'filled','MarkerFaceAlpha',1,'MarkerFaceColor',myColors(1))
  legend([h(2),hs(1),h(1)]); clear k l h hs

  %% polar plot of corrected unit data for significance test
  makeFigure('signif',"Significance test, p: "+string(stats.(type).p),polar=true);
  half_bin_width = (bins_phase(2)-bins_phase(1)) / 2;
  bin_edges = [bins_phase.'-half_bin_width,bins_phase(end)+half_bin_width];
  l = 1;
  for k = 1 : 5 : opt.shuffle
    hs(l) = polarhistogram('BinEdges',bin_edges,'BinCounts',movmean(shuffled.distr(:,k),5),'Normalization','pdf','EdgeColor',myColors(2),'EdgeAlpha',0.5,'DisplayStyle','stairs','LineWidth',1.1,'DisplayName','shuffled'); l = l + 1; end
  for k = 1 : opt.shuffle
    polarscatter(shuffled.phi(k),shuffled.R(k),75,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor',myColors(2)); end
  h(1) = polarhistogram(phase(:,2),100,'Normalization','pdf','EdgeColor',[0.7,0.7,0.7],'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','rhythm');
  h(2) = polarhistogram('BinEdges',bin_edges,'BinCounts',movmean(distr.(type),5),'Normalization','pdf','EdgeColor',myColors(1),'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','corrected');
  polarscatter(stats.(type).phi,stats.(type).R,75,'filled','MarkerFaceAlpha',1,'MarkerFaceColor',myColors(1))
  legend([h(2),hs(1),h(1)]); clear k l h hs half_bin_width bin_edges

  %% plot ripple pdf (corrected and not)
  makeFigure('stats',"Phase pdf, p: "+string(stats.(type).p),polar=true);
  half_bin_width = (bins_phase(2)-bins_phase(1)) / 2;
  bin_edges = [bins_phase.'-half_bin_width,bins_phase(end)+half_bin_width];
  h(1) = polarhistogram(phase(:,2),100,'Normalization','pdf','EdgeColor',[0.7,0.7,0.7],'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','rhythm');
  h(2) = polarhistogram(phase_event,50,'Normalization','pdf','EdgeColor',myColors(1),'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','uncorrected');
  polarscatter(stats.(type).phi0,stats.(type).R0,75,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor',myColors(1))
  h(3) = polarhistogram('BinEdges',bin_edges,'BinCounts',movmean(distr.(type),5),'Normalization','pdf','EdgeColor',myColors(2),'DisplayStyle','stairs','LineWidth',1.3,'DisplayName','corrected');
  polarscatter(stats.(type).phi,stats.(type).R,75,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor',myColors(2))
  rticks([0.2,0.4]), rlim([0,0.5]), legend(h(3:-1:1))
  clear half_bin_width bin_edges

end