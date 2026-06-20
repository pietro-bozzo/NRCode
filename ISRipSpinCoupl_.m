function [phases,is_coupl,not_us,us_phase,phases_sh,is_coupl_sh,animal] = ISRipSpinCoupl_(session,window,opt)
% ISRipSpinCoupl_ Analize ripple-spindle coupling w.r.t. ISR

arguments
  session (1,1) string
  window (1,1) {mustBeNumeric}
  opt.event (:,1) string = []
  opt.state (1,1) string = "all"
  opt.duration (1,1) {mustBeNumeric,mustBeNonnegative} = 0
  opt.shuffle (1,1) {mustBeNumeric,mustBeInteger,mustBeNonnegative} = 0
end

% load slow-rhythm intervals
R = regions(session,'states',opt.state,'events',["InfraSlowRhythm/slownr","InfraSlowRhythm/slowavalnr","ripples","spindles"],'load_spikes',false,'verbose',false);
animal = extractBefore(R.basename,7);
isr_intervals = R.eventIntervals('slownr');
isr_avals = R.eventIntervals('slowavalnr');

% remove extra ON intervals
[on,off,isr_intervals] = cleanISROnOff(isr_avals,isr_intervals);
phase = interpolatePhase(0.005,unique([on(:,1),off(:,2)]),0,unique([on(:,2),off(:,1)]),pi);
phase = Restrict(phase,isr_intervals);

[basepath,basename] = fileparts(session);
names = ["ripples","deltas","spindles"];
for n = names([1,3]) % load ripples and spindles
  try
    times.(n) = readmatrix(fullfile(basepath,basename+"."+n),'FileType','text','CommentStyle','%');
    times.(n) = times.(n)(:,2); % peak time
  catch
    times.(n) = [];
  end
end
try % load deltas
  load(fullfile(basepath,basename+".deltaWaves.events.mat"))
  times.deltas = deltaWaves.peaks;
catch
  times.deltas = [];
end

% identify coupled events, where couples are ripples preceding spindles by less than 'window' s   ADD DELTAS
[is_coupl.ripples,is_coupl.spindles] = isCoupled(times.ripples,times.spindles,window);
is_coupl.deltas = false(size(times.deltas));

% restrict in state and event (after identifying couples!)
restrict_intervals = R.eventIntervals(opt.state,opt.event,'regexp',true,'duration',opt.duration);
for n = names
  [times.(n),ind] = Restrict(times.(n),restrict_intervals,'verbose','off');
  is_coupl.(n) = is_coupl.(n)(ind); 
end

% analyze coupling
for n = names

  % restrict in US intervals
  [times.(n),ind] = Restrict(times.(n),isr_intervals,'verbose','off');
  not_us.(n)(1) = numel(is_coupl.(n)) - numel(ind); % not_us.(n)(1) is number of "n" events outside US rhythm
  coupl_not_us = is_coupl.(n);
  coupl_not_us(ind) = false;
  not_us.(n)(2) = sum(coupl_not_us); % not_us.(n)(2) is number of "n" events outside US rhythm coupled to "non-n" events
  is_coupl.(n) = is_coupl.(n)(ind); % keep only events inside US for the rest

  % phase of each event
  phases.(n) = interp1(phase(:,1),phase(:,2),times.(n));

  % declare variables for surrogate dta
  times_shifted.(n) = Restrict(times.(n),isr_intervals,'shift','on','verbose','off');
  phases_sh.(n) = zeros(numel(times.(n)),opt.shuffle);
  is_coupl_sh.(n) = false(numel(times.(n)),opt.shuffle);

end

% shuffle event times inside ISR to get phases and coupling
for i = 1 : opt.shuffle
  for n = names
    times_shuffled.(n) = Unshift(shuffleSpikes(times_shifted.(n)),isr_intervals);
    phases_sh.(n)(:,i) = interp1(phase(:,1),phase(:,2),times_shuffled.(n));
  end
  [is_coupl_sh.ripples(:,i),is_coupl_sh.spindles(:,i)] = isCoupled(times_shuffled.ripples,times_shuffled.spindles,window);
end

% return all ISR phase values to build overall distribution
us_phase = phase(:,2);

end

% --- helper functions ---

function [is_coupl_rip,is_coupl_spin] = isCoupled(ripple_t,spindle_t,window)
% identify coupled events, where couples are ripples preceding spindles by less than 'window' s

  delays = spindle_t - ripple_t.'; % each row is a spindle, each column a ripple
  delays(delays<0) = NaN;
  distan = min(delays,[],1);
  is_coupl_rip = (distan <= window).';
  distan = min(delays,[],2);
  is_coupl_spin = distan <= window;

  % code to find couples using a window centered on "n" event
  % for n = names
  % ind = knnsearch(times.(names(names~=n)),times.(n)); % ind(i) is index of closest "non-n" event for i-th "n" event
  % distan = abs(times.(n)-times.(names(names~=n))(ind));
  % is_coupl.(n) = distan <= window;
  % end

end

% --- Extra code to plot examples in debug mode ---

function plotInDebug()
  % this function is not meant to be called, rather its code can be executed in debug mode to produce plots

  phase = phase; % to avoid MATLAB compilation error due to unknown variable
  n_bins = 250; % high number to estimate pdf correction
  [dist_phase,bins_phase] = CircularDistribution(phase(:,2),'nBins',n_bins,'normalize','pdf'); % phase was linearly produced: bi-uniform distr

  %%
  R = regions(session,load_spikes=false,verbose=false);

  %% plot phase
  makeFigure('phase',"US phase, "+R.printBasename());
  PlotXY(phase)
  PlotIntervals(isr_intervals);
  PlotIntervals(isr_avals);

  %% plot events
  makeFigure('events',"Events, "+R.printBasename());
  PlotIntervals(isr_avals(:,1:2),'color',[1,1,0.7],'alpha',1);
  PlotIntervals(isr_intervals,'alpha',1);
  % A = [rip(:,[1,3]),nan(size(rip,1),1)].';
  % plot(A(:),0.55*ones(numel(A),1),'Color',myColors(1));
  A = [times.spindles-window,times.spindles+window,nan(size(times.spindles))].';
  plot(A(:),0.55*ones(numel(A),1),'Color',[0.7,0.7,0.7]);
  B = [raw(:,[1,3]),nan(size(raw,1),1)].';
  plot(B(:),0.5*ones(numel(B),1),'Color',myColors(2));
  RasterPlot([times.ripples,0.5*ones(size(times.ripples))],'Color',myColors(1)) % ripples are almost istantaneous in this context
  RasterPlot([times.spindles,0.5*ones(size(times.spindles))],'Color',myColors(2))

  %% plot events to check coupling detection
  makeFigure('events',"Events, "+R.printBasename());
  RasterPlot([times.ripples(is_coupl.ripples),0.5*ones(size(times.ripples(is_coupl.ripples)))],'Color',paperColors(10));
  RasterPlot([times.ripples(~is_coupl.ripples),0.5*ones(size(times.ripples(~is_coupl.ripples)))],'Color',0.4*[1,0,0]);
  RasterPlot([times.spindles(is_coupl.spindles),0.5*ones(size(times.spindles(is_coupl.spindles)))],'Color',paperColors(11));
  RasterPlot([times.spindles(~is_coupl.spindles),0.5*ones(size(times.spindles(~is_coupl.spindles)))],'Color',0.7*[1,1,1]);

  %% PETH
  makeFigure('peth');
  ripples = times.ripples;
  spindles = times.spindles;
  PETH(ripples,spindles,'Color',myColors(2))
  

  [~,ind_us_rip] = Restrict(times.ripples,isr_intervals);

  PETH(times.ripples(is_coupl_rip),times.spindles(is_coupl_spin))
  PETH(times.ripples(~is_coupl_rip),times.spindles(~is_coupl_spin))

  %% plot phase distribution for ripples
  makeFigure('rip',R.printBasename+", "+opt.type+" phase distribution");
  plot(bins_phase,dist_phase,'Color',[0.7,0.7,0.7],'DisplayName','US phase')
  plot(bins_phase,distr_orig,'Color',myColors(1),'DisplayName','ripple phase')
  plot(bins_phase,distr,'Color',myColors(2),'DisplayName','corrected ripple phase')
  xlim([0,2*pi]), xticks([0,pi/2,pi,3*pi/2,2*pi]), xticklabels(["0","π/2",'π','3π/2','2π']), xlabel('phase (rad)'), ylabel('p(phase)'), legend

end