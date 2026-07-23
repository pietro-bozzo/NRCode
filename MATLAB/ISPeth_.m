function [peth_on_off,peth_off_on,t] = ISPeth_(session,regs,events,opt)

arguments
  session (1,:) char
  regs (:,1) string
  events (:,1) string = []
  opt.shuffle (1,1) {mustBeNumeric,mustBeInteger,mustBeNonnegative} = 0 % IMPLEMENT
  opt.save (1,1) {mustBeLogical} = false
  opt.load (1,1) {mustBeLogical} = true
  opt.verbose (1,1) {mustBeLogical} = false
end

R = regions(session,'events',["InfraSlowRhythm/slownr","InfraSlowRhythm/slowavalnr"],'verbose',false);

% declare variables
n_bins = 201;
for r = regs'
  peth_on_off.(r) = nan(1,n_bins);
  peth_off_on.(r) = nan(1,n_bins);
end

% keep only regions found in data
[~,ind] = ismember(R.ids,regs);
ind = ind(ind~=0);
regs = regs(ind);

% slow-rhythm intervals
us_intervals = R.eventIntervals('slownr',events);
us_avals = R.eventIntervals('slowavalnr',events);
% transitions time stamps
[~,ind] = Restrict(us_avals(:,2)+1e-5,us_intervals,'verbose','off');
transition.on_off = us_avals(ind,2);
[~,ind] = Restrict(us_avals(:,1)-1e-5,us_intervals,'verbose','off');
transition.off_on = us_avals(ind,1);

% if no rhythm is detected
if all(isnan(us_intervals),'all') || isempty(us_intervals)
  t = nan(1,n_bins);
  return
end

for j = 1 : numel(regs)

  % region data
  spikes = R.spikes('all',regs(j));
  spikes = Restrict(spikes,us_intervals);
  reg_units = R.units(regs(j));
  unit_ind = ismember(reg_units,unique(spikes(:,2)));
  
  [~,t,P] = PETH(spikes(:,1),transition.on_off,'durations',[-5,5],'nBins',n_bins,'fast','on','group',spikes(:,2),'smooth',2);
  peth_on_off.(regs(j)) = nan(numel(reg_units),n_bins);
  peth_on_off.(regs(j))(unit_ind,:) = vertcat(P{:});

  [~,~,P] = PETH(spikes(:,1),transition.off_on,'durations',[-5,5],'nBins',n_bins,'fast','on','group',spikes(:,2),'smooth',2);
  peth_off_on.(regs(j)) = nan(numel(reg_units),n_bins);
  peth_off_on.(regs(j))(unit_ind,:) = vertcat(P{:});

  % overall statistics
  % stat.(regs(j)).h = holmBonferroni(stat.(regs(j)).p);
  % to_log = stat.(regs(j)).h(~isnan(stat.(regs(j)).h));

end

end

% --- Extra code to plot examples in debug mode ---

function plotInDebug()
  % this function is not meant to be called, rather its code can be executed in debug mode to produce plots

  % load entrainment info
  [filebase,basename] = fileparts(session);
  fname = fullfile(filebase,"InfraSlowRhythm",basename+"_unitEntrainment.mat");
  s = load(fname,'stats');

  %% PETH of a unit vs transitions
  i = 18;
  makeFigure('unit',"PETH, unit "+string(i)+", "+regs(j)+", h: "+s.stat.(regs(j)).h(i)+", s: "+s.stat.(regs(j)).phi(i));
  plot(t,peth_on_off.(regs(j))(i,:))
  hold on
  plot(t,peth_on_off.(regs(j))(i,:))

  %% average PETH by population PROBLEM WHEN ONE UNIT NEVER SPIKES AND THUS DISAPPEARS FROM group in PETH
  j = 2;
  color_ind = [2,3,1];
  makeFigure('unit',"Average PETH, "+regs(j));
  is_entrained = s.stat.(regs(j)).h;
  is_entrained(isnan(is_entrained)) = 0;
  is_entrained = logical(is_entrained);
  h(1) = plot(t,mean(zscore(peth_on_off.(regs(j))(is_entrained,:),0,2),1,'omitmissing'),'Color',paperColors(color_ind(j)),'DisplayName','ON -> OFF');
  plot(t,mean(zscore(peth_on_off.(regs(j))(~is_entrained,:),0,2),1,'omitmissing'),'Color',.7*[1,1,1])
  % ADD SAME FOR OTHER transition
  h(2) = plot(t,mean(zscore(peth_off_on.(regs(j))(is_entrained,:),0,2),1,'omitmissing'),'Color',mean([paperColors(color_ind(j));1,1,1]),'DisplayName','OFF -> ON');
  plot(t,mean(zscore(peth_off_on.(regs(j))(~is_entrained,:),0,2),1,'omitmissing'),'--','Color',.5*[1,1,1])
  xline(0,'--')
  legend(h);
  clearvars h is_entrained color_ind

end