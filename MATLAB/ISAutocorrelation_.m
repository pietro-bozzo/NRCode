function [firing_on,firing_off,peth] = ISAutocorrelation_(session,regs,events,opt)

arguments
  session (1,:) char
  regs (:,1) string
  events (:,1) string = "all"
  opt.window (1,1) {mustBeNumeric,mustBePositive} = 0.05
  opt.step (1,1) {mustBeNumeric,mustBeInteger,mustBePositive} = 5
  opt.smooth (1,1) {mustBeNumeric,mustBePositive} = 100
end

% load data
R = regions(session,'events',["InfraSlowRhythm/slownr","InfraSlowRhythm/slowavalnr"],'verbose',false);

% keep only regions found in data
is_reg = ismember(regs,R.ids);

% compute firing rate as for ISR detection
[FR,time] = R.firingRate('all',regs(is_reg),'window',opt.window,'smooth',opt.smooth,'step',opt.step,'mode','fr_norm');

r_count = 1;
for r = 1 : numel(regs)
  if is_reg(r)
    [autoc.(regs(r)),lags] = autocorr(FR(:,r_count));
  else
    [autoc.(regs(r)),lags] = deal(NaN);
  end
end