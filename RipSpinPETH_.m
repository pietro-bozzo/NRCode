function [peth,trans,time_bins] = RipSpinPETH_(session,opt)
% RipSpinPETH_ Ripple-spindle PETHs wrt ISR

arguments
  session (1,:) char
  opt.limits (2,2) {mustBeNumeric} = [-5,5;-5,5]
  opt.shuffle (1,1) {mustBeNumeric,mustBeInteger,mustBeNonnegative} = 0
end

% load slow-rhythm intervals
R = regions(session,'events',["InfraSlowRhythm/slownr","InfraSlowRhythm/slowavalnr","ripples","spindles"],'load_spikes',false,'verbose',false);
us_intervals = R.eventIntervals('slownr');
us_on = R.eventIntervals('slowavalnr');
us_off = SubtractIntervals(us_intervals,us_on(:,1:2));

% load ripples and spindles
ripples = R.eventInfo("ripples");
ripples = ripples(:,2); % peak time
spindles = R.eventInfo("spindles");
spindles = spindles(:,2); % peak time

% plot PETH, note how swapping first two args of PETH produces mirrored results
% makeFigure('peth');
% PETH(ripples,spindles,'nbins',n_bins,'duration',limits,'smooth',1.3,'Color',myColors(1))
% xline(0,'--')
% plot also at smaller scale, note nesting of ripples inside spindle cycles
% PETH(ripples,spindles,'nbins',n_bins,'duration',[-1,1],'smooth',1.5,'Color',myColors(1))
% xline(0,'--')

n_bins = round(diff(opt.limits,1,2) / 0.04); % 40 ms time bin
peth_f = @(x,y) PETH(x,y,'nbins',n_bins(1),'duration',opt.limits(1,:),'smooth',1.5);
peth_trans = @(x,y) PETH(x,y,'nbins',n_bins(2),'duration',opt.limits(2,:),'smooth',1.5);

% 1. reference PETHs
[~,time_bins.ref,peth.ref] = peth_f(ripples,spindles);
[~,~,peth.inverse] = peth_f(spindles,ripples);

% restrict to US ON intervals
[~,~,int_ind] = Restrict(us_on(:,1),us_intervals); % int_ind(i) is index of us interval ind for i-th ON
% ON which have OFF after
on_for_cycle = us_on(int_ind(1:end-1) >= int_ind(2:end),:);
on_off_transition = on_for_cycle(:,2);
on_off_half = [mean(on_for_cycle,2),mean(us_off,2)];
% ON which have OFF before
on_for_cycle = us_on([false;int_ind(1:end-1) >= int_ind(2:end)],:);
off_on_transition = on_for_cycle(:,1);
off_on_half = [mean(us_off,2),mean(on_for_cycle,2)];

% 2. PETH at transitions
[~,time_bins.trans,trans.rip_on_off] = peth_trans(ripples,on_off_transition);
[~,~,trans.rip_off_on] = peth_trans(ripples,off_on_transition);
[~,~,trans.spin_on_off] = peth_trans(spindles,on_off_transition);
[~,~,trans.spin_off_on] = peth_trans(spindles,off_on_transition);

% keep intervals longer than 2 s
us_on = us_on(diff(us_on(:,1:2),1,2) > 2,1:2);
us_off = us_off(diff(us_off(:,1:2),1,2) > 2,:);
on_off_half = on_off_half(diff(on_off_half(:,1:2),1,2) > 2,:);
off_on_half = off_on_half(diff(off_on_half(:,1:2),1,2) > 2,:);
% shorten them by 2 s
isr_short = us_intervals + [0.5,-0.5];
us_on_short = us_on + [0.5,-0.5];
us_off_short = us_off + [0.5,-0.5];
on_off_half_short = on_off_half + [0.5,-0.5];
off_on_half_short = off_on_half + [0.5,-0.5];
% keep spindles in them
[spin_isr,isr_ind] = Restrict(spindles,isr_short);
is_nisr = true(size(spindles));
is_nisr(isr_ind) = false;
spin_nisr = spindles(is_nisr);
spin_on = Restrict(spindles,us_on_short);
spin_off = Restrict(spindles,us_off_short);
spin_on_off = Restrict(spindles,on_off_half_short);
spin_off_on = Restrict(spindles,off_on_half_short);

% 3. PETH in specific portions of cycle
[~,~,peth.isr] = peth_f(ripples,spin_isr);
[~,~,peth.nisr] = peth_f(ripples,spin_nisr);
[~,~,peth.on] = peth_f(ripples,spin_on);
[~,~,peth.off] = peth_f(ripples,spin_off);
[~,~,peth.on_off] = peth_f(ripples,spin_on_off);
[~,~,peth.off_on] = peth_f(ripples,spin_off_on);