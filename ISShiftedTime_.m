function [IS_fraction,time] = ISShiftedTime_(session,states,target,events,opt)

arguments
  session (1,:) char
  states (:,1) string
  target (1,1) string % target state
  events (:,1) string
  opt.window (1,1) {mustBeNumeric} = 45
end

R = regions(session,'states',states,'events','InfraSlowRhythm/slownr','load_spikes',false);

% special case: regular expression
if isscalar(events) && contains(events,"*")
  ind = cellfun(@(x) ~isempty(x) && x==1, regexp(R.event_names,events));
  events = R.event_names(ind);
end

slow_intervals = R.eventIntervals('slownr');

% parameters
N = 200;
windows = opt.window * [0 : N-1; 1 : N].';
for i = 1 : numel(events)

  state_intervals = R.eventIntervals(target,events(i));
  isr_state = IntersectIntervals(slow_intervals,state_intervals);

  if isempty(isr_state)
    IS_fraction.(events(i)) = zeros(N,1);
  else
    % shift wrt state
    isr_state_shifted = [Restrict(isr_state(:,1),state_intervals,'shift','on'),Restrict(isr_state(:,2),state_intervals,'shift','on')];

    [intersection,ind] = IntersectIntervals(windows+[1,-1]*1e-7,isr_state_shifted); % shorten windows slightly to avoid their consolidation inside IntersectIntervals
    IS_fraction.(events(i)) = accumarray(ind,diff(intersection,1,2),[N,1]) / opt.window;
  end

end

time = mean(windows,2);

end

% --- Extra code to plot examples in debug mode ---

function plotInDebug()
  % this function is not meant to be called, rather its code can be executed in debug mode to produce plots
  R = regions(session,states=states);
  R.plotFiringRates([],[],regions=70,step=5,smooth=100);

  figure
  PlotIntervals(R.eventIntervals(events(1)),'bars')
  PlotIntervals(windows(1:2:end,:))
  PlotIntervals(slow_intervals,'color',[1,1,0])
  PlotIntervals(R.stateIntervals(target),'color',[1,0,0])
end