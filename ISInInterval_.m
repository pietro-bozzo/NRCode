function [IS_fraction,time] = ISInInterval_(session,states,target,events,opt)

arguments
  session (1,:) char
  states (:,1) string
  target (1,1) string % target state
  events (:,1) string
  opt.n_bins (1,1) {mustBeNumeric} = 50
end

R = regions(session,'states',states,'events','InfraSlowRhythm/slownr','load_spikes',false);

% special case: regular expression
if isscalar(events) && contains(events,"*")
  ind = cellfun(@(x) ~isempty(x) && x==1, regexp(R.event_names,events));
  events = R.event_names(ind);
end

slow_intervals = R.eventIntervals('slownr');

windows = linspace(0,1,opt.n_bins+1);
windows = [0,repelem(windows(2:end-1),1,2),1];
windows = reshape(windows,2,[]).';

for i = 1 : numel(events)

  % ISR in target state
  state_intervals = R.eventIntervals(target,events(i));
  [isr_state,state_ind] = IntersectIntervals(state_intervals,slow_intervals);

  IS_fraction.(events(i)) = zeros(size(state_intervals,1),opt.n_bins);
  for j = 1 : size(state_intervals,1)
    % rescale state interval length
    rescaled = (isr_state(state_ind==j,:)-state_intervals(j,1)) / (state_intervals(j,2)-state_intervals(j,1));
    [rescaled,win_ind] = IntersectIntervals(windows+[1,-1]*1e-7,rescaled); % shorten windows slightly to avoid their consolidation inside IntersectIntervals
    if ~isempty(rescaled)
      IS_fraction.(events(i))(j,:) = accumarray(win_ind,diff(rescaled,1,2),[opt.n_bins,1]) / windows(1,2);
    end
  end

end

time = mean(windows,2);

end