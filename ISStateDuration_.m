function [us_state_f,animal] = ISStateDuration_(session,states,events)
% ISStateDuration_ Get fraction of time spent in NR ISR for every interval of given states
%
% arguments:
%     session           string, path of session xml file
%     states            (n_states,1) string, brain states
%     events            (n_events,1) string = [], session epochs to restrict analysis to
%
% output:
%     us_state_f        struct, DESCRIBE (each row is [fraction,duration])

% Copyright (C) 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

arguments
  session (1,:) char
  states (:,1) string
  events (:,1) string = []
end

R = regions(session,'states',states,'events','InfraSlowRhythm/slownr','load_spikes',false);
animal = string(R.basename{1}(1:6));

% IS intervals in session epochs
isr_intervals = R.eventIntervals('slownr',events);

for j = 1 : numel(states)

  % state intervals in session epochs
  state_intervals = R.stateIntervals(states(j),events);

  % ratio of US time for every state interval
  state_dur = diff(state_intervals,1,2);
  [intersection,int_ind] = IntersectIntervals(state_intervals,isr_intervals); % int_ind(i) is index of state interval containing isr_intervals(i)
  if isempty(intersection)
    us_state_dur = zeros(size(state_dur));
  else
    us_state_dur = accumarray(int_ind,diff(intersection,1,2),size(state_dur));
  end
  us_state_f.(states(j)) = [us_state_dur./state_dur, state_dur];

end