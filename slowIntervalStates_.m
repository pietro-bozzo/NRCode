function [slow_f,cond_f,ref_f,intervals,animal] = slowIntervalStates_(session,states,events,opt)
% slowIntervals_ Load slow-rythm intervals and separate them by state and/or event
%
% arguments:
%     session           string, path of session xml file
%     states            (n_states,1) string, brain states
%     events            (n_events,1) string = [], experiment events
%
% name-value arguments:
%     event_restrict    (n_event_r,1) string = [], session epochs to restrict analysis to
%
% output:
%     slow_f            struct, fraction of slow-rythm time spent in a certain condition
%     cond_f            struct, fraction of condition spent in slow-rythm intervals
%     ref_f             struct, fraction of recording spent in condition

% Copyright (C) 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

arguments
  session (1,:) char
  states (:,1) string
  events (:,1) string = []
  opt.event_restrict (:,1) string = []
end

% regions handler
R = regions(session,'states',states,'events','InfraSlowRhythm/slownr','load_spikes',false);
animal = string(R.basename{1}(1:6));
session_duration = sum(diff(R.eventIntervals(opt.event_restrict),1,2));

% condition time spent in recording
for event = events'
  cond_int.(event) = R.eventIntervals(event,opt.event_restrict);
  ref_f.(event) = dur(cond_int.(event)) / session_duration;
end
for state = states'
  cond_int.(state) = R.eventIntervals(state,opt.event_restrict);
  ref_f.(state) = dur(cond_int.(state)) / session_duration;
end

try
  slow_intervals = R.eventIntervals('slownr',opt.event_restrict);
  slow_dur = dur(slow_intervals);
  if all(isnan(slow_intervals),'all')
    slow_intervals = [];
  end
  for condition = [events;states]'
    intervals.(condition) = IntersectIntervals(slow_intervals,cond_int.(condition));
    slow_f.(condition) = dur(intervals.(condition)) / slow_dur;
    cond_f.(condition) = dur(intervals.(condition)) / dur(cond_int.(condition));
  end
catch
  for condition = [events;states]'
    intervals.(condition) = [NaN,NaN];
    slow_f.(condition) = NaN;
    cond_f.(condition) = NaN;
  end
end

end

% helper function to get total duration from intervals
function dur = dur(intervals)
  dur = sum(diff(intervals,1,2));
  dur(isnan(dur)) = 0;
end