function [perf,data] = PacPerformance_(session,opt)
% PacPerformance_ Get rat performance in PacMaze task
%
% arguments:
%     session    string, path to session xml file
%
% output:
%     perf        struct, DESCRIBE
%     data        struct

% Copyright (C) 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

arguments
  session (1,:) char
  opt.week (1,1) {mustBeNumeric} = NaN
  opt.day (1,1) {mustBeNumeric} = NaN
end

R = regions(session,'events',["trials","successtrials"],'load_spikes',false); %,"successdirecttrials" % doesn't include return time

for task = ["tachem","tachea"]

  trials = R.eventIntervals('trials',task);
  success = R.eventIntervals('successtrials',task);

  [A,is_succ] = ismember(success,trials);

  % check raw data
  if ~all(A,'all') || any(is_succ(:,1) ~= is_succ(:,2)-size(trials,1))
    disp(session)
  end

  perf.(task) = Unfind(is_succ(:,1),size(trials,1));

end

data.sleep = ismember('sleepn',R.phase.names);
data.week = opt.week;
data.day = opt.day;