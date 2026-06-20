function [intervals,slow_avals] = slowIntervals_(session,regs,labels,opt)
% slowIntervals_ Find and save slow-oscillation intervals in spiking activity
%
% arguments:
%     session       string, 
%     regs          (n_regions,1) double, region ids
%     labels        (n_regions,1) string, region labels
%
% name-value arguments:
%     window        double = 0.05, time bin in s for avalanche computation
%     smooth        double = 25, gaussian kernel std in number of samples, default is no smoothing
%     threshold     double = 0.025, population-firing-rate threshold for avalanche computation, in #units / (s * #total_units)
%     save          logical = false, if true, save slow intervals
%     load          logical = true, if true, load slow intervals and corresponding avalanches
%
% output:
%     intervals     (n_intervals,2) double, each row is [start,stop] of a slow-rythm interval
%     slow_avals    (n_avals,2) double, each row is [start,stop] of an avalanche within a slow-rythm interval

% Copyright (C) 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

arguments
  session (1,:) char
  regs (:,1) {mustBeNumeric,mustBeInteger,mustBePositive}
  labels (:,1) string
  opt.window (1,1) {mustBeNumeric,mustBePositive} = 0.05
  opt.smooth (1,1) {mustBeNumeric,mustBePositive} = 25
  opt.threshold (1,1) {mustBeNumeric,mustBePositive} = 0.025
  opt.save (1,1) {mustBeLogical} = false
  opt.load (1,1) {mustBeLogical} = true
end

if numel(regs) ~= numel(labels)
  error('slowIntervals_:labelsSize','Arguments ''regs'' and ''labels'' must have the same number of elements')
end

% try loading intervals
name_parts = split(session,'.');
file_root = strjoin(name_parts(1:end-1),'.');
extensions = join([repmat(".slow",size(labels)),labels],'');
extensions_aval = join([repmat(".slowaval",size(labels)),labels],'');
found = false;
for i = 1 : numel(labels)
  if opt.load
    try
      intervals.(labels(i)) = readmatrix(file_root+extensions(i),FileType='text',CommentStyle='%');
      slow_avals.(labels(i)) = readmatrix(file_root+extensions_aval(i),FileType='text',CommentStyle='%');
      found = true;
    catch
      intervals.(labels(i)) = [NaN,NaN];
      slow_avals.(labels(i)) = [NaN,NaN];
    end
  else
    slow_avals.(labels(i)) = [NaN,NaN];
  end
end
% if something was loaded, return
if found
  return
end

% instantiate handler, load spikes
R = regions(session,verbose=false);

% keep only regions found in data
labels = labels(ismember(regs,R.ids));
regs = regs(ismember(regs,R.ids));
extensions = join([repmat(".slow",size(labels)),labels],'');

% avalanches
R = R.computeAvalanches(opt.window,opt.smooth,opt.threshold,perc=false,mode="ratio");

% parameters
defrag_time = 0.25; % s
dur_thresh = 20; % s
interval_stop = 10; % s

for i = 1 : numel(regs)
  aval_intervals = R.avalIntervals('all',regs(i));

  % 1. defragment avals
  frag_ind = aval_intervals(2:end,1) - aval_intervals(1:end-1,2) >= defrag_time;
  aval_intervals1 = [aval_intervals([true;frag_ind],1),aval_intervals([frag_ind;true],2)];

  % 2. duration criterium
  dur_ind = diff(aval_intervals1,1,2) < dur_thresh;
  aval_intervals2 = aval_intervals1(dur_ind,:);

  % 3. slow-rythm intervals
  if isempty(aval_intervals2)
    intervals.(labels(i)) = [NaN,NaN];
    slow_avals.(labels(i)) = [NaN,NaN];
  else
    deltas = aval_intervals2(2:end,1) - aval_intervals2(1:end-1,2);
    int_ind = deltas > interval_stop;
    n_avals = accumarray(cumsum([true;int_ind]),1); % n of avalanches in every slow-rythm interval
    intervals.(labels(i)) = [aval_intervals2([true;int_ind],1),aval_intervals2([int_ind;true],2)];
    % keep only intervals with at least two avalanches
    intervals.(labels(i)) = intervals.(labels(i))(n_avals>1,:);
    % keep only avalanches inside slow intervals
    slow_avals.(labels(i)) = Restrict(aval_intervals2,intervals.(labels(i)));
  end

  % save slow-rythm intervals
  if opt.save
    saveMatrix(intervals.(labels(i)),file_root + extensions(i),'beginnig','end of slow-rythm interval')
  end
end

return

% -- Extra code to plot examples in debug mode --

i = 3;
slow_dur.(labels(i)) = sum(intervals.(labels(i))(:,2) - intervals.(labels(i))(:,1));

R.plotFiringRates(0,0,opt.window,regions=regs(i),smooth=opt.smooth,mode="ratio");
h = findobj(gca,'Type','Line');
y_max = max(get(h(1), 'YData'));
title("Slow-rythm identification, "+R.printBasename()+', '+labels(i)+' (n: '+num2str(R.nNeurons(regs(i)))+'), w: '+num2str(opt.window)+' s, s: '+num2str(opt.smooth)+ ...
  ', t: '+num2str(opt.threshold)+', T: '+num2str(slow_dur.(labels(i)))+' s')
adjustAxes(gca,'YTickMode','auto','YTickLabelMode','auto','YLim',[0,y_max*1.2]);

% plot slow-rythm intervals
%PlotIntervals(aval_intervals,'color',[0,0,0],'legend','off','alpha',0.75)
%PlotIntervals(aval_intervals1,'color',[0.5,0.5,1],'legend','off','alpha',0.8)
PlotIntervals(aval_intervals2,'color',[0,0.6,0.6],'legend','off','alpha',0.2)
%yline(opt.threshold,HandleVisibility='off')
PlotIntervals(intervals.(labels(i)),'legend','off')