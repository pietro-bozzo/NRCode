%% prepare data for Brain Hack 01/2026

batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_units.batch';
sessions = readBatchFile(batch_file);
sessions = sessions([12,27])';
regs = ["pfc","hpc","nr"];

for s = 1 : numel(sessions)

  data = struct;
  R = regions(sessions(s),'states',["sws","rem"],'events',["ripples","spindles","InfraSlowRhythm/slownr","InfraSlowRhythm/slowavalnr"]);

  for r = regs
    data.("spikes_"+r) = R.spikes('all',r);
  end

  data.ripples = R.eventInfo('ripples');
  data.ripples = data.ripples(:,2);
  data.spindles = R.eventInfo('spindles');
  data.spindles = data.spindles(:,2);
  data.nrem = R.eventIntervals('sws');
  data.rem = R.eventIntervals('rem');
  data.wake = R.eventIntervals('other');
  data.protocol_times = R.eventIntervals();
  data.protocol_names = string(R.phase.names);
  data.nr_isr = R.eventIntervals('slownr');
  data.nr_isa = R.eventIntervals('slowavalnr');

  save("/mnt/hubel-data-103/Pietro/Data/BrainHack/session"+s,"data");

end

clearvars batch_file sessions regs data R s r