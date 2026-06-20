function [activ_rate,weights,neurons] = reuAsmb_(session,opt)
%reuAsmb_ Summary of this function goes here

arguments
  session (1,:) char
  opt.event (:,1) string = "tachem"    % allows to require specific events, e.g., "fcBOX"
  opt.state (1,1) string = "all"
  opt.regions (:,1) double = []
  opt.ICA_window (1,1) double {mustBePositive} = 0.03
  opt.plot (1,1) {mustBeLogical} = false
  opt.verbose (1,1) {mustBeLogical} = false
end

% load protocol phases
[event_names,stamps] = loadEvents(session);
s1_ind = find(event_names=='sleepm');
task_ind = find(event_names==opt.event);
s2_ind = find(event_names=='sleepn');
if numel(task_ind) ~= 1
  error(append('Multiple' ,opt.event,' found.'))
end

% instantiate handler
R = regions(session,events=stamps([s1_ind;task_ind;s2_ind]),states=opt.state,regions=opt.regions);
R = R.loadSpikes();
spikes = R.spikes();

% compute ICA over all neurons
weights = ICAssemblies(spikes,opt.ICA_window,restrict=stamps{task_ind}); % detect assemblies over requested time
% compute assemblies' activations over all spikes
ICs_activations = ICActivations(spikes,weights,opt.ICA_window);

% get sws stamps in thw two sleep sessions
sws_stamps = R.state_stamps{1};
s1_stamps = [];
for ind = s1_ind.'
  s1_stamps = [s1_stamps;sws_stamps(sws_stamps(:,2) > stamps{ind}(1) & sws_stamps(:,2) < stamps{ind}(2),:)];
end
s2_stamps = [];
for ind = s2_ind.'
  s2_stamps = [s2_stamps;sws_stamps(sws_stamps(:,2) > stamps{ind}(1) & sws_stamps(:,2) < stamps{ind}(2),:)];
end
s1_ICs_activations = Restrict(ICs_activations,s1_stamps);
s2_ICs_activations = Restrict(ICs_activations,s2_stamps);

% plot rasters as sanity check
if opt.plot
  [~,basename] = fileparts(session);
  basename = strrep(basename,'_','-');
  ICs_activations = Restrict(ICs_activations,sws_stamps);
  start = 0; stop = 0;
  figure(Name=append('raster'),NumberTitle='off',Position=get(0,'Screensize')); hold on
  times = ICs_activations(:,1);
  if stop <= 0
    stop = -times(end) - 0.001;
  end
  ICs = ICs_activations(times > start & times < abs(stop),2);
  times = times(times > start & times < abs(stop));
  raster([times,ICs],'color',myColors(1))
  % adjust plot
  set(gca,TickDir='out',XLim=[start;abs(stop)],YLim=[0,size(weights,2)]); xlabel('time (s)',FontSize=14); ylabel('ICs',FontSize=14);
  xline([stamps{task_ind}(1,1),stamps{task_ind}(1,2)],'r')
  title(basename,FontSize=17,FontWeight='normal')
end

% count activations
[n_activ1_nnz,unique_ICs1] = groupcounts(s1_ICs_activations(:,2));
n_activ1_nnz = n_activ1_nnz / sum(s1_stamps(:,2)-s1_stamps(:,1));
[n_activ2_nnz,unique_ICs2] = groupcounts(s2_ICs_activations(:,2));
n_activ2_nnz = n_activ2_nnz / sum(s2_stamps(:,2)-s2_stamps(:,1));
% handle ICs which never activate
activ_rate = zeros(size(weights,2),2);
activ_rate(unique_ICs1,1) = n_activ1_nnz;
activ_rate(unique_ICs2,2) = n_activ2_nnz;

% get neurons per region
neurons = {};
for i = 1 : numel(R.regions_array)
  neurons{i,1} = R.regions_array(i).neurons;
end