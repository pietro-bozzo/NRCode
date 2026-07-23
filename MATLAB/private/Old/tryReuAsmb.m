session = {'/mnt/hubel-data-131/perceval/Rat003_20231212/Rat003_20231212.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231213/Rat003_20231213.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231214/Rat003_20231214.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231215/Rat003_20231215.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231217/Rat003_20231217.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231218/Rat003_20231218.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231219/Rat003_20231219.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231221/Rat003_20231221.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231222/Rat003_20231222.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231223/Rat003_20231223.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231224/Rat003_20231224.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231226/Rat003_20231226.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231227/Rat003_20231227.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231228/Rat003_20231228.xml';
           '/mnt/hubel-data-131/perceval/Rat003_20231229/Rat003_20231229.xml'};
% choose session
i = 4;
session = session{i};
[~,basename] = fileparts(session);
phase = "all";
states = "all";
regs = [12;32;71;0];
% load protocol phases
[event_names,stamps] = loadEvents(session);
event = 'tachem';
s1_ind = find(event_names=='sleepm');
task_ind = find(event_names==event);
s2_ind = find(event_names=='sleepn');
if numel(task_ind) ~= 1
  error(append('Multiple' ,event,' found.'))
end
% instantiate handler
R = regions(session,phases=stamps([s1_ind;task_ind;s2_ind]),states=states,regions=regs);
R = R.loadSpikes();
spikes = R.regions_array(end).spikes;
first_sleep = [11250,12000];

%% compute ICA over all neurons
ICA_window = 0.030;
[weights,~,ICs_activity,IC_time,ICs_activations] = getICActivity(spikes,ICA_window,restrict=stamps{task_ind});

%% get sws
Rs = regions(session,phases=stamps([s1_ind;task_ind;s2_ind]),states="sws",regions=regs);
sws_stamps = Rs.state_stamps{1};
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

%% plot rasters
ICs_activations = Restrict(ICs_activations,sws_stamps);
start = 0; stop = 0;
fig = figure(Name=append('raster'),NumberTitle='off',Position=get(0,'Screensize')); hold on
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

%% count activations
[n_activ1,unique_ICs1] = groupcounts(s1_ICs_activations(:,2));
n_activ1 = n_activ1 / sum(s1_stamps(:,2)-s1_stamps(:,1));
[n_activ2,unique_ICs2] = groupcounts(s2_ICs_activations(:,2));
n_activ2 = n_activ2 / sum(s2_stamps(:,2)-s2_stamps(:,1));
if numel(unique_ICs1) ~= numel(unique_ICs2) || any(unique_ICs1 ~= unique_ICs2)
  for j = 1 : min([numel(unique_ICs1),numel(unique_ICs2)])
    if unique_ICs1(j) > unique_ICs2(j)
      unique_ICs1 = [unique_ICs1(1:j-1);unique_ICs2(j);unique_ICs1(j:end)];
      n_activ1 = [n_activ1(1:j-1);0;n_activ1(j:end)];
    elseif unique_ICs1(j) < unique_ICs2(j)
      unique_ICs2 = [unique_ICs2(1:j-1);unique_ICs1(j);unique_ICs2(j:end)];
      n_activ2 = [n_activ2(1:j-1);0;n_activ2(j:end)];
    end
  end
end
%result{i} = n_activ2 > n_activ1;
%result_fraction(i) = mean(result{i});