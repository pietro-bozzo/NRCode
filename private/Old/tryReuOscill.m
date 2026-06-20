%% choose session
session = { % cell array to allow for easy selection among sessions
'/mnt/hubel-data-131/perceval/Rat003_20231213/Rat003_20231213.xml';
'/mnt/hubel-data-131/perceval/Rat003_20231215/Rat003_20231215.xml';
'/mnt/hubel-data-131/perceval/Rat003_20231221/Rat003_20231221.xml'; % no sleep
'/mnt/hubel-data-131/perceval/Rat003_20231222/Rat003_20231222.xml'; % no sleep
'/mnt/hubel-data-131/perceval/Rat003_20231223/Rat003_20231223.xml'; % no sleep
'/mnt/hubel-data-131/perceval/Rat003_20231224/Rat003_20231224.xml'; % no sleep
'/mnt/hubel-data-139/perceval/Rat003_20231226/Rat003_20231226.xml';
'/mnt/hubel-data-139/perceval/Rat003_20231227/Rat003_20231227.xml';
'/mnt/hubel-data-139/perceval/Rat003_20231228/Rat003_20231228.xml';
'/mnt/hubel-data-139/perceval/Rat003_20231229/Rat003_20231229.xml';
'/mnt/hubel-data-139/karadoc/Rat004_20240226/Rat004_20240226.xml';
'/mnt/hubel-data-139/karadoc/Rat004_20240227/Rat004_20240227.xml';
'/mnt/hubel-data-139/karadoc/Rat004_20240228/Rat004_20240228.xml'; % some neurons don't follow (are they in TH?), also most TH neurons follow, REWATCH!
'/mnt/hubel-data-139/karadoc/Rat004_20240303/Rat004_20240303.xml'; % some neurons don't follow (are they in TH?), also most TH neurons follow, REWATCH!
'/mnt/hubel-data-140/karadoc/Rat004_20240305/Rat004_20240305.xml'; % no sleep
'/mnt/hubel-data-140/karadoc/Rat004_20240306/Rat004_20240306.xml'; % no sleep
'/mnt/hubel-data-140/karadoc/Rat004_20240308/Rat004_20240308.xml'; % no sleep
'/mnt/hubel-data-140/karadoc/Rat004_20240309/Rat004_20240309.xml'; % no sleep, some neurons don't follow
'/mnt/hubel-data-140/karadoc/Rat004_20240310/Rat004_20240310.xml'; % no sleep
'/mnt/hubel-data-140/karadoc/Rat004_20240311/Rat004_20240311.xml'; % no sleep
'/mnt/hubel-data-140/karadoc/Rat004_20240313/Rat004_20240313.xml';
'/mnt/hubel-data-140/karadoc/Rat004_20240314/Rat004_20240314.xml';
'/mnt/hubel-data-140/karadoc/Rat004_20240315/Rat004_20240315.xml';
'/mnt/hubel-data-140/karadoc/Rat004_20240316/Rat004_20240316.xml';
'/mnt/hubel-data-140/karadoc/Rat004_20240317/Rat004_20240317.xml';
'/mnt/hubel-data-140/karadoc/Rat004_20240318/Rat004_20240318.xml';
'/mnt/hubel-data-140/karadoc/Rat004_20240319/Rat004_20240319.xml';
'/mnt/hubel-data-143/peter/Rat005_20240520/Rat005_20240520.xml'; %  weird CC and behavior!
'/mnt/hubel-data-143/peter/Rat005_20240521/Rat005_20240521.xml'; %  weird CC and behavior!
'/mnt/hubel-data-143/peter/Rat005_20240523/Rat005_20240523.xml';
'/mnt/hubel-data-145/peter/Rat005_20240601/Rat005_20240601.xml';
'/mnt/hubel-data-145/peter/Rat005_20240602/Rat005_20240602.xml';
};
session = session{end-1}; % convert to character array
disp(session)

%%
sessions = readBatchFile('/mnt/hubel-data-103/Pietro/Data/BatchFiles/sasuke.batch');
%%
session = sessions{1};
disp(session)

%% choose parameter values
states = ["sws","rem"];
regs = [12;32;70;80];
labels = ["pfc","hpc","nr","th"];
window = 0.05; % smaller than 50 ms catches little spikes, this threshold could be computed analytically from Poisson

%% instantiate handler, load spikes
%R = regions(session,states=states);
R = regions(session);

%% see raster
R.plotSpikeRaster(14250,16000);

%% see firing rate
R.plotFiringRates(window,regions=70,smooth=25);

%% compute avalanches
R = R.computeAvalanches(window,15,1,perc=false);

%% see avalanches
R.plotFiringRates(window,avals=true);


% IDEA PER SHUFFLE:
% PRENDI n NEURONI DALLE TRE STRUTTURE A CASO, MOSTRA CHE FENOMENO NON C'E



%% 1a. example figure
R.plotSpikeRaster(1500,3000,regions=[70,32,12],colors=[0.7,0.7,0.7]);
legend('off')
yline(cumsum(R.nNeurons(regs([3,2])))+0.5,Color=myColors(1),LineWidth=1.7)



%% 1b. avalanches allow to identify slow rhythm intervals
% choose parameter values
aval_window = 0.05; % s
step = 5;
smo = 100;
%defrag_time = 0.25;
dur_thresh = 20; % s
aval_thresh = 5; % percentile
max_thresh = 10; % 25; % percentile; seems like Blinky 25, Kara 10
int_stop = 15; % s
% MAYBE lower off_thresh for Blinky, to detect less false rhythm

% avalanches
%[FR,time] = R.firingRate('all',regs,window=aval_window,smooth=smo,mode="ratio");
%R = R.computeAvalanches(aval_window,smo,aval_thresh,perc=false,mode="ratio");
[FR,time] = R.firingRate('all',window=aval_window,smooth=smo,step=step);
R = R.computeAvalanches(aval_window,smo,aval_thresh,step=step);
max_thresh = prctile(FR,max_thresh);

%%
plot_flag = [false,true,true,false];
disp(' ')
for i = 1 : numel(regs)
  aval_intervals = R.avalIntervals('all',regs(i));

  % 1. duration criteria
  dur_ind = diff(aval_intervals,1,2) < dur_thresh;
  aval_intervals1 = aval_intervals(dur_ind,:);

  % 3. remove avals with small max
  if ~isempty(aval_intervals1)
    [~,time_ind,aval_ind] = Restrict(time,aval_intervals1);
    maxes = accumarray(aval_ind,FR(time_ind,i),[],@max);
    aval_intervals3.(labels(i)) = aval_intervals1(maxes > max_thresh(i),:); %& prctls > avrg_thresh
    % check maxes identification
    % AA = aval_intervals1.';
    % h = plot(AA(:),repelem(maxes,2,1));
  end

  aval_durs.(labels(i)) = diff(aval_intervals3.(labels(i)),1,2);
  
  % % 4. lengthen silences
  % ttt = prctile(FR(:,i),10);
  % t_up_cross = time(FR(1:end-1,i) < ttt & FR(2:end,i) >= ttt);
  % ind = discretize(t_up_cross,aval_intervals3(:,1));
  % % remove t_up_cross wheere ind is nan and remove ind nan
  % t_up_cross = t_up_cross(~isnan(ind));
  % ind = ind(~isnan(ind));
  % [~,a_ind] = unique(ind);
  % aval_intervals3(1:numel(a_ind),1) = t_up_cross(a_ind);
  % t_down_cross = time(FR(1:end-1,i) > ttt & FR(2:end,i) <= ttt);
  % ind = discretize(t_up_cross,aval_intervals3(:,1));

  % 4. slow-rhythm intervals
  if isempty(aval_intervals3)   
    n_avals.(labels(i)) = NaN;
    slow_intervals.(labels(i)) = [NaN,NaN];
    slow_dur.(labels(i)) = 0;
    slow_aval.(labels(i)) = [NaN,NaN];
  else
    % THEN CHECK as SUPPL THAT REMOVING > 10 s silences doesn't separate rhythm intervals, i,e, that many avals are in the middle of rhythm int

    silence_intervals = [aval_intervals3.(labels(i))(1:end-1,2),aval_intervals3.(labels(i))(2:end,1)];
    silence_durs.(labels(i)) = aval_intervals3.(labels(i))(2:end,1) - aval_intervals3.(labels(i))(1:end-1,2);

    %silence_durs.(labels(i)) = silence_durs.(labels(i))(silence_durs.(labels(i)) < int_stop);
    %deltas = aval_intervals3.(labels(i))(2:end,1) - aval_intervals3.(labels(i))(1:end-1,2);
    int_ind = silence_durs.(labels(i)) > int_stop; % true iff a silence ends a slow-rhythm interval
    aval_ind = cumsum([true;int_ind]); % aval_ind(j) is avalanche containing silence_durs(j)
    % n of avalanches in every slow-rhythm interval
    n_avals.(labels(i)) = accumarray(aval_ind,1);
    slow_intervals.(labels(i)) = [aval_intervals3.(labels(i))([true;int_ind],1),aval_intervals3.(labels(i))([int_ind;true],2)];
    % keep only silences inside slow_intervals

    [~,time_ind,silence_ind] = Restrict(time,silence_intervals);
    mins = accumarray(silence_ind,FR(time_ind,i),[],@min);
    % check min identification
    % AA = silence_intervals.';
    % h = plot(AA(:),repelem(mins,2,1));

    ok_silence = accumarray(aval_ind(~int_ind),mins(~int_ind),[],@(x) sum(x<max_thresh(i)/3)/numel(x));
    disp(string(sum(ok_silence>0.33))+' out of '+string(size(slow_intervals.(labels(i)),1)))

% STUFF I'm using to plot
% R.plotFiringRates([],[],aval_window,'regions',regs(3),"step",step,"smooth",smo);
% yline(prctile(FR(:,i),25)); yline(prctile(FR(:,i),5)); yline(prctile(FR(:,i),25)/3);



    % keep only intervals with at least two avalanches
    slow_intervals.(labels(i)) = slow_intervals.(labels(i))(n_avals.(labels(i))>1,:);
    % interval durations
    slow_dur.(labels(i)) = sum(slow_intervals.(labels(i))(:,2) - slow_intervals.(labels(i))(:,1));
    % keep only avalanches inside slow intervals
    slow_aval.(labels(i)) = Restrict(aval_intervals3.(labels(i)),slow_intervals.(labels(i)));

    

    silence_durs.(labels(i)) = slow_aval.(labels(i))(2:end,1) - slow_aval.(labels(i))(1:end-1,2);
    silence_durs.(labels(i)) = silence_durs.(labels(i))(silence_durs.(labels(i)) <= int_stop);
  end

  % NEW IDEA: identify putative rhythm intervals as 4., then remove those
  % where less than half of silences go below max_thresh/3
  % in suppl, a fig with distribution of log silences dur to claim that nr
  % has longer ones

  disp(labels(i)+' '+num2str(size(aval_intervals3,1))+' out of '+num2str(size(aval_intervals3,1))+', T: '+num2str(slow_dur.(labels(i))))
  if labels(i) == "nr"
    nr_avals = aval_intervals3;
  end

  % plot
  if plot_flag(i)
    % plot firing rate used for detection
    R.plotFiringRates(0,0,aval_window,regions=regs(i),smooth=smo,mode="ratio");
    title("Slow-rhythm identification, "+R.printBasename()+', '+labels(i)+' (n: '+num2str(R.nNeurons(regs(i)))+'), w: '+num2str(aval_window)+' s, s: '+num2str(smo)+ ...
      ', t: '+num2str(aval_thresh)+', T: '+num2str(slow_dur.(labels(i)))+' s')
    adjustAxes(gca,'YTickMode','auto','YTickLabelMode','auto'); %'YLim',[0,max(FR(:,i))*1.1],

    % plot slow-rhythm intervals
    %PlotIntervals(aval_intervals,'color',[0,0,0],'legend','off','alpha',0.75)
    %PlotIntervals(aval_intervals1,'color',[0.5,0.5,1],'legend','off','alpha',0.8)
    PlotIntervals(aval_intervals3,'color',[0,0.6,0.6],'legend','off','alpha',0.2)
    %iii = ~ismember(aval_intervals2(:,1),aval_intervals3(:,1));
    %PlotIntervals(aval_intervals2(iii,:),'color',[1,0.8,0.8],'legend','off','alpha',0.5)
    %yline(aval_thresh,HandleVisibility='off')
    %yline([prctile(FR(:,i),aval_thresh),max_thresh(i)],HandleVisibility='off')

    PlotIntervals(slow_intervals.(labels(i)),'legend','off')
  end
end
clear plot_flag i aval_intervals frag_ind aval_intervals1 dur_ind aval_intervals3 deltas int_ind



%% make figure
start = 1800; stop = 2300;
[fig,axs] = makeFigure('identif',"Slow-rhythm identification, "+R.printBasename()+', NR (n: '+num2str(R.nNeurons(70))+'), w: '+num2str(aval_window)+' s, s: '+num2str(smo)+ ...
  ', t: '+num2str(aval_thresh)+', T: '+num2str(slow_dur.nr)+' s',[2,1],TileSpacig='none');
% raster
R.plotSpikeRaster(start,stop,states='all',regions=70,colors=[0.7,0.7,0.7],ax=axs(1));
set(axs(1),'XTick',[],'YTick',[]); legend(axs(1),'off')
% avalanches
R.plotFiringRates(start,stop,window,states='all',regions=70,smooth=smo,mode="ratio",ax=axs(2));
PlotIntervals(nr_avals,'color',[1,0.89,0.82],'legend','off','alpha',1)
PlotIntervals(slow_intervals.nr,'legend','off')
set(axs(2),'YLim',[0,0.25],'YTickMode','auto','YTickLabelMode','auto'); legend(axs(2),'off')



%% 1b. if necessary, recompute avalanches inside slow intervals with different parameters to better identify silences



%% 2. ripples
% load ripples
rip = readmatrix(fileparts(session)+"/"+R.basename+'.ripples',FileType='text');
rip_peak_t = rip(:,2);

% compute avalanches with higher threshold
rip_thresh = 0.045;
R = R.computeAvalanches(aval_window,smo,rip_thresh,perc=false,mode="ratio");
aval_intervals = R.avalIntervals('all',70);
% keep only avals which fall inside slow intervals
[~,ind1] = Restrict(aval_intervals(:,1),slow_intervals.nr);
[~,ind2] = Restrict(aval_intervals(:,2),slow_intervals.nr);
aval_intervals = aval_intervals(intersect(ind1,ind2),:);

rip_up_f = numel(Restrict(rip_peak_t,aval_intervals)) / sum(diff(aval_intervals,1,2));
down_interval = [aval_intervals(1:end-1,2),aval_intervals(2:end,1)];
% keep only down intervals which fall inside slow intervals
[~,ind] = Restrict(mean(down_interval,2),slow_intervals.nr);
down_interval = down_interval(ind,:);
rip_down_f = numel(Restrict(rip_peak_t,down_interval)) / sum(diff(down_interval,1,2));

% plot
R.plotFiringRates(0,0,window,regions=70,smooth=smo,mode="ratio");
h = findobj(gca,'Type','Line');
y_max = max(get(h(1), 'YData'));
title("Ripples in slow rhythm, "+R.printBasename()+', '+labels(3)+' (n: '+num2str(R.nNeurons(70))+'), w: '+num2str(aval_window)+' s, s: '+num2str(smo)+ ...
  ', t: '+num2str(rip_thresh)+', f: ['+num2str(rip_up_f)+','+num2str(rip_down_f)+'] Hz')
adjustAxes(gca,'YTickMode','auto','YTickLabelMode','auto','YLim',[0,y_max*1.2]);
raster([rip_peak_t,y_max*1.1*ones(size(rip_peak_t))],'Color',myColors(4,'IBMcb'),'Displayname','ripples',height=0.01)
PlotIntervals(aval_intervals,'color',[0,0.6,0.6],'legend','off','alpha',0.2)
PlotIntervals(slow_intervals.nr,'legend','off')
clear h y_max



%% 1. μ and σ allow to identify slow rhythm intervals DEPRECATED
% compute firing rate with smoothing
[FR,time] = R.firingRate('all',window=window,smooth=smooth_low);
[FR_smoothed,time] = R.firingRate('all',window=window,smooth=smooth_high);

% choose parameter values
mean_window = 500; % time points
mean_thresh = 15;
std_window = 500; % time points
std_thresh = 60;
power_perc = 0.5; %0.99;
f_thresh = 1;
manMix_perc = 0.1; %0.25
ratio_thresh = 2; %1.7

plot_flag = true;
disp(' ')
disp("p "+num2str(power_perc)+" f "+num2str(f_thresh)+" m "+num2str(manMix_perc)+" r "+num2str(ratio_thresh))
for i = 1 : numel(regs)
  % detect slow-oscillations intervals
  mov_mean = movmean(FR(:,i),mean_window);
  mov_std = movstd(FR(:,i),std_window);
  slow_ind = mov_mean <= prctile(mov_mean,mean_thresh) & mov_std >= prctile(mov_std,std_thresh);

  % exclude intervals where max(FR)/min(FR) is less than ??
  interval_ind = cumsum([slow_ind(1);~slow_ind(1:end-1) & slow_ind(2:end)]);
  interval_ind = interval_ind(slow_ind);
  %interval_f.(labels(i)) = accumarray(interval_ind,FR_smoothed(slow_ind,i),[],@(x) fEstimate(x,100)/window);
  interval_f.(labels(i)) = accumarray(interval_ind,FR(slow_ind,i),[],@(x) fPowerBound(x,1/window,power_perc));
  interval_ratio.(labels(i)) = accumarray(interval_ind,FR_smoothed(slow_ind,i),[],@(x) maxMinRatio(x,manMix_perc));
  %interval_ratio.(labels(i)) = accumarray(interval_ind,FR_smoothed(slow_ind,i),[],@(x) max(x)/min(x));
  %interval_ratio.(labels(i)) = accumarray(interval_ind,FR_smoothed(slow_ind,i),[],@std);
  zero_ind = isnan(interval_f.(labels(i)));
  no_f_ind = interval_f.(labels(i)) > f_thresh;
  no_r_ind = ~zero_ind & ~no_f_ind & (interval_ratio.(labels(i)) < ratio_thresh);
  good_ind = ~zero_ind & ~no_f_ind & ~no_r_ind;
  

  % intervals and duration
  slow_intervals.(labels(i)) = [time([slow_ind(1) ~= 0; slow_ind(1:end-1) == 0 & slow_ind(2:end) == 1]), ...
    time([slow_ind(1:end-1) == 1 & slow_ind(2:end) == 0; slow_ind(end) ~= 0])];
  slow_dur.(labels(i)) = sum(slow_intervals.(labels(i))(good_ind,2) - slow_intervals.(labels(i))(good_ind,1));

  disp(labels(i)+' '+num2str(sum(good_ind))+' out of '+num2str(numel(good_ind))+', T: '+num2str(slow_dur.(labels(i))))
  
  EXAMPLES.(labels(i)) = {};
  TIMESLOW = time(slow_ind);
  FRSLOW = FR(slow_ind,i);
  %FRSLOW = FR_smoothed(slow_ind,i);
  [~,ind] = sort(slow_intervals.(labels(i))(:,2) - slow_intervals.(labels(i))(:,1));
  for k = ind(end-4:end).'
      EXAMPLES.(labels(i)){end+1} = [TIMESLOW(interval_ind==k),FRSLOW(interval_ind==k)];
  end

  % plot
  if plot_flag
    [~,axs] = makeFigure('slow',"Slow-rhythm identification, "+R.printBasename()+', '+labels(i)+' (n: '+num2str(R.nNeurons(regs(i)))+'), w: '+num2str(window)+' s, s: '+num2str(smooth_low)+ ...
      ', mw: '+num2str(mean_window)+', μt: '+num2str(mean_thresh)+', σt: '+num2str(std_thresh)+', T: '+num2str(slow_dur.(labels(i)))+' s',[2,1],TileSpacig='none');
    linkaxes(axs,'x')

    % plot detection metrics
    scale = 1; height = 0;
    plot(axs(1),time,scale*mov_mean+height,Color=myColors(5,'IBMcb'),DisplayName='μ');
    yline(axs(1),scale*prctile(mov_mean,mean_thresh)+height,Color=myColors(5,'IBMcb'),HandleVisibility='off')
    plot(axs(1),time,scale*mov_std+height,Color=myColors(4,'IBMcb'),DisplayName='σ');
    yline(axs(1),scale*prctile(mov_std,std_thresh)+height,Color=myColors(4,'IBMcb'),HandleVisibility='off')
    legend(axs(1)); adjustAxes(axs(1),'XColor','none'); ylabel(axs(1),'(a.u.)')

    % plot firing rate used for detection
    R.plotFiringRates(0,0,window,regions=regs(i),smooth=smooth_high,ax=axs(end));
    adjustAxes(axs(end),'YLim',[0,max(FR(:,i))*1.1],'YTickMode','auto','YTickLabelMode','auto');

    % plot slow-rhythm intervals
    PlotIntervals(slow_intervals.(labels(i))(good_ind,:),'legend','off')
    PlotIntervals(slow_intervals.(labels(i))(zero_ind,:),'color',[0,0,1],'legend','off','alpha',0.2)
    PlotIntervals(slow_intervals.(labels(i))(no_f_ind,:),'color',[1,0,0],'legend','off','alpha',0.2)
    PlotIntervals(slow_intervals.(labels(i))(no_r_ind,:),'color',[1,1,0],'legend','off','alpha',0.2)

    % add bar to separate subplots
    pause(1) % to avoid position problems
    axes(Position=axs(1).Position.*[1,1,0.1,1],Color='none',YColor='none',YTick=[],XTick=[],LineWidth=1.7)
  end
end
clear mov_mean mov_std slow_ind interval_ind good_ind axs scale height



%% 2. try PSD
% 2.1 try non-whitened data
target_f = 0.1; % Hz
fields = ["reu","pfc","hpc"];
% compute firing rate without smoothing
[FR,time] = R.firingRate('all',window=window);

%% choose intervals to analyze
psd_intervals = slow_intervals.reu + [-5,5];
[~,slow_ind] = Restrict(time,slow_intervals.reu);
[~,psd_ind] = Restrict(time,psd_intervals);
duration = sum(slow_intervals.reu(:,2)-slow_intervals.reu(:,1));

%% compute wavelet spectrograms
low_f_bound = 20 / duration; % frequencies smaller than 2 / duration are meaningless
high_f_bound = 1 / (2 * window);
[spectrogram.reu,t.reu,f.reu] = WaveletSpectrogram([time(psd_ind),FR(psd_ind,3)],'range',[low_f_bound,high_f_bound]);
[spectrogram.pfc,t.pfc,f.pfc] = WaveletSpectrogram([time(psd_ind),FR(psd_ind,1)],'range',[low_f_bound,high_f_bound]);
[spectrogram.hpc,t.hpc,f.hpc] = WaveletSpectrogram([time(psd_ind),FR(psd_ind,2)],'range',[low_f_bound,high_f_bound]);

% compute psd
for field = fields
  psd.(field{1}) = mean(spectrogram.(field{1})(:,ismember(psd_ind,slow_ind)),2);
end

% compute fraction of power below target frequency
for field = fields
  ind = f.(field{1}) < target_f;
  low_power.(field{1}) = trapz(f.(field{1})(ind),psd.(field{1})(ind)) / trapz(f.(field{1}),psd.(field{1}));
end

% pad spectrogram with NaNs outside sws
spec = nan(numel(f.reu),numel(time));
% ismember(psd_ind,slow_ind) to pass from an array over psd time to an array over slow rhythm time
spec(:,slow_ind) = spectrogram.reu(:,ismember(psd_ind,slow_ind));
spectrogram.reu = spec;
t.reu = time;
clear field ind spec

%% plot average spectrum, PSD
makeFigure('spectrum',"Wavelet power spectral density, non whitened, " + R.printBasename() + ', w: ' + num2str(window));
plot(f.reu,psd.reu,Color=myColors(1,'IBMcb'),LineWidth=1.8)
plot(f.pfc,psd.pfc,Color=myColors(2,'IBMcb'),LineWidth=1.8)
plot(f.hpc,psd.hpc,Color=myColors(3,'IBMcb'),LineWidth=1.8)
xline(target_f,Color=myColors(4,'IBMcb'),LineWidth=1.7)
adjustAxes(gca(),'XScale','log','YScale','log','XLim',[f.reu(1)*0.9,f.reu(end)*1.1])
lp_labels = string.empty;
for field = fields
  lp_labels = [lp_labels,num2str(low_power.(field)*100,3)+" %"];
end
legend([join([fields;lp_labels],1),'target f'])
xlabel('frequency (log(• / Hz))');
ylabel('psd (log( • / (W * s)))');

%% 2.2 try whitening data BEFORE OR AFTER RESTRICT?
FR_whitened = Whitening(FR,2000/window,false);
FR_whitened = FR_whitened(psd_ind,:);
% OR whiten data after restriction
%FR_whitened = Whitening(FR(psd_ind,:),2000/window,false);

%% compute wavelet spectrograms
[spectrogram.reu,t.reu,f.reu] = WaveletSpectrogram([time(psd_ind),FR_whitened(:,3)],'range',[low_f_bound,high_f_bound]);
[spectrogram.pfc,t.pfc,f.pfc] = WaveletSpectrogram([time(psd_ind),FR_whitened(:,1)],'range',[low_f_bound,high_f_bound]);
[spectrogram.hpc,t.hpc,f.hpc] = WaveletSpectrogram([time(psd_ind),FR_whitened(:,2)],'range',[low_f_bound,high_f_bound]);

% compute psd
for field = fields
  psd.(field{1}) = mean(spectrogram.(field{1})(:,ismember(psd_ind,slow_ind)),2);
end

% compute fraction of power below target frequency
for field = fields
  ind = f.(field{1}) < target_f;
  low_power.(field{1}) = trapz(f.(field{1})(ind),psd.(field{1})(ind)) / trapz(f.(field{1}),psd.(field{1}));
end

% pad spectrogram with NaNs outside sws
spec = nan(numel(f.reu),numel(time));
spec(:,slow_ind) = spectrogram.reu(:,ismember(psd_ind,slow_ind));
spectrogram.reu = spec;
t.reu = time;

%% plot average spectrum, PSD
makeFigure('spectrum',"Wavelet power spectral density, whitened, " + R.printBasename() + ', w: ' + num2str(window));
plot(f.reu,psd.reu,Color=myColors(1,'IBMcb'),LineWidth=1.8)
plot(f.pfc,psd.pfc,Color=myColors(2,'IBMcb'),LineWidth=1.8)
plot(f.hpc,psd.hpc,Color=myColors(3,'IBMcb'),LineWidth=1.8)
xline(target_f,Color=myColors(4,'IBMcb'),LineWidth=1.7)
adjustAxes(gca(),'XScale','log','YScale','log','XLim',[f.reu(1)*0.9,f.reu(end)*1.1])
lp_labels = string.empty;
for field = fields
  lp_labels = [lp_labels,num2str(low_power.(field)*100,3)+" %"];
end
legend([join([fields;lp_labels],1),'target f'])
xlabel('frequency (log(• / Hz))');
ylabel('psd (log( • / (W * s)))');
















%% 2. compute phase of slow rhythm
% filter signal DEPRECATED
%filtered_FR = Filter([time,FR(:,3)],'passband',[0,0.15],'nyquist',0.5/window); % Nyquist frequency is half the sampling rate
%filtered_FR(:,2) = smoothdata(filtered_FR(:,2),'gaussian',50);

% get smoothed firing rate
slow_smooth = 15;
[filtered_FR,filtered_FR_t] = R.firingRate('all',71,window=window,smooth=slow_smooth);
filtered_FR = [filtered_FR_t,filtered_FR]; clear filtered_FR_t

% get peaks
[~,peak_ind] = findpeaks(filtered_FR(:,2),'MinPeakProminence',10);
[~,thro_ind] = findpeaks(-filtered_FR(:,2),'MinPeakProminence',10);

% interpolate phase
reu_phase = [time,nan(size(filtered_FR,1),1)];
% then do interpolation MAKE F phaseLinearInterp
if peak_ind(1) < thro_ind
  for i = 1 : numel(thro_ind)
    reu_phase(peak_ind(i):thro_ind(i),2) = linspace(0,pi,thro_ind(i)-peak_ind(i)+1);
    if i ~= numel(peak_ind)
      reu_phase(thro_ind(i):peak_ind(i+1),2) = linspace(-pi,0,peak_ind(i+1)-thro_ind(i)+1);
    end
  end
else
  for i = 1 : numel(peak_ind)
    reu_phase(thro_ind(i):peak_ind(i),2) = linspace(-pi,0,peak_ind(i)-thro_ind(i)+1);
    if i ~= numel(thro_ind)
      reu_phase(peak_ind(i):thro_ind(i+1),2) = linspace(0,pi,thro_ind(i+1)-peak_ind(i)+1);
    end
  end
end
reu_phase = reu_phase(~isnan(reu_phase(:,2)),:);

% keep only peaks in slow rhythm
[~,valid_peak_ind] = myRestrict(time(peak_ind),slow_intervals.reu.');
valid_peak_ind = peak_ind(valid_peak_ind);
[~,valid_thro_ind] = myRestrict(time(thro_ind),slow_intervals.reu.');
valid_thro_ind = thro_ind(valid_thro_ind);



%% plot filtered rhythm and peaks against firing rate
R.plotFiringRates(0,0,window,regions=71,states='all',smooth=smo);

% plot filtered signal
plot(filtered_FR(:,1),filtered_FR(:,2),Color=myColors(1,'IBMcb'),DisplayName='filtered');

% plot phase
scale = 10; height = 300;
plot(reu_phase(:,1),scale*reu_phase(:,2)+height,Color=myColors(5,'IBMcb'),DisplayName='phase');

% plot peaks
yLim = ylim;
raster([time(valid_peak_ind),repmat((yLim(2)+yLim(1))/2,size(valid_peak_ind))],'Color',myColors(4,'IBMcb'),'LineWidth',1.5,'DisplayName','peak',height=(yLim(2)+yLim(1))/2)
raster([time(valid_thro_ind),repmat((yLim(2)+yLim(1))/2,size(valid_thro_ind))],'Color',myColors(2,'IBMcb'),'LineWidth',1.5,'DisplayName','through',height=(yLim(2)+yLim(1))/2)

% plot slow-rhythm intervals
PlotIntervals(slow_intervals.reu,'legend','off')

clear yLim

%% plot peaks and phase against raster
R.plotSpikeRaster(2300,2700,regions=71);

% plot phase
scale = 1.3; height = 70;
plot(reu_phase(:,1),scale*reu_phase(:,2)+height,Color=myColors(5,'IBMcb'),DisplayName='phase');

% plot peaks
yLim = ylim;
raster([time(valid_peak_ind),repmat((yLim(2)+yLim(1))/2,size(valid_peak_ind))],'Color','k','LineWidth',1.5,'DisplayName','peak',height=(yLim(2)+yLim(1))/2)
raster([time(valid_thro_ind),repmat((yLim(2)+yLim(1))/2,size(valid_thro_ind))],'Color','r','LineWidth',1.5,'DisplayName','through',height=(yLim(2)+yLim(1))/2)

% adjust zoom
yLim = ylim; ylim([yLim(1),yLim(2)*1.3])





%% 3. plot ripples
rip = readmatrix(fileparts(session)+"/"+R.basename+'.ripples',FileType='text');
rip_peak_t = rip(:,2);
raster([rip_peak_t,ones(size(rip_peak_t))],'Color',myColors(4,'IBMcb'),'Displayname','ripples',height=3)
%rip_f = Frequency(rip_peak_t,'binSize',window,'smooth',smo);
%plot(rip_f(:,1),rip_f(:,2),Color=myColors(4,'IBMcb'),DisplayName='ripples f');

%% check portions against raster
R.plotSpikeRaster(1160,1700,regions=71);



%% 2. check PETH, I CAN'T REPRODUCE THE RESULTS, DON'T KNOW WHY
rip_ind = false(size(rip_peak_t));
for i = 1 : numel(rip_peak_t)
  rip_ind(i) = any(rip_peak_t(i) >= slow_intervals(:,1) & rip_peak_t(i) <= slow_intervals(:,2));
end

[peth,bins,peth_avrg] = PETH(reu_phase(slow_ind,:),rip_peak_t(rip_ind),'durations',[-2,2],'nBins',501);

%makeFigure('peth',"Average PETH for " + R.printBasename() + ', w: ' + num2str(window) + ' s, b: [0.05,0.2' + num2str('') + ']');
%[peth,bins,peth_avrg] = PETH(reu_phase(slow_ind,:),rip_peak_t(rip_ind),'durations',[-2,2],'nBins',501,'show','on');
%adjustAxes(gca())

%makeFigure('peth',"PETH for " + R.printBasename() + ', w: ' + num2str(window) + ' s, b: [0.05,0.2' + num2str('') + ']');
%PlotColorMap(peth,'x',bins,'bar','phase')
%adjustAxes(gca())

makeFigure('phase',"Ripple phase distribution for " + R.printBasename() + ', w: ' + num2str(window) + ' s, b: [0.05,0.2' + num2str('') + ']');
histogram(peth(:,bins==0),60,Normalization='pdf')
adjustAxes(gca())

makeFigure('phase',"Ripple phase distribution for " + R.printBasename() + ', w: ' + num2str(window) + ' s, b: [0.05,0.2' + num2str('') + ']');
histogram(peth(:,300),60,Normalization='pdf')
adjustAxes(gca())



%% 3. try distribution of Reu firing rate
% show that it is more bimodal than other regions



%% 4. try psd
% whiten data BEFORE OR AFTER RESTRICT?
FR_whitened = Whitening(FR,2000/window,false);
%FR = FR_whitened;

%% choose intervals to analyze
sws_stamps = R.state_stamps{R.states=='sws'};
sws_stamps = [0,10000]; % TRY A BIG CHUNK OF DATA
ind_sws = false(size(time));
ind = false(size(time));
for interval = sws_stamps.'
  ind_sws = ind_sws | (time > interval(1) & time < interval(2));
  ind = ind | (time > interval(1)-5 & time < interval(2)+5); % extend by 5 s to avoid border effects
end; clear interval
duration = sum(sws_stamps(:,2)-sws_stamps(:,1));

%% compute multi-taper spectrum
% target period is 10 s, i.e., 0.1 Hz
MT_window = 100; % window should be > 4 times target period
[MTspectrum.reu,MTf.reu,MTerr.reu] = MTSpectrum([time(ind),FR(ind,3)],'frequency',1/window,'window',MT_window);
[MTspectrum.pfc,MTf.pfc,MTerr.pfc] = MTSpectrum([time(ind),FR(ind,1)],'frequency',1/window,'window',MT_window);
[MTspectrum.hpc,MTf.hpc,MTerr.hpc] = MTSpectrum([time(ind),FR(ind,2)],'frequency',1/window,'window',MT_window);

%% plot multi-taper spectrum
makeFigure('spectrum',"MT power spectral density for " + R.printBasename() + ', w: ' + num2str(window) + ' s, s: ' + num2str(smo));
plot(MTf.reu,MTspectrum.reu,Color=myColors(1,'IBMcb'),LineWidth=1.8)
plot(MTf.pfc,MTspectrum.pfc,Color=myColors(2,'IBMcb'),LineWidth=1.8)
plot(MTf.hpc,MTspectrum.hpc,Color=myColors(3,'IBMcb'),LineWidth=1.8)
xline(2/duration,'k',LineWidth=1.5)
xline(0.1,Color=myColors(4,'IBMcb'),LineWidth=1.5)
adjustAxes(gca(),'XScale','log','YScale','log') % ,'XLim',[0,5])
legend([fieldnames(MTf);'2 / duration';'target f'])
xlabel('frequency (Hz)');
ylabel('psd (??)');

%% compute wavelet spectrograms
low_f_bound = 20 / duration; % frequencies smaller than 2 / duration are meaningless
high_f_bound = 1 / (2 * window);
[spectrogram.reu,t.reu,f.reu] = WaveletSpectrogram([time(ind),FR(ind,3)],'range',[low_f_bound,high_f_bound]);
[spectrogram.pfc,t.pfc,f.pfc] = WaveletSpectrogram([time(ind),FR(ind,1)],'range',[low_f_bound,high_f_bound]);
[spectrogram.hpc,t.hpc,f.hpc] = WaveletSpectrogram([time(ind),FR(ind,2)],'range',[low_f_bound,high_f_bound]);

% pad spectrogram with NaNs outside sws

spec = nan(numel(f.reu),numel(time));
spec(:,ind_sws) = spectrogram.reu(:,ind_sws(ind));
spectrogram.reu = spec;
t.reu = time;

%% plot average spectrum, PSD
makeFigure('spectrum',"Wavelet power spectral density for " + R.printBasename() + ', w: ' + num2str(window) + ' s, s: ' + num2str(smo));
plot(f.reu,mean(spectrogram.reu,2),Color=myColors(1,'IBMcb'),LineWidth=1.8)
plot(f.pfc,mean(spectrogram.pfc,2),Color=myColors(2,'IBMcb'),LineWidth=1.8)
plot(f.hpc,mean(spectrogram.hpc,2),Color=myColors(3,'IBMcb'),LineWidth=1.8)
xline(2/duration,'k',LineWidth=1.5)
xline(0.1,Color=myColors(4,'IBMcb'),LineWidth=1.5)
adjustAxes(gca(),'XScale','log','YScale','log') % ,'XLim',[0,5])
legend([fieldnames(f);'2 / duration';'target f'])
xlabel('frequency (Hz)');
ylabel('psd');

%% plot Reu spectrogram
reg = 'reu';
makeFigure('spectrogram',"Reu, " + R.printBasename() + ', w: ' + num2str(window) + ' s, s: ' + num2str(smo));
f_min = 0.015; f_max = 0.3; f_ind = f.(reg)>f_min & f.(reg)<f_max;
f_range = f.(reg)(f_ind);
PlotColorMap(log(spectrogram.(reg)(f_ind,:)),'x',t.(reg),'y',f_range,'cutoffs',[5,13],'piecewise','off','bar','psd (log( • / ??))');
% make yticks
%ticks = expIndToLin(f_range);
ticks = (log(f_range)-min(log(f_range))) ./ (max(log(f_range))-min(log(f_range))) * (max(f_range)-min(f_range)) + min(f_range);
labels = f_range;
adjustAxes(gca(),'XLim',[t.(reg)(1),t.(reg)(end)],'YLim',[f_range(1),f_range(end)],'YTick',ticks(1:10:end),'YTickLabel',labels(1:10:end))
xlabel('time (s)');
ylabel('frequency (Hz)');

%% plot mPFC spectrogram
reg = 'pfc';
fig = makeFigure('spectrogram');
f_min = 0.015; f_max = 1.01; f_ind = f.(reg)>f_min & f.(reg)<f_max;
f_range = f.(reg)(f_ind);
PlotColorMap(log(spectrogram.(reg)(f_ind,:)),'x',t.(reg),'y',f_range,'cutoffs',[10,16],'piecewise','off','bar','psd (log( • / ??))');
% make yticks
%ticks = expIndToLin(f_range);
ticks = (log(f_range)-min(log(f_range))) ./ (max(log(f_range))-min(log(f_range))) * (max(f_range)-min(f_range)) + min(f_range);
labels = f_range;
adjustAxes(gca(),'XLim',[t.(reg)(1),t.(reg)(end)],'YLim',[f_range(1),f_range(end)],'YTick',ticks(1:10:end),'YTickLabel',labels(1:10:end))
xlabel('time (s)');
ylabel('frequency (Hz)');
title('mPFC');

%% define as a f maybe
function lin = expIndToLin(exp_ind)
  lin = (log(exp_ind)-min(log(exp_ind))) ./ (max(log(exp_ind))-min(log(exp_ind))) * (max(exp_ind)-min(exp_ind)) + min(exp_ind);
end