function [dFC_peth,shuffled,t,animal] = ISRPSD_(session,regs,opt)

arguments
  session (1,:) char
  regs (:,1) string
  opt.events (:,1) string = 'all'
  opt.shuffle {mustBeInteger,mustBeNonnegative} = 0
end

R = regions(session,'phases',opt.events,'events',["InfraSlowRhythm/slownr","InfraSlowRhythm/slowavalnr"],'verbose',false);
animal = string(R.printBasename{1}(1:6));

% keep only regions found in data
[~,ind] = ismember(R.ids,regs);
ind = ind(ind~=0);
regs = regs(ind);

% show that PSD of NR FR in ISR has low f component using whitening, compare with normal regions + nISR NR
window = 0.05; % 50 ms
[fr,t] = R.firingRate('all',regs,'window',window);

is_isr = Restrict(t,R.eventIntervals('slownr'));

fr_nr_isr = fr(is_isr,regs=='nr');
fr_nr_nisr = fr(~is_isr,regs=='nr');
fr = [fr_nr_isr,fr_nr_nisr,fr(regs~='nr')];

% compute PSD for every trace
for i = 1 : size(fr,2)
  psd = WaveletSpectrogram(fr(:,i));
end


% declare variables
n_bins = 201;
for ind = nchoosek(1:numel(regs),2).'
  dFC_peth.(regs(ind(1))+"_"+regs(ind(2))) = nan(1,n_bins);
  shuffled.(regs(ind(1))+"_"+regs(ind(2))) = nan(opt.shuffle,n_bins);
end
t = nan(1,n_bins);



% firing rates of units which follow the ISR

[filebase,basename] = fileparts(session);
%load(fullfile(filebase,"InfraSlowRhythm",basename+"_ISRUnitEntrainment.mat"),'stat');
for i = 1 : numel(regs)
  [firing_rates.(regs(i)),time.fr] = R.unitFiringRates('all',regs(i),'window',window);
  %firing_rates.(regs(i)) = firing_rates.(regs(i))(:,stat.(regs(i)).h==1);
end

% make windows to compute dFC around transitions
width = 5; % time points
on = R.eventIntervals('slowavalnr');
if all(isnan(on),'all')
  return
end
[~,~,isr,on_off,off_on] = cleanISROnOff(on,R.eventIntervals('slownr')); % remove extra ON intervals
[time.dFC,time_ind] = Restrict(time.fr,off_on + [-5,5]);
time.dFC_peth = time.dFC(ceil(width/2):end-ceil(width/2)); % adjust time for PETH

% make windows for shuffled transitions
transition_sh = zeros(numel(off_on),opt.shuffle);
[time_ind_sh, time.dFC_sh, time.dFC_peth_sh] = deal(cell(opt.shuffle,1));
for k = 1 : opt.shuffle
  transition_sh(:,k) = Unshift(shuffleSpikes(Restrict(off_on,isr,'Shift','on')),isr);
  [time.dFC_sh{k},time_ind_sh{k}] = Restrict(time.fr,transition_sh(:,k) + [-5,5]);
  time.dFC_peth_sh{k} = time.dFC_sh{k}(ceil(width/2):end-ceil(width/2));
end

limits = [-0.5,0.5];
for ind = nchoosek(1:numel(regs),2).'

  field = regs(ind(1)) + "_" + regs(ind(2));

  % dFC for a region pair
  dFC.(field) = nan(numel(time_ind)-width,1);
  for t = 1 : numel(time_ind)
    if t < numel(time_ind)-width && time.dFC(t) == time.dFC(t+width)-window*width && ~isempty(firing_rates.(regs(ind(1)))) && ~isempty(firing_rates.(regs(ind(2))))
      t_ind = time_ind(t) : time_ind(t)+width-1;
      dFC.(field)(t,1) = dFC_helper(firing_rates.(regs(ind(1)))(t_ind,:), firing_rates.(regs(ind(2)))(t_ind,:));
    end
  end

  % PETH of dFC around transitions
  [~,t,dFC_peth.(field{1})] = PETH([time.dFC_peth,dFC.(field{1})],off_on,'durations',limits,'nBins',n_bins,'smooth',0,'DisplayName',strrep(field{1},'_','-'));

  % repeat for shuffled transitions
  for k = 1 : opt.shuffle
    dFC_sh.(field) = nan(numel(time_ind_sh{k})-width,1);
    for t = 1 : numel(time_ind_sh{k})
      if t < numel(time_ind_sh{k})-width && time.dFC_sh{k}(t) == time.dFC_sh{k}(t+width)-window*width && ~isempty(firing_rates.(regs(ind(1)))) && ~isempty(firing_rates.(regs(ind(2))))
        t_ind = time_ind_sh{k}(t) : time_ind_sh{k}(t)+width-1;
        dFC_sh.(field)(t,1) = dFC_helper(firing_rates.(regs(ind(1)))(t_ind,:), firing_rates.(regs(ind(2)))(t_ind,:));
      end
    end
    [~,t,shuffled.(field{1})(k,:)] = PETH([time.dFC_peth_sh{k},dFC_sh.(field{1})],transition_sh(:,k),'nBins',n_bins,'durations',limits,'smooth',0,'DisplayName',strrep(field{1},'_','-'));
  end
  
end

end

% --- Helper functions ---

function dfc = dFC_helper(a,b)
  dfc = mean(abs(corr(a,b,'Type','Spearman')),'all','omitmissing');
end

% --- Extra code to plot examples in debug mode ---

function plotInDebug
[nr_popul_fr,time.nr] = R.firingRate('all','nr','window',window,'smooth',25);

figure, hold on
plot(time.nr,nr_popul_fr,'Color',paperColors(3),'LineWidth',1.3)
plot(time.fr,firing_rates.hpc(:,[7,15,25]))
plot(time.fr,firing_rates.pfc(:,[9,15,25]))

plot(time.dFC,200*squeeze(dFC.(field)(10,14,:)))
plot(time.dFC,2000*mean(dFC,2,'omitmissing'),'Color',paperColors(1),'LineWidth',1.3)
plot(time.dFC,smoothdata(200*mean(dFC.hpc_pfc,2,'omitmissing'),'gaussian',20),'Color',paperColors(1),'LineWidth',1.3)
PlotIntervals(R.eventIntervals('slownr'))
PlotIntervals(R.eventIntervals('slowavalnr'))

%%
figure, hold on
fields = string(fieldnames(dFC_peth));
for i = [1,2,4] % fields to actually plot
  plot(t,dFC_peth.(fields{i}),'DisplayName',strrep(fields{i},'_','-'))
  plot(t,shuffled.(fields{i}),'Color',.7/i*[1,1,1],'DisplayName',strrep(fields{i},'_','-'))
end
legend
clearvars fields i

end