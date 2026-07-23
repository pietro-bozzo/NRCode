function [psd,power_ratio,duration] = slowPSD_(session,regs,labels,opt)
% slowPSD_ Compute power spectral density of population firing rate on slow-oscillation intervals

arguments
  session (1,:) char
  regs (:,1) {mustBeNumeric,mustBeInteger,mustBePositive}
  labels (:,1) string
  opt.window (1,1) {mustBeNumeric,mustBePositive} = 0.05 % s
  opt.target_f (1,1) {mustBeNumeric,mustBePositive} = 0.1; % Hz
  opt.f_min (1,1) {mustBeNumeric} = -1; % Hz
end

if numel(regs) ~= numel(labels)
  error('slowIntervals_:labelsSize','Arguments ''regs'' and ''labels'' must have the same number of elements')
end

% instantiate handler, load spikes
R = regions(session);
R = R.loadSpikes();

% keep only regions found in data
labels = labels(ismember(regs,R.ids));
regs = regs(ismember(regs,R.ids));

% firing rate without smoothing
[FR,time] = R.firingRate('all',window=opt.window);

% load slow-rythm intervals of Reu
name_parts = split(session,'.');
file_root = strjoin(name_parts(1:end-1),'.');
slow_intervals = readmatrix(file_root+".slownr",FileType='text');

% PSD on slightly larger intervals to avoid border effects
psd_intervals = slow_intervals + [-5,5];
[~,slow_ind] = Restrict(time,slow_intervals);
[~,psd_ind] = Restrict(time,psd_intervals);
duration = sum(slow_intervals(:,2)-slow_intervals(:,1));

% set default value, frequencies smaller than 2 / duration are meaningless
if opt.f_min < 0
  opt.f_min = 4 / duration;
end
high_f_bound = 1 / (2 * opt.window);
for i = 1 : numel(regs)
  % wavelet spectrogram
  [spectrogram,~,psd.f] = WaveletSpectrogram([time(psd_ind),FR(psd_ind,i)],'range',[opt.f_min,high_f_bound]);
  % psd, ismember(psd_ind,slow_ind) to pass from an array over psd time to an array over slow rythm time
  psd.(labels(i)) = mean(spectrogram(:,ismember(psd_ind,slow_ind)),2);
  % fraction of power below target frequency
  ind = psd.f < opt.target_f;
  power_ratio.(labels(i)) = trapz(psd.f(ind),psd.(labels(i))(ind)) / trapz(psd.f,psd.(labels(i)));
end

% extra code to plot and check results
% makeFigure('spectrum',"Wavelet power spectral density, non whitened, " + R.printBasename() + ', w: ' + num2str(opt.window));
% plot(f,psd.nr,Color=myColors(1,'IBMcb'),LineWidth=1.8)
% plot(f,psd.pfc,Color=myColors(2,'IBMcb'),LineWidth=1.8)
% plot(f,psd.hpc,Color=myColors(3,'IBMcb'),LineWidth=1.8)
% xline(opt.target_f,Color=myColors(4,'IBMcb'),LineWidth=1.7)
% adjustAxes(gca(),'XScale','log','YScale','log','XLim',[f(1)*0.9,f(end)*1.1])
% lp_labels = string.empty;
% for field = labels([3,1,2]).'
%   lp_labels = [lp_labels,num2str(power_ratio.(field)*100,3)+" %"];
% end
% legend([join([labels([3,1,2]).';lp_labels],1),'target f'])
% xlabel('frequency (log(• / Hz))');
% ylabel('psd (log( • / (W * s)))');