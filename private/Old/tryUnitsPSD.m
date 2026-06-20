% PSD of neurons, maybe differences are hidden if we don't look inside us rythm cos it's not enough time? 
% maybe change windows?

nr_spikes = R.spikes('sws',70);
pfc_spikes = R.spikes('sws',12);

sleep1 = R.event_stamps{1};
nr_spikes = Restrict(nr_spikes,sleep1);
pfc_spikes = Restrict(pfc_spikes,sleep1);

unit = 173;
f1 = Frequency(nr_spikes(nr_spikes(:,2)==unit,1),'limits',sleep1,'step',4,'binSize',0.15);

unit = 191;
f2 = Frequency(nr_spikes(nr_spikes(:,2)==unit,1),'limits',sleep1,'step',4,'binSize',0.15);

unit = 22;
f3 = Frequency(pfc_spikes(pfc_spikes(:,2)==unit,1),'limits',sleep1,'step',4,'binSize',0.15);

unit = 45;
f4 = Frequency(pfc_spikes(pfc_spikes(:,2)==unit,1),'limits',sleep1,'step',4,'binSize',0.15);

PlotXY(f1); hold on; PlotXY(f2);
PlotXY(f3); hold on; PlotXY(f4);

figure
[spectrogram,~,wf] = WaveletSpectrogram(f1); %,'range',[opt.f_min,high_f_bound]);
plot(wf,mean(spectrogram,2));
set(gca,'XScale','log','YScale','log'); hold on
[spectrogram,~,wf] = WaveletSpectrogram(f2); %,'range',[opt.f_min,high_f_bound]);
plot(wf,mean(spectrogram,2));

[spectrogram,~,wf] = WaveletSpectrogram(f3); %,'range',[opt.f_min,high_f_bound]);
plot(wf,mean(spectrogram,2));
[spectrogram,~,wf] = WaveletSpectrogram(f4); %,'range',[opt.f_min,high_f_bound]);
plot(wf,mean(spectrogram,2));