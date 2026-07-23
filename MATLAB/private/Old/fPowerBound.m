function f = fPowerBound(x,fs,power)
% fPowerBound

% % fast Fourier transform
% N = length(x);
% xdft = fft(x);
% % keep only positive frequencies as x is Real
% xdft = xdft(1:N/2+1);
% % psd
% psd = (1/(fs*N)) * abs(xdft).^2;
% psd(2:end-1) = 2 * psd(2:end-1);
% freq = 0 : fs/N : fs/2;

if numel(x) < 5
  f = NaN;
  return
end

% Lomb-Scargle psd
freq = (0 : 0.025: 10).';
psd_plomb = plomb(x,fs,freq);
% normalized cumulative power spectral density
cpsd = cumsum(psd_plomb);
cpsd = cpsd / cpsd(end);

% % psd wavelet
% [spectrogram,~,f_wavelet] = WaveletSpectrogram([1/fs*(0:numel(x)-1).',x]); %,'range',[0,2]
% psd_wave = mean(spectrogram,2);

% find target f
f_ind = find(cpsd>power,1);
if ~isscalar(f_ind)
  f_ind = numel(freq);
end
f = freq(f_ind);