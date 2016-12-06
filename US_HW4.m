%%Q1







dopplerData = importdata('homework4_2.mat');
fs = dopplerData.fs; 
signal = dopplerData.signal;
f0 = dopplerData.f0;

imagesc(data_delay(:,:,64));
colormap gray

baselineRF = squeeze(sum(data_delay,2));
baselineEnv = abs(hilbert(baeslineRF));
baselineEnv = baselinEnv./max(baselineEnv(:));
imagesc(20*log10(baselineEnv), [-60 -10]);

%%%welch windows
N = size(data_delay,2);
n = 0:(N-1);
welch = 1-((n-(N-1)/2)/((N-1)/2)).^2

%%%apply the windows
szDat = size(data_delay);
apodMask = repmat(welch,[szDat(1),1,szDat(3)]);

%Q2

signal = dopplerData.signal;
timeArray = timeArray(0:length(signal)-1)/signal;
inPhase = cos(2*pi*f0*timearray);
outPhase = sin(2*pi*f0*timearray);

sigHeterI = signal .* inPhase;
sigHeterQ = signal .* outPhase;

