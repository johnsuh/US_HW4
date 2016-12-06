%%Q1
%part a
size = 1000;

%rectangular window
window = -size/2:size/2;
for i = 1:length(window)
    rec(i) = 1;
end
recFFT = fft(rec);
figure
plot(rec)

%welch window
n = 0:(size-1);
welch = 1-((n-(size-1)/2)/((size-1)/2)).^2;
figure
plot(welch)

%Blackman-Harris Window
n = 0:(size-1);
a0 = 0.35875;
a1 = 0.48829;
a2 = 0.14128;
a3 = 0.01168;
BHwin = a0 - a1*cos(2*pi*n/(size-1)) + a2*cos(4*pi*n/(size-1)) - a3*cos(6*pi*n/(size-1));
figure
plot(BHwin)





% dopplerData = importdata('homework4_2.mat');
% fs = dopplerData.fs; 
% signal = dopplerData.signal;
% f0 = dopplerData.f0;
% 
% imagesc(data_delay(:,:,64));
% colormap gray
% 
% baselineRF = squeeze(sum(data_delay,2));
% baselineEnv = abs(hilbert(baeslineRF));
% baselineEnv = baselinEnv./max(baselineEnv(:));
% imagesc(20*log10(baselineEnv), [-60 -10]);
% 
% %%%welch windows
% N = size(data_delay,2);
% n = 0:(N-1);
% welch = 1-((n-(N-1)/2)/((N-1)/2)).^2
% 
% %%%apply the windows
% szDat = size(data_delay);
% apodMask = repmat(welch,[szDat(1),1,szDat(3)]);


%% Q2
%import data
signal = dopplerData.signal;
f0 = dopplerData.f0;
fs = dopplerData.fs; %%Hz?

%part a
t = 1:length(signal);
t = t./fs;
figure
plot( t(1:(length(t)/6)), signal(1:(length(t)/6)) )
xlabel('Time (s)')
ylabel('Volts (V)')
title('Part A')

%part b
S = fft(signal);
S = S(1:length(S)/2+1);
freq = fs * (0:length(signal)/2)/length(signal);
figure
plot(freq, abs(S/length(signal)))
xlabel('Frequency (Hz)')
ylabel('TBD')
title('Part B')

%%

timeArray = timeArray(0:length(signal)-1)/signal;
inPhase = cos(2*pi*f0*timearray);
outPhase = sin(2*pi*f0*timearray);

sigHeterI = signal .* inPhase;
sigHeterQ = signal .* outPhase;

