%%Q1

%part a
size = 128;

%rectangular window
rec = zeros(1024,1);
for i = length(rec)/2-64:length(rec)/2+64
    rec(i) = 1;
end
recFFT = abs(fftshift(fft(rec,1024)));
window = -length(recFFT)/2+1:length(recFFT)/2;
figure
plot(window, 20*log10(recFFT)-max(20*log10(recFFT)))
title('Part rec')


%welch window
welch = zeros(1,1024);
n = 0:(size-1);
welchS = 1-((n-(size-1)/2)/((size-1)/2)).^2;
welch((length(welch)/2-64+1):(length(welch)/2+64)) = welchS;
figure
plot(welch)
welchFFT = abs(fftshift(fft(welch,1024)));
window = -length(welchFFT)/2+1:length(welchFFT)/2;
figure
plot(window, 20*log10(welchFFT)-max(20*log10(welchFFT)))
title('Part welch')


%Blackman-Harris Window
BHwin = zeros(1,1024);
n = 0:(size-1);
a0 = 0.35875;
a1 = 0.48829;
a2 = 0.14128;
a3 = 0.01168;
BHwinS = a0 - a1*cos(2*pi*n/(size-1)) + a2*cos(4*pi*n/(size-1)) - a3*cos(6*pi*n/(size-1));
BHwin((length(BHwin)/2-64+1):(length(BHwin)/2+64)) = BHwinS;
figure
plot(BHwin)
BHwinFFT = abs(fftshift(fft(BHwin,1024)));
window = -length(BHwinFFT)/2+1:length(BHwinFFT)/2;
figure
plot(window, 20*log10(BHwinFFT)-max(20*log10(BHwinFFT)))
title('Part bh')

%%
%part b
%parameters
cyst = importdata('anecoicCystData.mat');
center = 64; %center channel
s_rate = cyst.samplingRateMHz * 10^6;
freq = cyst.frequencyMHz * 10^6;
c = 1540;

%data split into 5 equal segments
seg = zeros(1,5);
for i = 1:4
    seg(i) = i * floor(size(cyst.data,1)/5);
end
seg(5) = size(cyst.data,1);

%segment delay profiles
for i = 1:5
    dis = sqrt((abs((1:128)-center) * 0.1953/1000).^2 + (seg(i)/s_rate * c ).^2);
    t01 = 1/c .* dis;
    t01 = (t01 - min(t01));
    t(i,1:128) = t01;
end
t = t';
sampleD = zeros(128,5);
for i = 1:5
    sampleD(:,i) = floor(t(:,i) .* s_rate);
end
sampleD = sampleD';


%create delayed signal
seg(2:6) = seg(1:5);
seg(1) = 1;
RFdata = zeros(size(cyst.data,1),size(cyst.data,2));
for i = 2:6
    for j = 1:size(cyst.data,3)
        for n = 1:size(cyst.data,2)
            for k = seg(i-1):seg(i)
                delayedSignal( (k-seg(i-1)+1)-sampleD(i-1,n)+max(sampleD(i-1,:)),n) = cyst.data(k,n,j);
            end
        end
        delayedSignal = delayedSignal(max(sampleD(i-1,:))+1:end,:) 
        RFdata(seg(i-1):seg(i),j) = sum(delayedSignal,2);
        clear delayedSignal
    end
end

xmax = 0.1953/1000 * 128;
ymax = seg(i)/s_rate * c / 2;
bMode = abs(hilbert(RFdata));
bMode = bMode./max(max(bMode));
bModeComp = 20*log10(bMode);
figure
imagesc([0 xmax],[0 ymax],bModeComp, [-50 0])
colormap gray
xlabel('lateral width (m)')
ylabel('depth (m)')
title('Part B: Anechoic Cyst')
axis image

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
S = abs(fftshift(fft(signal)));
freq = -fs/2:fs/length(signal):(fs/2 - fs/length(signal));
figure
plot(freq,S)
xlabel('Frequency (Hz)')
ylabel('Amplitude (V)')
title('Part B')
%%
%part c
t(2:length(signal)+1) = t(1:length(signal));
t(1) = 0;
t = t(1:length(signal));
inPhase = cos(2*pi*f0*t);
outPhase = sin(2*pi*f0*t);
sigHeterI = signal .* inPhase';
sigHeterQ = signal .* outPhase';

%filter
filtBwArr = [0.01:.1:2];
for i=1:length(filtBwArr);
  filtBW = filtBwArr(i);
  filtCutoff = [f0-f0*filtBW/2 f0+f0*filtBW/2]/(fs/2);

  filtCoefs = fir1(100,filtCutoff,'low');

  sigFilt = filtfilt(filtCoefs,1,sigHeterI);
  sigFiltStore(:,i) = sigFilt;
  
  PSigErr(i) = sum((sigHeterI-sigFilt).^2);
  
  snrStore2(i) = PSigErr(i);
end

[Y,bw2] = min(snrStore2);
filtCutoff = [f0-f0*filtBwArr(bw2)/2 f0+f0*filtBwArr(bw2)/2]/(fs/2);
filtCoefs = fir1(100,filtCutoff,'low');
sigFilt2 = filtfilt(filtCoefs,1,sigHeterI);

figure
plot(sigFilt2)

figure
plot(sigHeterI)
xlabel('Frequency (Hz)')
ylabel('Amplitude (V)')
title('Part C: I Component')
figure
plot(sigHeterQ)
xlabel('Frequency (Hz)')
ylabel('Amplitude (V)')
title('Part C: Q Component')

%% Q4
%Problem 7.10






