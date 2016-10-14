%% setup
M = 64;
mdct_fbank = mdct_filterbank(M)
% [violinSound, violinFs] = audioread('testsignal3.wav');
[violinSound, violinFs] = audioread('project_test1.wav');

%% a
for i = 1 : M
    figure(1);
    [freqResponse, calculatedFreq] = freqz(mdct_fbank(i,:));
%     subplot(2, 5, i); plot(calculatedFreq, freqResponse);
end

%% b
nwindow=1024;
noverlap=512;
nFFT=1024;
filteredViolin = zeros(M, length(violinSound));
figure(2);
for i = 1 : M
    bandFilter = mdct_fbank(i,:);
    filteredBand = filter(bandFilter, 1, violinSound);
    filteredViolin(i,:) = filteredBand;
%     subplot(2, 5, i);
    spectrogram(filteredBand, nwindow, noverlap, nFFT, violinFs);
end

%% c
figure(3);
downsampledViolin = zeros(M, ceil(length(violinSound)/M));
for i = 1 : M
    filteredBand = filteredViolin(i,:);
    downsampled = filteredBand(1:M:end);
    downsampledViolin(i,:) = downsampled;
%     subplot(2, 5, i);
    spectrogram(downsampled, nwindow, noverlap, nFFT, violinFs);
end

%% d
figure(4);
zeroFilled = zeros(1, length(violinSound));
interpolatedViolin = zeros(M, length(zeroFilled));
for i = 1 : M
    downsampled = downsampledViolin(i,:);
    zeroFilled(1:M:length(violinSound)) = downsampled;
    zeroFilled(numel(violinSound)) = 0;
    interpolatedViolin(i,:) = zeroFilled;
%     subplot(2, 5, i);
    spectrogram(zeroFilled, nwindow, noverlap, nFFT, violinFs);
end

% the spectral content is spread out, instead of being concentrated on a
% band

%% e
figure(5);
inverseFilteredViolin = zeros(M, length(zeroFilled));
for i = 1 : M
    bandFilter = mdct_fbank(i,:);
    inverseFilter = fliplr(bandFilter);
    inverseFilteredBand = filter(inverseFilter, 1, interpolatedViolin(i,:));
    inverseFilteredViolin(i,:) = inverseFilteredBand;
%     subplot(2, 5, i);
    spectrogram(inverseFilteredBand, nwindow, noverlap, nFFT, violinFs);
end

%% f
processedSignal = zeros(1, length(inverseFilteredBand));
for i = 1: M
    processedSignal = processedSignal + inverseFilteredViolin(i,:);
end
% subplot(2, 1, 1); plot(violinSound);
% subplot(2, 1, 2); plot(processedSignal);
% plot(residual);

phaseDiff = M*2;

origPlot = violinSound(1:end-phaseDiff+1);
origPlot = origPlot./max(origPlot);
processedPlot = processedSignal(phaseDiff:end);
processedPlot = processedPlot./max(processedPlot);
residual = processedPlot-origPlot';

figure(6);
% plot(residual);
plot(origPlot);
hold on;
plot(processedPlot./2);
hold off;










    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    