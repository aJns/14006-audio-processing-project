%% filereads
[audioSignal, audioFs] = audioread('make_my_day.wav');

%% Task 1
frameDurationS = 20/1000;
frameSamples = frameDurationS * audioFs;

nextExp = nextpow2(frameSamples);
hiNextPower = 2^nextExp;
loNextPower = 2^(nextExp-1);

if abs(hiNextPower-frameSamples) < abs(loNextPower-frameSamples)
    frameSamples = hiNextPower;
else
    frameSamples = loNextPower;
end

frameCount = (length(audioSignal) / (frameSamples/2)) + 1;
frameCount = round(frameCount);

audioFrames = zeros(frameSamples, frameCount);
audioFrames(:, 1) = audioSignal(1:frameSamples);
prevStart = 0;

for i = 2: frameCount
    start = prevStart + (frameSamples/2);
    cutoff = start + frameSamples - 1;
    
    if start > length(audioSignal)
        break
    elseif cutoff > length(audioSignal)
        break
    end
    
    audioFrames(:, i) = audioSignal(start:cutoff);
    prevStart = start;
end

hanningWindow = hann(frameSamples);

for i = 1: frameCount
    audioFrames(:, i) = audioFrames(:, i) .* hanningWindow;
end

%% Task 2
frameMatrixSize = size(audioFrames);
processedFrames = zeros(frameMatrixSize);

for i = 1: frameCount
    processedFrames(:, i) = audioFrames(:, i) .* hanningWindow;
end

processedSignal = zeros(size(audioSignal));
processedSignal(1:frameSamples) = processedFrames(:, 1);
prevEnd = frameSamples;
halfFrame = frameSamples/2;

for i = 2: frameCount
    partToSumStart = prevEnd - halfFrame;
    partToSumEnd = prevEnd-1;
    partToSum = processedSignal(partToSumStart:partToSumEnd);
    newFrame = processedFrames(:, i);
    
    processedSignal(partToSumStart:partToSumEnd) = partToSum + newFrame(1:halfFrame);
    processedSignal(partToSumEnd:partToSumEnd+halfFrame) = newFrame(halfFrame:frameSamples);
    prevEnd = partToSumEnd+halfFrame;
end

figure(1);
plot(audioSignal);
hold on;
plot(processedSignal);
hold off;
    

%% Task 3

audioFFT = zeros(size(audioFrames));
rAudioFFT = zeros(size(audioFrames));

for i = 1: length(audioFrames)
    audioFFT(:, i) = fft(audioFrames(:, i));
end

for frame_index = 1: length(audioFFT)
    if frame_index==1
        nMAG = abs(audioFFT(:, frame_index));
    else
        nMAG = nMAG*.95 + 0.05*abs(audioFFT(:, frame_index));
    end
    
    xMAG = abs(audioFFT(:, frame_index)); % mixture magnitude
    xPHA = angle(audioFFT(:, frame_index)); % mixture phase angle
    
    % Reconstruct signal from noise reduced magnitude and original phase.
    reFrame = max((xMAG - nMAG),0) .* exp(i*xPHA);
    rAudioFFT(:, frame_index) = reFrame;
end

rAudioFFT(~isfinite(rAudioFFT))=0;
rAudioFFT(isnan(rAudioFFT))=0;


%% Task 4
forRealProcessedFrames = zeros(size(rAudioFFT));
for i = 1: length(rAudioFFT)
    forRealProcessedFrames(:, i) = ifft(rAudioFFT(:, i));
end

frameMatrixSize = size(forRealProcessedFrames);
processedFrames = zeros(frameMatrixSize);

for i = 1: frameCount
    processedFrames(:, i) = forRealProcessedFrames(:, i) .* hanningWindow;
end

processedSignal = zeros(size(audioSignal));
processedSignal(1:frameSamples) = processedFrames(:, 1);
prevEnd = frameSamples;
halfFrame = frameSamples/2;

for i = 2: frameCount
    partToSumStart = prevEnd - halfFrame;
    partToSumEnd = prevEnd-1;
    partToSum = processedSignal(partToSumStart:partToSumEnd);
    newFrame = processedFrames(:, i);
    
    processedSignal(partToSumStart:partToSumEnd) = partToSum + newFrame(1:halfFrame);
    processedSignal(partToSumEnd:partToSumEnd+halfFrame) = newFrame(halfFrame:frameSamples);
    prevEnd = partToSumEnd+halfFrame;
end

figure(2);
plot(audioSignal/max(audioSignal));
hold on;
plot(real(processedSignal)/max(real(processedSignal)));
hold off;



















    
    
    
    
    
    
    
    