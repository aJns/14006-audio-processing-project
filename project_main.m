%% Project header
% 14006 Audio Processing project work 2016
% Jonas Nikula  240497
% Vili Saura    #opnum

%% Load project test sample
[testSample, testSampleRate] = audioread('project_test1.wav');

%% Create and test Analysis-Synthesis filterbank
mBandCount = 64;
asFilterBank = mdct_filterbank( mBandCount);

% filter into subbands
bandFilteredSignals = zeros(mBandCount, length(testSample));
for i = 1: length(mBandCount)
    subFilter = asFilterBank(i);
    bandFilteredSignals(i, :) = filter(subFilter, 1, testSample);
end

% downsample
downsampledBands = zeros(mBandCount, ceil(length(testSample)/mBandCount));
for i = 1 : mBandCount
    filteredBand = bandFilteredSignals(i,:);
    downsampled = filteredBand(1:mBandCount:end);
    downsampledBands(i,:) = downsampled;
end

%% inverse filter subbands and combine them to reconstruct the signal
reconstructedSignal = zeros(size(testSample));
for i = 1: length(mBandCount)
    subFilter = fliplr(asFilterBank(i));
    signalToFilter = bandFilteredSignals(i,:)';
    filterResult = filter(1, subFilter, signalToFilter);
    reconstructedSignal = reconstructedSignal + filterResult;
end

figure(1);
plot(testSample);
hold on;
plot(reconstructedSignal);
hold off;

%% Psychoacoustic model
% Wanted time frame length is 20ms
frameLengthS = 20/1000;
frameSamples = frameLengthS*testSampleRate;
powerOf2 = ceil(log2(frameSamples));
frameSamples = 2^powerOf2;


%% Based on project assigment example loop
NDFT = frameSamples;
M = mBandCount;
x_input_signal = testSample;


X_dft = zeros(frameSamples,floor(length(testSample)/frameSamples));
% We now have the sub­band signals for M bands, given by the analysis
% block
si=1;ei=NDFT; % Start and End indices of input signal for each frame

win=0;
% Main loop
while(ei < length(x_input_signal))
    win=win+1; % keep track of number of processed windows/frames
    % indices of datapoints of full­band original signal
    % for the psychoacoustic model:
    time_index_fullband = si:ei;
    
    X_dft(:,win) = fft(x_input_signal(time_index_fullband));
    
    % Subband time indices at 1/M sampling rate (for each band).
    % these sample indices need to be then quantized in this window
    time_index_subband = max((si-1)/M,1):ei/M;
    for i = 1: mBandCount
    end
    
    % The most important part: SPL, Masking threshold, SMR and
    % quantization here for samples, and count the number of bits
    % used, e.g. sum up how many sub­band samples were quantized with
    % nbits over all sub­bands and all processing frames.
    
    % SPL calculation
    freq1k = 1000;
    sinusoid1k = sin(2*pi*freq1k*linspace(0,1,NDFT));
    maxlevel1k = max(fft(sinusoid1k));
    SPL = 96 + 20*log10(abs(X_dft(1:NDFT/2,win)/maxlevel1k));
    
    % visualize spectrum and spl
%     if mod(win, 50) == 0
%         figure();
%         dftFrame = X_dft(:,win);
%         dftFrame = dftFrame.*conj(dftFrame);
%         dftFrame = (dftFrame/max(dftFrame)) * max(SPL);
%         plot(dftFrame);
%         hold on;
%         plot(SPL+abs(min(SPL)));
%         hold off;
%     end
    
    si=si+NDFT; % advance 1 full frame in input data (no overlap)
    ei=si+NDFT-1;
end
% Proceed to reconstruct the signal from the quantized subband
% signals as in exercise 4.







%%






















