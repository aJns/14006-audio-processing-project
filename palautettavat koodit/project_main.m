%% Project header
% 14006 Audio Processing project work 2016
% Jonas Nikula  240497
% Vili Saura    240264

%% Load project test sample
[testSample, testSampleRate] = audioread('project_test1.wav');

%% Create and test Analysis-Synthesis filterbank
mBandCount = 64;
asFilterBank = mdct_filterbank(mBandCount);

downsampledBands = divide_into_bands(testSample, mBandCount, asFilterBank);

%% interpolate and inverse filter the downsampled subbands
%  and combine them to reconstruct the signal

reconstructedSignal = reconstruct_from_bands(downsampledBands, length(testSample), asFilterBank);

figure(1);
plot(testSample, 'DisplayName', 'Original signal');
hold on;
plot(reconstructedSignal, 'DisplayName', 'Reconstructed signal');
hold off;
legend('show');

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
MASK_dB = 16;

quantizedBands = zeros(size(downsampledBands));
frameCount = floor(length(testSample)/frameSamples);
X_dft = zeros(frameSamples,frameCount);
% We now have the sub band signals for M bands, given by the analysis
% block
si=1;ei=NDFT; % Start and End indices of input signal for each frame

% figure handles
figSPLAndSpectrum = figure(2);
figBitAllocationAndSMR = figure(3);
figMaskAndTIQ = figure(4);
figJoinedThresholds = figure(5);
figSMR = figure(6);

visualizeWindowIndex = 317; %randi(500);
disp(['Visualized window index ' num2str(visualizeWindowIndex)]);
maxFrequency = testSampleRate/2;

win=0;
% Main loop
th_mask = zeros(M,1);
SMR = zeros(M, frameCount);
matrix_of_bits = zeros(M, frameCount);
while(ei < length(x_input_signal))
    win=win+1; % keep track of number of processed windows/frames
    % indices of datapoints of full band original signal
    % for the psychoacoustic model:
    time_index_fullband = si:ei;
    
    X_dft(:,win) = fft(x_input_signal(time_index_fullband));
    
    % Subband time indices at 1/M sampling rate (for each band).
    % these sample indices need to be then quantized in this window
    time_index_subband = max((si-1)/M,1):ei/M;
    
    % The most important part: SPL, Masking threshold, SMR and
    % quantization here for samples, and count the number of bits
    % used, e.g. sum up how many sub band samples were quantized with
    % nbits over all sub bands and all processing frames.
    
    % normalized SPL calculation
    freq1k = 1000;
    sinusoid1k = sin(2*pi*freq1k*linspace(0,1,NDFT));
    maxlevel1k = max(fft(sinusoid1k));
    SPL = 96 + 20*log10(abs(X_dft(1:NDFT/2,win)/maxlevel1k));
    
    % visualize spectrum and spl
    if any(win == visualizeWindowIndex)
         figure(figSPLAndSpectrum);
         clf;
         dftFrame = X_dft(:,win);
         dftFrame = dftFrame.*conj(dftFrame);
         dftFrame = (dftFrame/max(dftFrame)) * max(SPL);
         spectrumData = dftFrame;
         spectrumXAxis = linspace(1, maxFrequency, length(spectrumData));
         plot(spectrumXAxis, spectrumData, 'DisplayName', 'Magnitude Spectrum');
         hold on;
         yyaxis right;
         SPLData = SPL;
         SPLXAxis = linspace(1, maxFrequency, length(SPLData));
         plot(SPLXAxis, SPLData, 'DisplayName', 'SPL (dB)');
         hold off;
         legend('show');
         xlabel('Frequency (Hz)');
    end
    
    % max SPL per band (per frame)
    dftIndexIncrease = (NDFT/2)/M;
    dftIndices = zeros(dftIndexIncrease, M);
    for i = 1: M
        startIndex = 1 + (i-1)*dftIndexIncrease;
        endIndex = startIndex + dftIndexIncrease - 1;
        dftIndices(:,i) = startIndex:endIndex;
    end
    SPL_band = zeros(M, length(x_input_signal));
    for bandIndex = 1: M
        SPL_band(bandIndex,win) = max(0, max(SPL(dftIndices(:,bandIndex))));
    end
    
    % threshold in quiet for center frequency f
    tiqDbForBands = zeros(M, M);
%     maskThrs = zeros(M, M);
    maxFrequency = testSampleRate/2;
    centerFreqs = linspace(1, maxFrequency/10^3, M);
    for i = 1: M
        f = centerFreqs(i);
        tiq = 3.64 * f.^(-0.8) - 6.5*exp(-0.6*(f-3.3).^2)+(10^-3)*f.^4;
        tiqDbForBands(i, win) = tiq;
    end
    % computing masking threshold from the SPL levels
    maskThr_current = max( conv(SPL_band(:,win),[0.05 0.6 0.3 0.05],'same') - MASK_dB , tiqDbForBands(:,win) );
%     maskThrs(:, win) = maskThr_current;
    th_mask_old = th_mask;
    th_mask = (max(th_mask*0.9,maskThr_current));
    
    % computing signal to mask ratio SMR
    SMR(:, win) = SPL_band(:,win) - th_mask;
    
    % Bit allocation
    dividerConstant = 6.02; % dB
    bits = 1; % Where is this value supposed to come from???
    samples_in_subband = frameSamples/M;
    
    bits_for_subband = zeros(1, M);
    for ind_subband = 1: M
        bits_for_current_subband = ceil(SMR(ind_subband, win)/(dividerConstant/bits));
        if bits_for_current_subband < 0
            bits_for_current_subband = 0;
        end
        bits_for_subband(ind_subband) = bits_for_current_subband;
%         sumOfBits = matrix_of_bits(ind_subband,win) + bits_for_current_subband;
        matrix_of_bits(ind_subband,win) = bits_for_current_subband*samples_in_subband;
    end
    
    % Quantization
    for i = 1: M
        quantizedBands(i, time_index_subband) = myquantize(downsampledBands(i, time_index_subband), bits_for_subband(i));
    end
    
    % Decoder
    
    if any(win == visualizeWindowIndex)
         figure(figBitAllocationAndSMR);
         clf;
         yyaxis left;
         plot(SMR(:,win), 'DisplayName', 'SMR (dB)');
         hold on;
         yyaxis right;
         plot(matrix_of_bits(:,win), 'DisplayName', 'Bits Allocated');
         hold off;
         legend('show');
         xlabel('Band Index');
         
         
         figure(figMaskAndTIQ);
         clf;
         hold on;
         plot(SPL_band(:,win), 'DisplayName', 'SPL(band)');
         plot(tiqDbForBands(:, win), 'DisplayName', 'Threshold in Quiet');
         plot(maskThr_current, 'DisplayName', 'Masking Threshold');
         hold off;
         ylabel('(in dB)');
         xlabel('Band Index');
         legend('show');
         
         
         figure(figJoinedThresholds);
         clf;
         plot(maskThr_current, 'o');
         hold on;
         plot(th_mask_old*0.9, '+');
         plot(th_mask);
         legend('Masking Threshold of current window', 'Mask Threshold of old window * 0.9', 'Joined mask');
         ylabel('Amplitude (dB)');
         xlabel('sub-band index');
         hold off;
    end
    
    si=si+NDFT; % advance 1 full frame in input data (no overlap)
    ei=si+NDFT-1;
end

%% SMR visualization
figure(figSMR);
sampleLengthSeconds = length(testSample)/testSampleRate;

[smrYElems, smrXElems] = size(SMR);

SMRXAxis = linspace(0, sampleLengthSeconds, smrXElems);
SMRYAxis = linspace(0, maxFrequency, smrYElems);

clf;
SMR(SMR<0) = 0;
imagesc(SMRXAxis, SMRYAxis, SMR); colormap jet; colorbar; set(gca,'YDir','normal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');


%% interpolate and inverse filter the downsampled subbands
%  and combine them to reconstruct the signal
quantizedBands(isnan(quantizedBands))=0;

reconstructedSignal = reconstruct_from_bands(quantizedBands, length(testSample), asFilterBank);

figure(7);
plot(testSample, 'DisplayName', 'Original signal');
hold on;
plot(reconstructedSignal, 'DisplayName', 'Reconstructed signal');
hold off;
legend('show');

% avg bits/sample (bitrate?)
averageBits = sum(matrix_of_bits(:))/length(testSample);
disp(['Average bitrate: ' num2str(averageBits) ' bits per sample']);


%% save figures
oldFolder = cd('figures');
print(figBitAllocationAndSMR, 'BitAllocationAndSMR', '-dpng')
print(figJoinedThresholds, 'JoinedThresholds', '-dpng')
print(figMaskAndTIQ, 'MaskAndTIQ', '-dpng')
print(figSMR, 'SMR', '-dpng')
print(figSPLAndSpectrum, 'SPLAndSpectrum', '-dpng')
cd(oldFolder);
















