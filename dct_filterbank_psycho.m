function [y, bitrate] = dct_filterbank_psycho(x_input, fs, M, MASK_dB)
%% Create and test Analysis-Synthesis filterbank
mBandCount = M;
asFilterBank = mdct_filterbank(mBandCount);

downsampledBands = divide_into_bands(x_input, mBandCount, asFilterBank);

%% Psychoacoustic model
% Wanted time frame length is 20ms
frameLengthS = 20/1000;
frameSamples = frameLengthS*fs;
powerOf2 = ceil(log2(frameSamples));
frameSamples = 2^powerOf2;


%% Based on project assigment example loop
NDFT = frameSamples;
x_input_signal = x_input;

quantizedBands = zeros(size(downsampledBands));
frameCount = floor(length(x_input)/frameSamples);
X_dft = zeros(frameSamples,frameCount);
% We now have the sub band signals for M bands, given by the analysis
% block
si=1;ei=NDFT; % Start and End indices of input signal for each frame

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
    maxFrequency = fs/2;
    centerFreqs = linspace(1, maxFrequency/10^3, M);
    for i = 1: M
        f = centerFreqs(i);
        tiq = 3.64 * f.^(-0.8) - 6.5*exp(-0.6*(f-3.3).^2)+(10^-3)*f.^4;
        tiqDbForBands(i, win) = tiq;
    end
    % computing masking threshold from the SPL levels
    maskThr_current = max( conv(SPL_band(:,win),[0.05 0.6 0.3 0.05],'same') - MASK_dB , tiqDbForBands(:,win) );
%     maskThrs(:, win) = maskThr_current;
    th_mask = (max(th_mask*0.9,maskThr_current));
    
    % computing signal to mask ratio SMR
    tempSMR = SPL_band(:,win) - th_mask;
    tempSMR(tempSMR<0) = 0;
    SMR(:, win) = tempSMR;
    
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
    
    si=si+NDFT; % advance 1 full frame in input data (no overlap)
    ei=si+NDFT-1;
end


%% interpolate and inverse filter the downsampled subbands
%  and combine them to reconstruct the signal
quantizedBands(isnan(quantizedBands))=0;

reconstructedSignal = reconstruct_from_bands(quantizedBands, length(x_input), asFilterBank);


y = reconstructedSignal;

% avg bits/sample (bitrate?)
averageBits = sum(matrix_of_bits(:))/length(x_input);
bitrate = averageBits;

end