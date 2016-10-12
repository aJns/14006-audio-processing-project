%% Project header
% 14006 Audio Processing project work 2016
% Jonas Nikula  240497
% Vili Saura    #opnum

%% Load project test sample
[testSample, testSampleRate] = audioread('project_test1.wav');

%% Create and test Analysis-Synthesis filterbank
mBandCount = 64;
asFilterBank = mdct_filterbank( mBandCount);

%% Psychoacoustic model
% Wanted time frame length is 20ms
frameLengthS = 20/1000;
frameSamples = frameLengthS*testSampleRate;
powerOf2 = ceil(log2(frameSamples));
frameSamples = 2^powerOf2;

NDFT = frameSamples;
M = mBandCount;
x_input_signal = testSample;

%% Based on project assigment example loop
% We now have the sub­band signals for M bands, given by the analysis
% block
si=0;ei=NDFT; %Start and End indices of input signal for each frame

win=0;
% Main loop
while(ei < length(x_input_signal))
    win=win+1; % keep track of number of processed windows/frames
    % indices of datapoints of full­band original signal
    % for the psychoacoustic model:
    time_index_fullband = si:ei;
    
    % Subband time indices at 1/M sampling rate (for each band).
    % these sample indices need to be then quantized in this window
    time_index_subband = max((si-1)/M,1):ei/M;
    
    % The most important part: SPL, Masking threshold, SMR and
    % quantization here for samples, and count the number of bits
    % used, e.g. sum up how many sub­band samples were quantized with
    % nbits over all sub­bands and all processing frames.
    
    si=si+NDFT; % advance 1 full frame in input data (no overlap)
    ei=si+NDFT-1;
end
% Proceed to reconstruct the signal from the quantized subband
% signals as in exercise 4.
%%