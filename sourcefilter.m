nsamples = 8192; % Duration of the glottal signal in samples
f0 = 100; % Pitch in Hz
fs=8192;
nharm = floor((fs/2)/f0); % Number of harmonics to generate band-limited impulse train

% Generate the glottal signal (source) here and visualize
%
glottalSignal = glottal(f0, fs, nharm, nsamples);
figure(1);
plot(glottalSignal);
%
% Normalize the glottal signal by dividing every sample by the maximum
% value. 

normGlottal = glottalSignal./(max(glottalSignal(:)));


%% Vocal tract Filter model
F =  [850 1610 760]; % Specify the formant frequencies for vowel /a/ (Hz) 
BW = [130,  70,  160];  % Formant bandwidths (Hz)

fs = 8192;              % Sampling rate (Hz)

R = exp(-pi*BW*(1/fs));     % Pole radii
theta = ((2*pi)/fs)*F;      % Pole angles
poles = R.*exp(theta); % Generate complex poles using R and theta. Don't forget the symmetry in the pole-zero plot! 
B = 1;  % Numerator coefficients of the transfer function of the filter
A = real(poly([poles,conj(poles)])); % Denominator coefficients of the transfer function

% FILL IN !! Visualize frequency response and pole-zero plot of the vocal tract filter:
figure(2);
zplane(B./A);
% Now synthesize the vowel and listen to it:
% vowel = ; % FILL IN !!
