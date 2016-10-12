% what about going to regents park
% and having a look at the spring flowers
% before they all get blown away
% by this horrible wind
% they say the tulips are magnificent this year

%% T1
f0 = 100;
fs = 8192;
nharm = floor((fs/2)/f0);
nsamples = 8192;
blImpulseTrain = glottal(f0, fs, nharm, nsamples);
figure(1);
plot(blImpulseTrain);