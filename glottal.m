function blImpulseTrain = glottal(f0, fs, nharm, nsamples)

%sum k=0 to nharm, cos(2 pi k f0 n Ts) n = length of sound

blImpulseTrain = 0;
Ts = 1/fs;
n = (0:nsamples);

for k = 0 : nharm
    kElem = cos(2*pi*k*f0*n*Ts);
    blImpulseTrain = blImpulseTrain + kElem;
end

end