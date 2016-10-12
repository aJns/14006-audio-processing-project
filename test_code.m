lol = 64;
tiqDbForBands = zeros(1, M);
maskThrs = zeros(lol, M);
centerFreqs = linspace(1, testSampleRate/10^3, M);
for i = 1: lol
    f = centerFreqs(i);
    tiq = 3.64 * f.^(-0.8) - 6.5*exp(-0.6*(f-3.3).^2)+(10^-3)*f.^4;
    tiqDbForBands(i) = tiq;
    % computing masking threshold from the SPL levels
    maskThr_current = max( conv(SPL_band(:,win),[0.05 0.6 0.3 0.05],'same') - MASK_dB , tiq' );
    maskThrs(i, :) = maskThr_current;
end


figure(1);
plot(SPL_band(:,win)./max(SPL_band(:,win)));
hold on;
plot(tiqDbForBands./max(tiqDbForBands));
plot(maskThr_current./max(maskThr_current));
hold off;