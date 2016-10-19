%% dct_filterbank_psycho test
[testSample, testSampleRate] = audioread('project_test1.wav');
bands = 64;
maskDb = 16;

[reconstructedSignal, averageBits] = dct_filterbank_psycho(testSample, testSampleRate, bands, maskDb);

figure(1);
plot(testSample, 'DisplayName', 'Original signal');
hold on;
plot(reconstructedSignal, 'DisplayName', 'Reconstructed signal');
hold off;
legend('show');

disp(['Average bitrate: ' num2str(averageBits) ' bits per sample']);
