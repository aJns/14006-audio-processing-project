%% Project header
% 14006 Audio Processing project work 2016
% Jonas Nikula  240497
% Vili Saura    240264

function reconstruction = reconstruct_from_bands(bands, signal_length, filter_bank, account_for_phase_diff)

[bandCount, ~] = size(bands);
reconstructedSignal = zeros(1, signal_length);

interpolated = zeros(1, signal_length);
for i = 1 : bandCount
    downsampled = bands(i,:);
    interpolated(1:bandCount:signal_length) = downsampled;
    
    bandFilter = filter_bank(i,:);
    inverseFilter = fliplr(bandFilter);
    inverseFilteredBand = filter(inverseFilter, 1, interpolated);
    
    reconstructedSignal = reconstructedSignal + inverseFilteredBand;
end

if nargin < 4
    account_for_phase_diff = true;
end

if account_for_phase_diff
    phaseDiff = 2*bandCount;
    reconstruction = reconstructedSignal(phaseDiff:end);
else
    reconstruction = reconstructedSignal;
end

end