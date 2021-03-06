%% Project header
% 14006 Audio Processing project work 2016
% Jonas Nikula  240497
% Vili Saura    240264

function bands = divide_into_bands(signal, band_count, filter_bank)

downsampledBands = zeros(band_count, ceil(length(signal)/band_count));
for i = 1: band_count
    bandFilter = filter_bank(i,:);
    filteredBand = filter(bandFilter, 1, signal);
    downsampled = filteredBand(1:band_count:end);
    downsampledBands(i,:) = downsampled;
end

bands = downsampledBands;

end