function [ y ] = myquantize( x, bits )
scale = max(abs(x));
x = x./scale;
levels = 1:1:2^bits;
levels = levels - mean(levels);
levels = levels / max(levels);
for in = 1:length(x)
    [~,ind] = min(abs(x(in) - levels));
    y(in) = levels(ind)*scale;
end
end
