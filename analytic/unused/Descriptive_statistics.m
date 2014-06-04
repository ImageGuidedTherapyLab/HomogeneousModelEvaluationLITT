% This function quickly does descriptive statistics for a 1-D distribution

function [stats] = Descriptive_statistics ( dist );

stats.mean = mean(dist);
stats.median = median(dist);
stats.std = std(dist);
stats.skew = skewness(dist);
stats.kurt = kurtosis(dist);
stats.n = length(dist);

end