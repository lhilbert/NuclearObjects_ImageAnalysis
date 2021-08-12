function [threshold,metric] = ...
    otsuLimit(intSupport,intCounts,intLimits)
% [threshold,metric] = otsuLimit(intSupport,intCounts,intLimits)
%
% Finds Otsu threshold based on intensity values (intSupport) and their
% respective counts of occurrence (intCounts), only checks for intensity
% values in interval intLimits (give as two scalar vector)
%
% Returns optimal threshold and Otsu metric for optimal threshold

intLimits = sort(intLimits);

[intSupport,sortInd] = sort(intSupport);
intCounts = intCounts(sortInd);

minInd = find(intSupport>=intLimits(1),1,'first');
maxInd = find(intSupport<=intLimits(2),1,'last');

intSupport = intSupport(minInd:maxInd);
intCounts = intCounts(minInd:maxInd);

numSupport = numel(intSupport);

pp = intCounts / sum(intCounts);
sigma_b = zeros(1,numSupport);

for tt = 1:numSupport
   q_L = sum(pp(1:tt));
   q_H = sum(pp(tt+1:end));
   mu_L = sum(pp(1:tt) .* intSupport(1:tt)) ./ q_L;
   mu_H = sum(pp(tt+1:end) .* intSupport(tt+1:numSupport)) ./ q_H;
   sigma_b(tt) = q_L.*q_H.*(mu_L-mu_H)^2;
end

[metric,maxInd] = max(sigma_b);

threshold = intSupport(maxInd);

end