function r = weightedVectorStrength(theta,r2)
% WEIGHTEDVECTORSTRENGTH(theta) calculates the weighted vector strength of angular data (in
% radians) in the column-vector theta weighted by respective column vector of r2 values.

r = sqrt(sum(cos(theta).*r2).^2 + sum(sin(theta).*r2).^2) / sum(r2);