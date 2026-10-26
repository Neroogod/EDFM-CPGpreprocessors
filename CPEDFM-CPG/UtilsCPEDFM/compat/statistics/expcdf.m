function p = expcdf(x, mu)
%EXPCDF Lightweight exponential CDF fallback used by bundled ADFNE.
% Supports the two-argument form required by ADFNE when the Statistics and
% Machine Learning Toolbox is unavailable.
if nargin < 2, mu = 1; end
assert(all(mu(:) > 0), 'Exponential mean MU must be positive.');
p = 1 - exp(-max(x, 0)./mu);
end
