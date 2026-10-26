function x = expinv(p, mu)
%EXPINV Lightweight exponential inverse CDF fallback used by ADFNE.
% Supports the two-argument form required by ADFNE when the Statistics and
% Machine Learning Toolbox is unavailable.
if nargin < 2, mu = 1; end
assert(all(mu(:) > 0), 'Exponential mean MU must be positive.');
assert(all(p(:) >= 0 & p(:) <= 1), 'Probabilities must lie in [0, 1].');
x = -mu.*log1p(-p);
end
