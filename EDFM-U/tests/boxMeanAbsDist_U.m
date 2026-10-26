function m = boxMeanAbsDist_U(lo, hi, n, p0)
% Exact mean of |n.(x - p0)| over an axis-aligned box. Test reference only.
%
% SYNOPSIS:
%   m = boxMeanAbsDist_U(lo, hi, n, p0)
%
% PARAMETERS:
%   lo, hi - 1x3 lower and upper box corners.
%   n      - 1x3 plane normal (need not be normalised).
%   p0     - 1x3 point on the plane.
%
% RETURNS:
%   m - (1/V) * INTEGRAL over the box of |n.(x-p0)| dV.
%
% This is the quantity calcDavg_U computes, evaluated by a route that
% shares no code with it: the integral along the dominant axis is done in
% closed form and the remaining two directions by adaptive quadrature. For
% a Cartesian cell the box IS the cell, so this is the exact <d> and can be
% used to grade both calcDavg_U and the sampling-based hfm/calcdavg.

n  = n / norm(n);
L  = hi - lo;
p0 = p0 - lo;                       % work in box-local coordinates

[~, k] = max(abs(n));
other  = setdiff(1:3, k);

f = @(u,v) lineIntegral(u, v, L, n, p0, k, other);

% The 2D integrand is only piecewise smooth (it kinks along the fracture
% trace), so pushing below ~1e-11 merely exhausts the evaluation budget.
m = integral2(f, 0, L(other(1)), 0, L(other(2)), ...
    'AbsTol', 1e-11, 'RelTol', 1e-11) / prod(L);
end

% -------------------------------------------------------------------------
function I = lineIntegral(u, v, L, n, p0, k, other)
% Integral over t in [0, L(k)] of |a + b t|, vectorised over (u,v).
%
% g(t) = (a+bt)|a+bt| / (2b) satisfies g'(t) = |a+bt| everywhere, including
% at the kink, so no case split on the root is needed.
sz = size(u);
x  = zeros(numel(u), 3);
x(:, other(1)) = u(:);
x(:, other(2)) = v(:);

a  = (x - p0) * n(:);               % value at t = 0 (x(:,k) is zero)
b  = n(k);
Lk = L(k);

if abs(b) < eps
    I = abs(a) * Lk;
else
    g = @(t) (a + b*t) .* abs(a + b*t) / (2*b);
    I = g(Lk) - g(0);
end
I = reshape(I, sz);
end
