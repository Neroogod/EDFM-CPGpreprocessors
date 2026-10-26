function P = orderPlanarLoop_U(P, n)
% Order coplanar points into a convex loop by angle about their mean.
%
% SYNOPSIS:
%   P = orderPlanarLoop_U(P, n)
%
% PARAMETERS:
%   P - Mx3 coplanar points (unordered).
%   n - 1x3 normal of their common plane.
%
% RETURNS:
%   P - the same points, permuted into a closed non-self-intersecting loop.
%
% Valid whenever the points are the vertices of a CONVEX polygon, which is
% the case for the cross-section of a convex polyhedron - the only use
% here. For a non-convex vertex set an angular sort does not recover the
% boundary.

n = n / norm(n);

% orthonormal in-plane basis
[~, k] = min(abs(n));
e = zeros(1,3); e(k) = 1;
u = cross(n, e); u = u / norm(u);
v = cross(n, u);

c = mean(P, 1);
D = P - c;
[~, order] = sort(atan2(D*v', D*u'));
P = P(order, :);
end
