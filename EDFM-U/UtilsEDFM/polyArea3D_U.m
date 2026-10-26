function a = polyArea3D_U(P)
% Area of a planar polygon in 3D, by fan triangulation about its centroid.
%
% SYNOPSIS:
%   a = polyArea3D_U(P)
%
% Unlike hfm/utils/polyArea3D this routine does not assert coplanarity and
% never returns a negative value, so it is safe on polygons produced by
% floating-point clipping (which are coplanar only to within roundoff).
% For an exactly planar polygon the two agree to machine precision.

if size(P,1) < 3
    a = 0;
    return
end

c = mean(P,1);
A = P - c;
B = A([2:end, 1], :);
a = 0.5 * sum(sqrt(sum(cross(A, B, 2).^2, 2)));
end
