function davg = calcDavg_U(G, cells, planenormal, planepoint, tol)
% Average normal fracture-matrix distance <d>, exact on unstructured cells.
%
% SYNOPSIS:
%   davg = calcDavg_U(G, cells, planenormal, planepoint, tol)
%
% PARAMETERS:
%   G           - MRST grid with computed geometry (the matrix grid).
%   cells       - Nx1 matrix cell indices.
%   planenormal - 1x3 fracture plane normal.
%   planepoint  - 1x3 point on the fracture plane.
%   tol         - geometric tolerance (length units).
%
% RETURNS:
%   davg - Nx1 average normal distance from the cell to the fracture plane,
%
%              <d> = (1/V) * INTEGRAL over the cell of |n.(x - x0)| dV,
%
%          the quantity entering the Moinfar et al. (2013) fracture-matrix
%          connectivity index CI = A / <d>.
%
% WHY THIS EXISTS:
%   The Cartesian implementation (hfm/.../calcdavg.m) evaluates <d> by
%   laying a 10x10x10 sample grid over the AXIS-ALIGNED BOUNDING BOX of the
%   cell nodes. For a Cartesian cell the bounding box IS the cell, so that
%   is right. For a PEBI, tetrahedral or prismatic cell the bounding box
%   can be several times the cell volume and has a different centroid, so
%   both the geometry and the sampling error corrupt <d> - and therefore
%   every fracture-matrix transmissibility.
%
% METHOD (exact, no sampling):
%   The signed distance d(x) = n.(x - x0) is affine, so over any region on
%   a single side of the plane the integral of |d| is exactly V*|d(C)| with
%   C the region's centroid. Splitting the cell at the plane into the
%   negative part (V-, C-) and the positive part (V+, C+),
%
%       <d> = ( V- * |d(C-)| + V+ * |d(C+)| ) / V .
%
%   One half-space clip (clipCellByPlane_U) supplies (V-, C-); the positive
%   part follows by complement from the grid's own volume and centroid, so
%   the two parts reproduce G.cells.volumes(c) and G.cells.centroids(c,:)
%   by construction. Cells the plane misses reduce to <d> = |d(centroid)|.
%
% NOTE:
%   Exactness relies on the cell being convex with planar faces - true for
%   Cartesian, PEBI/Voronoi, tetrahedral and extruded-prism grids. For
%   corner-point cells with warped faces the result is an approximation of
%   the same order as MRST's own geometry, and still far closer than the
%   bounding-box estimate.
%
% SEE ALSO:
%   clipCellByPlane_U, fracMatNNC_U

if nargin < 5 || isempty(tol)
    tol = 1e-5;
end

n  = planenormal / norm(planenormal);
p0 = planepoint;

cells = cells(:);
davg  = zeros(numel(cells), 1);

for i = 1:numel(cells)
    c  = cells(i);
    Vc = G.cells.volumes(c);
    Cc = G.cells.centroids(c,:);

    [Vm, Cm] = clipCellByPlane_U(G, c, n, p0, tol);
    Vm = min(max(Vm, 0), Vc);
    Vp = Vc - Vm;

    if Vm <= 0 || Vp <= tol^3
        % the plane does not cut this cell (or grazes it): <d> is the
        % distance of the cell centroid
        davg(i) = abs((Cc - p0) * n');
    else
        Cp = (Vc*Cc - Vm*Cm) / Vp;
        davg(i) = (Vm*abs((Cm - p0)*n') + Vp*abs((Cp - p0)*n')) / Vc;
    end

    if ~(davg(i) > 0)
        % degenerate (zero-volume or exactly symmetric to roundoff): fall
        % back on a length scale so CI stays finite
        davg(i) = max(davg(i), tol);
    end
end
end
