function [Vin, Cin] = clipCellByPlane_U(G, c, planenormal, planepoint, tol)
% Volume and centroid of the part of a convex cell on one side of a plane.
%
% SYNOPSIS:
%   [Vin, Cin] = clipCellByPlane_U(G, c, planenormal, planepoint, tol)
%
% PARAMETERS:
%   G           - MRST grid with computed geometry.
%   c           - cell index. The cell is treated as a convex polyhedron
%                 bounded by its faces.
%   planenormal - 1x3 plane normal (need not be normalised).
%   planepoint  - 1x3 point on the plane.
%   tol         - geometric tolerance (length units).
%
% RETURNS:
%   Vin - volume of the sub-region with signed distance
%         (x - planepoint)*n' <= 0. Zero if the cell lies wholly on the
%         positive side.
%   Cin - 1x3 centroid of that sub-region (NaN when Vin == 0).
%
% METHOD:
%   1. Each face polygon of the cell is clipped by the half-space, giving
%      the part of the face that survives.
%   2. The vertices left lying in the plane (|d| <= tol) are collected and
%      ordered by angle about their mean in an in-plane basis. Because the
%      cell is convex the cut is a single convex polygon, so this ordering
%      closes the sub-polyhedron with its "cap" face.
%   3. Volume and centroid follow from a tetrahedral fan: every face is
%      fanned about its own centroid and every resulting triangle is joined
%      to an interior reference point. Convexity guarantees the tetrahedra
%      tile the sub-polyhedron without overlap, so absolute tetrahedron
%      volumes may be summed and face winding is irrelevant.
%
% SEE ALSO:
%   calcDavg_U, clipPolyByCell_U

n  = planenormal / norm(planenormal);
p0 = planepoint;

fpos  = G.cells.facePos(c):(G.cells.facePos(c+1)-1);
faces = G.cells.faces(fpos, 1);
nF    = numel(faces);

polys  = cell(nF, 1);
cutPts = zeros(0, 3);
nKept  = 0;

for q = 1:nF
    f  = faces(q);
    np = G.faces.nodePos(f):(G.faces.nodePos(f+1)-1);
    P  = G.nodes.coords(G.faces.nodes(np), :);      % ordered around the face

    d = (P - p0) * n';
    if all(d > tol)
        continue                                    % face fully outside
    end
    if any(d > tol)
        P = clipPolyHalfspace(P, d, tol);
        P = dedupPoly_U(P, tol);
        d = (P - p0) * n';
    end
    if size(P,1) < 3
        continue
    end

    nKept       = nKept + 1;
    polys{nKept} = P;
    cutPts       = [cutPts; P(abs(d) <= tol, :)];   %#ok<AGROW>
end
polys = polys(1:nKept);

if nKept == 0
    Vin = 0;
    Cin = nan(1,3);
    return
end

% Cap face: the vertices lying in the cutting plane, ordered into a loop.
if size(cutPts,1) >= 3
    cutPts = uniquetol(cutPts, tol, 'ByRows', true, 'DataScale', 1);
    if size(cutPts,1) >= 3
        polys{end+1} = orderPlanarLoop_U(cutPts, n);
    end
end

[Vin, Cin] = polyhedronVolumeCentroid(polys);

if Vin <= tol^3
    Vin = 0;
    Cin = nan(1,3);
end
end

% -------------------------------------------------------------------------
function P = clipPolyHalfspace(P, d, tol)
% Sutherland-Hodgman on a 3D polygon: keep the d <= tol side.
nV   = size(P,1);
keep = zeros(0,3);
for i = 1:nV
    j  = mod(i, nV) + 1;
    di = d(i); dj = d(j);
    if di <= tol
        keep(end+1,:) = P(i,:); %#ok<AGROW>
    end
    if (di <= tol && dj > tol) || (di > tol && dj <= tol)
        t = max(min(-di / (dj - di), 1), 0);
        keep(end+1,:) = P(i,:) + t*(P(j,:) - P(i,:)); %#ok<AGROW>
    end
end
P = keep;
end

% -------------------------------------------------------------------------
function [V, C] = polyhedronVolumeCentroid(polys)
% Volume and centroid of a convex polyhedron given as a list of face
% polygons (winding irrelevant), by an interior-point tetrahedral fan.
allV = vertcat(polys{:});
ref  = mean(allV, 1);

V = 0;
M = zeros(1,3);
for q = 1:numel(polys)
    P  = polys{q};
    fc = mean(P, 1);
    A  = P - ref;
    B  = P([2:end, 1], :) - ref;
    r  = fc - ref;
    % |det[a b r]| / 6 for each fan triangle
    vt = abs(sum(cross(A, B, 2) .* r, 2)) / 6;
    ct = (ref + fc + P + P([2:end, 1], :)) / 4;
    V  = V + sum(vt);
    M  = M + sum(vt .* ct, 1);
end

if V > 0
    C = M / V;
else
    C = nan(1,3);
end
end
