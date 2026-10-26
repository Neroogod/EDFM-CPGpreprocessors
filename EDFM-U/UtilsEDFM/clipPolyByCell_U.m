function [xarea, polyOut, onFace] = clipPolyByCell_U(poly, G, c, tol)
% Intersection of a planar fracture polygon with a convex grid cell.
%
% SYNOPSIS:
%   [xarea, polyOut, onFace] = clipPolyByCell_U(poly, G, c, tol)
%
% PARAMETERS:
%   poly - Nx3 ordered vertices of the planar fracture polygon.
%   G    - any MRST grid with computed geometry. The cell is treated as
%          the intersection of its face half-spaces, which is exact for
%          convex cells (Cartesian, corner-point hexahedra, PEBI/Voronoi,
%          tetrahedra, extruded PEBI prisms).
%   c    - cell index.
%   tol  - geometric tolerance (length units).
%
% RETURNS:
%   xarea   - area of the clipped polygon (0 if no intersection).
%   polyOut - Mx3 ordered vertices of the clipped polygon (empty if none).
%   onFace  - global face index if the fracture polygon is coplanar with
%             one of the cell's faces (the 'fracmat boundary' situation),
%             otherwise 0.
%
% METHOD: 3D Sutherland-Hodgman clipping of the polygon against each face
% plane of the cell, keeping the inside (outward signed distance <= 0)
% part. Face planes are taken as (face centroid, outward unit normal);
% MRST face normals point from neighbors(:,1) to neighbors(:,2), so the
% outward normal for cell c flips sign when c is the second neighbour.
%
% When the polygon is coplanar with one of the cell's faces every vertex
% has |d| < tol for that face, so clipping against it is a no-op; the loop
% simply records the face index and carries on, which leaves polyOut equal
% to (fracture polygon) INTERSECT (that face) - the correct boundary
% footprint.
%
% SEE ALSO:
%   EDFMgrid_U, clipCellByPlane_U

fpos  = G.cells.facePos(c):(G.cells.facePos(c+1)-1);
faces = G.cells.faces(fpos, 1);

polyOut = poly;
onFace  = 0;

for q = 1:numel(faces)
    f = faces(q);
    n = G.faces.normals(f,:);
    n = n / norm(n);
    if G.faces.neighbors(f,1) ~= c
        n = -n;
    end

    d = (polyOut - G.faces.centroids(f,:)) * n';   % signed outward distance

    if onFace == 0 && all(abs(d) < tol)
        onFace = f;                                % coplanar: clip is a no-op
        continue
    end

    polyOut = clipHalfspace(polyOut, d, tol);
    if size(polyOut,1) < 3
        polyOut = zeros(0,3);
        xarea   = 0;
        onFace  = 0;
        return
    end
end

polyOut = dedupPoly_U(polyOut, tol);
xarea   = polyArea3D_U(polyOut);

if size(polyOut,1) < 3 || xarea < tol^2
    polyOut = zeros(0,3);
    xarea   = 0;
    onFace  = 0;
end
end

% -------------------------------------------------------------------------
function P = clipHalfspace(P, d, tol)
% Sutherland-Hodgman: keep the d <= tol side of the polygon P.
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
