function aUnion = projectedAreaPEDFM_CPG(G, fCells, face, tol)
%PROJECTEDAREAPEDFM_CPG Union area of fracture cells projected onto a face plane.
%
% SYNOPSIS:
%   aUnion = projectedAreaPEDFM_CPG(G, fCells, face, tol)
%
% Implements Eq. 29 of Olorode, Wang & Rashid (2020), SPE-201243-PA:
% when several fracture segments are colocated within one host matrix cell,
% the total projected area on a given side is the UNION of the individual
% projections, not their sum.
%
% PARAMETERS:
%   G      - global EDFM grid.
%   fCells - fracture cell indices (global numbering) whose mid-plane
%            polygons are to be projected.
%   face   - global face index defining the projection plane.
%   tol    - geometric tolerance.
%
% RETURNS:
%   aUnion - area of the union of the projected polygons, or NaN if the
%            projection could not be formed (degenerate fracture cell, or
%            polyshape failure). Callers must treat NaN as "keep the summed
%            area" and count it as a fallback.
%
% CORNER-POINT GENERALIZATION
%   The Cartesian reference (getUnionProjectArea, a subfunction of
%   pMatFracNNCs3D.m) projects onto the global YZ, XZ or XY plane selected
%   by an axis flag. That is only correct when the receiving face is
%   axis-aligned. Here the projection plane is the ACTUAL face plane, with
%   an orthonormal in-plane basis built from the face normal, so tilted and
%   sheared corner-point faces get the right area.
%   See pedfm_implementation.md section 2.4.

    nf = G.faces.normals(face,:);
    nrm = norm(nf);
    if nrm <= eps
        aUnion = NaN;
        return
    end
    nf = nf / nrm;

    % Orthonormal in-plane basis (u,v) for the face plane.
    [~, kmin] = min(abs(nf));
    e = zeros(1,3); e(kmin) = 1;
    u = cross(nf, e); u = u / norm(u);
    v = cross(nf, u);

    polys = [];
    for i = 1:numel(fCells)
        % Fracture cells are thin extruded boxes; averaging the two node
        % layers recovers the mid-plane polygon (same construction as the
        % Cartesian release).
        [cn, ~] = gridCellNodes(G, fCells(i));
        xyz = G.nodes.coords(cn,:);
        nV = size(xyz,1);
        if mod(nV,2) ~= 0 || nV < 6
            aUnion = NaN;
            return
        end
        mid = (xyz(1:nV/2,:) + xyz(nV/2+1:end,:)) / 2;
        mid = orderCoplanarCCW(mid);

        P2 = [mid*u', mid*v'];
        try
            ps = polyshape(P2(:,1), P2(:,2), 'Simplify', true);
        catch
            aUnion = NaN;
            return
        end
        if ps.NumRegions >= 1 && area(ps) > tol^2
            if isempty(polys)
                polys = ps;
            else
                polys(end+1) = ps; %#ok<AGROW>
            end
        end
    end

    if isempty(polys)
        aUnion = NaN;
        return
    end
    aUnion = area(union(polys));
end


function xyzOrder = orderCoplanarCCW(xyz)
% Order coplanar 3D points along the polygon edge (SVD trick, as in the
% Cartesian release).
    xyzc = mean(xyz, 1);
    P = xyz - xyzc;
    [~, ~, V] = svd(P, 0);
    [~, is] = sort(atan2(P*V(:,1), P*V(:,2)));
    xyzOrder = xyz(is, :);
end
