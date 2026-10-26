function [cells, CI, area, type, davg] = fracMatNNC_U(F, Gm, tol)
% Fracture-matrix NNCs for one fracture plane on an unstructured grid.
%
% SYNOPSIS:
%   [cells, CI, area, type, davg] = fracMatNNC_U(F, Gm, tol)
%
% PARAMETERS:
%   F   - one fracture grid G.FracGrid.FracN produced by EDFMgrid_U.
%   Gm  - the matrix grid (G.Matrix).
%   tol - geometric tolerance.
%
% RETURNS:
%   cells - Nx2 global cell indices [matrix cell, fracture cell].
%   CI    - Nx1 connectivity index, CI = A_surface / <d> (Moinfar et al.,
%           2013), with A_surface the fracture SURFACE area seen by the
%           matrix cell: twice the intersection area for a fracture inside
%           a cell, and the intersection area itself for a fracture lying
%           on the face between two cells (Tene et al., 2017).
%   area  - Nx1 intersection area of each connection.
%   type  - Nx1 cell array, 'fracmat boundary' or 'fracmat interior'.
%   davg  - Nx1 average normal distance used for each connection.
%
% DIFFERENCE FROM THE CARTESIAN hfm/fracmatnnc:
%   * <d> comes from calcDavg_U, integrated over the true cell polyhedron,
%     not over the cell's axis-aligned bounding box.
%   * the surface-area factor is taken from the recorded intersection type
%     instead of being inferred by counting how many rows share a fracture
%     cell, so a fracture that is face-coplanar in part of the domain and
%     interior elsewhere is handled correctly.
%
% SEE ALSO:
%   calcDavg_U, EDFMgrid_U, fracturematrixNNC_U

mc      = F.matrix_connection;
mcells  = mc.cells(:);
fcells  = mc.fraccells(:) + F.cells.start - 1;
area    = mc.area(:);
typenum = mc.type(:);

cells = [mcells, fcells];

planenormal = mc.normal(1,:);
planepoint  = mc.planepoint(1,:);

davg = calcDavg_U(Gm, mcells, planenormal, planepoint, tol);

% Fracture surface area seen by the matrix cell: an interior fracture
% exposes both of its sides to the same cell, a face-coplanar one exposes
% one side to each of the two neighbours.
surfarea = area;
surfarea(typenum == 2) = 2 * surfarea(typenum == 2);

CI = surfarea ./ davg;

type = cell(numel(typenum), 1);
type(typenum == 1) = {'fracmat boundary'};
type(typenum == 2) = {'fracmat interior'};
assert(all(~cellfun('isempty', type)), ...
    'fracMatNNC_U: unrecognised intersection type code.');
end
