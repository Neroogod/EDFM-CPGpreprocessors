function G = pMatFracNNCsCPEDFM_U(G, tol, varargin)
% Continuous projection-based EDFM (CPEDFM) NNCs on unstructured grids.
%
% SYNOPSIS:
%   G = pMatFracNNCsCPEDFM_U(G, tol)
%
% Call AFTER EDFMgrid_U (or any gridding producing the same structures)
% and fracturematrixNNC_U. Follow with setupCPEDFMOpsTPFA_CPG, which is
% grid-agnostic and reused from the CPEDFM-CPG module.
%
% Use fracturematrixNNC_U rather than the Cartesian fracturematrixShaleNNC3D:
% the latter builds CI = A/<d> from hfm/calcdavg, which integrates <d> over
% the cell's axis-aligned BOUNDING BOX. That is exact only for Cartesian
% cells; on the unstructured grids this module targets it biases every
% fracture-matrix transmissibility (~16% median on PEBI, ~44% on tetrahedra).
%
% ALGORITHM - this function loops the fracture planes and delegates the
% per-plane work to calculateProjectionCPEDFM_U, which seals every face of
% a fracture-hosting cell whose two neighboring cell centroids straddle
% the fracture surface (tip-guarded by an in-polygon test). See that file
% for the full description; in short, the structured-grid machinery of the
% CPEDFM paper (orientation branch table, IJK local neighborhoods,
% first/last-member collapse) has no meaning on unstructured grids and is
% replaced by the sign-flip cut, which is provably continuous for any cell
% shape.
%
% For each sealed face one projected-matrix<->fracture NNC is added with
% the Eq. 10 transmissibility
%   T = k_harm * A / dist(projCellCentroid, fracCellCentroid)
% where A is the complete face area ('fullface', default, per Figs. 10
% and 13 of the paper) or the clipped projection of the fracture-cell
% polygon onto the face ('partialsum'), and k_harm harmonically averages
% the normal-projected permeabilities n'*K*n of the fracture and
% projection cells.
%
% Appends to G.nnc: cells, T, type ('fracpmat'), area, and the
% face-indexed projection record pMMfaces / pMMarea / pMMneighs / pDir /
% pFracPlane / nncFaces consumed by transmultpEDFM_CPG and the plotting
% helpers (pDir holds the dominant axis of the face normal, for plot
% coloring only).
%
% FACE-COPLANAR FRACTURES are handled separately, by sealCoplanarFaces
% below. Such a fracture needs no projection at all - the face to close is
% the one it lies in - but it does need closing: left open, a sealing
% fracture on a face plane is simply bypassed by matrix-matrix flow across
% that same face, which is the leak CPEDFM exists to remove. Those faces
% are added to the pMM* record without new NNC rows, because the fracture's
% own 'fracmat boundary' NNCs to BOTH neighbours already carry the flow.
% (Earlier versions returned early when no interior NNC existed, so a
% purely coplanar fracture had nothing sealed.)

if nargin < 2, tol = 1e-5; end
opt = struct('areaMode', 'fullface');
opt = merge_options(opt, varargin{:});
assert(any(strcmp(opt.areaMode, {'fullface', 'partialsum'})), ...
    'areaMode must be ''fullface'' or ''partialsum''.');
useFullFace = strcmp(opt.areaMode, 'fullface');

assert(isfield(G, 'Matrix') && isfield(G, 'FracGrid'), ...
    'G must be a global EDFM grid (run EDFMgrid_U first).');
assert(isfield(G.nnc, 'normal'), ...
    'G.nnc.normal missing (run fracturematrixNNC_U first).');

interiorIdx = find(strcmp(G.nnc.type, 'fracmat interior'));

fracNames = fieldnames(G.FracGrid);
nFrac = numel(fracNames);
cstart = zeros(nFrac,1); cend = zeros(nFrac,1);
for fn = 1:nFrac
    Fg = G.FracGrid.(fracNames{fn});
    cstart(fn) = Fg.cells.start;
    cend(fn)   = Fg.cells.start + Fg.cells.num - 1;
end

E.MFcells  = zeros(0,2); E.MMneighs = zeros(0,2);
E.face = zeros(0,1); E.area = zeros(0,1); E.T = zeros(0,1);
E.dir  = zeros(0,1); E.plane = zeros(0,1);

for fn = 1:nFrac
    rows = interiorIdx(G.nnc.cells(interiorIdx,2) >= cstart(fn) & ...
                       G.nnc.cells(interiorIdx,2) <= cend(fn));
    if isempty(rows), continue; end

    hostCells = G.nnc.cells(rows,1);
    fCells    = G.nnc.cells(rows,2);

    % Plane normal, plane point and polygon: single source of truth is
    % what the gridding stage stored. All rows of one plane agree.
    plane = struct('normal', G.nnc.normal(rows(1),:), ...
                   'point',  G.nnc.planepoint(rows(1),:), ...
                   'poly',   G.FracGrid.(fracNames{fn}).points, ...
                   'number', fn);

    % Flip-cut projection for this plane (see calculateProjectionCPEDFM_U).
    Ef = calculateProjectionCPEDFM_U(G, plane, hostCells, fCells, ...
                                     tol, useFullFace);

    E.MFcells  = [E.MFcells;  Ef.MFcells];
    E.MMneighs = [E.MMneighs; Ef.MMneighs];
    E.face     = [E.face;     Ef.face];
    E.area     = [E.area;     Ef.area];
    E.T        = [E.T;        Ef.T];
    E.dir      = [E.dir;      Ef.dir];
    E.plane    = [E.plane;    Ef.plane];
end

%% Face-coplanar fractures: seal the face the fracture occupies
% These need no projection - the face to close is the one the fracture
% lies in. Sealing it forces all cross-face flow through the fracture's
% own 'fracmat boundary' NNCs, which fracturematrixNNC_U created for BOTH
% neighbouring cells. No new NNC row is required, only a face record for
% the multiplier stage, so these entries are kept apart from E.
S = sealCoplanarFaces(G, cstart, cend, tol);

if isempty(E.face) && isempty(S.face)
    warning('pMatFracNNCsCPEDFM_U:noProjections', ...
        'No projected or sealed faces were generated.');
    return
end

nNew = numel(E.face);
if nNew > 0
    G.nnc.cells = [G.nnc.cells; E.MFcells];
    G.nnc.T     = [G.nnc.T;     E.T];
    G.nnc.type  = [G.nnc.type;  repmat({'fracpmat'}, nNew, 1)];
    G.nnc.area  = [G.nnc.area;  E.area];

    % keep every per-connection column the same height as G.nnc.cells: the
    % projected rows have no fracture plane of their own, and the diagnostic
    % columns fracturematrixNNC_U adds carry no meaning for them
    G.nnc.normal     = [G.nnc.normal;     zeros(nNew,3)];
    G.nnc.planepoint = [G.nnc.planepoint; zeros(nNew,3)];
    if isfield(G.nnc, 'CI'),   G.nnc.CI   = [G.nnc.CI;   nan(nNew,1)]; end
    if isfield(G.nnc, 'davg'), G.nnc.davg = [G.nnc.davg; nan(nNew,1)]; end
end

G.nnc.pMMneighs  = [E.MMneighs; S.neighs];
G.nnc.pMMfaces   = [E.face;     S.face];
G.nnc.pMMarea    = [E.area;     S.area];
G.nnc.pDir       = [E.dir;      S.dir];
G.nnc.pFracPlane = [E.plane;    S.plane];
G.nnc.nncFaces   = [G.nnc.pMMfaces, G.nnc.pDir, G.nnc.pFracPlane];
end

% =========================================================================
function S = sealCoplanarFaces(G, cstart, cend, tol)
% Face records closing the matrix-matrix faces that face-coplanar
% ('fracmat boundary') fracture cells lie in.
%
% Each such fracture cell has one NNC row per neighbouring matrix cell, so
% the face is registered ONCE per fracture cell rather than once per row -
% accumulating both rows would double the projected area and over-seal
% neighbouring faces. The sealed area is the fracture's clipped footprint
% on that face, so a fracture covering only part of a face closes only
% that part, consistent with how transmultpEDFM_CPG accumulates areas.

S.face = zeros(0,1); S.area = zeros(0,1); S.neighs = zeros(0,2);
S.dir  = zeros(0,1); S.plane = zeros(0,1);

bndIdx = find(strcmp(G.nnc.type, 'fracmat boundary'));
if isempty(bndIdx), return; end

Gm = G.Matrix;
[fCells, ~, grp] = unique(G.nnc.cells(bndIdx, 2));

for k = 1:numel(fCells)
    rows  = bndIdx(grp == k);
    hosts = G.nnc.cells(rows, 1);
    hosts = hosts(hosts <= Gm.cells.num);
    if numel(hosts) < 2
        continue    % fracture on a domain boundary face: nothing to seal
    end

    f = sharedFace(Gm, hosts(1), hosts(2));
    if f == 0
        continue
    end

    nrm = Gm.faces.normals(f,:) / norm(Gm.faces.normals(f,:));
    [~, domAx] = max(abs(nrm));

    fp = find(fCells(k) >= cstart & fCells(k) <= cend, 1);

    S.face  (end+1,1) = f;
    S.area  (end+1,1) = min(G.nnc.area(rows(1)), Gm.faces.areas(f));
    S.neighs(end+1,:) = [hosts(1), hosts(2)];
    S.dir   (end+1,1) = domAx;
    S.plane (end+1,1) = fp;
end

% a face reached from more than one fracture cell must still be listed once
if ~isempty(S.face)
    [S.face, ia] = unique(S.face);
    S.area   = S.area(ia);
    S.neighs = S.neighs(ia,:);
    S.dir    = S.dir(ia);
    S.plane  = S.plane(ia);
end

assert(all(S.area > tol^2 | S.area == 0), ...
    'sealCoplanarFaces: degenerate sealed area.');
end

% -------------------------------------------------------------------------
function f = sharedFace(Gm, c1, c2)
% The single face separating two neighbouring cells (0 if they share none).
f1 = Gm.cells.faces(Gm.cells.facePos(c1):(Gm.cells.facePos(c1+1)-1), 1);
f2 = Gm.cells.faces(Gm.cells.facePos(c2):(Gm.cells.facePos(c2+1)-1), 1);
shared = intersect(f1, f2);
if isempty(shared)
    f = 0;
else
    f = shared(1);
end
end
