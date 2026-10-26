function G = pMatFracNNCsCPEDFM_CPG(G, tol, varargin)
% Continuous projection-based EDFM (CPEDFM) NNCs on corner-point grids.
%
% SYNOPSIS:
%   G = pMatFracNNCsCPEDFM_CPG(G, tol)
%   G = pMatFracNNCsCPEDFM_CPG(G, tol, 'areaMode', 'fullface')
%
% Corner-point generalization of pMatFracNNCs3DCPEDFM.m (Rashid & Olorode,
% SPE-217469-PA, Figs. 8-13). Call AFTER EDFMGridCPG (or EDFMshalegrid)
% and fracturematrixShaleNNC3D. Follow with setupCPEDFMOpsTPFA_CPG.
%
% Algorithm (per fracture plane):
%   1. Classify the plane orientation (dip, dip direction, normal-z) and
%      select one face tag per logical direction (fracOrientation, a
%      local function of calculateProjectionCPEDFM_CPG); branch table
%      identical to the Cartesian release.
%   2. For every host matrix cell, resolve the selected tags into actual
%      face sets via getLogicalFaces (fault-split faces give sets with
%      several members; pinch-outs give empty sets) and compute the
%      partial projected area  A_part = A_fm * |cos(faceNormal, fracNormal)|.
%   3. Group host cells into local neighborhoods by logical IJK indices
%      (buildNeighborhoods) - the corner-point-safe replacement for
%      centroid-coordinate matching.
%   4. Seal, per direction d, the union of two face sets, and build one
%      projected-matrix<->fracture NNC per sealed face with Eq. 10
%      transmissibility
%      T = k_harm * A / dist(projCellCentroid, fracSegmentCentroid):
%
%      (i) ENVELOPE faces - the paper's collapse rule: in every local
%          neighborhood keep the selected-side face set of the EXTREME
%          member (first member when the selected side is the odd tag,
%          last when even). Reproduces the Cartesian release for the
%          1-2-cell-thick host columns of its validation cases.
%      (ii) FLIP faces - any direction-d face of a host cell whose two
%          neighboring cells lie on opposite sides of the fracture plane
%          (signed-distance sign change). This set is, by construction,
%          a cut of the matrix connectivity graph separating the two
%          sides of the fracture wherever the fracture extends, so the
%          barrier is watertight even where the host slab is several
%          cells thick (steep/coarse columns, fault throw) - a case in
%          which the envelope alone provably leaks.
%
%      The area A depends on 'areaMode':
%      'fullface' (default) - the COMPLETE area of the kept face, per the
%          paper (Fig. 10: "the complete or continuous area projected
%          using CPEDFM"; Fig. 13: "the area ... equals the area of the
%          blue-colored interface between Cells 9 and 10"). This is what
%          seals the stair-step: the matrix-matrix multiplier on every
%          kept face becomes exactly zero, and the kept faces of all
%          neighborhoods form a closed envelope around the fracture.
%          For fractures parallel to a grid axis the summed partial areas
%          equal the full face area anyway, so both modes coincide there
%          (which is why the Cartesian release, which sums partials,
%          seals its vertical validation cases correctly).
%      'partialsum' - the sum of the members' partial projected areas,
%          reproducing the Cartesian release bit-for-bit. Under-seals
%          fractures that are inclined to ALL grid axes (the partial sum
%          covers only a fraction of each face; no face multiplier
%          reaches zero and the barrier leaks).
%   5. ('partialsum' only) Across fracture planes, faces receiving
%      projections from several planes get their areas rescaled so the
%      total equals the union of the projected polygons
%      (getUnionProjectArea). Under 'fullface' the multiplier clamp
%      makes this correction moot.
%
% Appends to G.nnc:
%   cells, T, type ('fracpmat'), area   - one row per projected NNC,
%                                         consumed by setupEDFMOperatorsTPFA.
%   pMMneighs  - [hostCell projCell] pairs (ordered as in the Cartesian
%                release: [proj host] for negative sides, [host proj] for
%                positive sides).
%   pMMfaces   - global matrix face carrying each projection. THE field
%                transmultpEDFM_CPG consumes.
%   pMMarea    - projected area per face entry (union-corrected).
%   pDir       - logical direction (1,2,3) of each projection.
%   pFracPlane - fracture-grid number of the projecting plane.
%   nncFaces   - [faceID dir fracPlane] rows, for projection plotting.
%
% Differences from the Cartesian release, by design:
%   (a) Neighborhoods from IJK indices, not centroid coordinates.
%   (b) The kept face is always the extreme member's SELECTED face. The
%       Cartesian release zeroes faces with a fixed pattern that, for
%       positive-side selections in the x/y loops, retains no face and
%       returns zero merged transmissibility while still counting the
%       area. The rule used here is the geometrically consistent reading
%       of the paper ("project onto the first or the last matrix face in
%       the neighborhood, depending on the orientation angles").
%   (c) Directional permeability column (min(d, ncols)) instead of always
%       column 1. Identical for isotropic rock.

if nargin < 2, tol = 1e-5; end
opt = struct('areaMode', 'fullface');
opt = merge_options(opt, varargin{:});
assert(any(strcmp(opt.areaMode, {'fullface', 'partialsum'})), ...
    'areaMode must be ''fullface'' or ''partialsum''.');
useFullFace = strcmp(opt.areaMode, 'fullface');

assert(isfield(G, 'Matrix') && isfield(G, 'FracGrid'), ...
    'G must be a global EDFM grid (run EDFMGridCPG / EDFMshalegrid first).');
assert(isfield(G.nnc, 'normal'), ...
    'G.nnc.normal missing (run fracturematrixShaleNNC3D first).');

% ---- interior fracture-matrix NNCs, grouped by fracture plane ----------
interiorIdx = find(strcmp(G.nnc.type, 'fracmat interior'));
if isempty(interiorIdx)
    warning('pMatFracNNCsCPEDFM_CPG:noInterior', ...
        'No interior fracture-matrix NNCs found; nothing to project.');
    return
end

fracNames = fieldnames(G.FracGrid);
nFrac = numel(fracNames);
cstart = zeros(nFrac,1); cend = zeros(nFrac,1);
for fn = 1:nFrac
    Fg = G.FracGrid.(fracNames{fn});
    cstart(fn) = Fg.cells.start;
    cend(fn)   = Fg.cells.start + Fg.cells.num - 1;
end

% ---- accumulators for projected NNC entries -----------------------------
E.MFcells  = zeros(0,2); % [projMatrixCell, fractureCell]
E.MMneighs = zeros(0,2); % [hostCell projCell] ordered by side convention
E.face     = zeros(0,1); % global face carrying the projection
E.area     = zeros(0,1); % projected area assigned to that face
E.T        = zeros(0,1); % pM-F transmissibility
E.dir      = zeros(0,1); % logical direction 1..3
E.plane    = zeros(0,1); % fracture-grid number

for fn = 1:nFrac
    rows = interiorIdx(G.nnc.cells(interiorIdx,2) >= cstart(fn) & ...
                       G.nnc.cells(interiorIdx,2) <= cend(fn));
    if isempty(rows)
        warning('pMatFracNNCsCPEDFM_CPG:allBoundary', ...
            '%s has no interior NNCs (fracture lies on cell faces); skipped.', ...
            fracNames{fn});
        continue
    end

    hostCells = G.nnc.cells(rows,1);
    fCells    = G.nnc.cells(rows,2);
    fmAreas   = G.nnc.area(rows);

    % Plane normal and plane point: single source of truth is what the
    % gridding stage stored. All rows of one plane carry the same point.
    plane = struct('normal', G.nnc.normal(rows(1),:), ...
                   'point',  G.nnc.planepoint(rows(1),:), ...
                   'number', fn);

    % Steps 1-4 for this plane (see calculateProjectionCPEDFM_CPG).
    Ef = calculateProjectionCPEDFM_CPG(G, plane, hostCells, fCells, ...
                                       fmAreas, tol, useFullFace);

    E.MFcells  = [E.MFcells;  Ef.MFcells];
    E.MMneighs = [E.MMneighs; Ef.MMneighs];
    E.face     = [E.face;     Ef.face];
    E.area     = [E.area;     Ef.area];
    E.T        = [E.T;        Ef.T];
    E.dir      = [E.dir;      Ef.dir];
    E.plane    = [E.plane;    Ef.plane];
end

if isempty(E.face)
    warning('pMatFracNNCsCPEDFM_CPG:noProjections', ...
        'No projected NNCs were generated.');
    return
end

% ---- Step 5 ('partialsum' only): union correction for shared faces -----
if ~useFullFace
    h = calculateProjectionCPEDFM_CPG('helpers');
    [uf, ~, ic] = unique(E.face);
    for u = 1:numel(uf)
        entries = find(ic == u);
        planes = unique(E.plane(entries));
        if numel(planes) < 2, continue; end
        % Union of the projected fracture-cell polygons on this face plane
        fcList = unique(E.MFcells(entries,2));
        aUnion = h.getUnionProjectArea(G, fcList, uf(u), tol);
        aSum = sum(E.area(entries));
        if isfinite(aUnion) && aUnion > 0 && aSum > aUnion
            scale = aUnion / aSum;
            E.area(entries) = E.area(entries) * scale;
            E.T(entries)    = E.T(entries)    * scale;
        end
    end
end

% ---- append to G.nnc -----------------------------------------------------
nNew = numel(E.face);
G.nnc.cells = [G.nnc.cells; E.MFcells];
G.nnc.T     = [G.nnc.T;     E.T];
G.nnc.type  = [G.nnc.type;  repmat({'fracpmat'}, nNew, 1)];
G.nnc.area  = [G.nnc.area;  E.area];

G.nnc.pMMneighs  = E.MMneighs;
G.nnc.pMMfaces   = E.face;
G.nnc.pMMarea    = E.area;
G.nnc.pDir       = E.dir;
G.nnc.pFracPlane = E.plane;
G.nnc.nncFaces   = [E.face, E.dir, E.plane];
end
