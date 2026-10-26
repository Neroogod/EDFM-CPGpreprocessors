function G = pMatFracNNCsPEDFM_U(G, tol, varargin)
%PMATFRACNNCSPEDFM_U 3D pEDFM projections on unstructured grids.
%
% SYNOPSIS:
%   G = pMatFracNNCsPEDFM_U(G, tol)
%   G = pMatFracNNCsPEDFM_U(G, tol, 'unionArea', false)
%
% Unstructured implementation of the 3D projection-based EDFM of
% Olorode, Wang & Rashid (2020), SPE-201243-PA. The paper's algorithm is
% restricted to structured Cartesian grids and defers everything else
% (p. 2147); this function covers general unstructured grids whose cells
% are convex polyhedra: PEBI/Voronoi, tetrahedra, extruded prisms.
%
% Call after EDFMGridU (from CPEDFM-U) and fracturematrixShaleNNC3D.
% Follow with setupPEDFMOpsTPFA3D_CPG from PEDFM-CPG, which is
% grid-agnostic.
%
% PARAMETERS:
%   G   - global EDFM grid with G.Matrix, G.FracGrid, G.nnc and G.rock.
%   tol - absolute geometric tolerance (default 1e-5).
%
% OPTIONAL ARGUMENTS:
%   'relTol'    - relative tolerance for the equidistant-pair test, as a
%                 fraction of the host cell's characteristic length
%                 V^(1/3) (default 1e-6).
%   'unionArea' - apply the Eq. 29 union correction when several fracture
%                 cells project onto the same face (default true).
%   'defaultFace' - 'bottom' (default) | 'top' | 'left' | 'right' |
%                 'front' | 'back'. Seeds the paper's Case 3, where every
%                 antipodal pair is equidistant and the first face is
%                 chosen arbitrarily. The unstructured reading is "the face
%                 whose outward normal is most aligned with that
%                 direction"; 'bottom' matches the Cartesian release's
%                 default so the two agree on Cartesian grids.
%
% RETURNS:
%   G - input grid with G.nnc extended:
%         .cells      appended [projCell, fracCell] pM-F pairs
%         .T          appended pM-F transmissibilities (Eqs. 24-25)
%         .type       appended 'fracpmat'
%         .area       appended pM-F contact areas
%         .pMMneighs  [hostCell, projCell]
%         .pMMfaces   global face index carrying each projection
%         .pMMarea    projected area on that face (union-corrected)
%         .pDir       dominant axis of the face normal (plot colouring)
%         .pFracPlane fracture-grid number
%         .pStats     nHosts, nTies, nProjections, nUnionFallback
%
% METHOD
%   The paper selects, for each host cell, the nearer face of each of the
%   three axis pairs. Unstructured cells have no logical directions, no
%   face tags and any number of faces, so that step is replaced by the
%   antipodal-pair rule described in selectProjectionFacesPEDFM_U: pair
%   each face with the face whose outward normal is most anti-parallel to
%   it, and keep the one the fracture is closer to. On a hexahedron the
%   antipodal pairs ARE the axis pairs and the rule reproduces the
%   published selection exactly (tests/test20_regression_cartesian_U.m
%   matches pMatFracNNCs3D on all 101 projection pairs).
%
%   A tip guard then drops selected faces the fracture polygon does not
%   actually reach, since the selection rule sees only the infinite plane.
%
% RELATION TO CPEDFM-U
%   This is not the continuous variant. CPEDFM-U seals every face whose two
%   neighbouring cell centroids straddle the fracture, which is by
%   construction a cut of the matrix connectivity graph and therefore
%   watertight. pEDFM selects faces independently per host cell, so the
%   sealed set need not be connected and a sealing fracture throttles flow
%   without necessarily compartmentalizing the domain. That difference is
%   the reason the continuous variant exists; use CPEDFM-U when guaranteed
%   compartmentalization is required.
%
% SEE ALSO: calculateProjectionPEDFM_U, selectProjectionFacesPEDFM_U,
%           setupPEDFMOpsTPFA3D_CPG, transmultPEDFM3D_CPG.

    if nargin < 2 || isempty(tol), tol = 1e-5; end

    opt = struct('relTol', 1e-6, 'unionArea', true, ...
                 'defaultFace', 'bottom', 'tol', tol);
    opt = merge_options(opt, varargin{:});
    opt.tol = tol;

    assert(isfield(G, 'Matrix') && isfield(G, 'FracGrid'), ...
        'pMatFracNNCsPEDFM_U:notGridded', ...
        'G must be a global EDFM grid (run EDFMGridU first).');
    assert(isfield(G, 'nnc') && isfield(G.nnc, 'normal'), ...
        'pMatFracNNCsPEDFM_U:noNNC', ...
        'G.nnc.normal missing (run fracturematrixShaleNNC3D first).');

    stats = struct('nHosts', 0, 'nTies', 0, 'nFallback', 0, ...
                   'nProjections', 0, 'nUnionFallback', 0);

    interiorIdx = find(strcmp(G.nnc.type, 'fracmat interior'));
    if isempty(interiorIdx)
        warning('pMatFracNNCsPEDFM_U:noInterior', ...
            'No interior fracture-matrix NNCs found; nothing to project.');
        G.nnc.pStats = stats;
        return
    end

    fracNames = fieldnames(G.FracGrid);
    nFrac = numel(fracNames);

    E.MFcells = zeros(0,2); E.MMneighs = zeros(0,2);
    E.face = zeros(0,1); E.area = zeros(0,1); E.T = zeros(0,1);
    E.dir  = zeros(0,1); E.plane = zeros(0,1);

    for fn = 1:nFrac
        Fg = G.FracGrid.(fracNames{fn});
        cstart = Fg.cells.start;
        cend   = cstart + Fg.cells.num - 1;

        rows = interiorIdx(G.nnc.cells(interiorIdx,2) >= cstart & ...
                           G.nnc.cells(interiorIdx,2) <= cend);
        if isempty(rows), continue; end

        plane = struct('normal', G.nnc.normal(rows(1),:), ...
                       'point',  G.nnc.planepoint(rows(1),:), ...
                       'poly',   Fg.points, ...
                       'number', fn);

        Ef = calculateProjectionPEDFM_U(G, plane, G.nnc.cells(rows,1), ...
            G.nnc.cells(rows,2), G.nnc.area(rows), tol, opt);

        stats.nHosts    = stats.nHosts + numel(rows);
        stats.nTies     = stats.nTies + Ef.nTies;
        stats.nFallback = stats.nFallback + Ef.nFallback;

        E.MFcells  = [E.MFcells;  Ef.MFcells];
        E.MMneighs = [E.MMneighs; Ef.MMneighs];
        E.face     = [E.face;     Ef.face];
        E.area     = [E.area;     Ef.area];
        E.T        = [E.T;        Ef.T];
        E.dir      = [E.dir;      Ef.dir];
        E.plane    = [E.plane;    Ef.plane];
    end

    if isempty(E.face)
        warning('pMatFracNNCsPEDFM_U:noProjections', ...
            'No projected NNCs were generated.');
        G.nnc.pStats = stats;
        return
    end

    % Collapse duplicate (face, projCell, fracCell) records.
    [~, keep, ic] = unique([E.face, E.MFcells], 'rows', 'stable');
    if numel(keep) < numel(E.face)
        areaS  = accumarray(ic, E.area);
        transS = accumarray(ic, E.T);
        E.MFcells = E.MFcells(keep,:); E.MMneighs = E.MMneighs(keep,:);
        E.face = E.face(keep); E.dir = E.dir(keep); E.plane = E.plane(keep);
        E.area = areaS; E.T = transS;
    end

    % Eq. 29: several fracture cells reaching one face contribute their
    % UNION, not their sum. Only the per-face area driving the Eq. 27
    % multiplier is corrected; the pM-F transmissibilities stay per
    % fracture cell.
    pArea = E.area;
    if opt.unionArea
        [uf, ~, fidx] = unique(E.face);
        for k = 1:numel(uf)
            rows = find(fidx == k);
            if numel(rows) < 2, continue; end
            aU = projectedAreaPEDFM_CPG(G, E.MFcells(rows,2), uf(k), tol);
            if isnan(aU)
                stats.nUnionFallback = stats.nUnionFallback + 1;
                continue
            end
            s = sum(E.area(rows));
            if s > 0
                pArea(rows) = E.area(rows) * ...
                    (min(aU, G.faces.areas(uf(k))) / s);
            end
        end
    end

    nNew = numel(E.face);
    stats.nProjections = nNew;

    G.nnc.cells = [G.nnc.cells; E.MFcells];
    G.nnc.T     = [G.nnc.T;     E.T];
    G.nnc.type  = [G.nnc.type;  repmat({'fracpmat'}, nNew, 1)];
    G.nnc.area  = [G.nnc.area;  E.area];

    G.nnc.pMMneighs  = E.MMneighs;
    G.nnc.pMMfaces   = E.face;
    G.nnc.pMMarea    = pArea;
    G.nnc.pDir       = E.dir;
    G.nnc.pFracPlane = E.plane;
    G.nnc.nncFaces   = [E.face, E.dir, E.plane];
    G.nnc.pStats     = stats;
end
