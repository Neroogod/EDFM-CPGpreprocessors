function G = pMatFracNNCsPEDFM3D_CPG(G, tol, varargin)
%PMATFRACNNCSPEDFM3D_CPG 3D pEDFM projections on corner-point grids.
%
% SYNOPSIS:
%   G = pMatFracNNCsPEDFM3D_CPG(G, tol)
%   G = pMatFracNNCsPEDFM3D_CPG(G, tol, 'defaultFace', 'left')
%
% Corner-point implementation of the 3D projection-based EDFM of
% Olorode, Wang & Rashid (2020), SPE-201243-PA. The paper's algorithm is
% restricted to structured Cartesian grids ("Although our algorithm could
% be adapted for simple corner-point grids with six faces, we leave this
% ... to a future study", p. 2147); this function is that extension.
%
% Call after the grid has been through EDFMGridCPG and
% fracturematrixShaleNNC3D. It appends the projection matrix-fracture NNCs
% (Eqs. 24-25) to G.nnc and records the per-face projected areas that
% setupPEDFMOpsTPFA3D_CPG / transmultPEDFM3D_CPG turn into the
% matrix-matrix transmissibility multipliers of Eq. 27.
%
% PARAMETERS:
%   G   - global EDFM grid with G.Matrix, G.FracGrid, G.nnc and G.rock.
%   tol - absolute geometric tolerance (default 1e-5).
%
% OPTIONAL ARGUMENTS:
%   'defaultFace'  - 'left' | 'front' | 'bottom' (default 'bottom'). The
%                    arbitrary face pick for Case 3, where the fracture
%                    centroid is equidistant from all three pairs of faces.
%                    Matches the Cartesian reference's default.
%   'relTol'       - relative tolerance for the equidistance test, as a
%                    fraction of the host cell's characteristic length
%                    V^(1/3) (default 1e-6).
%   'distanceMode' - 'centroid' (default, the paper's rule) or 'plane'.
%   'unionArea'    - apply the Eq. 29 union correction when several
%                    fracture cells project onto the same face
%                    (default true).
%
% RETURNS:
%   G - input grid with G.nnc extended:
%         .cells      appended [projCell, fracCell] pM-F pairs
%         .T          appended pM-F transmissibilities (Eq. 24)
%         .type       appended 'fracpmat'
%         .area       appended pM-F contact areas
%         .pMMneighs  [hostCell projCell], negative side first
%         .pMMfaces   global face index carrying each projection
%         .pMMarea    projected area on that face (union-corrected)
%         .pDir       logical direction 1..3
%         .pFracPlane fracture-grid number
%         .pStats     nCase0..nCase3, nFallback, nUnionFallback, nDropped
%
% DIFFERENCES FROM THE CARTESIAN REFERENCE (pMatFracNNCs3D.m)
%   The five Cartesian assumptions that fail on CPG, and their fixes, are
%   tabulated in pedfm_implementation.md section 2. In summary: faces are
%   addressed by tag set rather than position; projection cells are
%   resolved as "the neighbour that is not the host" rather than by
%   neighbors-column; the projection direction comes from the tag rather
%   than from find() on the face normal; the carrying face index is
%   recorded so fault-split interfaces work; and areas use the real face
%   normals and face planes.
%
%   One behavioural difference is deliberate. When several fracture cells
%   sit in one host cell, the reference collapses them to a single NNC per
%   matrix-cell pair, keeping only the first fracture cell's connection.
%   This function keeps one pM-F NNC per (projection cell, fracture cell)
%   pair - every fracture cell stays coupled to the projection cell - and
%   applies the Eq. 29 union only to the per-face area that drives the
%   Eq. 27 multiplier, which is where over-sealing would otherwise occur.
%   With one fracture cell per host cell the two are identical, which is
%   what tests/test00_regression_cartesian.m pins down.
%
% SEE ALSO: selectProjectionFacesPEDFM_CPG, transmultPEDFM3D_CPG,
%           setupPEDFMOpsTPFA3D_CPG.

    if nargin < 2 || isempty(tol), tol = 1e-5; end

    opt = struct('defaultFace',  'bottom', ...
                 'relTol',       1e-6, ...
                 'distanceMode', 'centroid', ...
                 'unionArea',    true, ...
                 'tol',          tol);
    opt = merge_options(opt, varargin{:});
    opt.tol = tol;

    assert(isfield(G, 'Matrix') && isfield(G, 'FracGrid'), ...
        'pMatFracNNCsPEDFM3D_CPG:notGridded', ...
        'Run EDFMGridCPG and fracturematrixShaleNNC3D first.');
    assert(isfield(G, 'nnc') && isfield(G.nnc, 'type'), ...
        'pMatFracNNCsPEDFM3D_CPG:noNNC', ...
        'G.nnc missing; run fracturematrixShaleNNC3D first.');

    nMat     = G.Matrix.cells.num;
    perm     = G.rock.perm;
    nPermCol = size(perm, 2);

    emptyStats = struct('nCase0', 0, 'nCase1', 0, 'nCase2', 0, 'nCase3', 0, ...
                        'nFallback', 0, 'nUnionFallback', 0, 'nDropped', 0);

    interior = find(strcmp(G.nnc.type, 'fracmat interior'));
    if isempty(interior)
        warning('pMatFracNNCsPEDFM3D_CPG:noInterior', ...
            'No interior fracture-matrix NNCs found; nothing to project.');
        G.nnc.pStats = emptyStats;   % keep the field present for callers
        return
    end

    % Map every fracture cell to the fracture grid (plane) it belongs to,
    % so projections can be attributed to a plane.
    fracNames = fieldnames(G.FracGrid);
    planeOf = zeros(G.cells.num, 1);
    for fn = 1:numel(fracNames)
        F = G.FracGrid.(fracNames{fn});
        planeOf(F.cells.start:(F.cells.start + F.cells.num - 1)) = fn;
    end

    nInt = numel(interior);
    maxRows = nInt * 3 * 4;      % 3 directions, headroom for split sides
    MF    = zeros(maxRows, 2);
    MM    = zeros(maxRows, 2);
    face  = zeros(maxRows, 1);
    area  = zeros(maxRows, 1);
    trans = zeros(maxRows, 1);
    dirv  = zeros(maxRows, 1);
    plane = zeros(maxRows, 1);
    nRec  = 0;

    stats = emptyStats;

    for ii = 1:nInt
        row  = interior(ii);
        host = G.nnc.cells(row, 1);
        fc   = G.nnc.cells(row, 2);
        fmArea = G.nnc.area(row);

        if G.cells.volumes(host) <= tol^3
            stats.nDropped = stats.nDropped + 1;   % degenerate host cell
            continue
        end

        nrm = G.nnc.normal(row, :);
        if norm(nrm) <= eps
            stats.nDropped = stats.nDropped + 1;
            continue
        end
        nrm = nrm / norm(nrm);
        ppt = G.nnc.planepoint(row, :);
        fracCentroid = G.cells.centroids(fc, :);

        % --- Fig. 2: which face is selected in each direction -----------
        [selTags, info] = selectProjectionFacesPEDFM_CPG(G, host, ...
            fracCentroid, nrm, ppt, opt);

        stats.(sprintf('nCase%d', info.caseId)) = ...
            stats.(sprintf('nCase%d', info.caseId)) + 1;
        stats.nFallback = stats.nFallback + info.nFallback;

        faceSets = getLogicalFacesPEDFM_CPG(G, host);

        for d = 1:3
            tag = selTags(d);
            if tag == 0, continue; end
            faces = faceSets{tag};
            if isempty(faces), continue; end

            isNeg = mod(tag, 2) == 1;
            kCol  = min(d, nPermCol);

            % Eq. 24 numerator: area of the fracture projected onto this
            % side, using the real (area-weighted) face normal.
            nf = sum(G.faces.normals(faces,:), 1);
            if norm(nf) <= eps, continue; end
            projTotal = fmArea * abs(dot(nf, nrm)) / norm(nf);
            if projTotal <= tol^2, continue; end

            % Split over the members of a fault-split side by face area.
            w = G.faces.areas(faces);
            w = w / sum(w);

            for q = 1:numel(faces)
                f  = faces(q);
                nb = G.faces.neighbors(f,:);
                projCell = nb(nb ~= host & nb ~= 0);
                if numel(projCell) ~= 1 || projCell > nMat
                    % Domain boundary or a fracture cell on the far side:
                    % no projection matrix cell exists. Matches
                    % 'ProjMatrix outside' in the Cartesian reference.
                    stats.nDropped = stats.nDropped + 1;
                    continue
                end

                a = min(projTotal * w(q), G.faces.areas(f));

                kf = perm(fc, kCol);
                kp = perm(projCell, kCol);
                if (kf + kp) > 0
                    kHarm = 2 * kf * kp / (kf + kp);   % Eq. 25
                else
                    kHarm = 0;   % both sealing: keep the seal, no exchange
                end
                dist = max(norm(G.cells.centroids(projCell,:) - fracCentroid), tol);

                if isNeg
                    mmPair = [projCell, host];
                else
                    mmPair = [host, projCell];
                end

                nRec = nRec + 1;
                MF(nRec,:)   = [projCell, fc];
                MM(nRec,:)   = mmPair;
                face(nRec)   = f;
                area(nRec)   = a;
                trans(nRec)  = kHarm * a / dist;      % Eq. 24
                dirv(nRec)   = d;
                plane(nRec)  = planeOf(fc);
            end
        end
    end

    if nRec == 0
        warning('pMatFracNNCsPEDFM3D_CPG:noProjections', ...
            'No pEDFM projections were generated.');
        G.nnc.pStats = stats;
        return
    end

    MF = MF(1:nRec,:);  MM = MM(1:nRec,:);  face = face(1:nRec);
    area = area(1:nRec); trans = trans(1:nRec);
    dirv = dirv(1:nRec); plane = plane(1:nRec);

    % --- Collapse duplicate (face, projCell, fracCell) records ----------
    % A fracture cell can reach the same face twice only through repeated
    % members of a split side; sum those contributions.
    [~, keep, ic] = unique([face, MF], 'rows', 'stable');
    if numel(keep) < nRec
        areaS  = accumarray(ic, area);
        transS = accumarray(ic, trans);
        MF = MF(keep,:); MM = MM(keep,:); face = face(keep);
        dirv = dirv(keep); plane = plane(keep);
        area = areaS; trans = transS;
    end

    % --- Eq. 29: union of colocated projections, per face ---------------
    % The pM-F transmissibilities stay per fracture cell; only the per-face
    % area that drives the Eq. 27 multiplier is union-corrected, since that
    % is where summing would over-seal the face.
    pArea = area;
    if opt.unionArea
        [uf, ~, fidx] = unique(face);
        for k = 1:numel(uf)
            rows = find(fidx == k);
            if numel(rows) < 2, continue; end
            aU = projectedAreaPEDFM_CPG(G, MF(rows,2), uf(k), tol);
            if isnan(aU)
                stats.nUnionFallback = stats.nUnionFallback + 1;
                continue                      % keep the summed areas
            end
            % Distribute the union area over the contributing records in
            % proportion to their individual projections, so the per-face
            % total is the union while each record keeps its share.
            s = sum(area(rows));
            if s > 0
                pArea(rows) = area(rows) * (min(aU, G.faces.areas(uf(k))) / s);
            end
        end
    end

    nNew = numel(face);
    G.nnc.cells = [G.nnc.cells; MF];
    G.nnc.T     = [G.nnc.T;     trans];
    G.nnc.type  = [G.nnc.type;  repmat({'fracpmat'}, nNew, 1)];
    G.nnc.area  = [G.nnc.area;  area];

    G.nnc.pMMneighs  = MM;
    G.nnc.pMMfaces   = face;
    G.nnc.pMMarea    = pArea;
    G.nnc.pDir       = dirv;
    G.nnc.pFracPlane = plane;
    G.nnc.pStats     = stats;
end
