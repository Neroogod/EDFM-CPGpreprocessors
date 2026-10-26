function [selTags, info] = selectProjectionFacesPEDFM_CPG(G, host, fracCentroid, ...
                                                          planeNormal, planePoint, opt)
%SELECTPROJECTIONFACESPEDFM_CPG Choose the three pEDFM projection faces of a host cell.
%
% SYNOPSIS:
%   [selTags, info] = selectProjectionFacesPEDFM_CPG(G, host, fracCentroid, ...
%                                                    planeNormal, planePoint, opt)
%
% This is the flow chart of Fig. 2 in Olorode, Wang & Rashid (2020),
% SPE-201243-PA, generalized to corner-point geometry. For one host matrix
% cell and one fracture cell inside it, it returns the logical face tag
% selected in each of the three grid directions.
%
% PARAMETERS:
%   G            - global EDFM grid.
%   host         - host matrix cell index (scalar).
%   fracCentroid - 1x3 centroid of the fracture cell inside `host`.
%   planeNormal  - 1x3 fracture plane normal (need not be unit length).
%   planePoint   - 1x3 point on the fracture plane.
%   opt          - options struct with fields:
%                    .tol         absolute geometric tolerance
%                    .relTol      relative tolerance for the equidistance
%                                 test (fraction of the cell's
%                                 characteristic length)
%                    .defaultFace 'left' | 'front' | 'bottom', the
%                                 arbitrary pick for Case 3
%                    .distanceMode 'centroid' (paper) | 'plane'
%
% RETURNS:
%   selTags - 1x3. selTags(d) is the selected tag for direction d, i.e. one
%             of {2d-1, 2d}, or 0 if direction d is degenerate (no face on
%             either side, e.g. a pinched-out cell).
%   info    - struct: .caseId    number of equidistant pairs (0..3), as in
%                               the paper's Cases 0-3
%                     .nFallback times the node-containment rule failed and
%                               the centroid fallback was used
%
% ALGORITHM
%   Step 1  Group the host cell's faces by tag (never positionally; a
%           faulted CPG cell has more than six faces and they are not
%           sorted). Each direction d owns the tag pair (2d-1, 2d).
%   Step 2  Distance test. The representative point of a tag set is its
%           AREA-WEIGHTED centroid, which is the true centroid of a
%           fault-split interface and reduces to the single face centroid
%           on Cartesian grids. Compare the two sides; equal within
%           tolerance means an equidistant pair.
%   Step 3  Tie-break for Cases 1-3. Table 1 of the paper maps (direction
%           of the selected face, index of its node farthest from the
%           fracture plane) to the two remaining faces. That table is a
%           hexahedron incidence table: in every row the two listed faces
%           are exactly the other two faces of the cell that CONTAIN the
%           chosen node. This code applies that rule directly, which needs
%           no node ordering, no fixed face count, and works on fault-split
%           sides. See pedfm_implementation.md section 2.3.
%
% SEE ALSO: pMatFracNNCsPEDFM3D_CPG, getLogicalFacesPEDFM_CPG.

    faceSets = getLogicalFacesPEDFM_CPG(G, host);

    n = planeNormal(:)' / norm(planeNormal);
    selTags = zeros(1, 3);
    info = struct('caseId', 0, 'nFallback', 0);

    % Characteristic cell length for the relative equidistance test. An
    % absolute tolerance (as in the Cartesian reference) is meaningless on
    % a grid whose cells span orders of magnitude in size.
    Lc = G.cells.volumes(host)^(1/3);
    tolEq = max(opt.relTol * Lc, opt.tol);

    % --- Step 2: distance test, one direction at a time -----------------
    isTied = false(1, 3);
    for d = 1:3
        tagNeg = 2*d - 1;
        tagPos = 2*d;
        hasNeg = ~isempty(faceSets{tagNeg});
        hasPos = ~isempty(faceSets{tagPos});

        if ~hasNeg && ~hasPos
            selTags(d) = 0;              % degenerate side: nothing to seal
            continue
        elseif ~hasNeg
            selTags(d) = tagPos;         % only one side exists: no tie
            continue
        elseif ~hasPos
            selTags(d) = tagNeg;
            continue
        end

        dNeg = tagDistance(G, faceSets{tagNeg}, fracCentroid, n, opt);
        dPos = tagDistance(G, faceSets{tagPos}, fracCentroid, n, opt);

        if abs(dNeg - dPos) <= tolEq
            isTied(d) = true;            % equidistant pair
        elseif dNeg < dPos
            selTags(d) = tagNeg;
        else
            selTags(d) = tagPos;
        end
    end

    nEqui = sum(isTied);
    info.caseId = nEqui;
    if nEqui == 0
        return                            % Case 0: nothing further to do
    end

    % --- Step 3: tie-break ---------------------------------------------
    resolved = find(selTags > 0 & ~isTied);

    switch nEqui
        case 1
            % Case 1: two faces already selected. The candidate nodes are
            % those shared by the two selected sides (an edge on a
            % hexahedron); pick the one farthest from the fracture plane.
            nodesA = tagNodes(G, faceSets{selTags(resolved(1))});
            nodesB = tagNodes(G, faceSets{selTags(resolved(2))});
            cand = intersect(nodesA, nodesB);
            if isempty(cand)
                % Fault-split sides that do not share an edge. Fall back to
                % all nodes of the first selected side.
                cand = nodesA;
            end
            nStar = farthestNode(G, cand, n, planePoint);
            tied = find(isTied);
            [selTags(tied), used] = resolveByNode(G, faceSets, tied, nStar, ...
                                                  n, planePoint);
            info.nFallback = info.nFallback + double(~used);

        case 2
            % Case 2: one face selected. Use its own nodes.
            nStar = farthestNode(G, tagNodes(G, faceSets{selTags(resolved(1))}), ...
                                 n, planePoint);
            tiedDirs = find(isTied);
            for k = 1:numel(tiedDirs)
                [selTags(tiedDirs(k)), used] = resolveByNode(G, faceSets, ...
                    tiedDirs(k), nStar, n, planePoint);
                info.nFallback = info.nFallback + double(~used);
            end

        case 3
            % Case 3: no face selected anywhere. The paper picks one at
            % random; expose it as an option instead of hiding it.
            switch lower(opt.defaultFace)
                case 'left',   dDef = 1; tDef = 1;
                case 'front',  dDef = 2; tDef = 3;
                case 'bottom', dDef = 3; tDef = 5;
                otherwise
                    error('selectProjectionFacesPEDFM_CPG:defaultFace', ...
                        'Unknown defaultFace ''%s'' (use left, front or bottom).', ...
                        opt.defaultFace);
            end
            if isempty(faceSets{tDef})
                tDef = tDef + 1;                     % opposite side
                if isempty(faceSets{tDef})
                    selTags(:) = 0;                  % fully degenerate cell
                    return
                end
            end
            selTags(dDef) = tDef;
            isTied(dDef) = false;

            nStar = farthestNode(G, tagNodes(G, faceSets{tDef}), n, planePoint);
            tiedDirs = find(isTied);
            for k = 1:numel(tiedDirs)
                [selTags(tiedDirs(k)), used] = resolveByNode(G, faceSets, ...
                    tiedDirs(k), nStar, n, planePoint);
                info.nFallback = info.nFallback + double(~used);
            end
    end
end


% ========================================================================
%                          LOCAL HELPER FUNCTIONS
% ========================================================================

function dist = tagDistance(G, faces, fracCentroid, n, opt)
% Distance from the fracture cell centroid to a tag set's representative
% point. 'centroid' is the paper's rule (Euclidean distance to the face
% centroid); 'plane' measures along the fracture normal instead, which is
% better behaved on strongly sheared cells where a face centroid can be
% far away laterally while the face itself sits right against the
% fracture.
    a = G.faces.areas(faces);
    x = sum(G.faces.centroids(faces,:) .* a, 1) / sum(a);   % area-weighted
    switch lower(opt.distanceMode)
        case 'centroid'
            dist = norm(x - fracCentroid);
        case 'plane'
            dist = abs(dot(x - fracCentroid, n));
        otherwise
            error('selectProjectionFacesPEDFM_CPG:distanceMode', ...
                'Unknown distanceMode ''%s'' (use centroid or plane).', ...
                opt.distanceMode);
    end
end


function nodes = tagNodes(G, faces)
% Union of the node indices of every face in a tag set.
    nodes = [];
    for q = 1:numel(faces)
        f = faces(q);
        nodes = [nodes; G.faces.nodes(G.faces.nodePos(f):(G.faces.nodePos(f+1)-1))]; %#ok<AGROW>
    end
    nodes = unique(nodes);
end


function nStar = farthestNode(G, cand, n, planePoint)
% The candidate node farthest from the fracture plane (absolute normal
% distance, matching pointplanedistance used by the Cartesian reference).
    if isempty(cand)
        nStar = 0;
        return
    end
    d = abs((G.nodes.coords(cand,:) - planePoint) * n');
    [~, k] = max(d);
    nStar = cand(k);
end


function [tag, usedNode] = resolveByNode(G, faceSets, d, nStar, n, planePoint)
% Select the tag of direction d whose face set contains node nStar. This is
% Table 1 of the paper expressed as hexahedron incidence.
%
% Returns usedNode = false when the containment test could not decide and
% the centroid fallback was used, so the caller can count fallbacks.

    tagNeg = 2*d - 1;
    tagPos = 2*d;
    hasNeg = ~isempty(faceSets{tagNeg});
    hasPos = ~isempty(faceSets{tagPos});

    if hasNeg && ~hasPos
        tag = tagNeg; usedNode = true; return
    elseif hasPos && ~hasNeg
        tag = tagPos; usedNode = true; return
    elseif ~hasNeg && ~hasPos
        tag = 0; usedNode = true; return
    end

    if nStar > 0
        inNeg = ismember(nStar, tagNodes(G, faceSets{tagNeg}));
        inPos = ismember(nStar, tagNodes(G, faceSets{tagPos}));
        if inNeg && ~inPos
            tag = tagNeg; usedNode = true; return
        elseif inPos && ~inNeg
            tag = tagPos; usedNode = true; return
        end
    end

    % Fallback: the node test was inconclusive (a fault-split side that
    % does not touch nStar, or a degenerate cell where both sides do).
    % Keep the intent of Table 1 - push away from the fracture - by taking
    % the side whose area-weighted centroid is farther from the fracture
    % plane.
    dNeg = tagPlaneDistance(G, faceSets{tagNeg}, n, planePoint);
    dPos = tagPlaneDistance(G, faceSets{tagPos}, n, planePoint);
    if dNeg >= dPos
        tag = tagNeg;
    else
        tag = tagPos;
    end
    usedNode = false;
end


function dist = tagPlaneDistance(G, faces, n, planePoint)
% Absolute normal distance from a tag set's area-weighted centroid to the
% fracture plane.
    a = G.faces.areas(faces);
    x = sum(G.faces.centroids(faces,:) .* a, 1) / sum(a);
    dist = abs(dot(x - planePoint, n));
end
