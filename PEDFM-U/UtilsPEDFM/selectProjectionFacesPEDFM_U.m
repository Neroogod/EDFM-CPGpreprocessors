function [selected, info] = selectProjectionFacesPEDFM_U(G, host, fracCentroid, ...
                                                         planeNormal, planePoint, opt)
%SELECTPROJECTIONFACESPEDFM_U pEDFM projection-face selection on unstructured cells.
%
% SYNOPSIS:
%   [selected, info] = selectProjectionFacesPEDFM_U(G, host, fracCentroid, ...
%                                                   planeNormal, planePoint, opt)
%
% Generalizes the face-selection step of the 3D pEDFM of Olorode, Wang &
% Rashid (2020), SPE-201243-PA, to cells that are arbitrary convex
% polyhedra: PEBI/Voronoi, tetrahedra, extruded prisms.
%
% PARAMETERS:
%   G            - global EDFM grid.
%   host         - host matrix cell index (scalar).
%   fracCentroid - 1x3 centroid of the fracture cell inside `host`.
%   planeNormal  - 1x3 fracture plane normal (need not be unit length).
%   planePoint   - 1x3 point on the fracture plane.
%   opt          - options struct with fields .tol and .relTol.
%
% RETURNS:
%   selected - column vector of global face indices selected for projection.
%   info     - struct: .nFaces number of faces of the host cell
%                      .nTies  equidistant antipodal pairs encountered
%                      .nFallback  ties the node rule could not decide
%
% THE PROBLEM
%   The paper's algorithm assumes the host cell has exactly six faces
%   forming three opposing pairs, one per grid direction, and picks the
%   nearer face of each pair (Fig. 2, Case 0). Unstructured grids have no
%   logical directions, no face tags and any number of faces per cell, so
%   "the pair of faces in the X-direction" is undefined.
%
% RULE 1 - ANTIPODAL PAIRS (replaces the three axis pairs)
%   The three axis pairs of a hexahedron are just its three antipodal face
%   pairs, and that notion survives on any convex polyhedron:
%
%     For each face f, let its partner g be the face whose outward unit
%     normal is most anti-parallel to f's. Select f when the fracture is
%     closer to f's plane than to g's, measuring perpendicular distance
%     from the fracture cell centroid.
%
%   On a hexahedron the partner of "left" is "right" (dot product -1), and
%   the perpendicular-distance comparison is identical to the paper's
%   centroid-distance comparison: both faces of an axis pair share the
%   cell's mid-plane in the other two coordinates, so those terms cancel
%   from the comparison. The rule therefore reproduces the published
%   selection exactly - pinned down by tests/test20_regression_cartesian_U.m.
%
%   On a PEBI cell with n faces this yields roughly n/2 selections, one per
%   antipodal direction, which is the natural reading of "project the
%   fracture once per direction" for that cell shape.
%
% RULE 2 - EQUIDISTANT PAIRS (the paper's Cases 1-3)
%   When a pair is equidistant the paper consults Table 1, which maps
%   (direction of the selected face, index of its node farthest from the
%   fracture plane) to the remaining faces. As shown in
%   PEDFM-CPG/pedfm_implementation.md section 2.3, every row of that table
%   names exactly the other faces of the cell that CONTAIN the chosen node
%   - it is a hexahedron incidence table written longhand. That reading
%   needs no node ordering and no fixed face count, so it transfers here
%   unchanged:
%
%     Take the node n* of the already-selected faces that is farthest from
%     the fracture plane. For a tied pair, keep the face containing n*.
%
%   Ties are counted in info.nTies, and cases the node rule cannot decide
%   (it falls back to the lower global face index, which is antisymmetric
%   so exactly one face of the pair is always kept) in info.nFallback.
%
% NOTES
%   Pairing need not be mutual on a general polyhedron: f's partner may be
%   g while g's partner is h, so both f and g can be selected. This is
%   harmless - each selected face carries at most its own projected area
%   and the Eq. 27 multiplier is clamped at zero - and cannot occur on a
%   hexahedron.
%
% SEE ALSO: calculateProjectionPEDFM_U, pMatFracNNCsPEDFM_U.

    fpos  = G.cells.facePos(host):(G.cells.facePos(host+1)-1);
    faces = G.cells.faces(fpos, 1);
    nF    = numel(faces);

    info = struct('nFaces', nF, 'nTies', 0, 'nFallback', 0);
    selected = zeros(0,1);
    if nF == 0, return; end

    n = planeNormal(:)' / norm(planeNormal);

    % Outward unit normals and perpendicular distances from the fracture
    % cell centroid to each face plane.
    m = zeros(nF, 3);
    d = inf(nF, 1);
    for q = 1:nF
        f = faces(q);
        nrm = G.faces.normals(f,:);
        len = norm(nrm);
        if len <= eps, continue; end        % degenerate face: never nearer
        nrm = nrm / len;
        if G.faces.neighbors(f,1) ~= host
            nrm = -nrm;                     % orient outward from the host
        end
        m(q,:) = nrm;
        d(q) = abs(dot(G.faces.centroids(f,:) - fracCentroid, nrm));
    end

    Lc = G.cells.volumes(host)^(1/3);
    tolEq = max(opt.relTol * Lc, opt.tol);

    % --- Pass 1: antipodal pairing and the distance test ----------------
    partner = zeros(nF,1);
    keep = false(nF,1);
    tied = false(nF,1);
    for q = 1:nF
        if ~isfinite(d(q)), continue; end
        dots = m * m(q,:)';
        dots(q) = inf;                      % never pair a face with itself
        [~, g] = min(dots);
        partner(q) = g;
        if ~isfinite(dots(g)) || ~isfinite(d(g))
            keep(q) = true;                 % no usable partner: keep
        elseif d(q) < d(g) - tolEq
            keep(q) = true;
        elseif d(q) > d(g) + tolEq
            keep(q) = false;
        else
            tied(q) = true;
        end
    end

    % --- Pass 2: resolve equidistant pairs by node containment ----------
    if any(tied)
        info.nTies = sum(tied);

        if ~any(keep)
            % Every pair is equidistant - the paper's Case 3, where the
            % first face is chosen arbitrarily ("we randomly select one out
            % of the cell faces"). The Cartesian release exposes that as
            % 'defaultFace' with default 'bottom'; the unstructured
            % equivalent is the face whose outward normal is most aligned
            % with the corresponding direction. Seeding it here means the
            % node rule below then reproduces Table 1 exactly, so Case 3
            % agrees with the reference instead of making its own arbitrary
            % choice.
            seed = defaultSeed(opt);
            cand = find(tied);
            [~, kk] = max(m(cand,:) * seed(:));
            q0 = cand(kk);
            keep(q0) = true;
            tied(q0) = false;
            if partner(q0) > 0
                tied(partner(q0)) = false;
                keep(partner(q0)) = false;
            end
        end

        nStar = farthestNode(G, faces, keep, n, planePoint);
        for q = 1:nF
            if ~tied(q), continue; end
            g = partner(q);
            [keep(q), used] = tieWins(G, faces, q, g, nStar);
            if ~used, info.nFallback = info.nFallback + 1; end
        end
    end

    selected = faces(keep);
end


% ========================================================================
%                          LOCAL HELPER FUNCTIONS
% ========================================================================

function seed = defaultSeed(opt)
% Direction seeding the paper's Case 3, the unstructured analogue of the
% Cartesian release's 'defaultFace' option.
    if isfield(opt, 'defaultFace') && ~isempty(opt.defaultFace)
        name = lower(opt.defaultFace);
    else
        name = 'bottom';
    end
    switch name
        case 'bottom', seed = [0 0 -1];
        case 'top',    seed = [0 0  1];
        case 'left',   seed = [-1 0 0];
        case 'right',  seed = [ 1 0 0];
        case 'front',  seed = [0 -1 0];
        case 'back',   seed = [0  1 0];
        otherwise
            error('selectProjectionFacesPEDFM_U:defaultFace', ...
                ['Unknown defaultFace ''%s'' (use bottom, top, left, ', ...
                 'right, front or back).'], name);
    end
end


function nStar = farthestNode(G, faces, keep, n, planePoint)
% The node, among the already-selected faces, farthest from the fracture
% plane. Where two or more faces are already selected the candidates are
% the nodes they share (the corner the paper's Table 1 keys on); if they
% share none, all their nodes are candidates. With nothing selected yet,
% every node of the cell's faces is a candidate.
    sel = find(keep);
    if isempty(sel), sel = 1:numel(faces); end

    cand = faceNodes(G, faces(sel(1)));
    for k = 2:numel(sel)
        nxt = faceNodes(G, faces(sel(k)));
        shared = intersect(cand, nxt);
        if isempty(shared)
            cand = union(cand, nxt);
        else
            cand = shared;
        end
    end
    if isempty(cand), nStar = 0; return; end

    dd = abs((G.nodes.coords(cand,:) - planePoint) * n');
    [~, k] = max(dd);
    nStar = cand(k);
end


function [win, usedNode] = tieWins(G, faces, q, g, nStar)
% Keep the face of an equidistant pair that contains nStar (the paper's
% Table 1, read as hexahedron incidence). Falls back to the lower global
% face index, which is antisymmetric, so exactly one face of the pair is
% always kept.
    if nStar > 0
        inQ = any(faceNodes(G, faces(q)) == nStar);
        inG = any(faceNodes(G, faces(g)) == nStar);
        if inQ && ~inG
            win = true;  usedNode = true;  return
        elseif inG && ~inQ
            win = false; usedNode = true;  return
        end
    end
    win = faces(q) < faces(g);
    usedNode = false;
end


function nd = faceNodes(G, f)
    nd = G.faces.nodes(G.faces.nodePos(f):(G.faces.nodePos(f+1)-1));
end
