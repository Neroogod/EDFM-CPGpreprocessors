function E = calculateProjectionCPEDFM_CPG(G, plane, hostCells, fCells, fmAreas, tol, useFullFace)
% Continuous projection of ONE fracture plane onto corner-point cell faces.
%
% SYNOPSIS:
%   E = calculateProjectionCPEDFM_CPG(G, plane, hostCells, fCells, ...
%                                     fmAreas, tol, useFullFace)
%   h = calculateProjectionCPEDFM_CPG('helpers')
%
% Corner-point counterpart of calculateProjectionCPEDFM.m in the Cartesian
% release. Called once per fracture plane by pMatFracNNCsCPEDFM_CPG, which
% concatenates the returned entries and appends them to G.nnc.
%
% PARAMETERS:
%   G           - global EDFM grid (after EDFMGridCPG / EDFMshalegrid and
%                 fracturematrixShaleNNC3D).
%   plane       - struct describing the fracture plane:
%                   .normal - unit normal (sign convention as stored by
%                             the gridding stage in G.nnc.normal)
%                   .point  - a point on the plane (G.nnc.planepoint)
%                   .number - fracture-grid number, recorded in E.plane
%   hostCells   - column vector of matrix cells hosting this plane's
%                 fracture cells (G.nnc.cells(rows,1)).
%   fCells      - matching fracture cells (G.nnc.cells(rows,2)).
%   fmAreas     - matching fracture-matrix contact areas (G.nnc.area(rows)).
%   tol         - geometric tolerance.
%   useFullFace - true for 'fullface' area mode (paper rule), false for
%                 'partialsum' (Cartesian-release reproduction).
%
% RETURNS:
%   E - projection entries for this plane, one row per sealed face:
%         .MFcells  [projMatrixCell, fractureCell]
%         .MMneighs [hostCell projCell] in the side-order convention
%         .face     global matrix face carrying the projection
%         .area     projected area assigned to that face
%         .T        pM-F transmissibility (Eq. 10)
%         .dir      logical direction 1..3
%         .plane    fracture-grid number (= plane.number)
%
% ACCESSING THE HELPERS:
%   The four projection helpers are local subfunctions of this file, as in
%   the Cartesian release. Callers that need them directly - the
%   conventional pEDFM variant (pMatFracNNCsPEDFM_CPG), the cross-plane
%   union correction in pMatFracNNCsCPEDFM_CPG, and the unit tests - obtain
%   handles with
%       h = calculateProjectionCPEDFM_CPG('helpers');
%       tags = h.fracOrientation(nrm);
%   The returned struct has fields fracOrientation, getLogicalFaces,
%   buildNeighborhoods and getUnionProjectArea.
%
% Algorithm (Rashid & Olorode, SPE-217469-PA, Figs. 8-13); see the header
% of pMatFracNNCsCPEDFM_CPG for the full description and the list of
% deliberate deviations from the Cartesian release.

% --- helper-handle accessor (keeps the local functions reachable) --------
if nargin == 1 && (ischar(G) || isstring(G)) && strcmp(G, 'helpers')
    E = struct('fracOrientation',     @fracOrientation, ...
               'getLogicalFaces',     @getLogicalFaces, ...
               'buildNeighborhoods',  @buildNeighborhoods, ...
               'getUnionProjectArea', @getUnionProjectArea);
    return
end

nMat     = G.Matrix.cells.num;
permMat  = G.rock.perm;
nPermCol = size(permMat, 2);
nHost    = numel(hostCells);

nrm = plane.normal(:)';
nrm = nrm / norm(nrm);

E.MFcells  = zeros(0,2);
E.MMneighs = zeros(0,2);
E.face     = zeros(0,1);
E.area     = zeros(0,1);
E.T        = zeros(0,1);
E.dir      = zeros(0,1);
E.plane    = zeros(0,1);

% Step 1: face-tag selection per direction (Fig. 8 branch table)
selTags = fracOrientation(nrm);

% Step 2: per-host-cell selected face sets and partial projected areas
faceSel   = cell(nHost, 3);
pAreaPart = zeros(nHost, 3);
for r = 1:nHost
    fs = getLogicalFaces(G, hostCells(r));
    for d = 1:3
        faces = fs{selTags(d)};
        faceSel{r,d} = faces;
        if isempty(faces)
            % Pinched/degenerate side: no face, but the fracture
            % segment still shades this direction. Fall back to the
            % axis-aligned cosine so the neighborhood area stays
            % consistent.
            pAreaPart(r,d) = fmAreas(r) * abs(nrm(d));
        else
            % Area-weighted mean normal of the (near-coplanar) set.
            % G.faces.normals are area-scaled in MRST, so the plain
            % sum is the area weighting.
            nf = sum(G.faces.normals(faces,:), 1);
            pAreaPart(r,d) = fmAreas(r) * abs(dot(nf, nrm)) / norm(nf);
        end
    end
end

% Step 3: local neighborhoods from logical IJK
nb = buildNeighborhoods(G.Matrix, hostCells);

% Signed distance of matrix cell centroids to the fracture plane,
% used by the flip-face rule.
sdMat = (G.cells.centroids(1:nMat,:) - plane.point) * nrm';
isHost = false(nMat,1); isHost(hostCells) = true;
hostFcell = zeros(nMat,1); hostFcell(hostCells) = fCells;

% Step 4: seal envelope faces + flip faces, one NNC per sealed face.
% sealedFace bookkeeping avoids duplicate entries when the envelope
% and flip rules (or two host cells sharing a flip face) pick the
% same face.
planeSealed = false(G.faces.num, 1);

for d = 1:3
    isNegSide = mod(selTags(d), 2) == 1;
    kCol = min(d, nPermCol);

    % ---- (i) envelope faces (paper's collapse rule) ----
    for g = 1:numel(nb{d})
        members = nb{d}{g};
        if isNegSide
            rep = members(1);
        else
            rep = members(end);
        end
        keptFaces = faceSel{rep,d};
        mergedArea = sum(pAreaPart(members, d));
        if isempty(keptFaces) || mergedArea <= tol^2
            % Nothing to seal (degenerate side or negligible
            % projection). The Cartesian release likewise drops these.
            continue
        end

        % Fracture segment centroid: midpoint of the first and last
        % member's fracture-cell centroids (as in the original).
        segCntrd = (G.cells.centroids(fCells(members(1)),:) + ...
                    G.cells.centroids(fCells(members(end)),:)) / 2;
        fcRep = fCells(rep);

        keptAreas = G.faces.areas(keptFaces);
        totKept   = sum(keptAreas);
        for q = 1:numel(keptFaces)
            f = keptFaces(q);
            if planeSealed(f), continue; end
            nbrs = G.faces.neighbors(f,:);
            projCell = nbrs(nbrs ~= hostCells(rep));
            if numel(projCell) ~= 1 || projCell == 0 || projCell > nMat
                % Domain boundary: no projection cell outside the
                % host -> drop ('ProjMatrix outside' in the original).
                continue
            end

            if useFullFace
                aShare = keptAreas(q); % paper rule (Figs. 10, 13)
            else
                aShare = mergedArea * keptAreas(q) / totKept;
            end
            dist  = max(norm(G.cells.centroids(projCell,:) - segCntrd), tol);
            kSum  = permMat(fcRep,kCol) + permMat(projCell,kCol);
            if kSum > 0
                kHarm = 2 * (permMat(fcRep,kCol) * permMat(projCell,kCol)) / kSum;
            else
                kHarm = 0; % degenerate (zero-perm) cells: keep the
                           % seal, no fracture exchange through it
            end

            if isNegSide
                mmPair = [projCell, hostCells(rep)];
            else
                mmPair = [hostCells(rep), projCell];
            end

            planeSealed(f) = true;
            E.MFcells (end+1,:) = [projCell, fcRep];       %#ok<AGROW>
            E.MMneighs(end+1,:) = mmPair;                  %#ok<AGROW>
            E.face    (end+1,1) = f;                       %#ok<AGROW>
            E.area    (end+1,1) = aShare;                  %#ok<AGROW>
            E.T       (end+1,1) = kHarm*aShare/dist;       %#ok<AGROW>
            E.dir     (end+1,1) = d;                       %#ok<AGROW>
            E.plane   (end+1,1) = plane.number;            %#ok<AGROW>
        end
    end

    % ---- (ii) flip faces (watertight-cut rule) ----
    % Any direction-d face of a host cell whose neighbors straddle
    % the fracture plane. Checked on both tags of direction d so the
    % cut is complete regardless of the branch side.
    for r = 1:nHost
        hc = hostCells(r);
        for tag = [2*d-1, 2*d]
            if tag == selTags(d)
                faces = faceSel{r,d};
            else
                fsAll = getLogicalFaces(G, hc);
                faces = fsAll{tag};
            end
            for q = 1:numel(faces)
                f = faces(q);
                if planeSealed(f), continue; end
                nbrs = G.faces.neighbors(f,:);
                other = nbrs(nbrs ~= hc);
                if numel(other) ~= 1 || other == 0 || other > nMat
                    continue
                end
                if sdMat(hc) * sdMat(other) > 0
                    continue % same side: fracture does not cross here
                end

                % Fig. 13 pairing: the projected NNC connects the
                % fracture segment to the cell whose direct
                % matrix-matrix path this face seal severs.
                if ~isHost(other)
                    % one host, one plain cell: the plain cell is the
                    % projection cell, paired with the host's segment
                    pM = other; hostM = hc;
                else
                    % both cells host the fracture: take the cell on
                    % the branch-selected grid side as the projection
                    % cell, paired with the other side's segment
                    cmp = G.cells.centroids(other,d) - G.cells.centroids(hc,d);
                    if (cmp > 0) == ~isNegSide
                        pM = other; hostM = hc;
                    else
                        pM = hc; hostM = other;
                    end
                end
                fcM = hostFcell(hostM);

                if useFullFace
                    aShare = G.faces.areas(f);
                else
                    aShare = pAreaPart(r,d);
                end
                dist  = max(norm(G.cells.centroids(pM,:) - ...
                             G.cells.centroids(fcM,:)), tol);
                kSum  = permMat(fcM,kCol) + permMat(pM,kCol);
                if kSum > 0
                    kHarm = 2 * (permMat(fcM,kCol) * permMat(pM,kCol)) / kSum;
                else
                    kHarm = 0;
                end

                if isNegSide
                    mmPair = [pM, hostM];
                else
                    mmPair = [hostM, pM];
                end

                planeSealed(f) = true;
                E.MFcells (end+1,:) = [pM, fcM];           %#ok<AGROW>
                E.MMneighs(end+1,:) = mmPair;              %#ok<AGROW>
                E.face    (end+1,1) = f;                   %#ok<AGROW>
                E.area    (end+1,1) = aShare;              %#ok<AGROW>
                E.T       (end+1,1) = kHarm*aShare/dist;   %#ok<AGROW>
                E.dir     (end+1,1) = d;                   %#ok<AGROW>
                E.plane   (end+1,1) = plane.number;        %#ok<AGROW>
            end
        end
    end
end
end


% ========================================================================
%                          LOCAL HELPER FUNCTIONS
% ========================================================================

function [selTags, dip, dipdir, azm] = fracOrientation(planeNormal, tolAng)
% Orientation-based selection of the three projection face tags (Fig. 8 of
% Rashid & Olorode, SPE-217469-PA), one per logical grid direction.
%
% PARAMETERS:
%   planeNormal - fracture plane normal (need not be unit length). Use the
%                 normal stored by the gridding stage (G.nnc.normal /
%                 fracplanes.normal) so the sign convention matches the
%                 rest of the pipeline.
%   tolAng      - angular tolerance in degrees for the trivial
%                 horizontal/vertical branches (default 1e-8).
%
% RETURNS:
%   selTags - 1x3 vector of selected face tags, selTags(d) in {2d-1, 2d}:
%             selTags(1) in {1,2} (left/right),
%             selTags(2) in {3,4} (front/back),
%             selTags(3) in {5,6} (bottom/top).
%   dip     - dip angle in [0,90] degrees.
%   dipdir  - dip direction in [0,360) degrees.
%   azm     - atand(nz), the signed z-inclination used by the flow chart.
%
% The branch table below is ported verbatim from
% calculateProjectionCPEDFM.m (lines 44-89) so that the CPG implementation
% selects exactly the faces the Cartesian CPEDFM release selects.

if nargin < 2, tolAng = 1e-8; end

n = planeNormal(:)'/norm(planeNormal);
dipdir = mod(atan2d(n(2), n(1)), 360);
dip    = 90 - asind(abs(n(3)));
azm    = atand(n(3));

if abs(dip) < tolAng && abs(dipdir) < tolAng
    % horizontal fracture (trivial case)
    selTags = [1 3 5];
elseif abs(dip-90) < tolAng && abs(dipdir) < tolAng
    % vertical fracture (trivial case)
    selTags = [1 3 5];
elseif dipdir >= 0 && dipdir <= 90 && azm > 0        % branch 1
    selTags = [1 3 5];
elseif dipdir >= 0 && dipdir <= 90 && azm <= 0       % branch 2
    selTags = [1 3 6];
elseif dipdir > 90 && dipdir <= 180 && azm > 0       % branch 3
    selTags = [2 3 5];
elseif dipdir > 90 && dipdir <= 180 && azm <= 0      % branch 4
    selTags = [1 4 5];
elseif dipdir > 180 && dipdir <= 270 && azm > 0      % branch 5
    selTags = [2 4 5];
elseif dipdir > 180 && dipdir <= 270 && azm <= 0     % branch 6
    selTags = [1 3 5];
elseif dipdir > 270 && dipdir < 360 && azm > 0       % branch 7
    selTags = [1 4 5];
elseif dipdir > 270 && dipdir < 360 && azm <= 0      % branch 8
    selTags = [2 3 5];
else
    error('calculateProjectionCPEDFM_CPG:orientation', ...
        'Could not classify fracture orientation (dip=%g, dipdir=%g, azm=%g).', ...
        dip, dipdir, azm);
end
end


function faceSets = getLogicalFaces(G, c)
% Group the faces of cell c by their logical direction tag.
%
% PARAMETERS:
%   G - grid with G.cells.faces carrying tags in column 2 (1=left/West,
%       2=right/East, 3=front/South, 4=back/North, 5=bottom, 6=top), as
%       produced by cartGrid, tensorGrid and processGRDECL.
%   c - single cell index.
%
% RETURNS:
%   faceSets - 1x6 cell array. faceSets{t} holds the global face indices
%              of cell c with tag t. On a Cartesian grid every set is a
%              singleton. On a faulted corner-point grid a set can hold
%              several faces (fault-split interface) and can be empty
%              (pinch-out / degenerate side).
%
% NOTE: This function exists so that no caller ever indexes cell faces
% positionally. Corner-point grids do not keep cell faces sorted by tag,
% and faulted cells have more than six faces.

fpos = G.cells.facePos(c):(G.cells.facePos(c+1)-1);
f = G.cells.faces(fpos, 1);
assert(size(G.cells.faces,2) >= 2, ...
    'G.cells.faces has no tag column; grid type not supported.');
t = G.cells.faces(fpos, 2);

faceSets = cell(1,6);
for tag = 1:6
    faceSets{tag} = f(t == tag);
end
end


function nb = buildNeighborhoods(Gm, hostCells)
% Build the CPEDFM "local neighborhoods" (Fig. 9 of SPE-217469-PA) from
% logical IJK indices instead of centroid coordinates.
%
% PARAMETERS:
%   Gm        - matrix grid (G.Matrix) with cartDims and cells.indexMap.
%   hostCells - column vector of matrix cell indices hosting fracture
%               cells of ONE fracture plane.
%
% RETURNS:
%   nb - 1x3 cell array. nb{d} is a cell array of neighborhoods for
%        logical direction d. Each neighborhood is a column vector of
%        indices INTO hostCells (not cell ids), sorted by the direction-d
%        logical index (ascending).
%
%        nb{1}: cells sharing (J,K), varying I  -> left/right faces
%        nb{2}: cells sharing (I,K), varying J  -> front/back faces
%        nb{3}: cells sharing (I,J), varying K  -> bottom/top faces
%
% The Cartesian release groups host cells by exact equality of centroid
% coordinates (calculateProjectionCPEDFM.m lines 127-139). That test is
% only meaningful on Cartesian/tensor grids; on a distorted corner-point
% grid every group degenerates to a single cell and the method silently
% loses its continuity property. Logical IJK grouping is equivalent on
% Cartesian grids and remains valid under any geometric distortion,
% including across fault throws.

hostCells = hostCells(:);
[I, J, K] = ind2sub(reshape(Gm.cartDims, 1, []), Gm.cells.indexMap(hostCells));
IJK = [I(:), J(:), K(:)];

nb = cell(1,3);
for d = 1:3
    others = setdiff(1:3, d);
    [~, ~, gid] = unique(IJK(:, others), 'rows');
    groups = accumarray(gid, (1:numel(hostCells))', [], @(x){x});
    for g = 1:numel(groups)
        [~, order] = sort(IJK(groups{g}, d));
        groups{g} = groups{g}(order);
    end
    nb{d} = groups;
end
end


function aUnion = getUnionProjectArea(G, fCells, face, tol)
% Union area of fracture-cell polygons projected onto a given face plane.
%
% PARAMETERS:
%   G      - global EDFM grid.
%   fCells - fracture cell indices (global numbering) whose mid-plane
%            polygons are to be projected.
%   face   - global face index defining the projection plane (its normal
%            is used; on corner-point grids this is NOT axis-aligned).
%   tol    - geometric tolerance.
%
% RETURNS:
%   aUnion - area of the union of the projected polygons (NaN on failure;
%            callers should treat NaN as "keep the summed area").
%
% Corner-point generalization of getUnionProjectArea (subfunction of
% pMatFracNNCs3DCPEDFM.m): the original projects onto the global YZ/XZ/XY
% plane chosen by an axis flag; here the projection plane is the actual
% face plane, with an orthonormal basis built from the face normal.

nf = G.faces.normals(face,:);
nf = nf / norm(nf);

% Orthonormal in-plane basis (u,v) for the face plane
[~, kmin] = min(abs(nf));
e = zeros(1,3); e(kmin) = 1;
u = cross(nf, e); u = u / norm(u);
v = cross(nf, u);

polys = [];
for i = 1:numel(fCells)
    % Extract the fracture cell's mid-plane polygon. Fracture cells are
    % thin extruded boxes; averaging the two node layers recovers the
    % mid-plane polygon (same construction as the Cartesian release).
    [cn, ~] = gridCellNodes(G, fCells(i));
    xyz = G.nodes.coords(cn,:);
    nV = size(xyz,1);
    if mod(nV,2) ~= 0 || nV < 6
        aUnion = NaN;
        return
    end
    mid = (xyz(1:nV/2,:) + xyz(nV/2+1:end,:)) / 2;
    mid = poly3Dccw_local(mid);

    % Project onto the face plane basis
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


function xyz_order = poly3Dccw_local(xyz)
% Order coplanar 3D points along the polygon edge (SVD trick, as in the
% Cartesian release).
xyzc = mean(xyz, 1);
P = xyz - xyzc;
[~, ~, V] = svd(P, 0);
[~, is] = sort(atan2(P*V(:,1), P*V(:,2)));
xyz_order = xyz(is, :);
end
