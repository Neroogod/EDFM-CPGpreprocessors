function E = calculateProjectionCPEDFM_U(G, plane, hostCells, fCells, tol, useFullFace)
% Flip-cut projection of ONE fracture plane onto unstructured cell faces.
%
% SYNOPSIS:
%   E = calculateProjectionCPEDFM_U(G, plane, hostCells, fCells, tol, ...
%                                   useFullFace)
%
% Unstructured counterpart of calculateProjectionCPEDFM.m (Cartesian) and
% calculateProjectionCPEDFM_CPG.m (corner-point). Called once per fracture
% plane by pMatFracNNCsCPEDFM_U, which concatenates the returned entries
% and appends them to G.nnc.
%
% PARAMETERS:
%   G           - global EDFM grid (after EDFMGridU and
%                 fracturematrixShaleNNC3D).
%   plane       - struct describing the fracture plane:
%                   .normal - unit normal as stored by the gridding stage
%                             (G.nnc.normal), keeping one sign convention
%                   .point  - a point on the plane (G.nnc.planepoint)
%                   .poly   - fracture polygon vertices (FracGrid points),
%                             used by the tip guard
%                   .number - fracture-grid number, recorded in E.plane
%   hostCells   - matrix cells hosting this plane's fracture cells.
%   fCells      - matching fracture cells.
%   tol         - geometric tolerance.
%   useFullFace - true for 'fullface' area mode (paper rule), false for
%                 'partialsum' (clipped projection).
%
% RETURNS:
%   E - projection entries for this plane, one row per sealed face:
%         .MFcells  [projMatrixCell, fractureCell]
%         .MMneighs [hostCell projCell]
%         .face     global matrix face carrying the projection
%         .area     projected area assigned to that face
%         .T        pM-F transmissibility (Eq. 10)
%         .dir      dominant axis of the face normal (plot colouring only)
%         .plane    fracture-grid number (= plane.number)
%
% ALGORITHM - on unstructured grids the structured-grid machinery of the
% CPEDFM paper (orientation branch table, IJK local neighborhoods,
% first/last-member collapse) has no meaning: cells have arbitrarily many
% faces with arbitrary normals and there are no logical directions. All of
% that machinery approximates one thing - a continuous set of sealed faces
% - which is computed here directly:
%
%   SEAL every face of a fracture-hosting cell whose two neighboring cell
%   centroids lie on opposite sides of the fracture surface (the sign of
%   the signed distance flips), provided the face lies within the
%   fracture's lateral extent (tip guard).
%
% The sign-flip set is by construction a cut of the matrix connectivity
% graph separating the two sides of the fracture, so the projection is
% continuous - provably, for any cell shape - which is exactly the property
% CPEDFM constructs on Cartesian grids.
%
% TIP GUARD - a face is only sealed if its centroid, projected onto the
% fracture plane, lies inside the fracture polygon (with a buffer of half
% the face's characteristic length), so a fracture spanning part of the
% domain does not seal beyond its extent.

nMat    = G.Matrix.cells.num;
permMat = G.rock.perm;

nrm        = plane.normal(:)';  nrm = nrm / norm(nrm);
planePoint = plane.point(:)';
fracPoly   = plane.poly;

E.MFcells  = zeros(0,2); E.MMneighs = zeros(0,2);
E.face = zeros(0,1); E.area = zeros(0,1); E.T = zeros(0,1);
E.dir  = zeros(0,1); E.plane = zeros(0,1);

% in-plane orthonormal basis for the tip-guard polygon test
[~, kmin] = min(abs(nrm));
e = zeros(1,3); e(kmin) = 1;
u = cross(nrm, e); u = u / norm(u);
v = cross(nrm, u);
polyU = (fracPoly - planePoint) * u';
polyV = (fracPoly - planePoint) * v';

% signed distance of matrix centroids to the fracture plane
sdMat = (G.cells.centroids(1:nMat,:) - planePoint) * nrm';

isHost = false(nMat,1); isHost(hostCells) = true;
hostFcell = zeros(nMat,1); hostFcell(hostCells) = fCells;

planeSealed = false(G.faces.num, 1);

for r = 1:numel(hostCells)
    hc = hostCells(r);
    fpos = G.cells.facePos(hc):(G.cells.facePos(hc+1)-1);
    faces = G.cells.faces(fpos, 1);

    for q = 1:numel(faces)
        f = faces(q);
        if planeSealed(f), continue; end
        nbrs = G.faces.neighbors(f,:);
        other = nbrs(nbrs ~= hc);
        if numel(other) ~= 1 || other == 0 || other > nMat
            continue % boundary face or fracture-grid face
        end
        if sdMat(hc) * sdMat(other) > 0
            continue % both centroids on the same side: no crossing
        end

        % tip guard: face centroid projected onto the fracture plane
        % must lie inside the fracture polygon, with a buffer of half
        % the face's characteristic length
        fc = G.faces.centroids(f,:);
        pu = (fc - planePoint) * u';
        pv = (fc - planePoint) * v';
        buf = 0.5 * sqrt(G.faces.areas(f));
        if ~inpolygon(pu, pv, polyU, polyV) && ...
           dist2poly(pu, pv, polyU, polyV) > buf
            continue
        end

        % Fig. 13 pairing: pM = cell across the fracture from the
        % segment; F = fracture cell of the other side's host.
        if ~isHost(other)
            pM = other; hostM = hc;
        else
            % both host: deterministic global convention - the cell
            % on the positive-normal side is the projection cell
            if sdMat(other) >= sdMat(hc)
                pM = other; hostM = hc;
            else
                pM = hc; hostM = other;
            end
        end
        fcM = hostFcell(hostM);

        fnrm = G.faces.normals(f,:); fnrm = fnrm / norm(fnrm);
        if useFullFace
            aShare = G.faces.areas(f);
        else
            aShare = clippedProjArea(G, fcM, f, fnrm, tol);
            if aShare <= tol^2, continue; end
        end

        dist  = max(norm(G.cells.centroids(pM,:) - G.cells.centroids(fcM,:)), tol);
        kF = permNormal(permMat(fcM,:), fnrm);
        kM = permNormal(permMat(pM,:),  fnrm);
        if kF + kM > 0
            kHarm = 2*kF*kM / (kF + kM);
        else
            kHarm = 0; % degenerate (zero-perm) cells: keep the seal,
                       % no fracture exchange through it
        end

        [~, domAx] = max(abs(fnrm));

        planeSealed(f) = true;
        E.MFcells (end+1,:) = [pM, fcM];             %#ok<AGROW>
        E.MMneighs(end+1,:) = [hostM, pM];           %#ok<AGROW>
        E.face    (end+1,1) = f;                     %#ok<AGROW>
        E.area    (end+1,1) = aShare;                %#ok<AGROW>
        E.T       (end+1,1) = kHarm*aShare/dist;     %#ok<AGROW>
        E.dir     (end+1,1) = domAx;                 %#ok<AGROW>
        E.plane   (end+1,1) = plane.number;          %#ok<AGROW>
    end
end
end


% ========================================================================
%                          LOCAL HELPER FUNCTIONS
% ========================================================================

function k = permNormal(permRow, n)
% normal-projected permeability n'*K*n for 1-, 3- or 6-component perm
switch numel(permRow)
    case 1
        k = permRow;
    case 3
        k = permRow(1)*n(1)^2 + permRow(2)*n(2)^2 + permRow(3)*n(3)^2;
    case 6
        K = [permRow(1) permRow(2) permRow(3);
             permRow(2) permRow(4) permRow(5);
             permRow(3) permRow(5) permRow(6)];
        k = n*K*n';
    otherwise
        error('Unsupported permeability format (%d components).', numel(permRow));
end
end


function d = dist2poly(pu, pv, polyU, polyV)
% distance from point (pu,pv) to a 2D polygon boundary
n = numel(polyU);
d = inf;
for i = 1:n
    j = mod(i, n) + 1;
    a = [polyU(i), polyV(i)]; b = [polyU(j), polyV(j)];
    ab = b - a; p = [pu, pv] - a;
    t = max(min(dot(p, ab)/max(dot(ab, ab), eps), 1), 0);
    d = min(d, norm(p - t*ab));
end
end


function a = clippedProjArea(G, fCell, f, fnrm, tol)
% area of the fracture cell's mid-plane polygon projected onto face f's
% plane and clipped against the face polygon
[cn, ~] = gridCellNodes(G, fCell);
xyz = G.nodes.coords(cn,:);
nV = size(xyz,1);
if mod(nV,2) ~= 0 || nV < 6, a = 0; return; end
mid = (xyz(1:nV/2,:) + xyz(nV/2+1:end,:)) / 2;

% order and project both polygons into face-plane coordinates
[~, kmin] = min(abs(fnrm));
e = zeros(1,3); e(kmin) = 1;
u = cross(fnrm, e); u = u / norm(u);
v = cross(fnrm, u);
fc = G.faces.centroids(f,:);

mid = orderCCW(mid);
P2 = [(mid - fc)*u', (mid - fc)*v'];

fnodes = G.faces.nodes(G.faces.nodePos(f):G.faces.nodePos(f+1)-1);
fx = G.nodes.coords(fnodes,:);
F2 = [(fx - fc)*u', (fx - fc)*v'];

try
    psF = polyshape(F2(:,1), F2(:,2), 'Simplify', true);
    psP = polyshape(P2(:,1), P2(:,2), 'Simplify', true);
    a = area(intersect(psF, psP));
catch
    a = 0;
end
if a < tol^2, a = 0; end
end


function xyz = orderCCW(xyz)
c = mean(xyz, 1);
P = xyz - c;
[~, ~, V] = svd(P, 0);
[~, is] = sort(atan2(P*V(:,1), P*V(:,2)));
xyz = xyz(is, :);
end
