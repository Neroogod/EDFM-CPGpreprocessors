function E = calculateProjectionPEDFM_U(G, plane, hostCells, fCells, fmAreas, tol, opt)
%CALCULATEPROJECTIONPEDFM_U pEDFM projection of ONE fracture plane, unstructured grid.
%
% SYNOPSIS:
%   E = calculateProjectionPEDFM_U(G, plane, hostCells, fCells, fmAreas, tol, opt)
%
% Unstructured counterpart of the per-plane projection stage of the 3D
% pEDFM of Olorode, Wang & Rashid (2020), SPE-201243-PA. Called once per
% fracture plane by pMatFracNNCsPEDFM_U.
%
% PARAMETERS:
%   G         - global EDFM grid (after EDFMGridU and fracturematrixShaleNNC3D).
%   plane     - struct with fields:
%                 .normal - unit normal as stored by the gridding stage
%                 .point  - a point on the plane
%                 .poly   - fracture polygon vertices, used by the tip guard
%                 .number - fracture-grid number, recorded in E.plane
%   hostCells - matrix cells hosting this plane's fracture cells.
%   fCells    - matching fracture cells.
%   fmAreas   - matching fracture-matrix contact areas (G.nnc.area).
%   tol       - geometric tolerance.
%   opt       - options struct (see pMatFracNNCsPEDFM_U).
%
% RETURNS:
%   E - projection entries, one row per (selected face, fracture cell):
%         .MFcells  [projMatrixCell, fractureCell]
%         .MMneighs [hostCell, projCell]
%         .face     global matrix face carrying the projection
%         .area     projected area (Eq. 24 numerator)
%         .T        pM-F transmissibility (Eqs. 24-25)
%         .dir      dominant axis of the face normal (plot colouring only)
%         .plane    fracture-grid number
%         .nTies    equidistant antipodal pairs encountered
%
% ALGORITHM
%   1. Face selection - the antipodal-pair rule of
%      selectProjectionFacesPEDFM_U, which reproduces the paper's Cartesian
%      selection exactly on hexahedra and generalizes to any convex
%      polyhedron.
%   2. Tip guard - a selected face only carries a projection if its
%      centroid, projected onto the fracture plane, lies inside the
%      fracture polygon with a buffer of half the face's characteristic
%      length. Without it a fracture spanning part of the domain would
%      project beyond its own extent, because the selection rule sees only
%      the infinite plane. Same construction as CPEDFM-U.
%   3. Areas and transmissibilities - Eq. 24 with the real face normal,
%      Eq. 25 with normal-projected permeabilities n'*K*n so that
%      anisotropic rock is handled on faces of arbitrary orientation.
%
% NOTE ON pDir - on a structured grid this is the logical direction. There
% is no such thing here, so it holds the dominant axis of the face normal
% and is used only for plot colouring.

    nMat    = G.Matrix.cells.num;
    permMat = G.rock.perm;

    nrm        = plane.normal(:)';  nrm = nrm / norm(nrm);
    planePoint = plane.point(:)';
    fracPoly   = plane.poly;

    E.MFcells = zeros(0,2); E.MMneighs = zeros(0,2);
    E.face = zeros(0,1); E.area = zeros(0,1); E.T = zeros(0,1);
    E.dir  = zeros(0,1); E.plane = zeros(0,1);
    E.nTies = 0; E.nFallback = 0;

    % In-plane orthonormal basis for the tip-guard polygon test.
    [~, kmin] = min(abs(nrm));
    e = zeros(1,3); e(kmin) = 1;
    u = cross(nrm, e); u = u / norm(u);
    v = cross(nrm, u);
    polyU = (fracPoly - planePoint) * u';
    polyV = (fracPoly - planePoint) * v';

    for r = 1:numel(hostCells)
        hc = hostCells(r);
        fc = fCells(r);
        fmArea = fmAreas(r);

        if G.cells.volumes(hc) <= tol^3, continue; end

        [selected, info] = selectProjectionFacesPEDFM_U(G, hc, ...
            G.cells.centroids(fc,:), nrm, planePoint, opt);
        E.nTies = E.nTies + info.nTies;
        E.nFallback = E.nFallback + info.nFallback;

        for q = 1:numel(selected)
            f  = selected(q);
            nb = G.faces.neighbors(f,:);
            projCell = nb(nb ~= hc & nb ~= 0);
            if numel(projCell) ~= 1 || projCell > nMat
                continue    % domain boundary or a fracture-grid face
            end

            fnrm = G.faces.normals(f,:);
            if norm(fnrm) <= eps, continue; end
            fnrm = fnrm / norm(fnrm);

            % Eq. 24 numerator, using the actual face normal.
            aProj = fmArea * abs(dot(fnrm, nrm));
            if aProj <= tol^2, continue; end
            aProj = min(aProj, G.faces.areas(f));

            % Tip guard: does the fracture actually reach this face?
            fcx = G.faces.centroids(f,:);
            pu = (fcx - planePoint) * u';
            pv = (fcx - planePoint) * v';
            buf = 0.5 * sqrt(G.faces.areas(f));
            if ~inpolygon(pu, pv, polyU, polyV) && ...
                    dist2poly(pu, pv, polyU, polyV) > buf
                continue
            end

            % Eq. 25 with normal-projected permeabilities.
            kF = permNormal(permMat(fc,:), fnrm);
            kM = permNormal(permMat(projCell,:), fnrm);
            if (kF + kM) > 0
                kHarm = 2*kF*kM / (kF + kM);
            else
                kHarm = 0;  % both sealing: keep the seal, no exchange
            end
            dist = max(norm(G.cells.centroids(projCell,:) - ...
                            G.cells.centroids(fc,:)), tol);

            [~, domAx] = max(abs(fnrm));

            E.MFcells (end+1,:) = [projCell, fc];
            E.MMneighs(end+1,:) = [hc, projCell];
            E.face    (end+1,1) = f;
            E.area    (end+1,1) = aProj;
            E.T       (end+1,1) = kHarm*aProj/dist;
            E.dir     (end+1,1) = domAx;
            E.plane   (end+1,1) = plane.number;
        end
    end
end


% ========================================================================
%                          LOCAL HELPER FUNCTIONS
% ========================================================================

function k = permNormal(permRow, n)
% Normal-projected permeability n'*K*n for 1-, 3- or 6-component perm.
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
            error('calculateProjectionPEDFM_U:perm', ...
                'Unsupported permeability format (%d components).', numel(permRow));
    end
end


function d = dist2poly(pu, pv, polyU, polyV)
% Distance from a point to a 2D polygon boundary.
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
