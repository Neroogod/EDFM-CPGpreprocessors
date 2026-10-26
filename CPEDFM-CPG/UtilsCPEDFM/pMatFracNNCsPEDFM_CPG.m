function G = pMatFracNNCsPEDFM_CPG(G, tol)
%PMATFRACNNCSPEDFM_CPG Conventional pEDFM projection on CPG geometry.
%
% This is the corner-point counterpart of pMatFracNNCs3D. It deliberately
% retains pEDFM's independent per-host-cell face selection (no CPEDFM
% neighborhood collapse and no continuity/flip-face completion), while
% resolving logical sides to actual CPG face IDs. The resulting pMMfaces
% records are consumed by setupCPEDFMOpsTPFA_CPG/transmultpEDFM_CPG.

if nargin < 2, tol = 1e-5; end
assert(isfield(G, 'Matrix') && isfield(G, 'FracGrid'), ...
    'Run EDFMGridCPG and fracturematrixShaleNNC3D first.');

nMat = G.Matrix.cells.num;
interior = find(strcmp(G.nnc.type, 'fracmat interior'));
if isempty(interior)
    warning('pMatFracNNCsPEDFM_CPG:noInterior', ...
        'No interior fracture-matrix NNCs found; nothing to project.');
    return
end

fracNames = fieldnames(G.FracGrid);
perm = G.rock.perm;
nPermCol = size(perm, 2);

% Orientation and face-tag helpers live as local functions of
% calculateProjectionCPEDFM_CPG; fetch handles once.
h = calculateProjectionCPEDFM_CPG('helpers');

MF = zeros(0,2); MM = zeros(0,2); face = zeros(0,1);
area = zeros(0,1); trans = zeros(0,1); dir = zeros(0,1);
plane = zeros(0,1);

for fn = 1:numel(fracNames)
    F = G.FracGrid.(fracNames{fn});
    fFirst = F.cells.start;
    fLast = fFirst + F.cells.num - 1;
    rows = interior(G.nnc.cells(interior,2) >= fFirst & ...
                    G.nnc.cells(interior,2) <= fLast);
    if isempty(rows), continue; end

    nrm = G.nnc.normal(rows(1),:);
    nrm = nrm/norm(nrm);
    tags = h.fracOrientation(nrm);

    for rr = reshape(rows, 1, [])
        host = G.nnc.cells(rr,1);
        fc = G.nnc.cells(rr,2);
        fmArea = G.nnc.area(rr);
        fs = h.getLogicalFaces(G, host);

        for d = 1:3
            selected = fs{tags(d)};
            if isempty(selected), continue; end
            isNeg = mod(tags(d),2) == 1;

            % Project the local fracture polygon onto the selected CPG
            % side. Split its area over fault-split members according to
            % their actual face areas.
            nf = sum(G.faces.normals(selected,:), 1);
            if norm(nf) <= eps, continue; end
            projected = fmArea*abs(dot(nf, nrm))/norm(nf);
            if projected <= tol^2, continue; end
            weights = G.faces.areas(selected);
            weights = weights/sum(weights);

            for q = 1:numel(selected)
                f = selected(q);
                nb = G.faces.neighbors(f,:);
                proj = nb(nb ~= host & nb ~= 0);
                if numel(proj) ~= 1 || proj > nMat, continue; end

                a = projected*weights(q);
                kCol = min(d, nPermCol);
                kSum = perm(fc,kCol) + perm(proj,kCol);
                if kSum > 0
                    kh = 2*perm(fc,kCol)*perm(proj,kCol)/kSum;
                else
                    kh = 0;
                end
                dist = max(norm(G.cells.centroids(proj,:) - ...
                                G.cells.centroids(fc,:)), tol);
                if isNeg
                    pair = [proj, host];
                else
                    pair = [host, proj];
                end

                MF(end+1,:) = [proj, fc]; %#ok<AGROW>
                MM(end+1,:) = pair;       %#ok<AGROW>
                face(end+1,1) = f;        %#ok<AGROW>
                area(end+1,1) = a;        %#ok<AGROW>
                trans(end+1,1) = kh*a/dist; %#ok<AGROW>
                dir(end+1,1) = d;         %#ok<AGROW>
                plane(end+1,1) = fn;      %#ok<AGROW>
            end
        end
    end
end

if isempty(face)
    warning('pMatFracNNCsPEDFM_CPG:noProjections', ...
        'No pEDFM projections were generated.');
    return
end

% Sum duplicate contributions to the same face/fracture cell projection.
key = [face, MF];
[~, ~, ic] = unique(key, 'rows', 'stable');
if numel(unique(ic)) < numel(ic)
    keep = accumarray(ic, (1:numel(ic))', [], @min);
    areaSum = accumarray(ic, area);
    transSum = accumarray(ic, trans);
    MF = MF(keep,:); MM = MM(keep,:); face = face(keep);
    dir = dir(keep); plane = plane(keep);
    area = areaSum; trans = transSum;
end

nNew = numel(face);
G.nnc.cells = [G.nnc.cells; MF];
G.nnc.T = [G.nnc.T; trans];
G.nnc.type = [G.nnc.type; repmat({'fracpmat'}, nNew, 1)];
G.nnc.area = [G.nnc.area; area];
G.nnc.pMMneighs = MM;
G.nnc.pMMfaces = face;
G.nnc.pMMarea = area;
G.nnc.pDir = dir;
G.nnc.pFracPlane = plane;
G.nnc.nncFaces = [face, dir, plane];
end
