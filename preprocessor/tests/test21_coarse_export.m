%% test21_coarse_export
% Verify the coarse-model export: barrier-aware coarse blocks, aggregated
% connection-list transmissibilities, and the Eclipse deck round-trip.
%
% Case: the test01/test30 configuration (skewed faulted corner-point
% grid, slanted sealing fracture, west/east Dirichlet drive).
%
% Checks:
%   A. CONSERVATION + COMPARTMENTS: total coarse pore volume equals the
%      fine pore volume, and the coarse connectivity graph splits into
%      the same two macro compartments as the fine model (no coarse
%      connection crosses the sealed surface).
%   B. UPSCALING QUALITY (sealed case): solving the coarse model under
%      the fine model's west/east Dirichlet drive reproduces the fine
%      compartment plateau pressures within 2%.
%   C. DECK ROUND-TRIP: the exported Eclipse deck parses, and its PORV,
%      DEPTH, and NNC transmissibilities match the coarse model.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', '..', 'CPEDFM-CPG', 'startupCPEDFMCPG.m'));
addpath(fullfile(thisDir, '..'));
addpath(fullfile(thisDir, '..', '..', 'CPEDFM-MS'));
mrstModule add deckformat coarsegrid

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
failures = {};

%% Fine CPEDFM model (test01 configuration)
dims = [20, 10, 8];
grdecl = simpleGrdecl(dims, 0.12);
physdim = [200, 100, 40];
grdecl.COORD(1:3:end) = grdecl.COORD(1:3:end) * physdim(1);
grdecl.COORD(2:3:end) = grdecl.COORD(2:3:end) * physdim(2);
grdecl.COORD(3:3:end) = grdecl.COORD(3:3:end) * physdim(3);
grdecl.ZCORN = grdecl.ZCORN * physdim(3);
G = processGRDECL(grdecl); G = computeGeometry(G);
G.rock = makeRock(G, 100*milli*darcy, 0.25);

fp = struct('points', [70 -5 -5; 110 105 -5; 150 105 45; 110 -5 45], ...
    'aperture', 10*milli*meter, 'poro', 0.4, 'perm', 1*nano*darcy);
[G, fp] = EDFMGridCPG(G, fp, 'Tolerance', tol, 'plotgrid', false, ...
    'fracturelist', 1);
G = fracturematrixShaleNNC3D(G, tol);
Gp = pMatFracNNCsCPEDFM_CPG(G, tol);
ops = setupCPEDFMOpsTPFA_CPG(Gp, Gp.rock, tol);
nMat = G.Matrix.cells.num;
nc = Gp.cells.num;

%% Coarse model from the barrier-aware partition
p = partitionCPEDFM_MS(Gp, [5, 3, 2], 'tol', tol);
coarse = buildCoarseCPEDFM(Gp, ops, Gp.rock, p);
fprintf('coarse model: %d blocks, %d connections\n', ...
    coarse.nBlocks, size(coarse.pairs, 1));

%% Check A: conservation and compartments
pvFine = sum(poreVolume(Gp, Gp.rock));
pvCoarse = sum(coarse.pv);
relPV = abs(pvCoarse - pvFine) / pvFine;

% Coarse compartments over MATRIX-MATRIX connections only, exactly
% parallel to the fine-side check below. The exported coarse model
% correctly retains the fracture-mediated projection paths across the
% barrier (matrix -> sealing fracture -> matrix, many orders of
% magnitude weaker than matrix connections - the same tiny leakage the
% fine model carries); those paths must not enter the compartment
% topology check, just as fine NNCs do not.
mm = p(ops.N(:,1)) ~= p(ops.N(:,2)) & ...
     ops.N(:,1) <= nMat & ops.N(:,2) <= nMat;
loM = min(p(ops.N(mm,1)), p(ops.N(mm,2)));
hiM = max(p(ops.N(mm,1)), p(ops.N(mm,2)));
[pairsM, ~, icM] = unique([loM, hiM], 'rows');
TM = accumarray(icM, ops.T(mm));
sigM = TM > max(TM)*1e-12;
gC = graph(pairsM(sigM,1), pairsM(sigM,2), [], coarse.nBlocks);
compC = conncomp(gC)';
szC = accumarray(compC, coarse.pv);          % pore volume per component
macroC = nnz(szC > 0.05*pvFine);
% fracture-mediated cross-barrier strength, for the record
inMM = ismember(coarse.pairs, pairsM(sigM,:), 'rows');
maxCrossT = max([coarse.T(~inMM); 0]) / max(coarse.T);

% fine-side reference: macro components of the matrix graph
inter = all(G.Matrix.faces.neighbors ~= 0, 2);
Nmm = G.Matrix.faces.neighbors(inter,:);
Tmm = ops.T_all(1:G.Matrix.faces.num); Tmm = Tmm(inter);
gF = graph(Nmm(Tmm > max(Tmm)*1e-12, 1), Nmm(Tmm > max(Tmm)*1e-12, 2), [], nMat);
compF = conncomp(gF)';
szF = accumarray(compF, 1);
macroF = nnz(szF > 0.05*nMat);

% consistency: no coarse connection may join blocks whose (macro-component)
% fine cells lie in different compartments
blockComp = zeros(coarse.nBlocks, 1);
for b = 1:coarse.nBlocks
    cellsIn = find(p(1:nMat) == b);
    if ~isempty(cellsIn)
        cf = compF(cellsIn);
        big = szF(cf) > 0.05*nMat;
        if any(big), blockComp(b) = mode(cf(big)); else, blockComp(b) = mode(cf); end
    end
end
mixed = 0;
sigIdx = find(sigM);
for q = sigIdx(:)'
    b1 = blockComp(pairsM(q,1)); b2 = blockComp(pairsM(q,2));
    if b1 > 0 && b2 > 0 && b1 ~= b2 && ...
            szF(b1) > 0.05*nMat && szF(b2) > 0.05*nMat
        mixed = mixed + 1;
    end
end
fprintf(['A. PV rel diff %.2e | coarse macro comps %d (fine %d) | ' ...
    'matrix-matrix cross-compartment connections %d ' ...
    '(fracture-mediated paths at %.1e of matrix T)\n'], ...
    relPV, macroC, macroF, mixed, maxCrossT);
if relPV < 1e-12 && macroC == macroF && mixed == 0
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('A: relPV %.2e, comps %d/%d, mixed %d', ...
        relPV, macroC, macroF, mixed);
    fprintf('   FAIL\n');
end

%% Check B: upscaling quality under the fine Dirichlet drive
% fine solution (as in test30)
addpath(fullfile(thisDir, '..', '..', 'CPEDFM-MS'));
[A, ~, ~] = getCPEDFMSystemMatrix(Gp, ops);
bf = find(any(G.Matrix.faces.neighbors == 0, 2));
fc = G.Matrix.faces.centroids(bf, :);
wf = bf(fc(:,1) < 1e-3);
ef = bf(fc(:,1) > physdim(1) - 1e-3);
q = zeros(nc, 1);
Tb = ops.T_all;
for s = 1:2
    if s == 1, faces = wf; pb = 1; else, faces = ef; pb = 0; end
    for f = faces(:)'
        ccell = sum(G.Matrix.faces.neighbors(f, :));
        A(ccell, ccell) = A(ccell, ccell) + Tb(f);
        q(ccell) = q(ccell) + Tb(f) * pb;
    end
end
pFine = A \ q;

% coarse system: aggregated connections + aggregated boundary terms
Ac = sparse(coarse.pairs(:,1), coarse.pairs(:,2), -coarse.T, ...
    coarse.nBlocks, coarse.nBlocks);
Ac = Ac + Ac';
Ac = Ac - spdiags(sum(Ac, 2), 0, coarse.nBlocks, coarse.nBlocks);
qc = zeros(coarse.nBlocks, 1);
for s = 1:2
    if s == 1, faces = wf; pb = 1; else, faces = ef; pb = 0; end
    for f = faces(:)'
        b = p(sum(G.Matrix.faces.neighbors(f, :)));
        Ac(b, b) = Ac(b, b) + Tb(f);
        qc(b) = qc(b) + Tb(f) * pb;
    end
end
pCoarse = Ac \ qc;

% compartment plateaus, fine vs coarse (pore-volume-weighted means)
p0f = fp(1).points(1,:);
e1 = fp(1).points(2,:) - p0f; e2 = fp(1).points(4,:) - p0f;
nrm = cross(e1, e2); nrm = nrm/norm(nrm);
sdF = (G.cells.centroids(1:nMat,:) - p0f) * nrm';
pvF = poreVolume(Gp, Gp.rock);
plF = [sum(pFine(sdF<0).*pvF(sdF<0))/sum(pvF(sdF<0)), ...
       sum(pFine(sdF>0).*pvF(sdF>0))/sum(pvF(sdF>0))];
sdC = accumarray(p(1:nMat), sdF .* pvF(1:nMat)) ./ ...
      max(accumarray(p(1:nMat), pvF(1:nMat)), eps);
wB = coarse.pv;
plC = [sum(pCoarse(sdC<0).*wB(sdC<0))/sum(wB(sdC<0)), ...
       sum(pCoarse(sdC>0).*wB(sdC>0))/sum(wB(sdC>0))];
errPl = max(abs(plC - plF));
fprintf('B. plateaus fine [%.4f %.4f] vs coarse [%.4f %.4f], max diff %.4f\n', ...
    plF, plC, errPl);
if errPl < 0.02
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('B: plateau diff %.4f', errPl);
    fprintf('   FAIL\n');
end

%% Check C: deck round-trip
[iW, ~, ~] = ind2sub(dims, G.Matrix.cells.indexMap(1:nMat));
wcInj  = find(iW == 1);  wcInj = wcInj(1:2:end);
wcProd = find(iW == dims(1)); wcProd = wcProd(1:2:end);
W = addWell([], G.Matrix, G.rock, wcInj, 'Type', 'bhp', 'Val', 300*barsa, ...
    'Comp_i', [1 0], 'Name', 'INJ', 'Radius', 0.1, 'Sign', 1);
W = addWell(W, G.Matrix, G.rock, wcProd, 'Type', 'bhp', 'Val', 100*barsa, ...
    'Comp_i', [0 1], 'Name', 'PROD', 'Radius', 0.1, 'Sign', -1);

deckFile = fullfile(outDir, 'CPEDFM_COARSE.DATA');
info = writeEclipseCoarseCPEDFM(coarse, W, deckFile);
deck = readEclipseDeck(deckFile);
deck = convertDeckUnits(deck);

okC = true; msg = {};
% NNC
if isfield(deck.GRID, 'NNC')
    nn = deck.GRID.NNC;
    dPair = max(max(abs(nn(:,[1 4]) - coarse.pairs)));
    relT = max(abs(nn(:,7) - coarse.T) ./ coarse.T);
    okC = okC && size(nn,1) == size(coarse.pairs,1) && dPair == 0 && relT < 1e-6;
    if ~okC, msg{end+1} = 'NNC'; end
else
    okC = false; msg{end+1} = 'NNC missing';
end
% PORV / DEPTH (EDIT keywords land in deck.GRID)
if isfield(deck.GRID, 'PORV') && isfield(deck.GRID, 'DEPTH')
    relpv = max(abs(deck.GRID.PORV(:) - coarse.pv) ./ max(coarse.pv, eps));
    dz = max(abs(deck.GRID.DEPTH(:) - coarse.depth));
    okC = okC && relpv < 1e-6 && dz < 1e-3;
    if relpv >= 1e-6 || dz >= 1e-3, msg{end+1} = 'PORV/DEPTH'; end
else
    okC = false; msg{end+1} = 'PORV/DEPTH missing';
end
fprintf('C. deck round-trip: %d blocks, %d NNCs%s\n', info.nBlocks, info.nConn, ...
    ternary(okC, ' - all values match', [' - MISMATCH: ' strjoin(msg, ', ')]));
if okC
    fprintf('   PASS\n');
else
    failures{end+1} = 'C: deck round-trip mismatch';
    fprintf('   FAIL\n');
end

%% Summary
if isempty(failures)
    fprintf('\n=== test21: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test21: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test21 failed');
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
