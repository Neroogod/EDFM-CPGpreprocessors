%% test10_regression_cartesian_U
% Chain the unstructured implementation back to the validated stack on a
% Cartesian grid (a cartGrid is a valid unstructured grid whose cells
% are convex, so both pipelines apply).
%
% Case: the Fig. 11a slanted sealing fracture (same as CPEDFM-CPG
% test00, which itself matches the original Cartesian CPEDFM release).
%
% Checks:
%   A. GRIDDING: EDFMGridU host-cell set and intersection areas match
%      EDFMshalegrid (the existing, validated gridding) exactly.
%   B. SEAL: the unstructured flip-cut compartmentalizes the matrix into
%      2 components, as CPEDFM-CPG does; plain EDFM stays at 1.
%   C. CUT MINIMALITY: every sealed face separates cells on opposite
%      sides of the fracture plane, and the sealed set contains no face
%      outside the fracture's lateral extent (tip guard on a
%      full-domain fracture must not discard interior faces).

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupCPEDFMU.m'));

tol = 1e-5;
failures = {};

%% Grid and fracture (Fig. 11a case)
celldim = [30, 30, 5];
physdim = [60, 60, 10];
G0 = cartGrid(celldim, physdim);
G0 = computeGeometry(G0);
G0.rock = makeRock(G0, 1*milli*darcy, 0.1);

fpDef = struct('points', [10 0 0; 50 60 0; 50 60 10; 10 0 10], ...
    'aperture', 10*milli*meter, 'poro', 0.4, 'perm', 1*nano*darcy);

%% Check A: gridding regression (EDFMGridU vs EDFMshalegrid)
fpA = fpDef;
[GA, ~] = EDFMshalegrid(G0, fpA, 'Tolerance', tol, 'plotgrid', false, ...
    'fracturelist', 1);
fpB = fpDef;
[GB, ~] = EDFMGridU(G0, fpB, 'Tolerance', tol, 'fracturelist', 1);

mcA = GA.FracGrid.Frac1.matrix_connection.cells;
mcB = GB.FracGrid.Frac1.matrix_connection.cells;
arA = GA.FracGrid.Frac1.matrix_connection.area;
arB = GB.FracGrid.Frac1.matrix_connection.area;

sameCells = isequal(sort(mcA), sort(mcB));
if sameCells
    [~, ia] = sort(mcA); [~, ib] = sort(mcB);
    dA = max(abs(arA(ia) - arB(ib)));
    rel = dA / max(arA);
else
    rel = inf;
end
fprintf('A. gridding: cells %d vs %d, same set: %d, max area rel diff: %.3e\n', ...
    numel(mcA), numel(mcB), sameCells, rel);
if sameCells && rel < 1e-6
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('A: gridding mismatch (rel %.3e)', rel);
    fprintf('   FAIL\n');
end

%% Check B: seal on the unstructured path
GB = fracturematrixNNC_U(GB, tol);
GBp = pMatFracNNCsCPEDFM_U(GB, tol);
opsU = setupCPEDFMOpsTPFA_CPG(GBp, GBp.rock, tol);
opsE = setupShaleEDFMOpsTPFA(GB, GB.rock, tol);

nMat = GB.Matrix.cells.num;
inter = all(GB.Matrix.faces.neighbors ~= 0, 2);
Nmm = GB.Matrix.faces.neighbors(inter,:);
TmU = opsU.T_all(1:GB.Matrix.faces.num); TmU = TmU(inter);
TmE = opsE.T_all(1:GB.Matrix.faces.num); TmE = TmE(inter);
gU = graph(Nmm(TmU > max(TmU)*1e-12, 1), Nmm(TmU > max(TmU)*1e-12, 2), [], nMat);
gE = graph(Nmm(TmE > max(TmE)*1e-12, 1), Nmm(TmE > max(TmE)*1e-12, 2), [], nMat);
nU = max(conncomp(gU)); nE = max(conncomp(gE));
fprintf('B. components: CPEDFM-U %d | EDFM %d (%d sealed faces)\n', ...
    nU, nE, numel(unique(GBp.nnc.pMMfaces)));
if nU >= 2 && nE == 1
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('B: components U=%d EDFM=%d', nU, nE);
    fprintf('   FAIL\n');
end

%% Check C: every sealed face is a genuine plane-crossing face
p0 = fpDef.points(1,:);
e1 = fpDef.points(2,:) - p0; e2 = fpDef.points(4,:) - p0;
nrm = cross(e1, e2); nrm = nrm/norm(nrm);
sd = (GB.cells.centroids(1:nMat,:) - p0) * nrm';
okC = true;
for q = 1:numel(GBp.nnc.pMMfaces)
    nb = GB.Matrix.faces.neighbors(GBp.nnc.pMMfaces(q),:);
    if any(nb == 0) || sd(nb(1))*sd(nb(2)) > 0
        okC = false; break
    end
end
fprintf('C. sealed faces all straddle the fracture plane: %d\n', okC);
if okC
    fprintf('   PASS\n');
else
    failures{end+1} = 'C: non-crossing face sealed';
    fprintf('   FAIL\n');
end

%% Summary
if isempty(failures)
    fprintf('\n=== test10: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test10: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test10 failed');
end
