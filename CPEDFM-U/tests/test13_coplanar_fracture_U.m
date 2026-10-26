%% test13_coplanar_fracture_U
% Regression test for the face-coplanar ('fracmat boundary') path, which
% no existing CPEDFM-U test exercised - test10/11/12 all use fractures
% slanted relative to the grid, so the bugs this covers went unnoticed.
%
% A fracture placed exactly in the plane of an interior grid face used to
% fail three ways:
%
%   A. fracPolyCellIntersect abandoned its clip loop on detecting the
%      coplanar face, returning the WHOLE fracture polygon as the
%      intersection with every cell along the plane (a factor-of-12 area
%      overestimate here). The resulting duplicate polygons then made the
%      fracture grid degenerate and EDFMGridU threw outright.
%   B. EDFMGridU registered such a fracture cell against only one of the
%      two matrix cells bordering it, leaving the far side disconnected,
%      after which the row-counting surface-area rule doubled the
%      surviving transmissibility as if the fracture were interior.
%   C. the fracture-matrix CI was built from a bounding-box <d>.
%
% Checks
%   A. Clipped areas: each cell gets its own face footprint, and the total
%      over both sides is exactly 2 x the fracture polygon area.
%   B. Registration: 2 connections per fracture cell, one host on each
%      side of the fracture plane; no matrix cell along the plane is left
%      unconnected.
%   C. Surface-area rule: a boundary fracture cell must NOT get the
%      interior x2 factor, so CI = A/<d> exactly (not 2A/<d>).
%   D. The full CPEDFM-U pipeline runs on a coplanar fracture and seals
%      it: a sealing coplanar fracture must split the matrix into 2
%      connected components, where plain EDFM leaves 1.
%   E. Mixed case: one coplanar and one slanted fracture in the same
%      model, which the Cartesian implementation cannot represent at all
%      (fracmatnnc asserts Nm == Nf || Nm == 2*Nf).

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupCPEDFMU.m'));

tol = 1e-7;
failures = {};
mrstVerbose off

celldim = [4 4 3];
physdim = [40 40 30];
poly = [20 2 2; 20 38 2; 20 38 28; 20 2 28];   % exactly on the x = 20 faces
Apoly = polyArea3D_U(poly);

fpDef = struct('points', poly, 'aperture', 10*milli*meter, ...
               'poro', 0.4, 'perm', 1*nano*darcy);

G0 = computeGeometry(cartGrid(celldim, physdim));
G0.rock = makeRock(G0, 1*milli*darcy, 0.1);

%% A. clipped areas
onPlane = find(abs(G0.cells.centroids(:,1) - 15) < 1e-9 | ...
               abs(G0.cells.centroids(:,1) - 25) < 1e-9);
Asum = 0; nOnFace = 0; maxA = 0;
for c = onPlane'
    [a, ~, onFace] = fracPolyCellIntersect(poly, G0, c, tol);
    Asum = Asum + a;
    nOnFace = nOnFace + (onFace > 0);
    maxA = max(maxA, a);
end
relA = abs(Asum - 2*Apoly) / (2*Apoly);
fprintf('A. clipped area over %d cells: %.1f m^2 vs 2 x %.1f = %.1f (rel %.2e)\n', ...
    numel(onPlane), Asum, Apoly, 2*Apoly, relA);
fprintf('   coplanar face detected in %d/%d cells; largest single clip %.1f m^2\n', ...
    nOnFace, numel(onPlane), maxA);
if relA < 1e-12 && nOnFace == numel(onPlane) && maxA < Apoly
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('A: area sum rel %.2e, maxA %.1f', relA, maxA);
    fprintf('   FAIL\n');
end

%% B. registration: both neighbours, one per side
G = EDFMGridU(G0, fpDef, 'Tolerance', tol, 'Verbose', false);
mc  = G.FracGrid.Frac1.matrix_connection;
nfc = G.FracGrid.Frac1.cells.num;
xh  = G.Matrix.cells.centroids(mc.cells, 1);
share = accumarray(mc.fraccells, 1, [nfc 1]);

fprintf('B. %d fracture cells, %d connections; hosts west/east of x=20: %d/%d\n', ...
    nfc, numel(mc.cells), sum(xh < 20), sum(xh > 20));
okB = all(mc.type == 1) && all(share == 2) && ...
      sum(xh < 20) == nfc && sum(xh > 20) == nfc && ...
      isequal(sort(unique(mc.cells)), sort(onPlane));
if okB
    fprintf('   PASS (every fracture cell shared by exactly 2 hosts)\n');
else
    failures{end+1} = 'B: coplanar registration wrong';
    fprintf('   FAIL\n');
end

%% C. surface-area rule: boundary must not get the interior x2
G = fracturematrixNNC_U(G, tol, 'Verbose', false);
[~, CI, area, type, davg] = fracMatNNC_U(G.FracGrid.Frac1, G.Matrix, tol);
isBnd = strcmp(type, 'fracmat boundary');
ratio = CI(isBnd) .* davg(isBnd) ./ area(isBnd);   % 1 if correct, 2 if doubled
fprintf('C. boundary CI*<d>/A = %.6f (1 = correct, 2 = interior rule misapplied)\n', ...
    max(ratio));
if all(isBnd) && max(abs(ratio - 1)) < 1e-12
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('C: CI factor %.4f', max(ratio));
    fprintf('   FAIL\n');
end

%% D. the pipeline runs and a sealing coplanar fracture compartmentalizes
% Uses a fracture that spans the FULL cross-section (the one above stops
% short of the domain edges, so fluid would simply flow around it and no
% method could split the grid). Each x = 20 face is then fully covered, so
% the seal must be complete.
fpFull = fpDef;
fpFull.points = [20 -5 -5; 20 45 -5; 20 45 35; 20 -5 35];

Gd = EDFMGridU(G0, fpFull, 'Tolerance', tol, 'Verbose', false);
Gd = fracturematrixNNC_U(Gd, tol, 'Verbose', false);
Gp = pMatFracNNCsCPEDFM_U(Gd, tol);

opsU = setupCPEDFMOpsTPFA_CPG(Gp, Gp.rock, tol);
opsE = setupShaleEDFMOpsTPFA(Gd, Gd.rock, tol);

nMat  = Gd.Matrix.cells.num;
inter = all(Gd.Matrix.faces.neighbors ~= 0, 2);
Nmm   = Gd.Matrix.faces.neighbors(inter,:);
nComp = @(ops) countComponents(ops, Gd, inter, Nmm, nMat);
nU = nComp(opsU); nE = nComp(opsE);

fprintf('D. matrix components: CPEDFM-U %d | EDFM %d (%d sealed faces)\n', ...
    nU, nE, numel(unique(Gp.nnc.pMMfaces)));
% every G.nnc column must still line up after the projection stage
nRows = size(Gp.nnc.cells,1);
aligned = numel(Gp.nnc.T) == nRows && numel(Gp.nnc.type) == nRows && ...
          numel(Gp.nnc.area) == nRows && size(Gp.nnc.normal,1) == nRows && ...
          numel(Gp.nnc.CI) == nRows && numel(Gp.nnc.davg) == nRows;
fprintf('   G.nnc columns aligned after projection: %d\n', aligned);
if nU >= 2 && nE == 1 && aligned
    fprintf('   PASS (EDFM leaves the sealing fracture bypassed; CPEDFM-U closes it)\n');
else
    failures{end+1} = sprintf('D: components U=%d E=%d aligned=%d', nU, nE, aligned);
    fprintf('   FAIL\n');
end

%% E. mixed coplanar + slanted fractures in one model
fpMix(1) = fpDef;
fpMix(2) = struct('points', [4 4 5; 36 30 5; 36 30 25; 4 4 25], ...
                  'aperture', 10*milli*meter, 'poro', 0.4, 'perm', 1*darcy);
Gm = EDFMGridU(G0, fpMix, 'Tolerance', tol, 'Verbose', false);
Gm = fracturematrixNNC_U(Gm, tol, 'Verbose', false);

t1 = Gm.FracGrid.Frac1.matrix_connection.type;
t2 = Gm.FracGrid.Frac2.matrix_connection.type;
fprintf('E. mixed model: Frac1 types %s (%d rows), Frac2 types %s (%d rows)\n', ...
    mat2str(unique(t1)'), numel(t1), mat2str(unique(t2)'), numel(t2));
okE = all(t1 == 1) && all(t2 == 2) && ...
      numel(t1) == 2*Gm.FracGrid.Frac1.cells.num && ...
      numel(t2) ==   Gm.FracGrid.Frac2.cells.num;
try
    setupShaleEDFMOpsTPFA(Gm, Gm.rock, tol);
catch me
    okE = false;
    fprintf('   operator setup failed: %s\n', me.message);
end
if okE
    fprintf('   PASS (boundary and interior fractures coexist)\n');
else
    failures{end+1} = 'E: mixed coplanar/slanted model wrong';
    fprintf('   FAIL\n');
end

%% Summary
if isempty(failures)
    fprintf('\n=== test13: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test13: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test13 failed');
end

%% ------------------------------ helpers ---------------------------------
function n = countComponents(ops, G, inter, Nmm, nMat)
% Connected components of the matrix-matrix graph, keeping only faces that
% still carry transmissibility after the multipliers are applied.
T = ops.T_all(1:G.Matrix.faces.num);
T = T(inter);
open = T > max(T)*1e-12;
n = max(conncomp(graph(Nmm(open,1), Nmm(open,2), [], nMat)));
end
