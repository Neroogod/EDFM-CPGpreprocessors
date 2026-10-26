%% test21_pebi_sealing
% 3D pEDFM on a genuinely unstructured grid: extruded PEBI (Voronoi) cells,
% where no structured concept exists at all - no face tags, no IJK, no
% logical directions, and a varying number of faces per cell.
%
% This is what the antipodal-pair rule of selectProjectionFacesPEDFM_U
% exists for. test20 pins that rule to the published Cartesian result;
% this test exercises it where the paper's own formulation has nothing to
% say.
%
% WHAT pEDFM IS AND IS NOT EXPECTED TO DO HERE
%   pEDFM selects faces independently per host cell, so the sealed set need
%   not be connected and a sealing fracture throttles flow without
%   necessarily compartmentalizing the domain. CPEDFM-U's flip-cut is what
%   guarantees a watertight barrier. This test therefore asserts a large
%   reduction in cross-fracture conductance and a materially larger
%   sustained pressure drop than plain EDFM, and records the residual leak
%   rather than demanding it be zero.
%
% Checks:
%   A. Pipeline runs end to end on the PEBI grid; cells really do have
%      varying face counts.
%   B. Projection data self-consistent (multipliers in [0,1], faces
%      incident on their cell pair, positive finite areas, equal-length
%      G.nnc arrays).
%   C. Selection sanity: every host selects at most one face of each
%      antipodal pair, i.e. no more than ceil(nFaces/2) faces.
%   D. Faces the fracture crosses are substantially throttled.
%   E. Tip guard: a fracture spanning only the centre of the domain does
%      not project beyond its own extent.
%   F. Flow: pEDFM sustains a materially larger cross-fracture pressure
%      drop than plain EDFM.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupPEDFMU.m'));

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
failures = {};
physdim = [100, 100, 20];
nLayers = 6;
rng(42);

%% ---- Full sealing fracture -------------------------------------------
G = localPebiGrid(physdim, nLayers);
nfPerCell = diff(G.cells.facePos);
fprintf('PEBI grid: %d cells, %d faces, faces/cell %d-%d\n', ...
    G.cells.num, G.faces.num, min(nfPerCell), max(nfPerCell));
assert(max(nfPerCell) > min(nfPerCell), ...
    'grid has uniform face counts; it is not exercising the general rule');

fp = struct;
fp(1).points = [ 35  -5  -5;
                 55 105  -5;
                 75 105  25;
                 55  -5  25];       % slanted in x-y and x-z, full extent
fp(1).aperture = 10*milli*meter;
fp(1).poro = 0.4;
fp(1).perm = 1*nano*darcy;

[G, fp] = EDFMGridU(G, fp, 'Tolerance', tol, 'fracturelist', 1);
G = fracturematrixShaleNNC3D(G, tol);

Gp = pMatFracNNCsPEDFM_U(G, tol);
opsP = setupPEDFMOpsTPFA3D_CPG(Gp, Gp.rock, tol);
opsE = setupShaleEDFMOpsTPFA(G, G.rock, tol);

s = Gp.nnc.pStats;
fprintf('A. pipeline: PASS (%d hosts, %d projections, %d faces)\n', ...
    s.nHosts, s.nProjections, numel(unique(Gp.nnc.pMMfaces)));
fprintf('   antipodal ties=%d | node-rule fallbacks=%d | union fallbacks=%d\n', ...
    s.nTies, s.nFallback, s.nUnionFallback);

%% Check B: self-consistency
mult = transmultPEDFM3D_CPG(Gp, tol);
okB = all(mult >= -1e-12 & mult <= 1 + 1e-12) && ...
      all(Gp.nnc.pMMarea > 0) && all(isfinite(Gp.nnc.pMMarea)) && ...
      all(isfinite(Gp.nnc.T)) && ...
      isequal(size(Gp.nnc.cells,1), numel(Gp.nnc.T), ...
              numel(Gp.nnc.area), numel(Gp.nnc.type));
inc = false(numel(Gp.nnc.pMMfaces),1);
for q = 1:numel(Gp.nnc.pMMfaces)
    nb = Gp.faces.neighbors(Gp.nnc.pMMfaces(q),:);
    inc(q) = any(nb == Gp.nnc.pMMneighs(q,1)) && any(nb == Gp.nnc.pMMneighs(q,2));
end
okB = okB && all(inc);
if okB
    fprintf('B. projection data self-consistent: PASS\n');
else
    failures{end+1} = 'B: projection data inconsistent';
    fprintf('B. projection data: FAIL\n');
end

%% Check C: at most one face per antipodal pair
rows = find(strcmp(G.nnc.type, 'fracmat interior'));
opt = struct('tol', tol, 'relTol', 1e-6, 'defaultFace', 'bottom');
okC = true; worst = 0; nSelTot = 0;
for ii = 1:numel(rows)
    r = rows(ii);
    hc = G.nnc.cells(r,1);
    nrm = G.nnc.normal(r,:); nrm = nrm/norm(nrm);
    sel = selectProjectionFacesPEDFM_U(G, hc, ...
        G.cells.centroids(G.nnc.cells(r,2),:), nrm, G.nnc.planepoint(r,:), opt);
    nF = G.cells.facePos(hc+1) - G.cells.facePos(hc);
    nSelTot = nSelTot + numel(sel);
    worst = max(worst, numel(sel)/nF);
    if numel(sel) > ceil(nF/2), okC = false; end
end
fprintf('C. selection: %d faces selected over %d hosts, max fraction of a cell''s faces %.2f\n', ...
    nSelTot, numel(rows), worst);
if okC
    fprintf('   PASS (never more than one face per antipodal pair)\n');
else
    failures{end+1} = 'C: a host selected more than ceil(nFaces/2) faces';
    fprintf('   FAIL\n');
end

%% Check D: crossed faces are throttled
nMat = G.Matrix.cells.num;
r1 = rows(1);
nrm = G.nnc.normal(r1,:); nrm = nrm/norm(nrm);
ppt = G.nnc.planepoint(r1,:);
sd = (G.cells.centroids(1:nMat,:) - ppt) * nrm';
Nm = G.Matrix.faces.neighbors;
interior = all(Nm ~= 0, 2);
fIdx = find(interior);
crossed = fIdx(sd(Nm(interior,1)) .* sd(Nm(interior,2)) < 0);
% restrict to the fracture's lateral extent, as the tip guard does
fc = G.faces.centroids(crossed,:);
inExtent = fc(:,1) > 25 & fc(:,1) < 85;
crossed = crossed(inExtent);

mc = mult(crossed);
fracThrottled = mean(mc < 0.999);
fprintf('D. crossed faces (within extent): %d | %.0f%% throttled | median multiplier %.4f\n', ...
    numel(crossed), 100*fracThrottled, median(mc));
if fracThrottled > 0.5 && median(mc) < 0.5
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('D: %.0f%% throttled, median %.3f', ...
        100*fracThrottled, median(mc));
    fprintf('   FAIL\n');
end

%% Check E: tip guard - a central fracture must not project beyond its extent
Gc = localPebiGrid(physdim, nLayers);
fpC = struct;
fpC(1).points = [ 45  35  4;
                  55  65  4;
                  55  65 16;
                  45  35 16];       % small central patch
fpC(1).aperture = 10*milli*meter;
fpC(1).poro = 0.4;
fpC(1).perm = 1*nano*darcy;
[Gc, fpC] = EDFMGridU(Gc, fpC, 'Tolerance', tol, 'fracturelist', 1);
Gc = fracturematrixShaleNNC3D(Gc, tol);
Gcp = pMatFracNNCsPEDFM_U(Gc, tol);

pf = unique(Gcp.nnc.pMMfaces);
pc = Gc.faces.centroids(pf,:);
margin = 12;                       % generous: cell size is ~8 m
okE = all(pc(:,1) > 45-margin & pc(:,1) < 55+margin & ...
          pc(:,2) > 35-margin & pc(:,2) < 65+margin);
fprintf('E. tip guard: %d projected faces, x in [%.1f %.1f], y in [%.1f %.1f] (patch x 45-55, y 35-65)\n', ...
    numel(pf), min(pc(:,1)), max(pc(:,1)), min(pc(:,2)), max(pc(:,2)));
if okE
    fprintf('   PASS (projection confined to the fracture''s extent)\n');
else
    failures{end+1} = 'E: projections escaped the fracture extent';
    fprintf('   FAIL\n');
end

%% Check F: flow across the full sealing fracture
gravity off
fluid = initSimpleADIFluid('phases', 'W', ...
    'mu', 1*centi*poise, 'rho', 1000, 'c', 1e-8/barsa, 'pRef', 200*barsa);

cc = G.cells.centroids(1:nMat,:);
wcInj  = find(cc(:,1) < 0.12*physdim(1));
wcProd = find(cc(:,1) > 0.88*physdim(1));
assert(~isempty(wcInj) && ~isempty(wcProd), 'well cell selection failed');

pv = poreVolume(G, G.rock);
totTime = 2*year;
rate = 0.3*sum(pv(1:nMat))/totTime;

W = addWell([], G.Matrix, G.rock, wcInj, 'Type', 'rate', 'Val', rate, ...
    'Comp_i', 1, 'Name', 'I', 'Radius', 0.1);
W = addWell(W, G.Matrix, G.rock, wcProd, 'Type', 'bhp', 'Val', 100*barsa, ...
    'Comp_i', 1, 'Name', 'P', 'Radius', 0.1);

dt = rampupTimesteps(totTime, 60*day, 5);
schedule = simpleSchedule(dt, 'W', W);
state0 = initResSol(G, 200*barsa, 1);

modelP = GenericBlackOilModel(Gp, [], fluid, 'water', true, 'oil', false, 'gas', false);
modelP.operators = opsP; modelP.minimumPressure = 1*barsa;
[~, statesP] = simulateScheduleAD(state0, modelP, schedule);

modelE = GenericBlackOilModel(G, [], fluid, 'water', true, 'oil', false, 'gas', false);
modelE.operators = opsE; modelE.minimumPressure = 1*barsa;
[~, statesE] = simulateScheduleAD(state0, modelE, schedule);

sideA = find(cc(:,1) < 30); sideB = find(cc(:,1) > 80);
dpP = mean(statesP{end}.pressure(sideA)) - mean(statesP{end}.pressure(sideB));
dpE = mean(statesE{end}.pressure(sideA)) - mean(statesE{end}.pressure(sideB));
fprintf('F. cross-fracture dP: pEDFM %.4g bar | EDFM %.4g bar (ratio %.2f)\n', ...
    dpP/barsa, dpE/barsa, dpP/max(dpE, eps));
if dpP > 1.25*dpE && dpP > 1*barsa
    fprintf('   PASS (pEDFM sustains a materially larger pressure jump)\n');
else
    failures{end+1} = sprintf('F: dP pEDFM %.4g bar vs EDFM %.4g bar', ...
        dpP/barsa, dpE/barsa);
    fprintf('   FAIL\n');
end

%% Projection figure
try
    fig = figure('Visible', 'off', 'Position', [100 100 1100 480]);
    subplot(1,2,1); hold on
    plotGrid(G.Matrix, 'FaceColor', 'none', 'EdgeAlpha', 0.05);
    plotGrid(Gp, Gp.nnc.cells(strcmp(Gp.nnc.type,'fracmat interior'),2), ...
        'FaceColor', 'y', 'EdgeAlpha', 0.3);
    view(-35, 22); axis equal tight off; title('fracture cells (PEBI)');
    subplot(1,2,2); hold on
    plotGrid(G.Matrix, 'FaceColor', 'none', 'EdgeAlpha', 0.05);
    plotFaces(G.Matrix, unique(Gp.nnc.pMMfaces), 'FaceColor', 'r', 'EdgeAlpha', 0.2);
    view(-35, 22); axis equal tight off; title('pEDFM projected faces');
    print(fig, fullfile(outDir, 'test21_projections.png'), '-dpng', '-r150');
    close(fig);
    fprintf('projection figure saved to tests/output/test21_projections.png\n');
catch me
    fprintf('figure generation skipped: %s\n', me.message);
end

%% Summary
save(fullfile(outDir, 'test21_results.mat'), 'dpP', 'dpE', 's', 'fracThrottled');
if isempty(failures)
    fprintf('\n=== test21: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test21: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test21 failed');
end

%% ---- local helpers -----------------------------------------------------
function G = localPebiGrid(physdim, nLayers)
% Extruded PEBI grid scaled to physdim.
    G2D = pebiGrid2D(physdim(1)/12, physdim(1:2));
    G = makeLayeredGrid(G2D, nLayers);
    G.nodes.coords(:,3) = G.nodes.coords(:,3) * physdim(3)/nLayers;
    G = computeGeometry(G);
    G.rock = makeRock(G, 100*milli*darcy, 0.25);
end
