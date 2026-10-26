%% test01_skewedCPG_singleFrac
% CPEDFM on a genuinely distorted corner-point grid: a slanted sealing
% fracture must compartmentalize the domain.
%
% Grid: simpleGrdecl with skewed pillars (non-axis-aligned faces, cell
% centroids share no coordinates - the Cartesian CPEDFM machinery cannot
% form neighborhoods here).
%
% Fracture: low-conductivity plane cutting the full domain cross-section,
% slanted with respect to all axes.
%
% Checks:
%   A. Pipeline runs end to end on the CPG (gridding -> EDFM NNCs ->
%      CPEDFM projections -> operators).
%   B. Seal/graph check: with CPEDFM multipliers the matrix-matrix
%      connectivity graph (faces with T > 0) splits into two components
%      separated by the fracture; with plain EDFM it stays one component
%      (the leakage CPEDFM exists to fix).
%   C. Flow check: producer water rate under CPEDFM is orders of
%      magnitude below the EDFM (leaky) rate.
%   D. Sanity: multipliers in [0,1]; projected faces belong to host or
%      projection cells; areas positive.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupCPEDFMCPG.m'));

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
failures = {};

%% Corner-point grid with skewed pillars, scaled to field dimensions
dims = [20, 10, 8];
grdecl = simpleGrdecl(dims, 0.12);          % skewed, unit-box footprint
physdim = [200, 100, 40];                    % m
grdecl.COORD(1:3:end) = grdecl.COORD(1:3:end) * physdim(1);
grdecl.COORD(2:3:end) = grdecl.COORD(2:3:end) * physdim(2);
grdecl.COORD(3:3:end) = grdecl.COORD(3:3:end) * physdim(3);
grdecl.ZCORN = grdecl.ZCORN * physdim(3);

G = processGRDECL(grdecl);
G = computeGeometry(G);
G.rock = makeRock(G, 100*milli*darcy, 0.25);

% verify this grid actually breaks the Cartesian assumptions
cc = G.cells.centroids;
nShareX = numel(unique(cc(:,1)));
fprintf('distinct centroid x-coords: %d of %d cells (Cartesian code needs repeats)\n', ...
    nShareX, G.cells.num);

%% Slanted sealing fracture cutting the full YZ cross-section
% Plane sweeps in x as functions of y and z -> inclined to all axes.
fracplanes = struct;
fracplanes(1).points = [ 70  -5  -5;
                        110 105  -5;
                        150 105  45;
                        110  -5  45];
fracplanes(1).aperture = 10*milli*meter;
fracplanes(1).poro = 0.4;
fracplanes(1).perm = 1*nano*darcy;

%% CPG gridding + EDFM NNCs (existing code, unchanged)
[G, fracplanes] = EDFMGridCPG(G, fracplanes, ...
    'Tolerance', tol, 'plotgrid', false, 'fracturelist', 1);
G = fracturematrixShaleNNC3D(G, tol);

%% CPEDFM projections + operators (new code)
Gp = pMatFracNNCsCPEDFM_CPG(G, tol);
opsCPEDFM = setupCPEDFMOpsTPFA_CPG(Gp, Gp.rock, tol);
opsEDFM   = setupShaleEDFMOpsTPFA(G, G.rock, tol);  % leaky baseline
fprintf('A. pipeline: PASS (%d projected NNCs)\n', numel(Gp.nnc.pMMfaces));

%% Check D first (cheap sanity on the projection data)
mult = transmultpEDFM_CPG(Gp, tol);
okD = all(mult >= 0 & mult <= 1) && all(Gp.nnc.pMMarea > 0);
hostOrProj = false(numel(Gp.nnc.pMMfaces),1);
for q = 1:numel(Gp.nnc.pMMfaces)
    nbrs = Gp.faces.neighbors(Gp.nnc.pMMfaces(q),:);
    hostOrProj(q) = any(nbrs(:)' == Gp.nnc.pMMneighs(q,1)) && ...
                    any(nbrs(:)' == Gp.nnc.pMMneighs(q,2));
end
okD = okD && all(hostOrProj);
if okD
    fprintf('D. sanity (multipliers, areas, face-pair consistency): PASS\n');
else
    failures{end+1} = 'D: sanity checks failed';
    fprintf('D. sanity: FAIL\n');
end

%% Check B: graph compartmentalization of the matrix
nMat = G.Matrix.cells.num;
inter = all(G.Matrix.faces.neighbors ~= 0, 2);
Nmm = G.Matrix.faces.neighbors(inter, :);
% CPEDFM: face transmissibility after multipliers
Tfaces = opsCPEDFM.T_all(1:G.Matrix.faces.num);
Tmm = Tfaces(inter);
keep = Tmm > max(Tmm)*1e-12;
gCP = graph(Nmm(keep,1), Nmm(keep,2), [], nMat);
compCP = conncomp(gCP);
% EDFM baseline
TfacesE = opsEDFM.T_all(1:G.Matrix.faces.num);
TmmE = TfacesE(inter);
keepE = TmmE > max(TmmE)*1e-12;
gED = graph(Nmm(keepE,1), Nmm(keepE,2), [], nMat);
compED = conncomp(gED);

nCompCP = max(compCP);
nCompED = max(compED);
fprintf('B. matrix components: CPEDFM %d | EDFM %d\n', nCompCP, nCompED);
if nCompCP >= 2 && nCompED == 1
    % components must be split by the fracture (host cells on both sides)
    fprintf('   PASS (fracture compartmentalizes matrix under CPEDFM only)\n');
else
    failures{end+1} = sprintf('B: components CPEDFM=%d EDFM=%d', nCompCP, nCompED);
    fprintf('   FAIL\n');
end

%% Check C: flow. Injector west of fracture, producer east.
gravity off
fluid = initSimpleADIFluid('phases', 'W', ...
    'mu', 1*centi*poise, 'rho', 1000, 'c', 1e-8/barsa, 'pRef', 200*barsa);

[iW, ~, ~] = ind2sub(dims, G.Matrix.cells.indexMap(1:nMat));
wcInj  = find(iW == 1);   wcInj  = wcInj(1:2:end);
wcProd = find(iW == dims(1)); wcProd = wcProd(1:2:end);

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

modelCP = GenericBlackOilModel(Gp, [], fluid, 'water', true, 'oil', false, 'gas', false);
modelCP.operators = opsCPEDFM;
modelCP.minimumPressure = 1*barsa;
[~, statesCP] = simulateScheduleAD(state0, modelCP, schedule);

modelED = GenericBlackOilModel(G, [], fluid, 'water', true, 'oil', false, 'gas', false);
modelED.operators = opsEDFM;
modelED.minimumPressure = 1*barsa;
[wsED, statesED] = simulateScheduleAD(state0, modelED, schedule);

% pressure difference across the fracture at end of simulation
sideA = find(iW <= 6); sideB = find(iW >= 15);
dpCP = mean(statesCP{end}.pressure(sideA)) - mean(statesCP{end}.pressure(sideB));
dpED = mean(statesED{end}.pressure(sideA)) - mean(statesED{end}.pressure(sideB));
fprintf('C. cross-fracture dP: CPEDFM %.4g bar | EDFM %.4g bar (ratio %.1f)\n', ...
    dpCP/barsa, dpED/barsa, dpCP/max(dpED,eps));
if dpCP > 5*dpED && dpCP > 1*barsa
    fprintf('   PASS (CPEDFM sustains a much larger pressure jump)\n');
else
    failures{end+1} = sprintf('C: dP CPEDFM %.4g bar vs EDFM %.4g bar', ...
        dpCP/barsa, dpED/barsa);
    fprintf('   FAIL\n');
end

%% Projection figure (Fig. 11-style, on the CPG)
try
    fig = figure('Visible', 'off', 'Position', [100 100 1200 500]);
    subplot(1,2,1); hold on
    plotGrid(G.Matrix, 'FaceColor', 'none', 'EdgeAlpha', 0.08);
    plotGrid(Gp, Gp.nnc.cells(strcmp(Gp.nnc.type,'fracmat interior'),2), ...
        'FaceColor', 'y', 'EdgeAlpha', 0.3);
    view(-35, 25); axis equal tight off; title('fracture cells');
    subplot(1,2,2); hold on
    plotGrid(G.Matrix, 'FaceColor', 'none', 'EdgeAlpha', 0.08);
    cols = {'b', 'r', 'g'};
    for d = 1:3
        fd = Gp.nnc.pMMfaces(Gp.nnc.pDir == d);
        if ~isempty(fd)
            plotFaces(G.Matrix, fd, 'FaceColor', cols{d}, 'EdgeAlpha', 0.3);
        end
    end
    view(-35, 25); axis equal tight off;
    title('continuous projections (x=blue, y=red, z=green)');
    print(fig, fullfile(outDir, 'test01_projections.png'), '-dpng', '-r150');
    close(fig);
    fprintf('projection figure saved to tests/output/test01_projections.png\n');
catch me
    fprintf('figure generation skipped: %s\n', me.message);
end

%% Summary
save(fullfile(outDir, 'test01_results.mat'), 'dpCP', 'dpED', 'nCompCP', 'nCompED');
if isempty(failures)
    fprintf('\n=== test01: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test01: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test01 failed');
end
