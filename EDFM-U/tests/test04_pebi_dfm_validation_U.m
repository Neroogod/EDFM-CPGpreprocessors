%% test04_pebi_dfm_validation_U
% Flow validation of EDFM-U on a NON-CONFORMING unstructured grid against
% an explicit fracture-conforming (DFM) reference.
%
%   REFERENCE  : fracture-conforming PEBI grid; the fracture is a band of
%                real grid cells whose conductivity k*w matches the EDFM
%                fracture's k*aperture.
%   EDFM-U     : non-conforming extruded PEBI grid + EDFMgrid_U +
%                fracturematrixNNC_U.
%   LEGACY <d> : identical to EDFM-U except that CI is rebuilt from the
%                bounding-box <d> of hfm/calcdavg (T is proportional to
%                1/<d>). This isolates the effect of the exact <d> at flow
%                level on the same grid and the same connections.
%
% Two flow regimes are run, because they are sensitive to different things.
%
%   CASE 1 - conductive throughflow. A conductive fracture lies along the
%   flow direction, internal to the domain, driven by west/east boundary
%   pressures. This validates the whole pipeline end to end, but it is
%   deliberately reported as INSENSITIVE to <d>: the fracture is close to
%   an equipotential, so matrix-fracture transfer is not the rate-limiting
%   step and a large error in CI barely moves the throughflow.
%
%   CASE 2 - transfer-limited drainage. The domain is closed and the
%   fracture is pinned at a low pressure, so the entire response is the
%   matrix draining into the fracture. Here CI is rate-limiting and the
%   accuracy of <d> shows up directly in the matrix pressure decline.
%
% Checks
%   A. Case 1 throughflow: EDFM-U within 10% of the DFM reference.
%   B. Case 1 probe-line pressure through the fracture: rel L2 < 10%.
%   C. Case 2 matrix decline curve: EDFM-U within 15% of the reference AND
%      closer to it than the legacy-<d> variant. See the note at check C
%      for why the absolute tolerance is looser here than in case 1.
%
% Writes tests/output/test04_validation.png.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupEDFMU.m'));
outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
failures = {};
gravity off
mrstVerbose off
rng(9);

physdim  = [100, 100, 20];
nLayers  = 4;
pWest    = 200*barsa;
pEast    = 100*barsa;
pInit    = 200*barsa;
pDrain   = 100*barsa;
matPerm  = 1*milli*darcy;   matPoro  = 0.20;
fracPerm = 10*darcy;        fracPoro = 0.40;
aperture = 5*milli*meter;
bandHalf = 6;               % matrix outside |y-50| > bandHalf is compared

lineIn = [20 50; 80 50];    % internal fracture along x at y = 50

fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, ...
    'rho', 1000, 'c', 1e-8/barsa, 'pRef', 150*barsa);

% Case 2 needs a resolvable transient. With the near-incompressible fluid
% above the hydraulic diffusivity k/(phi*mu*ct) is ~49 m^2/s, so a 100 m
% block equilibrates in about a minute and every model returns the same
% converged answer - the comparison would be vacuous. A water-like
% compressibility puts the diffusivity at ~1e-2 m^2/s, i.e. a few days
% across the block, which is the transfer-limited window CI controls.
fluidC = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, ...
    'rho', 1000, 'c', 5e-5/barsa, 'pRef', 150*barsa);

%% ---------------- DFM reference (fracture-conforming PEBI) --------------
% A narrow constraint band matters: the reference drains from y = 50 +- w/2
% whereas EDFM drains from the plane y = 50, so a wide band biases the
% reference toward faster early drainage in case 2.
G2Dref = pebiGrid2D(4, physdim(1:2), 'cellConstraints', {lineIn}, ...
    'CCRefinement', true, 'CCFactor', 0.05, 'CCEps', 1/60);
Gref = makeLayeredGrid(G2Dref, nLayers);
Gref.nodes.coords(:,3) = Gref.nodes.coords(:,3) * physdim(3)/nLayers;
Gref = computeGeometry(Gref);

fracCellsRef = find(repmat(G2Dref.cells.tag(:), nLayers, 1));
dLine = pointSegDist2D(Gref.cells.centroids(fracCellsRef,1:2), lineIn(1,:), lineIn(2,:));
assert(max(dLine) < 1.5, 'fracture-cell tagging check failed (max offset %.2f m)', max(dLine));

% match conductivity: k_ref * w_meas = fracPerm * aperture
lineLen  = norm(lineIn(2,:) - lineIn(1,:));
wMeas    = sum(Gref.cells.volumes(fracCellsRef)) / (lineLen * physdim(3));
kRefFrac = fracPerm * aperture / wMeas;

permRef = matPerm*ones(Gref.cells.num,1);  permRef(fracCellsRef) = kRefFrac;
poroRef = matPoro*ones(Gref.cells.num,1);  poroRef(fracCellsRef) = fracPoro;
Gref.rock = makeRock(Gref, permRef, poroRef);

fprintf('reference: %d cells (%d fracture cells), band width %.2f m -> k %.4g D\n', ...
    Gref.cells.num, numel(fracCellsRef), wMeas, kRefFrac/darcy);

%% ---------------- EDFM-U on a non-conforming PEBI grid ------------------
G2D = pebiGrid2D(physdim(1)/16, physdim(1:2));
G = makeLayeredGrid(G2D, nLayers);
G.nodes.coords(:,3) = G.nodes.coords(:,3) * physdim(3)/nLayers;
G = computeGeometry(G);
G.rock = makeRock(G, matPerm, matPoro);
nMat = G.cells.num;

fp = struct('points', [lineIn(1,1) lineIn(1,2) 0;
                       lineIn(2,1) lineIn(2,2) 0;
                       lineIn(2,1) lineIn(2,2) physdim(3);
                       lineIn(1,1) lineIn(1,2) physdim(3)], ...
            'aperture', aperture, 'poro', fracPoro, 'perm', fracPerm);

G = EDFMgrid_U(G, fp, 'Tolerance', tol, 'Verbose', false);
G = fracturematrixNNC_U(G, tol, 'Verbose', false);
nFrac = G.FracGrid.Frac1.cells.num;
fracCellsU = G.FracGrid.Frac1.cells.start + (0:nFrac-1)';
fprintf('EDFM-U: %d matrix cells, %d fracture cells, %d NNCs\n', ...
    nMat, nFrac, size(G.nnc.cells,1));
opsU = setupShaleEDFMOpsTPFA(G, G.rock, tol);

%% ---------------- Legacy bounding-box <d> variant ------------------------
% Same grid, same connections; only CI is rebuilt from the bounding-box
% <d>. Since T is proportional to CI and CI to 1/<d>, rescaling T by
% <d>_exact / <d>_bbox reproduces what the Cartesian estimator would have
% produced on this unstructured grid.
Gl = G;
[cn, cpos] = gridCellNodes(G.Matrix, 1:G.Matrix.cells.num);
mcells = G.nnc.cells(:,1);
nrm = G.nnc.normal(1,:); pp = G.nnc.planepoint(1,:);
dBB = zeros(numel(mcells),1);
for j = 1:numel(mcells)
    nd = G.Matrix.nodes.coords(cn(cpos(mcells(j)):cpos(mcells(j)+1)-1), :);
    dBB(j) = calcdavg(nd, nrm, pp, tol);
end
Gl.nnc.T = G.nnc.T .* (G.nnc.davg ./ dBB);
dErr = abs(dBB - G.nnc.davg)./G.nnc.davg;
fprintf('legacy <d>: |bbox-exact|/exact over %d NNCs, median %.1f%%, max %.1f%%\n', ...
    numel(dBB), 100*median(dErr), 100*max(dErr));
opsL = setupShaleEDFMOpsTPFA(Gl, Gl.rock, tol);

%% ================= CASE 1: conductive throughflow ========================
[pRef, qRef] = steadyFlow(Gref, Gref.rock, fluid, [],   pWest, pEast, physdim);
[pU,   qU]   = steadyFlow(G,    G.rock,    fluid, opsU, pWest, pEast, physdim);
[pL,   qL]   = steadyFlow(Gl,   Gl.rock,   fluid, opsL, pWest, pEast, physdim);

eU = abs(qU - qRef)/qRef;
eL = abs(qL - qRef)/qRef;
fprintf('\nA. throughflow [m^3/day]: ref %.5g | EDFM-U %.5g (%.2f%%) | legacy <d> %.5g (%.2f%%)\n', ...
    qRef*day, qU*day, 100*eU, qL*day, 100*eL);
if eU < 0.10
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('A: throughflow error %.2f%%', 100*eU);
    fprintf('   FAIL\n');
end

xq = linspace(3, physdim(1)-3, 80)';
yq = 50*ones(size(xq));  zq = physdim(3)/2*ones(size(xq));
matRef = setdiff((1:Gref.cells.num)', fracCellsRef);
Fr = scatteredInterpolant(Gref.cells.centroids(matRef,1), ...
    Gref.cells.centroids(matRef,2), Gref.cells.centroids(matRef,3), ...
    pRef(matRef), 'linear', 'nearest');
Fu = matrixInterp(G, pU, nMat);
Fl = matrixInterp(G, pL, nMat);
pr = Fr(xq,yq,zq);  pu = Fu(xq,yq,zq);  pl = Fl(xq,yq,zq);
den  = norm(pr - mean(pr));
errU = norm(pu - pr)/den;
errL = norm(pl - pr)/den;
fprintf('B. probe-line rel L2 error: EDFM-U %.4f | legacy <d> %.4f\n', errU, errL);
fprintf('   (case 1 is transport-limited, so both <d> variants land within %.1f%% ', ...
    100*max(eU,eL));
fprintf('of the reference despite a %.0f%% median <d> error - case 2 is the CI test)\n', ...
    100*median(dErr));
if errU < 0.10
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('B: probe error %.4f', errU);
    fprintf('   FAIL\n');
end

%% ================= CASE 2: transfer-limited drainage =====================
% Closed domain; the fracture is pinned at pDrain by a very high
% productivity completion, so the response is pure matrix -> fracture
% transfer. The comparison metric is the volume-weighted mean pressure of
% the matrix OUTSIDE the reference's fracture band, a region both models
% represent identically.
% Log-spaced steps spanning well past equilibration, so the transient is
% captured wherever it falls; the comparison below then masks to the part
% of the curve that is actually still transient.
dtD = rampupTimesteps(6*day, 0.25*day, 9);

selRef = matRef(abs(Gref.cells.centroids(matRef,2) - 50) > bandHalf);
selU   = find(abs(G.Matrix.cells.centroids(:,2) - 50) > bandHalf);

[tD, mRef] = drainCurve(Gref, Gref.rock, fluidC, [],   fracCellsRef, ...
    pInit, pDrain, dtD, selRef, Gref.cells.volumes(selRef));
[~,  mU]   = drainCurve(G,    G.rock,    fluidC, opsU, fracCellsU, ...
    pInit, pDrain, dtD, selU,   G.Matrix.cells.volumes(selU));
[~,  mL]   = drainCurve(Gl,   Gl.rock,   fluidC, opsL, fracCellsU, ...
    pInit, pDrain, dtD, selU,   Gl.Matrix.cells.volumes(selU));

% Compare only where the reference is genuinely mid-transient. Once it has
% equilibrated at pDrain every model returns the same answer regardless of
% CI, so including that tail would dilute the test into a vacuous pass.
dRef = pInit - mRef;  dU = pInit - mU;  dL = pInit - mL;
frac = dRef / (pInit - pDrain);
mask = frac > 0.05 & frac < 0.90;
fprintf('   transient window: %d of %d samples, t = %.2f to %.2f days, decline %.0f%%-%.0f%%\n', ...
    sum(mask), numel(mask), min(tD(mask))/day, max(tD(mask))/day, ...
    100*min(frac(mask)), 100*max(frac(mask)));
assert(sum(mask) >= 5, ...
    ['case 2 has only %d samples in the transfer-limited transient; ', ...
     'refine the time stepping.'], sum(mask));

cU = norm(dU(mask) - dRef(mask))/norm(dRef(mask));
cL = norm(dL(mask) - dRef(mask))/norm(dRef(mask));
fprintf('\nC. matrix depletion over the transient window (rel L2 vs reference):\n');
fprintf('   EDFM-U %.4f | legacy <d> %.4f\n', cU, cL);
kk = find(mask, 1, 'last');
fprintf('   mean matrix pressure at t = %.2f days [bar]: ref %.2f | EDFM-U %.2f | legacy <d> %.2f\n', ...
    tD(kk)/day, mRef(kk)/barsa, mU(kk)/barsa, mL(kk)/barsa);
% The absolute tolerance here is looser than in case 1 on purpose. A
% transfer-limited regime is the hardest thing to ask of EDFM: the method
% assumes a single average distance <d> per host cell, so with host cells
% several metres across it resolves the near-fracture pressure gradient
% only to first order, and it drains measurably more slowly than a
% conforming reference. That residual is NOT a resolution artefact -
% refining the EDFM grid from 608 to ~1200 cells and narrowing the
% reference band from 2.5 m to ~1.5 m moved it only from 10.8% to 10.4%.
% The discriminating part of this check is therefore the RELATIVE one:
% the exact <d> must land closer to the reference than the bounding-box
% <d>, which overestimates <d> here by 29% (median) and so under-predicts
% transfer even further.
if cU < 0.15 && cU < cL
    fprintf('   PASS (exact <d> is %.1f%% closer to the reference than bbox <d>)\n', ...
        100*(1 - cU/cL));
else
    failures{end+1} = sprintf('C: cU %.4f cL %.4f', cU, cL);
    fprintf('   FAIL\n');
end

%% ---------------- Figure --------------------------------------------------
try
    fig = figure('Visible','off','Position',[50 50 1500 800]);
    topRef = find(Gref.cells.centroids(:,3) < physdim(3)/nLayers);
    topU   = find(G.Matrix.cells.centroids(:,3) < physdim(3)/nLayers);

    subplot(2,3,1);
    plotCellData(Gref, pRef/barsa, topRef, 'EdgeAlpha', 0.05);
    view(0,90); axis equal tight off; caxis([pEast pWest]/barsa); colorbar
    title(sprintf('case 1: DFM reference (%d cells)', Gref.cells.num));

    subplot(2,3,2);
    plotCellData(G.Matrix, pU(1:nMat)/barsa, topU, 'EdgeAlpha', 0.05);
    view(0,90); axis equal tight off; caxis([pEast pWest]/barsa); colorbar
    title(sprintf('case 1: EDFM-U (%d cells)', nMat));

    subplot(2,3,3); hold on
    plot(xq, pr/barsa, 'k-',  'LineWidth', 1.8);
    plot(xq, pu/barsa, 'b--', 'LineWidth', 1.5);
    plot(xq, pl/barsa, 'r:',  'LineWidth', 1.5);
    xlabel('x [m]'); ylabel('pressure [bar]'); grid on
    legend('DFM reference', 'EDFM-U', 'legacy bbox <d>', 'Location','northeast');
    title('case 1: probe line through the fracture');

    subplot(2,3,4); hold on
    plotGrid(G.Matrix, 'FaceColor','none','EdgeAlpha',0.06);
    plotGrid(G, fracCellsU, 'FaceColor','r','EdgeAlpha',0.3);
    view(-25,50); axis equal tight off
    title('EDFM-U fracture cells on the PEBI grid');

    subplot(2,3,5); hold on
    plot(tD/day, mRef/barsa, 'k-',  'LineWidth', 1.8);
    plot(tD/day, mU/barsa,   'b--', 'LineWidth', 1.5);
    plot(tD/day, mL/barsa,   'r:',  'LineWidth', 1.5);
    plot(tD(mask)/day, mRef(mask)/barsa, 'ko', 'MarkerSize', 4);
    set(gca, 'XScale', 'log');
    xlabel('time [days]'); ylabel('mean matrix pressure [bar]'); grid on
    legend('DFM reference', 'EDFM-U', 'legacy bbox <d>', 'compared window', ...
        'Location','southwest');
    title('case 2: transfer-limited decline');

    subplot(2,3,6);
    histogram(100*dErr, 20);
    xlabel('|bbox <d> - exact <d>| / exact  [%]'); ylabel('NNCs'); grid on
    title('per-NNC error of the bounding-box <d>');

    print(fig, fullfile(outDir,'test04_validation.png'), '-dpng', '-r150');
    close(fig);
    fprintf('figure saved to tests/output/test04_validation.png\n');
catch me
    fprintf('figure generation skipped: %s\n', me.message);
end

%% Summary
if isempty(failures)
    fprintf('\n=== test04: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test04: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test04 failed');
end

%% ------------------------------ helpers ---------------------------------
function F = matrixInterp(G, p, nMat)
F = scatteredInterpolant(G.cells.centroids(1:nMat,1), ...
    G.cells.centroids(1:nMat,2), G.cells.centroids(1:nMat,3), ...
    p(1:nMat), 'linear', 'nearest');
end

function model = makeModel(G, rock, fluid, ops)
model = GenericBlackOilModel(G, [], fluid, ...
    'water', true, 'oil', false, 'gas', false);
if isempty(ops)
    model.operators = setupOperatorsTPFA(G, rock);
else
    model.operators = ops;
end
model.minimumPressure = 1*barsa;
end

function [p, q] = steadyFlow(G, rock, fluid, ops, pW, pE, physdim)
% single-phase flow to steady state under fixed west/east boundary
% pressures; returns the pressure field and the west-boundary influx
bf = find(any(G.faces.neighbors == 0, 2));
if isfield(G, 'Matrix')
    bf = bf(bf <= G.Matrix.faces.num);
end
fc = G.faces.centroids(bf,:);
wf = bf(fc(:,1) < 1e-2);
ef = bf(fc(:,1) > physdim(1) - 1e-2);
assert(~isempty(wf) && ~isempty(ef), 'boundary faces not found');

bc = addBC([], wf, 'pressure', pW, 'sat', 1);
bc = addBC(bc, ef, 'pressure', pE, 'sat', 1);

state0   = initResSol(G, (pW+pE)/2, 1);
schedule = simpleSchedule(rampupTimesteps(5*year, 200*day, 6), 'bc', bc);
[~, states] = simulateScheduleAD(state0, makeModel(G, rock, fluid, ops), ...
    schedule, 'Verbose', false);

p = states{end}.pressure;
q = abs(sum(states{end}.flux(wf, 1)));
end

function [t, meanP] = drainCurve(G, rock, fluid, ops, drainCells, ...
                                 pInit, pDrain, dt, selCells, selVol)
% Closed domain; drainCells are pinned at pDrain by a completion whose
% productivity index is made large enough to be non-limiting. Returns the
% volume-weighted mean pressure of selCells over time.
model = makeModel(G, rock, fluid, ops);

W = addWell([], G, rock, drainCells, 'Type', 'bhp', 'Val', pDrain, ...
    'Comp_i', 1, 'Name', 'DRAIN', 'Radius', 0.1, 'Sign', -1, ...
    'refDepth', 0);
% pin rather than throttle: PI far above every connected transmissibility
Tscale = max([model.operators.T; 1e-12]);
W(1).WI = 1e4 * Tscale * ones(numel(W(1).cells), 1);

state0   = initResSol(G, pInit, 1);
schedule = simpleSchedule(dt, 'W', W);
[~, states] = simulateScheduleAD(state0, model, schedule, 'Verbose', false);

t     = cumsum(dt);
meanP = cellfun(@(s) sum(s.pressure(selCells).*selVol)/sum(selVol), states);
end

function d = pointSegDist2D(P, a, b)
ab = b - a;
t  = max(min(((P - a)*ab')/(ab*ab'), 1), 0);
d  = sqrt(sum((P - (a + t.*ab)).^2, 2));
end
