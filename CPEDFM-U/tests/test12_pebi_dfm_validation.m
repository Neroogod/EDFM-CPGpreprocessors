%% test12_pebi_dfm_validation
% Quantitative validation of CPEDFM-U against an explicit (DFM)
% reference, following the validation pattern of the CPEDFM paper
% (Figs. 14-20): west-to-east single-phase flow across a sealing
% fracture, driven by fixed boundary pressures.
%
%   REFERENCE : fracture-conforming PEBI grid (pebiGrid2D with the
%               fracture line as a cell constraint, extruded); the
%               fracture cells are actual grid cells with 1 nD perm.
%   CPEDFM-U  : NON-conforming extruded PEBI grid + EDFMGridU +
%               flip-cut projection.
%   EDFM      : same non-conforming grid, no projection (leaky control).
%
% Fracture: vertical, slanted 16 degrees from the y-axis in plan view,
% cutting the full domain cross-section.
%
% Checks:
%   A. Pressure profile along the mid-domain probe line: relative L2
%      error of CPEDFM-U vs the DFM reference is small (< 8%) and at
%      least 2x smaller than the EDFM error.
%   B. Cross-fracture pressure jump: CPEDFM-U reproduces the reference
%      jump; EDFM underestimates it badly.
%   C. West-boundary influx: CPEDFM-U within 25% of the reference;
%      EDFM overestimates (reported).
%
% Writes tests/output/test12_validation.png (pressure fields + probe).

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupCPEDFMU.m'));
outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
failures = {};

physdim = [100, 100, 20];
nLayers = 6;
pWest = 200*barsa; pEast = 100*barsa;
matPerm = 100*milli*darcy; matPoro = 0.25;
fracPerm = 1*nano*darcy;  fracPoro = 0.4;

% fracture line in plan view (x as function of y), full z extent
% endpoints inside the domain (for the conforming constraint) and
% extended beyond (for the EDFM fracture plane)
lineIn = [35 + 30*(0+2)/104,   0;
          35 + 30*(100+2)/104, 100];

gravity off
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, ...
    'rho', 1000, 'c', 1e-8/barsa, 'pRef', 150*barsa);

rng(7);

%% ---------------- DFM reference (fracture-conforming PEBI) -------------
G2Dref = pebiGrid2D(5, physdim(1:2), ...
    'cellConstraints', {lineIn}, ...
    'CCRefinement', true, 'CCFactor', 0.15, 'CCEps', 1/25);
nc2D = G2Dref.cells.num;
Gref = makeLayeredGrid(G2Dref, nLayers);
Gref.nodes.coords(:,3) = Gref.nodes.coords(:,3) * physdim(3)/nLayers;
Gref = computeGeometry(Gref);

% fracture cells: replicate the 2D constraint tag across layers, then
% verify geometrically
tag2D = G2Dref.cells.tag(:);
tag3D = repmat(tag2D, nLayers, 1);
fracCellsRef = find(tag3D);
dLine = pointLineDist2D(Gref.cells.centroids(fracCellsRef,1:2), ...
    lineIn(1,:), lineIn(2,:));
assert(max(dLine) < 2.0, ...
    'layer/tag mapping check failed (max frac-cell offset %.2f m)', max(dLine));

permRef = matPerm*ones(Gref.cells.num, 1);
poroRef = matPoro*ones(Gref.cells.num, 1);
permRef(fracCellsRef) = fracPerm;
poroRef(fracCellsRef) = fracPoro;
Gref.rock = makeRock(Gref, permRef, poroRef);

fprintf('reference grid: %d cells (%d fracture cells)\n', ...
    Gref.cells.num, numel(fracCellsRef));

[pRef, fluxRefW] = solveSteady(Gref, Gref.rock, fluid, [], ...
    pWest, pEast, physdim, tol);

%% ---------------- CPEDFM-U on a non-conforming PEBI grid ----------------
G2D = pebiGrid2D(physdim(1)/12, physdim(1:2));
G = makeLayeredGrid(G2D, nLayers);
G.nodes.coords(:,3) = G.nodes.coords(:,3) * physdim(3)/nLayers;
G = computeGeometry(G);
G.rock = makeRock(G, matPerm, matPoro);

fp = struct;
dirv = (lineIn(2,:) - lineIn(1,:)); dirv = dirv/norm(dirv);
pA = lineIn(1,:) - 5*dirv;  pB = lineIn(2,:) + 5*dirv;
fp(1).points = [pA(1) pA(2) -5;
                pB(1) pB(2) -5;
                pB(1) pB(2) physdim(3)+5;
                pA(1) pA(2) physdim(3)+5];
fp(1).aperture = 10*milli*meter;
fp(1).poro = fracPoro;
fp(1).perm = fracPerm;

[G, fp] = EDFMGridU(G, fp, 'Tolerance', tol, 'fracturelist', 1);
G = fracturematrixNNC_U(G, tol);

Gcp = pMatFracNNCsCPEDFM_U(G, tol);
opsCP = setupCPEDFMOpsTPFA_CPG(Gcp, Gcp.rock, tol);
opsED = setupShaleEDFMOpsTPFA(G, G.rock, tol);

nMat = G.Matrix.cells.num;
fprintf('non-conforming grid: %d matrix cells, %d sealed faces\n', ...
    nMat, numel(unique(Gcp.nnc.pMMfaces)));

[pCP, fluxCPW] = solveSteady(Gcp, Gcp.rock, fluid, opsCP, ...
    pWest, pEast, physdim, tol);
[pED, fluxEDW] = solveSteady(G, G.rock, fluid, opsED, ...
    pWest, pEast, physdim, tol);

%% ---------------- Comparison --------------------------------------------
% probe line at y = 50, z = 10 (mid-domain), excluding the barrier zone
xq = linspace(2, physdim(1)-2, 60)';
yq = 50*ones(size(xq)); zq = 10*ones(size(xq));

matRef = setdiff((1:Gref.cells.num)', fracCellsRef);
Fr = scatteredInterpolant(Gref.cells.centroids(matRef,1), ...
    Gref.cells.centroids(matRef,2), Gref.cells.centroids(matRef,3), ...
    pRef(matRef), 'linear', 'nearest');
Fc = scatteredInterpolant(G.cells.centroids(1:nMat,1), ...
    G.cells.centroids(1:nMat,2), G.cells.centroids(1:nMat,3), ...
    pCP(1:nMat), 'linear', 'nearest');
Fe = scatteredInterpolant(G.cells.centroids(1:nMat,1), ...
    G.cells.centroids(1:nMat,2), G.cells.centroids(1:nMat,3), ...
    pED(1:nMat), 'linear', 'nearest');

pr = Fr(xq, yq, zq); pc = Fc(xq, yq, zq); pe = Fe(xq, yq, zq);

% exclude a band of one coarse cell around the fracture trace where both
% discretizations smear the discontinuity
xFrac = lineIn(1,1) + (50 - lineIn(1,2))*(lineIn(2,1)-lineIn(1,1))/(lineIn(2,2)-lineIn(1,2));
band = abs(xq - xFrac) > physdim(1)/12;

errCP = norm(pc(band) - pr(band)) / norm(pr(band) - mean(pr(band)));
errED = norm(pe(band) - pr(band)) / norm(pr(band) - mean(pr(band)));
fprintf('A. probe-line rel L2 error: CPEDFM-U %.3f | EDFM %.3f\n', errCP, errED);
if errCP < 0.08 && errED > 2*errCP
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('A: errCP=%.3f errED=%.3f', errCP, errED);
    fprintf('   FAIL\n');
end

% cross-fracture pressure jump (mean plateau difference)
west = xq < xFrac - 10; east = xq > xFrac + 10;
dpRef = mean(pr(west)) - mean(pr(east));
dpCP  = mean(pc(west)) - mean(pc(east));
dpED  = mean(pe(west)) - mean(pe(east));
fprintf('B. cross-fracture dP [bar]: ref %.2f | CPEDFM-U %.2f | EDFM %.2f\n', ...
    dpRef/barsa, dpCP/barsa, dpED/barsa);
if abs(dpCP - dpRef)/dpRef < 0.15 && dpED < 0.7*dpRef
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('B: dP ref %.2f, CP %.2f, ED %.2f bar', ...
        dpRef/barsa, dpCP/barsa, dpED/barsa);
    fprintf('   FAIL\n');
end

% west-boundary influx. With a sealing fracture both the reference and
% CPEDFM-U influxes are (near) zero, so they must be compared on the
% scale of the EDFM leak, not relative to each other.
fprintf('C. west influx [m^3/d]: ref %.4g | CPEDFM-U %.4g | EDFM %.4g\n', ...
    fluxRefW*day, fluxCPW*day, fluxEDW*day);
if abs(fluxCPW - fluxRefW) < 1e-3*abs(fluxEDW)
    fprintf('   PASS (CPEDFM-U leak matches reference to <0.1%% of the EDFM leak)\n');
else
    failures{end+1} = sprintf('C: flux ref %.4g, CP %.4g, ED %.4g', ...
        fluxRefW*day, fluxCPW*day, fluxEDW*day);
    fprintf('   FAIL\n');
end

%% ---------------- Case 2: finite-conductivity barrier -------------------
% The sealing case above produces two plateaus - a binary outcome. A
% finite-permeability barrier makes the cross-fracture pressure jump a
% quantitative function of the projected NNC transmissibilities (Eq. 10),
% testing the full-face area rule, not just the seal.
%
% For a like-for-like comparison the EDFM fracture aperture is set to
% the measured mean width of the reference grid's fracture cells, so the
% cross-fracture resistance w/k is identical in both models.
fracPerm2 = 0.02*milli*darcy;

% measured mean width of the conforming fracture cells (from the 3D
% volumes: total fracture-cell volume = lineLen * width * height)
lineLen = norm(lineIn(2,:) - lineIn(1,:));
wMeas = sum(Gref.cells.volumes(fracCellsRef)) / (lineLen * physdim(3));
fprintf('case 2: finite-conductivity barrier (k = %.3g mD, w = %.2f m)\n', ...
    fracPerm2/(milli*darcy), wMeas);

% reference with finite-permeability fracture cells
permRef2 = permRef; permRef2(fracCellsRef) = fracPerm2;
Gref.rock = makeRock(Gref, permRef2, poroRef);
[pRef2, ~] = solveSteady(Gref, Gref.rock, fluid, [], pWest, pEast, physdim, tol);

% CPEDFM-U with matched aperture and permeability
G2Db = pebiGrid2D(physdim(1)/12, physdim(1:2));
Gb = makeLayeredGrid(G2Db, nLayers);
Gb.nodes.coords(:,3) = Gb.nodes.coords(:,3) * physdim(3)/nLayers;
Gb = computeGeometry(Gb);
Gb.rock = makeRock(Gb, matPerm, matPoro);

% fresh fracplane struct (the earlier one carries processed=true from
% its own gridding call and would be skipped by EDFMGridU)
fp2 = struct('points', fp(1).points, 'aperture', wMeas, ...
    'poro', fracPoro, 'perm', fracPerm2);
[Gb, fp2] = EDFMGridU(Gb, fp2, 'Tolerance', tol, 'fracturelist', 1);
Gb = fracturematrixNNC_U(Gb, tol);
Gbp = pMatFracNNCsCPEDFM_U(Gb, tol);
opsCP2 = setupCPEDFMOpsTPFA_CPG(Gbp, Gbp.rock, tol);
[pCP2, ~] = solveSteady(Gbp, Gbp.rock, fluid, opsCP2, pWest, pEast, physdim, tol);

nMatB = Gb.Matrix.cells.num;
Fr2 = scatteredInterpolant(Gref.cells.centroids(matRef,1), ...
    Gref.cells.centroids(matRef,2), Gref.cells.centroids(matRef,3), ...
    pRef2(matRef), 'linear', 'nearest');
Fc2 = scatteredInterpolant(Gb.cells.centroids(1:nMatB,1), ...
    Gb.cells.centroids(1:nMatB,2), Gb.cells.centroids(1:nMatB,3), ...
    pCP2(1:nMatB), 'linear', 'nearest');
pr2 = Fr2(xq, yq, zq); pc2 = Fc2(xq, yq, zq);

errCP2 = norm(pc2(band) - pr2(band)) / norm(pr2(band) - mean(pr2(band)));
dpRef2 = mean(pr2(west)) - mean(pr2(east));
dpCP2  = mean(pc2(west)) - mean(pc2(east));
fprintf('D. finite-k probe rel L2 error: %.3f | dP ref %.2f bar, CPEDFM-U %.2f bar\n', ...
    errCP2, dpRef2/barsa, dpCP2/barsa);
if errCP2 < 0.15 && abs(dpCP2 - dpRef2)/dpRef2 < 0.30
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('D: err %.3f, dP ref %.2f vs CP %.2f bar', ...
        errCP2, dpRef2/barsa, dpCP2/barsa);
    fprintf('   FAIL\n');
end

%% ---------------- Figure -------------------------------------------------
try
    fig = figure('Visible', 'off', 'Position', [50 50 1500 420]);
    ttl = {'DFM reference (conforming PEBI)', 'CPEDFM-U (non-conforming)', 'EDFM (leaky)'};
    Gs = {Gref, G.Matrix, G.Matrix};
    Ps = {pRef, pCP(1:nMat), pED(1:nMat)};
    for s = 1:3
        subplot(1,4,s);
        plotCellData(Gs{s}, Ps{s}/barsa, 'EdgeAlpha', 0.05);
        view(0,90); axis equal tight off; caxis([pEast pWest]/barsa);
        title(ttl{s});
    end
    subplot(1,4,4); hold on
    plot(xq, pr/barsa, 'k-', 'LineWidth', 1.5);
    plot(xq, pc/barsa, 'b--', 'LineWidth', 1.5);
    plot(xq, pe/barsa, 'r:', 'LineWidth', 1.5);
    plot(xq, pr2/barsa, 'k-', 'LineWidth', 0.8);
    plot(xq, pc2/barsa, 'b--', 'LineWidth', 0.8);
    xlabel('x [m]'); ylabel('pressure [bar]');
    legend('DFM ref (sealing)', 'CPEDFM-U (sealing)', 'EDFM (sealing)', ...
        'DFM ref (finite k)', 'CPEDFM-U (finite k)', 'Location', 'southwest');
    title('probe line y=50, z=10'); grid on
    print(fig, fullfile(outDir, 'test12_validation.png'), '-dpng', '-r150');
    close(fig);
    fprintf('figure saved to tests/output/test12_validation.png\n');
catch me
    fprintf('figure generation skipped: %s\n', me.message);
end

%% Summary
save(fullfile(outDir, 'test12_results.mat'), 'errCP', 'errED', ...
    'dpRef', 'dpCP', 'dpED', 'fluxRefW', 'fluxCPW', 'fluxEDW', ...
    'errCP2', 'dpRef2', 'dpCP2', 'wMeas');
if isempty(failures)
    fprintf('\n=== test12: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test12: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test12 failed');
end

%% ---------------- local helpers ------------------------------------------
function [p, influxW] = solveSteady(G, rock, fluid, ops, pW, pE, physdim, tol)
% single-phase flow to (near) steady state under fixed west/east
% boundary pressures; returns final pressure and west-boundary influx
bfaces = find(any(G.faces.neighbors == 0, 2));
if isfield(G, 'Matrix')
    bfaces = bfaces(bfaces <= G.Matrix.faces.num);
end
fc = G.faces.centroids(bfaces,:);
wf = bfaces(fc(:,1) < tol*1e3);
ef = bfaces(fc(:,1) > physdim(1) - tol*1e3);
assert(~isempty(wf) && ~isempty(ef), 'boundary faces not found');

bc = addBC([], wf, 'pressure', pW, 'sat', 1);
bc = addBC(bc, ef, 'pressure', pE, 'sat', 1);

model = GenericBlackOilModel(G, [], fluid, ...
    'water', true, 'oil', false, 'gas', false);
if ~isempty(ops)
    model.operators = ops;
else
    model.operators = setupOperatorsTPFA(G, rock);
end
model.minimumPressure = 1*barsa;

state0 = initResSol(G, (pW + pE)/2, 1);
dt = rampupTimesteps(2*year, 100*day, 5);
schedule = simpleSchedule(dt, 'bc', bc);
[~, states] = simulateScheduleAD(state0, model, schedule);
p = states{end}.pressure;

% west influx from the final state's boundary fluxes
flux = states{end}.flux(:, 1);
sgn = 1 - 2*(G.faces.neighbors(wf,1) == 0); % positive = into the domain
influxW = sum(sgn .* flux(wf) .* -1);
influxW = abs(influxW);
end

function d = pointLineDist2D(P, a, b)
ab = b - a;
t = max(min(((P - a)*ab')/(ab*ab'), 1), 0);
proj = a + t.*ab;
d = sqrt(sum((P - proj).^2, 2));
end
