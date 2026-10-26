%% examplePebiFracturedWaterflood_U
% Waterflood of a naturally fractured reservoir on an unstructured
% (extruded PEBI) grid, modelled with EDFM-U.
%
% Setup
%   - 400 x 400 x 40 m reservoir, extruded PEBI grid (~2500 cells), 50 mD
%     matrix.
%   - A five-plane fracture network: a long conductive backbone with two
%     branches crossing it (so fracture-fracture NNCs are exercised), plus
%     two shorter isolated fractures.
%   - Water injector in the west, oil producer in the east; the network
%     lies between them and channels the flood.
%
% The model is run twice: once with EDFM-U's exact <d>, and once with the
% connectivity indices rebuilt from the bounding-box <d> of the Cartesian
% implementation. On THIS problem the two are indistinguishable, and that
% is the point worth taking away: a bhp-driven waterflood through a
% conductive network is transport-dominated, so even a 7% median error in
% <d> leaves cumulative oil unchanged. The exact <d> matters where
% matrix-to-fracture transfer is the rate-limiting step - tight matrix,
% low-conductivity or partially sealing features, early-time transients,
% and unstructured cells far from cubic. tests/test03 and tests/test04
% quantify both regimes.
%
% Runtime: a couple of minutes. Figures go to examples/output/.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupEDFMU.m'));
outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
gravity off
mrstVerbose off
rng(21);

%% Reservoir: extruded PEBI grid
physdim = [400, 400, 40];
nLayers = 5;
G2D = pebiGrid2D(physdim(1)/18, physdim(1:2));
G = makeLayeredGrid(G2D, nLayers);
G.nodes.coords(:,3) = G.nodes.coords(:,3) * physdim(3)/nLayers;
G = computeGeometry(G);
G.rock = makeRock(G, 50*milli*darcy, 0.22);
nMat = G.cells.num;
fprintf('PEBI reservoir: %d cells, %d..%d faces per cell\n', nMat, ...
    min(diff(G.cells.facePos)), max(diff(G.cells.facePos)));

%% Fracture network
zb = 2; zt = physdim(3) - 2;
fracplanes = struct;

% (1) conductive backbone, SW-NE across the middle of the domain
fracplanes(1).points = [110  90 zb; 300 300 zb; 300 300 zt; 110  90 zt];
% (2)-(3) branches crossing the backbone
fracplanes(2).points = [150 260 zb; 260 120 zb; 260 120 zt; 150 260 zt];
fracplanes(3).points = [200 320 zb; 310 170 zb; 310 170 zt; 200 320 zt];
% (4)-(5) isolated fractures nearer the wells
fracplanes(4).points = [ 70 250 zb; 130 320 zb; 130 320 zt;  70 250 zt];
fracplanes(5).points = [300  70 zb; 350 150 zb; 350 150 zt; 300  70 zt];

for i = 1:numel(fracplanes)
    fracplanes(i).aperture = 5*milli*meter;
    fracplanes(i).poro     = 0.5;
    fracplanes(i).perm     = 2*darcy;
end

%% EDFM-U pipeline
[G, fracplanes] = EDFMgrid_U(G, fracplanes, 'Tolerance', tol);
G = fracturematrixNNC_U(G, tol);
[G, fracplanes] = fracturefractureNNC_U(G, fracplanes, tol, 'Verbose', false);

isFF = strcmp(G.nnc.type, 'fracfrac');
fprintf('NNCs: %d fracture-matrix, %d fracture-fracture\n', sum(~isFF), sum(isFF));
fprintf('<d> range: %.3g to %.3g m\n', min(G.nnc.davg), max(G.nnc.davg));

ops = setupShaleEDFMOpsTPFA(G, G.rock, tol);

%% Legacy bounding-box <d> variant, for comparison
% Same grid, same connections; CI (and hence T) rebuilt from the
% bounding-box <d>. T is proportional to 1/<d>.
Gb = G;
[cn, cpos] = gridCellNodes(G.Matrix, 1:G.Matrix.cells.num);
fmRows = find(~isFF);
dBB = zeros(numel(fmRows), 1);
for j = 1:numel(fmRows)
    r  = fmRows(j);
    c  = G.nnc.cells(r,1);
    nd = G.Matrix.nodes.coords(cn(cpos(c):cpos(c+1)-1), :);
    dBB(j) = calcdavg(nd, G.nnc.normal(r,:), G.nnc.planepoint(r,:), tol);
end
Gb.nnc.T(fmRows) = G.nnc.T(fmRows) .* (G.nnc.davg(fmRows) ./ dBB);
dErr = abs(dBB - G.nnc.davg(fmRows)) ./ G.nnc.davg(fmRows);
fprintf('bounding-box <d> error on this grid: median %.1f%%, max %.1f%%\n', ...
    100*median(dErr), 100*max(dErr));
opsB = setupShaleEDFMOpsTPFA(Gb, Gb.rock, tol);

%% Fluid and wells
fluid = initSimpleADIFluid('phases', 'WO', ...
    'mu',  [1, 4]*centi*poise, ...
    'rho', [1000, 800], ...
    'n',   [2, 2], ...
    'c',   [1e-6, 1e-5]/barsa, ...
    'pRef', 200*barsa);

cc = G.Matrix.cells.centroids;
zL = unique(round(cc(:,3), 6));
wInj = zeros(nLayers,1); wPrd = zeros(nLayers,1);
for L = 1:nLayers
    inL = find(abs(cc(:,3) - zL(L)) < 1e-6);
    [~, ii] = min((cc(inL,1) -  60).^2 + (cc(inL,2) -  60).^2);
    [~, ip] = min((cc(inL,1) - 340).^2 + (cc(inL,2) - 340).^2);
    wInj(L) = inL(ii);
    wPrd(L) = inL(ip);
end

W = addWell([], G.Matrix, G.rock, wInj, 'Type', 'bhp', 'Val', 320*barsa, ...
    'Comp_i', [1 0], 'Name', 'INJ', 'Radius', 0.1, 'Sign', 1);
W = addWell(W, G.Matrix, G.rock, wPrd, 'Type', 'bhp', 'Val', 120*barsa, ...
    'Comp_i', [1 0], 'Name', 'PROD', 'Radius', 0.1, 'Sign', -1);

totTime  = 8*year;
schedule = simpleSchedule(rampupTimesteps(totTime, year/3, 6), 'W', W);
state0   = initResSol(G, 200*barsa, [0 1]);

%% Simulate
model = GenericBlackOilModel(G, [], fluid, 'water', true, 'oil', true, 'gas', false);
model.operators = ops;
[wsU, statesU] = simulateScheduleAD(state0, model, schedule);

modelB = GenericBlackOilModel(Gb, [], fluid, 'water', true, 'oil', true, 'gas', false);
modelB.operators = opsB;
[wsB, statesB] = simulateScheduleAD(state0, modelB, schedule);

%% Results
t   = cumsum(schedule.step.val);
qoU = -cellfun(@(w) w(2).qOs, wsU);
qoB = -cellfun(@(w) w(2).qOs, wsB);
qwU = -cellfun(@(w) w(2).qWs, wsU);
cumU = trapz(t, qoU);
cumB = trapz(t, qoB);

fprintf('\nCumulative oil after %d years [m^3]: EDFM-U %.4g | bounding-box <d> %.4g (%+.2f%%)\n', ...
    round(totTime/year), cumU, cumB, 100*(cumB-cumU)/cumU);
fprintf(['  A %.1f%% median error in <d> moves cumulative oil by %.2f%%: this ', ...
         'waterflood is\n  transport-dominated, so it is insensitive to CI. ', ...
         'See tests/test04 case 2 for a\n  transfer-limited case, where the ', ...
         'same <d> error is what sets the answer.\n'], ...
    100*median(dErr), abs(100*(cumB-cumU)/cumU));
bt = find(qwU > 0.01*max(qwU), 1);
if ~isempty(bt)
    fprintf('Water breakthrough at the producer: %.1f years\n', t(bt)/year);
end

%% Figures
fig = figure('Visible','off','Position',[50 50 1400 460]);
subplot(1,3,1); hold on
plotGrid(G.Matrix, 'FaceColor','none','EdgeAlpha',0.05);
cols = lines(numel(fieldnames(G.FracGrid)));
fn = fieldnames(G.FracGrid);
for q = 1:numel(fn)
    Fg = G.FracGrid.(fn{q});
    plotGrid(G, Fg.cells.start:(Fg.cells.start+Fg.cells.num-1), ...
        'FaceColor', cols(q,:), 'EdgeAlpha', 0.2);
end
plotGrid(G.Matrix, [wInj; wPrd], 'FaceColor','k');
view(-30, 45); axis equal tight off
title(sprintf('fracture network (%d planes) + wells', numel(fn)));

topCells = find(abs(cc(:,3) - zL(1)) < 1e-6);
subplot(1,3,2);
plotCellData(G.Matrix, statesU{end}.s(1:nMat,1), topCells, 'EdgeAlpha', 0.05);
view(0,90); axis equal tight off; caxis([0 1]); colorbar
title(sprintf('EDFM-U water saturation, %d years', round(totTime/year)));

subplot(1,3,3); hold on
plot(t/year, qoU*day, 'b-',  'LineWidth', 1.8);
plot(t/year, qoB*day, 'r--', 'LineWidth', 1.5);
xlabel('time [years]'); ylabel('oil rate [m^3/day]'); grid on
legend('EDFM-U (exact <d>)', 'bounding-box <d>', 'Location','northeast');
title('producer oil rate');
print(fig, fullfile(outDir,'example_edfmu.png'), '-dpng', '-r150');
close(fig);
fprintf('figure saved to examples/output/example_edfmu.png\n');

if usejava('desktop')
    figure; plotToolbar(G, statesU); view(-30, 45); colormap jet
    title('EDFM-U waterflood (interactive)');
end
