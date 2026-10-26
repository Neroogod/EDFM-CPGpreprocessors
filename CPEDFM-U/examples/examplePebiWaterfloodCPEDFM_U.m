%% examplePebiWaterfloodCPEDFM_U
% Waterflood in an extruded PEBI (Voronoi) reservoir with a sealing
% fault and conductive natural fractures, modeled with CPEDFM on the
% unstructured grid.
%
% Setup:
%   - 500 x 500 x 50 m reservoir, extruded PEBI grid (~2000 cells).
%   - A sealing fault (1 nD), slanted in plan view and dipping, cutting
%     the full reservoir between the wells.
%   - Two conductive natural fractures (10 D) near the injector,
%     channeling the flood on the injector side.
%   - Water injector (west), oil producer (east, across the fault).
%
% The same model is run twice: with CPEDFM (continuous projection, the
% fault seals) and with plain EDFM (the known leakage). The comparison
% shows the fault-block compartmentalization CPEDFM captures: with EDFM
% the producer receives pressure support and water through the "sealing"
% fault; with CPEDFM the producer block stays isolated.
%
% Runtime: a few minutes. Figures are saved to examples/output/ and,
% when a MATLAB desktop is available, plotToolbar opens for interactive
% inspection.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupCPEDFMU.m'));
outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
gravity off

%% Reservoir: extruded PEBI grid
physdim = [500, 500, 50];
nLayers = 8;
rng(11);
G2D = pebiGrid2D(physdim(1)/14, physdim(1:2));
G = makeLayeredGrid(G2D, nLayers);
G.nodes.coords(:,3) = G.nodes.coords(:,3) * physdim(3)/nLayers;
G = computeGeometry(G);
G.rock = makeRock(G, 100*milli*darcy, 0.25);
fprintf('PEBI reservoir: %d cells\n', G.cells.num);

%% Fractures
fracplanes = struct;

% (1) sealing fault: slanted 15 deg from N-S in plan, dipping in x-z,
% cutting the full cross-section between the wells
fracplanes(1).points = [230 -10 -10;
                        300 510 -10;
                        340 510  60;
                        270 -10  60];
fracplanes(1).aperture = 0.01;
fracplanes(1).poro = 0.4;
fracplanes(1).perm = 1*nano*darcy;

% (2)-(3) conductive natural fractures near the injector (disjoint, so
% no fracture-fracture NNCs are needed)
fracplanes(2).points = [ 60 120  5;
                        190 260  5;
                        190 260 45;
                         60 120 45];
fracplanes(2).aperture = 0.01;
fracplanes(2).poro = 0.5;
fracplanes(2).perm = 10*darcy;

fracplanes(3).points = [ 80 380  5;
                        200 290  5;
                        200 290 45;
                         80 380 45];
fracplanes(3).aperture = 0.01;
fracplanes(3).poro = 0.5;
fracplanes(3).perm = 10*darcy;

%% CPEDFM-U pipeline
[G, fracplanes] = EDFMGridU(G, fracplanes, 'Tolerance', tol, ...
    'fracturelist', 1:numel(fracplanes));
G = fracturematrixNNC_U(G, tol);

Gp = pMatFracNNCsCPEDFM_U(G, tol);
opsCP = setupCPEDFMOpsTPFA_CPG(Gp, Gp.rock, tol);
opsED = setupShaleEDFMOpsTPFA(G, G.rock, tol);
fprintf('projection: %d sealed faces on %d fracture planes\n', ...
    numel(unique(Gp.nnc.pMMfaces)), numel(unique(Gp.nnc.pFracPlane)));

%% Fluid: water-oil
fluid = initSimpleADIFluid('phases', 'WO', ...
    'mu',  [1, 5]*centi*poise, ...
    'rho', [1000, 800], ...
    'n',   [2, 2], ...
    'c',   [1e-8, 1e-6]/barsa, ...
    'pRef', 200*barsa);

%% Wells: injector west of the fault, producer east of it
% Three wells: injector and one producer INSIDE the western fault block
% (a proper intra-block waterflood), plus an observation producer EAST of
% the sealing fault. Under CPEDFM the eastern producer is isolated;
% under EDFM it receives spurious water and pressure support through the
% "sealing" fault.
nMat = G.Matrix.cells.num;
cc = G.cells.centroids(1:nMat,:);
zLayers = unique(round(cc(:,3), 6));
wcInj = zeros(nLayers,1); wcProdW = zeros(nLayers,1); wcProdE = zeros(nLayers,1);
for L = 1:nLayers
    inL = find(abs(cc(:,3) - zLayers(L)) < 1e-6);
    [~, ii] = min((cc(inL,1) -  60).^2 + (cc(inL,2) - 250).^2);
    [~, iw] = min((cc(inL,1) - 185).^2 + (cc(inL,2) - 250).^2);
    [~, ie] = min((cc(inL,1) - 440).^2 + (cc(inL,2) - 250).^2);
    wcInj(L)   = inL(ii);
    wcProdW(L) = inL(iw);
    wcProdE(L) = inL(ie);
end

totTime = 10*year;

W = addWell([], G.Matrix, G.rock, wcInj, 'Type', 'bhp', 'Val', 300*barsa, ...
    'Comp_i', [1 0], 'Name', 'INJ', 'Radius', 0.1, 'Sign', 1);
W = addWell(W, G.Matrix, G.rock, wcProdW, 'Type', 'bhp', 'Val', 100*barsa, ...
    'Comp_i', [1 0], 'Name', 'PROD-W', 'Radius', 0.1, 'Sign', -1);
W = addWell(W, G.Matrix, G.rock, wcProdE, 'Type', 'bhp', 'Val', 100*barsa, ...
    'Comp_i', [1 0], 'Name', 'PROD-E', 'Radius', 0.1, 'Sign', -1);

dt = rampupTimesteps(totTime, year/2, 6);
schedule = simpleSchedule(dt, 'W', W);
state0 = initResSol(G, 200*barsa, [0 1]);

%% Simulate: CPEDFM vs EDFM
modelCP = GenericBlackOilModel(Gp, [], fluid, ...
    'water', true, 'oil', true, 'gas', false);
modelCP.operators = opsCP;
[wsCP, statesCP] = simulateScheduleAD(state0, modelCP, schedule);

modelED = GenericBlackOilModel(G, [], fluid, ...
    'water', true, 'oil', true, 'gas', false);
modelED.operators = opsED;
[wsED, statesED] = simulateScheduleAD(state0, modelED, schedule);

%% Compare producer responses
t = cumsum(schedule.step.val)/year;
dtv = cumsum(schedule.step.val);
qoW_CP = -cellfun(@(w) w(2).qOs, wsCP)*day;   % western (in-block) producer
qoW_ED = -cellfun(@(w) w(2).qOs, wsED)*day;
qoE_CP = -cellfun(@(w) w(3).qOs, wsCP)*day;   % eastern (across-fault) producer
qoE_ED = -cellfun(@(w) w(3).qOs, wsED)*day;

cumW_CP = trapz(dtv, -cellfun(@(w) w(2).qOs, wsCP));
cumW_ED = trapz(dtv, -cellfun(@(w) w(2).qOs, wsED));
cumE_CP = trapz(dtv, -cellfun(@(w) w(3).qOs, wsCP));
cumE_ED = trapz(dtv, -cellfun(@(w) w(3).qOs, wsED));

fprintf('\nCumulative oil after %d years [m^3]:\n', round(totTime/year));
fprintf('  PROD-W (inside fault block) : CPEDFM-U %.3g | EDFM %.3g\n', cumW_CP, cumW_ED);
fprintf('  PROD-E (across the fault)   : CPEDFM-U %.3g | EDFM %.3g\n', cumE_CP, cumE_ED);
fprintf('  EDFM overestimates across-fault recovery by a factor of %.1f\n', ...
    cumE_ED/max(cumE_CP, eps));

pPrdE_CP = mean(statesCP{end}.pressure(wcProdE))/barsa;
pPrdE_ED = mean(statesED{end}.pressure(wcProdE))/barsa;
fprintf('Eastern block end pressure [bar]: CPEDFM-U %.0f (isolated) | EDFM %.0f (supported)\n', ...
    pPrdE_CP, pPrdE_ED);

%% Figures
% fracture system + sealed faces
fig = figure('Visible', 'off', 'Position', [50 50 1200 500]);
subplot(1,2,1); hold on
plotGrid(G.Matrix, 'FaceColor', 'none', 'EdgeAlpha', 0.04);
cols = {'r', 'c', 'c'};
fn = fieldnames(G.FracGrid);
for q = 1:numel(fn)
    Fg = G.FracGrid.(fn{q});
    plotGrid(G, Fg.cells.start:(Fg.cells.start+Fg.cells.num-1), ...
        'FaceColor', cols{q}, 'EdgeAlpha', 0.2);
end
plotGrid(G.Matrix, [wcInj; wcProdW; wcProdE], 'FaceColor', 'k');
view(-30, 40); axis equal tight off
title('sealing fault (red), conductive fractures (cyan), wells (black)');
subplot(1,2,2); hold on
plotGrid(G.Matrix, 'FaceColor', 'none', 'EdgeAlpha', 0.04);
plotFaces(G.Matrix, unique(Gp.nnc.pMMfaces), 'FaceColor', 'b', 'EdgeAlpha', 0.2);
view(-30, 40); axis equal tight off
title('continuous projection: sealed matrix faces');
print(fig, fullfile(outDir, 'example_fracsystem.png'), '-dpng', '-r150');
close(fig);

% water saturation, top layer, end of flood
topCells = find(abs(cc(:,3) - zLayers(1)) < 1e-6);
fig = figure('Visible', 'off', 'Position', [50 50 1200 500]);
subplot(1,2,1);
plotCellData(G.Matrix, statesCP{end}.s(1:nMat,1), topCells, 'EdgeAlpha', 0.05);
view(0,90); axis equal tight off; caxis([0 1]); colorbar
title(sprintf('CPEDFM-U: water saturation after %d years', round(totTime/year)));
subplot(1,2,2);
plotCellData(G.Matrix, statesED{end}.s(1:nMat,1), topCells, 'EdgeAlpha', 0.05);
view(0,90); axis equal tight off; caxis([0 1]); colorbar
title('EDFM: water leaks through the sealing fault');
print(fig, fullfile(outDir, 'example_saturation.png'), '-dpng', '-r150');
close(fig);

% producer oil rates
fig = figure('Visible', 'off', 'Position', [50 50 1100 430]);
subplot(1,2,1); hold on
plot(t, qoW_CP, 'b-', 'LineWidth', 1.6);
plot(t, qoW_ED, 'r--', 'LineWidth', 1.6);
xlabel('time [years]'); ylabel('oil rate [m^3/day]');
legend('CPEDFM-U', 'EDFM', 'Location', 'northeast'); grid on
title('PROD-W (inside the fault block)');
subplot(1,2,2); hold on
plot(t, qoE_CP, 'b-', 'LineWidth', 1.6);
plot(t, qoE_ED, 'r--', 'LineWidth', 1.6);
xlabel('time [years]'); ylabel('oil rate [m^3/day]');
legend('CPEDFM-U', 'EDFM', 'Location', 'northeast'); grid on
title('PROD-E (across the sealing fault)');
print(fig, fullfile(outDir, 'example_wellrates.png'), '-dpng', '-r150');
close(fig);
fprintf('figures saved to examples/output/\n');

if usejava('desktop')
    figure; plotToolbar(G, statesCP); view(-30, 40); colormap jet
    title('CPEDFM-U waterflood (interactive)');
end
