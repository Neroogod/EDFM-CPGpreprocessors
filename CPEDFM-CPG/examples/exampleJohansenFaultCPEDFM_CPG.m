%% exampleJohansenFaultCPEDFM_CPG
% CO2 injection in the Johansen formation with sealing faults represented
% by CPEDFM on the corner-point grid.
%
% This example mirrors the case in
%   CPEDFMcornerpointaug24/PaperReviewCases/Case3/JohansenTEDFMNoFaultCase.m
% but replaces the MANUAL fault handling (hand-picked south-face lists,
% `model.operators.T_all(...) = 0`) with the automatic continuous
% projection: the same four fault planes are handed to
% pMatFracNNCsCPEDFM_CPG, which finds and seals the fault-crossing faces
% by itself.
%
% NOTE: needs the co2lab module and the Johansen dataset (downloaded by
% MRST on first use). Full schedule (100 y injection + 1000 y migration)
% takes tens of minutes; set quickTest = false for the full run.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupCPEDFMCPG.m'));
mrstModule add co2lab;

quickTest = true;   % true: short schedule to exercise the pipeline
tol = 1e-5;

%% Johansen grid and rock (as in the TEDFM case)
[G, rock, bcIx] = makeJohansenVEgridFull();

topLayer_global = false(G.cartDims);
topLayer_global(1:end, 1:end, 1) = true;
topLayer = find(topLayer_global(G.cells.indexMap));
rock.perm(topLayer, 1:2) = rock.perm(topLayer, 1:2)*1e4;

G.rock = rock;

%% Fault planes (identical geometry to the TEDFM script)
[i, j, k] = ind2sub(G.cartDims, G.cells.indexMap);

frac = struct; nf = 0;
mkplane = @(sub, off1, off2, off3, off4) deal(sub);

% Fault 1 part 1
sub = extractSubgrid(G, i==43 & j>0 & j<22 & k==6);
fc = G.cells.centroids(sub.cells.global,:);
pt1 = fc(1,:); pt2 = fc(end,:);
nf = nf+1;
frac(nf).points = [pt1-[-500,0,2000]; pt2-[-500,0,2000]; ...
                   pt2+[ 500,0,2000]; pt1+[ 500,0,2000]];

% Fault 1 part 2
sub = extractSubgrid(G, i==43 & j>19 & j<48 & k==6);
fc = G.cells.centroids(sub.cells.global,:);
pt1 = fc(1,:); pt2 = fc(end,:);
nf = nf+1;
frac(nf).points = [pt1-[-500,0,2000]; pt2-[ 400,0,2000]; ...
                   pt2+[-400,0,2000]; pt1+[ 500,0,2000]];

% Fault 2 part 1
sub = extractSubgrid(G, i==45 & j>0 & j<22 & k==6);
fc = G.cells.centroids(sub.cells.global,:);
pt1 = fc(1,:); pt2 = fc(end,:);
nf = nf+1;
frac(nf).points = [pt1-[-500,0,2000]; pt2-[-500,0,2000]; ...
                   pt2+[ 500,0,2000]; pt1+[ 500,0,2000]];

% Fault 2 part 2
sub = extractSubgrid(G, i==45 & j>19 & j<48 & k==6);
fc = G.cells.centroids(sub.cells.global,:);
pt1 = fc(1,:); pt2 = fc(end,:);
nf = nf+1;
frac(nf).points = [pt1-[-500,0,2000]; pt2-[ 600,0,2000]; ...
                   pt2+[-600,0,2000]; pt1+[ 500,0,2000]];

for q = 1:nf
    frac(q).aperture = 1;
    frac(q).poro = 0.5;
    % low-permeability (sealing) fault filling; the TEDFM script instead
    % zeroed cross-fault transmissibilities by hand
    frac(q).perm = [1*micro*darcy, 1*micro*darcy, 1*micro*darcy];
end

%% CPEDFM pipeline
[G, frac] = EDFMGridCPG(G, frac, 'Tolerance', tol, 'plotgrid', false, ...
    'fracturelist', 1:nf);
G = fracturematrixShaleNNC3D(G, tol);
G = pMatFracNNCsCPEDFM_CPG(G, tol);
TPFAoperators = setupCPEDFMOpsTPFA_CPG(G, G.rock, tol);

% Johansen data quirk: pinched cells with zero pore volume make the
% EDFM NNC transmissibilities of fracturematrixShaleNNC3D evaluate to
% 0/0 = NaN. The original TEDFM script applies the same guard
% (JohansenTEDFMNoFaultCase.m line 354).
nBad = nnz(~isfinite(TPFAoperators.T));
if nBad > 0
    fprintf('Sanitizing %d non-finite NNC transmissibilities (zero-PV cells).\n', nBad);
    TPFAoperators.T(~isfinite(TPFAoperators.T)) = 1e-10;
    TPFAoperators.T_all(~isfinite(TPFAoperators.T_all)) = 1e-10;
end

%% Fluid (same as the TEDFM script)
gravity on;
g = gravity;
rhow = 1000;
initState.pressure = rhow * g(3) * G.cells.centroids(:,3);
initState.s = repmat([1, 0], G.cells.num, 1);
initState.sGmax = initState.s(:,2);

co2     = CO2props();
p_ref   = 30*mega*Pascal;
t_ref   = 94 + 273.15;
rhoc    = co2.rho(p_ref, t_ref);
cf_co2  = co2.rhoDP(p_ref, t_ref)/rhoc;
cf_wat  = 0;
cf_rock = 4.35e-5/barsa;
muw     = 8e-4*Pascal*second;
muco2   = co2.mu(p_ref, t_ref)*Pascal*second;

fluid = initSimpleADIFluid('phases', 'WG', ...
    'mu', [muw, muco2], 'rho', [rhow, rhoc], 'pRef', p_ref, ...
    'c', [cf_wat, cf_co2], 'cR', cf_rock, 'n', [2 2]);

p_sc = 1*atm; T_sc = 288.706;
fluid.rhoSinj = [rhow, co2.rho(p_sc, T_sc)];

srw = 0.27; src = 0.01;
fluid.krW = @(s) fluid.krW(max((s-srw)./(1-srw), 0));
fluid.krG = @(s) fluid.krG(max((s-src)./(1-src), 0));
pe = 5*kilo*Pascal;
pcWG = @(sw) pe*sw.^(-1/2);
fluid.pcWG = @(sg) pcWG(max((1-sg-srw)./(1-srw), 1e-5));

%% Well and schedule
wc_global = false(G.cartDims);
wc_global(48, 48, 6:10) = true;
wc = find(wc_global(G.cells.indexMap));
inj_rate = 3.5*mega*1e3/year/fluid.rhoSinj(end);

W = addWell([], G.Matrix, rock, wc, 'type', 'rate', 'val', inj_rate, ...
    'comp_i', [0 1], 'Sign', 1);

if quickTest
    stepsInj = repmat(year, 5, 1);
    stepsMig = repmat(10*year, 3, 1);
else
    stepsInj = repmat(year, 100, 1);
    stepsMig = repmat(10*year, 100, 1);
end
[W_inj, W_shut] = deal(W);
W_shut(1).status = false;
schedule.control = [struct('W', W_inj); struct('W', W_shut)];
schedule.step.val     = [stepsInj; stepsMig];
schedule.step.control = [ones(numel(stepsInj),1); 2*ones(numel(stepsMig),1)];

%% Model and simulation
model = TwoPhaseWaterGasModel(G, [], fluid, 0, 0);
model.operators = TPFAoperators;

[wellSol, states] = simulateScheduleAD(initState, model, schedule);

%% Plots
if usejava('desktop')
    % interactive inspector (requires the MATLAB desktop)
    figure;
    plotToolbar(G, states); view(-63, 68); colormap jet;
    title('CPEDFM-CPG Johansen fault case');
end

sat_end = states{end}.s(:,2);
figure;
plotGrid(G, 'facecolor', 'none', 'edgeAlpha', 0.1);
plotGrid(G, sat_end > 1e-4, 'facecolor', 'red');
axis off; view(-75, 30);
title('CO2 plume extent (CPEDFM-sealed faults)');
