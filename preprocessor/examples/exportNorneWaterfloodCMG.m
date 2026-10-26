%% exportNorneWaterfloodCMG
% End-to-end demonstration of the MRST -> CMG preprocessor on a real field
% case: a sector of the Norne corner-point model with a water injector and a
% producer.
%
%   1. build the case in MRST (faulted corner-point grid, real perm/poro)
%   2. solve it in MRST to have a reference
%   3. export the grid and the wells with writeCMG
%   4. the resulting .dat runs in CMG IMEX
%
% The sector [10 25; 35 55; 1 22] is the same cut the hand-written
% NorneOut2CMGOMO.m script produced the committed CMG include files from,
% so the exported grid is directly comparable with them.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
mrstModule add deckformat incomp

outDir = fullfile(thisDir, '..', 'tests', 'output_cmg');
if ~isfolder(outDir), mkdir(outDir); end

%% ---- 1. Build the MRST case ------------------------------------------
fprintf('--- building Norne sector in MRST ---\n');
grdecl = readGRDECL(fullfile(getDatasetPath('norne'), 'NORNE.GRDECL'));
grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));
grdecl = cutGrdecl(grdecl, [10 25; 35 55; 1 22]);

G = processGRDECL(grdecl);
G = computeGeometry(G(1));
rock = grdecl2Rock(grdecl, G.cells.indexMap);
rock.perm = max(rock.perm, 1e-4*milli*darcy);
rock.poro = max(rock.poro, 0.01);
G.rock = rock;

fprintf('grid   : %d active cells of %d (cartDims %s)\n', ...
    G.cells.num, prod(G.cartDims), mat2str(G.cartDims));
fprintf('perm   : %.3g .. %.3g mD\n', ...
    convertTo(min(rock.perm(:,1)), milli*darcy), ...
    convertTo(max(rock.perm(:,1)), milli*darcy));
fprintf('poro   : %.3f .. %.3f\n', min(rock.poro), max(rock.poro));
fprintf('pv     : %.4g m3\n', sum(poreVolume(G, rock)));

%% ---- wells: injector west, producer east ------------------------------
[ii, jj, ~] = ind2sub(G.cartDims, G.cells.indexMap);
injCol  = find(ii == 4  & jj == 6);
prodCol = find(ii == 13 & jj == 16);
assert(~isempty(injCol) && ~isempty(prodCol), 'Well columns are inactive.');

W = addWell([], G, rock, injCol, 'Type', 'rate', ...
    'Val', 2000*meter^3/day, 'Radius', 0.1, 'Name', 'INJ', ...
    'Comp_i', [1 0], 'sign', 1);
W = addWell(W, G, rock, prodCol, 'Type', 'bhp', ...
    'Val', 250*barsa, 'Radius', 0.1, 'Name', 'PROD', ...
    'Comp_i', [0 1], 'sign', -1);

fprintf('wells  : INJ %d perfs (WI %.3g..%.3g), PROD %d perfs\n', ...
    numel(W(1).cells), min(W(1).WI), max(W(1).WI), numel(W(2).cells));

%% ---- 2. MRST reference solve -----------------------------------------
fprintf('\n--- MRST reference solve (single-phase incompressible) ---\n');
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1000*kilogram/meter^3);
hT = computeTrans(G, rock);
state = initResSol(G, 300*barsa, 1);
state = incompTPFA(state, G, hT, fluid, 'wells', W);

pMRST = convertTo(state.pressure, barsa);
fprintf('pressure : %.2f .. %.2f bar (mean %.2f)\n', ...
    min(pMRST), max(pMRST), mean(pMRST));
fprintf('INJ  bhp : %.2f bar\n', convertTo(state.wellSol(1).pressure, barsa));
fprintf('PROD rate: %.2f m3/day\n', ...
    convertTo(sum(state.wellSol(2).flux), meter^3/day));

%% ---- 3. Export to CMG -------------------------------------------------
fprintf('\n--- exporting to CMG ---\n');
datFile = fullfile(outDir, 'norne_waterflood.dat');
info = writeCMG(G, datFile, ...
    'simulator',  'IMEX', ...
    'cmgVersion', '2024.30', ...
    'grdecl',     grdecl, ...
    'wells',      W, ...
    'timesteps',  repmat(30*day, 12, 1), ...
    'title',      'NORNE_SECTOR_WATERFLOOD', ...
    'hostDeck',   'runnable');   % placeholder fluid model, for the IMEX check

fprintf('grid family   : %s\n', info.gridFamily);
fprintf('fracture model: %s\n', info.fractureModel);
fprintf('blocks        : %d (%d inactive)\n', info.counts.blocks, info.counts.inactive);
fprintf('wells         : %d\n', info.nWells);
fprintf('dataset       : %s\n', info.dataset);
fprintf('include       : %s\n', info.include);

% Figures to reconcile against the CMG run, printed at full precision so the
% comparison is exact rather than eyeballed from 4 significant figures.
mdl = buildCMGModel(G, rock, 'grdecl', grdecl, 'wells', W);
fprintf('\n--- reconciliation targets for the CMG run ---\n');
fprintf('active blocks     : %d\n', G.cells.num);
fprintf('inactive blocks   : %d\n', info.counts.inactive);
fprintf('pore volume       : %.4f m3  ( = %.3f M rm3 in CMG''s report)\n', ...
    sum(poreVolume(G, rock)), sum(poreVolume(G, rock))/1e3);
fprintf('datum (grid top)  : %.3f m\n', mdl.grid.ztop);
fprintf('grid base         : %.3f m\n', mdl.grid.zbot);
fprintf('injection rate    : %.4f m3/day\n', convertTo(W(1).val, meter^3/day));
fprintf('total well layers : %d\n', numel(W(1).cells) + numel(W(2).cells));

%% ---- 4. Report for the CMG-side comparison ---------------------------
save(fullfile(outDir, 'norne_reference.mat'), 'pMRST', 'G', 'rock', 'W');
fprintf('\nMRST reference saved for comparison against the CMG run.\n');
fprintf('=== norne export complete ===\n');
