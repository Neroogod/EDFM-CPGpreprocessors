%% test52_fractured_full_deck
% The combination that matters: an embedded-fracture MRST case exported as a
% COMPLETE deck - grid, fracture cells, NNCs, PVT, rock-fluid, initial state,
% wells and schedule - and run in IMEX.
%
% modules/preprocessor already exports the fracture reservoir description and
% test20/test32 verify it. What is new here is that the deck runs: the
% fracture representation and the fluid model have to be consistent for IMEX
% to accept and advance the case.
%
% Checks:
%   A. The exported block grid carries the fracture padding, and the initial
%      state is scattered onto the addresses buildCMGModel actually used
%      (deckBlockArray), not appended to the end of the matrix array.
%   B. Fracture cell pressures land on their padded blocks, verified against
%      grid.fracIJK directly.
%   C. The dataset is written and IMEX runs it to Normal Termination.
%   D. IMEX's pore volume matches the MRST model's, so the fracture pore
%      volume survived the export.
%   E. IMEX's cumulative injection matches the schedule.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
mrstModule add preprocessor hfm ad-core ad-props ad-blackoil

cart = fullfile(fileparts(fileparts(thisDir)), 'CPEDFM', 'UpdateMay27_2023');
addpath(fullfile(cart, 'hfmFns'));
addpath(fullfile(cart, 'pEDFM_nnc'));
addpath(fullfile(cart, 'utils'));
addpath(genpath(fullfile(cart, 'ADFNE15')));
Globals;

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
failures = {};

%% ---- Embedded-fracture MRST case ---------------------------------------
gravity reset on
G = cartGrid([8, 8, 2], [80, 80, 20]);
G.nodes.coords(:,3) = G.nodes.coords(:,3) + 1800;
G = computeGeometry(G);
G.rock = makeRock(G, 100*milli*darcy, 0.20);

fp = struct('points', [20  5 1800; 60 75 1800; 60 75 1820; 20  5 1820], ...
            'aperture', 5*milli*meter, 'poro', 0.4, 'perm', 5000*milli*darcy);
[G, fp] = EDFMshalegrid(G, fp, 'Tolerance', tol, 'plotgrid', false, 'fracturelist', 1);
G = fracturematrixShaleNNC3D(G, tol);

nMat = G.Matrix.cells.num;
nFrac = G.cells.num - nMat;
fprintf('EDFM grid: %d cells (%d matrix + %d fracture), %d NNCs\n', ...
    G.cells.num, nMat, nFrac, size(G.nnc.cells, 1));
assert(nFrac > 0, 'no fracture cells were created');

fluid = initSimpleADIFluid('phases', 'WO', ...
    'mu',  [1, 4]*centi*poise, 'rho', [1000, 820], ...
    'c',   [4.4e-5, 9e-5]/barsa, 'n', [2, 2], 'pRef', 220*barsa);

model = GenericBlackOilModel(G, G.rock, fluid, ...
    'water', true, 'oil', true, 'gas', false);
model.operators = setupShaleEDFMOpsTPFA(G, G.rock, tol);

state0 = initResSol(G, 220*barsa, [0.2, 0.8]);
% Give the fracture cells a distinguishable pressure so check B can tell
% whether they landed on the right blocks.
pFrac = 233*barsa;
state0.pressure(nMat+1:end) = pFrac;

cInj  = find(G.cells.centroids(1:nMat,1) < 12);
cProd = find(G.cells.centroids(1:nMat,1) > 68);
rateInj = 40*meter^3/day;
W = addWell([], G.Matrix, G.rock, cInj,  'Type','rate','Val', rateInj, ...
    'Comp_i',[1 0], 'Name','INJ',  'refDepth', 1800);
W = addWell(W, G.Matrix, G.rock, cProd, 'Type','bhp','Val', 180*barsa, ...
    'Comp_i',[0 1], 'Name','PROD', 'refDepth', 1800);
totTime = 200*day;
schedule = simpleSchedule(rampupTimesteps(totTime, 50*day, 2), 'W', W);

%% ---- Build the deck ----------------------------------------------------
deck = buildFullDeck(model, state0, schedule);
res = deck.res;
fprintf('exported: %s grid %s, %d blocks (%d active), fracture model %s\n', ...
    deck.meta.gridFamily, mat2str(res.grid.dims), res.props.nTot, ...
    sum(res.props.active), deck.meta.fractureModel);

%% ---- Check A: padding present and mapping available --------------------
okA = res.props.nTot > G.cells.num && isfield(res.grid, 'fracIJK') && ...
      size(res.grid.fracIJK, 1) == nFrac && ...
      ~strcmpi(deck.meta.fractureModel, 'matrix');
fprintf('A. padding: %d blocks for %d cells, fracIJK rows %d (fracture cells %d)\n', ...
    res.props.nTot, G.cells.num, size(res.grid.fracIJK,1), nFrac);
if okA, fprintf('   PASS\n'); else
    failures{end+1} = 'A: fracture padding or fracIJK mapping missing'; fprintf('   FAIL\n');
end

%% ---- Check B: initial state lands on the right blocks ------------------
pBlocks = deckBlockArray(res, state0.pressure, deck.init.pDatum);
linF = sub2ind(reshape(double(res.grid.dims),1,[]), ...
    res.grid.fracIJK(:,1), res.grid.fracIJK(:,2), res.grid.fracIJK(:,3));
linM = sub2ind(reshape(double(res.grid.dims),1,[]), ...
    res.grid.matIJK(:,1), res.grid.matIJK(:,2), res.grid.matIJK(:,3));
dFrac = max(abs(pBlocks(linF) - pFrac));
dMat  = max(abs(pBlocks(linM) - state0.pressure(1:nMat)));
% How much does the mapped scatter differ from simply appending the fracture
% values to the end of the matrix array? Reported, not asserted: for a single
% fracture whose cells all fall in the first padding layer the two can
% coincide, and the point of the check is that the scatter is right, not that
% the naive form is always wrong.
naive = [state0.pressure; repmat(deck.init.pDatum, res.props.nTot-G.cells.num, 1)];
nDiff = sum(abs(naive - pBlocks) > 1e-6);
fprintf('B. fracture-block pressure err %.3e, matrix-block err %.3e\n', dFrac, dMat);
fprintf('   naive append differs from the mapped scatter on %d of %d blocks\n', ...
    nDiff, res.props.nTot);
okB = dFrac < 1e-9 && dMat < 1e-9;
if okB, fprintf('   PASS (state scattered by fracIJK/matIJK)\n'); else
    failures{end+1} = 'B: initial state did not land on the mapped blocks'; fprintf('   FAIL\n');
end

%% ---- Check C: write and run --------------------------------------------
datFile = fullfile(outDir, 'edfm_full.dat');
info = writeCMGDeck(deck, datFile, 'title', 'MRST_EDFM_FULL');
fprintf('C. dataset: %s (%.1f kB)\n', info.dataset, localBytes(info.dataset)/1024);

imex = localFindIMEX();
okC = false; logTxt = '';
if isempty(imex)
    fprintf('   SKIP: IMEX not found on this machine.\n');
else
    outFile = strrep(info.dataset, '.dat', '.out');
    if isfile(outFile), delete(outFile); end
    [st, ~] = system(sprintf('cd /d "%s" && "%s" -f "%s" -log', ...
        outDir, imex, info.dataset));
    if isfile(outFile), logTxt = fileread(outFile); end
    nErr = numel(regexp(logTxt, 'ERROR \(from subroutine', 'start'));
    normalEnd = contains(logTxt, 'Normal Termination');
    fprintf('   IMEX exit=%d errors=%d NormalTermination=%d\n', st, nErr, normalEnd);
    okC = normalEnd && nErr == 0;
    if okC
        fprintf('   PASS (fractured deck runs as it stands)\n');
    else
        failures{end+1} = sprintf('C: IMEX errors=%d normalEnd=%d', nErr, normalEnd);
        fprintf('   FAIL\n');
        e = regexp(logTxt, 'ERROR \(from subroutine[^\n]*\n[^\n]*\n[^\n]*', 'match');
        for q = 1:min(3, numel(e)), fprintf('     %s\n', e{q}); end
    end
end

%% ---- Checks D, E: pore volume and injection ----------------------------
if isempty(logTxt)
    fprintf('D-E. SKIP (no IMEX log).\n');
else
    pvMRST = sum(poreVolume(G, G.rock));
    % Anchor on the unit string: a looser pattern captures the '3' out of
    % 'M rm3' instead of the value.  IMEX prints 'M rm3' = 1e3 reservoir m3.
    t = regexp(logTxt, 'Total Pore Volume\.?\s+M rm3\s+([\d.eE+-]+)', 'tokens', 'once');
    pvIMEX = NaN;
    if ~isempty(t), pvIMEX = str2double(t{1}) * 1000; end   % M rm3 -> m3
    relPV = abs(pvIMEX - pvMRST)/max(pvMRST, eps);
    fprintf('D. pore volume IMEX %.3f m3 vs MRST %.3f m3 (rel %.2e)\n', ...
        pvIMEX, pvMRST, relPV);
    okD = isfinite(relPV) && relPV < 1e-3;
    if okD, fprintf('   PASS (fracture pore volume survived the export)\n'); else
        failures{end+1} = sprintf('D: pore volume rel diff %.3e', relPV);
        fprintf('   FAIL\n');
    end

    cumInj = localGrabInj(logTxt);
    expInj = convertTo(rateInj, meter^3/day) * convertTo(totTime, day) / 1000;
    relInj = abs(cumInj - expInj)/max(expInj, eps);
    fprintf('E. cum water injection IMEX %.4f vs schedule %.4f (rel %.2e)\n', ...
        cumInj, expInj, relInj);
    okE = relInj < 1e-6;
    if okE, fprintf('   PASS\n'); else
        failures{end+1} = sprintf('E: injection rel diff %.3e', relInj);
        fprintf('   FAIL\n');
    end
end

%% ---- Summary ------------------------------------------------------------
if isempty(failures)
    fprintf('\n=== test52: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test52: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test52 failed');
end

function n = localBytes(f), d = dir(f); n = d.bytes; end

function exe = localFindIMEX()
    exe = '';
    root = getenv('CMG_ROOT');
    if isempty(root), root = 'C:\Program Files (x86)\CMG'; end
    c = dir(fullfile(root, 'IMEX', '*', 'Win_x64', 'EXE', 'mx*.exe'));
    if ~isempty(c)
        [~, ix] = max([c.datenum]);
        exe = fullfile(c(ix).folder, c(ix).name);
    end
end

function v = localGrabInj(txt)
    v = NaN;
    t = regexp(txt, 'Cumulative Injection\s+(NA|[\d.eE+-]+)\s+(NA|[\d.eE+-]+)\s+(NA|[\d.eE+-]+)', ...
        'tokens', 'once');
    if ~isempty(t), v = str2double(t{3}); end
end
