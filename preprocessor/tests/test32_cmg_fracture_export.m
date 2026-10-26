%% test32_cmg_fracture_export
% Phase 3/5 acceptance test: export EDFM, pEDFM and CPEDFM models to CMG
% from one shared embedded-fracture grid, and check that the deck conserves
% pore volume and carries one connection per NNC.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
cpgDir = fullfile(thisDir, '..', '..', 'CPEDFM-CPG');
if isfolder(cpgDir)
    addpath(cpgDir);
    if exist('startupCPEDFMCPG', 'file') == 2, startupCPEDFMCPG(); end
end
% hfm_v2 is the canonical EDFM module for this workspace. It is not merely a
% newer copy of hfm: calcdavg, which sets the fracture-matrix connectivity
% index, subdivides the matrix cell differently (10x10x3 perturbed vs
% 10x10x10 uniform), so the two modules yield different NNC
% transmissibilities. Mixing them across tests would compare unlike things.
mrstModule add deckformat hfm_v2

outDir = fullfile(thisDir, 'output_cmg');
if ~isfolder(outDir), mkdir(outDir); end

failures = {};

%% ---- shared embedded-fracture grid ------------------------------------
tol = 1e-5;
dims = [12, 10, 2];
physdim = [120, 100, 20];

G0 = computeGeometry(cartGrid(dims, physdim));
G0.rock = makeRock(G0, 100*milli*darcy, 0.2);

frac = struct('points', [35 0 0; 85 100 0; 85 100 20; 35 0 20], ...
    'aperture', 5*milli*meter, 'poro', 0.4, 'perm', 1*nano*darcy);

[Gbase, frac] = EDFMshalegrid(G0, frac, 'Tolerance', tol, ...
    'plotgrid', false, 'fracturelist', 1); %#ok<ASGLU>
Gbase = fracturematrixShaleNNC3D(Gbase, tol);

cases = struct('name', {'EDFM', 'pEDFM', 'CPEDFM'}, 'G', {[], [], []});
cases(1).G = prepareEDFMExport(Gbase);
cases(2).G = preparePEDFMExport(pMatFracNNCs3D(Gbase, tol));
cases(3).G = pMatFracNNCsCPEDFM_CPG(Gbase, tol);

%% ---- export and check each --------------------------------------------
for q = 1:numel(cases)
    name = cases(q).name;
    Gq = cases(q).G;

    datFile = fullfile(outDir, sprintf('frac_%s.dat', lower(name)));
    info = writeCMG(Gq, datFile, 'rock', Gq.rock, 'simulator', 'IMEX', ...
        'cmgVersion', '2024.30', 'datum', 2000, 'tol', tol, ...
        'hostDeck', 'runnable', 'title', ['FRAC_' name]);

    model = buildCMGModel(Gq, Gq.rock, 'datum', 2000, 'tol', tol);
    c = model.meta.counts;

    % --- pore volume: *PORVOL carries it directly ------------------------
    pvCMG = sum(model.props.porvol(model.props.active));
    pvMRST = sum(poreVolume(Gq, Gq.rock));
    relPV = abs(pvCMG - pvMRST)/pvMRST;

    % --- connections: one per NNC, plus any irregular-pair seals ---------
    nNNC = size(Gq.nnc.cells, 1);
    nSealConn = c.connections - nNNC;

    okModel = strcmp(info.fractureModel, name);
    okCells = c.fractureCells == (Gq.cells.num - Gq.Matrix.cells.num);
    okPV = relPV < 1e-10;
    okConn = c.connections >= nNNC;
    ok = okModel && okCells && okPV && okConn;

    fprintf(['%s %-7s blocks=%5d frac=%4d pad=%d NNC=%4d conn=%4d ', ...
             'seals=%3d(+%d irr) pvRel=%.2e\n'], ...
        subsref({'FAIL','PASS'}, substruct('{}', {ok+1})), name, ...
        c.blocks, c.fractureCells, c.paddingLayers, nNNC, c.connections, ...
        c.multipliers, nSealConn, relPV);

    if ~ok
        failures{end+1} = name; %#ok<AGROW>
        fprintf('     model=%d cells=%d pv=%d conn=%d\n', ...
            okModel, okCells, okPV, okConn);
    end
end

%% ---- model-specific expectations --------------------------------------
% EDFM must produce no seals at all; the projected methods must produce some.
mEdfm  = buildCMGModel(cases(1).G, cases(1).G.rock, 'tol', tol);
mPedfm = buildCMGModel(cases(2).G, cases(2).G.rock, 'tol', tol);
mCpedfm= buildCMGModel(cases(3).G, cases(3).G.rock, 'tol', tol);

okSeals = mEdfm.meta.counts.multipliers == 0 && ...
    mPedfm.meta.counts.multipliers > 0 && ...
    mCpedfm.meta.counts.multipliers > 0;
fprintf('%s seals: EDFM %d, pEDFM %d, CPEDFM %d\n', ...
    subsref({'FAIL','PASS'}, substruct('{}', {okSeals+1})), ...
    mEdfm.meta.counts.multipliers, mPedfm.meta.counts.multipliers, ...
    mCpedfm.meta.counts.multipliers);
if ~okSeals, failures{end+1} = 'seals'; end

okGrow = mEdfm.meta.counts.connections < mPedfm.meta.counts.connections && ...
    mEdfm.meta.counts.connections < mCpedfm.meta.counts.connections;
fprintf('%s connections grow with projection: %d < %d, %d\n', ...
    subsref({'FAIL','PASS'}, substruct('{}', {okGrow+1})), ...
    mEdfm.meta.counts.connections, mPedfm.meta.counts.connections, ...
    mCpedfm.meta.counts.connections);
if ~okGrow, failures{end+1} = 'growth'; end

%% ---- summary -----------------------------------------------------------
fprintf('\n');
if isempty(failures)
    fprintf('=== test32: ALL CHECKS PASSED ===\n');
else
    fprintf('=== test32: %d FAILURES: %s ===\n', numel(failures), ...
        strjoin(failures, ', '));
end
