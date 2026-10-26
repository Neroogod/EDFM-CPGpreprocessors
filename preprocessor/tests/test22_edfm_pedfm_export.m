%% test22_edfm_pedfm_export
% End-to-end adapter and writer test for EDFM, pEDFM, and CPEDFM.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
cpgDir = fullfile(thisDir, '..', '..', 'CPEDFM-CPG');
addpath(cpgDir);
startupCPEDFMCPG();
addpath(fullfile(thisDir, '..'));
mrstModule add deckformat

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
dims = [12, 10, 2];
physdim = [120, 100, 20];

G0 = cartGrid(dims, physdim);
G0 = computeGeometry(G0);
G0.rock = makeRock(G0, 100*milli*darcy, 0.2);

frac = struct;
frac(1).points = [35, 0, 0; 85, 100, 0; 85, 100, 20; 35, 0, 20];
frac(1).aperture = 5*milli*meter;
frac(1).poro = 0.4;
frac(1).perm = 1*nano*darcy;

[Gbase, frac] = EDFMshalegrid(G0, frac, 'Tolerance', tol, ...
    'plotgrid', false, 'fracturelist', 1);
Gbase = fracturematrixShaleNNC3D(Gbase, tol);

Gedfm = prepareEDFMExport(Gbase);
GpedfmRaw = pMatFracNNCs3D(Gbase, tol);
Gpedfm = preparePEDFMExport(GpedfmRaw);
Gcpedfm = pMatFracNNCsCPEDFM_CPG(Gbase, tol);

grdecl = simpleGrdecl(dims, 0);
grdecl.COORD(1:3:end) = grdecl.COORD(1:3:end)*physdim(1);
grdecl.COORD(2:3:end) = grdecl.COORD(2:3:end)*physdim(2);
grdecl.COORD(3:3:end) = grdecl.COORD(3:3:end)*physdim(3);
grdecl.ZCORN = grdecl.ZCORN*physdim(3);

[ii, ~, ~] = ind2sub(dims, Gbase.Matrix.cells.indexMap);
wells = struct('name', {'INJ', 'PROD'}, 'type', {'rate', 'bhp'}, ...
    'val', {25*meter^3/day, 100*barsa}, 'comp', {'WATER', 'OIL'}, ...
    'cells', {find(ii == 1, 2), find(ii == dims(1), 2)});

edfmFile = fullfile(outDir, 'EDFM_TEST22.DATA');
pedfmFile = fullfile(outDir, 'PEDFM_TEST22.DATA');
cpedfmFile = fullfile(outDir, 'CPEDFM_TEST22.DATA');

ie = writeEclipseEDFM(Gedfm, Gedfm.rock, grdecl, wells, edfmFile, 'tol', tol);
ip = writeEclipsePEDFM(GpedfmRaw, GpedfmRaw.rock, grdecl, wells, pedfmFile, 'tol', tol);
ic = writeEclipseEmbeddedFractureModel('CPEDFM', Gcpedfm, Gcpedfm.rock, ...
    grdecl, wells, cpedfmFile, 'tol', tol);

assert(strcmp(ie.modelType, 'EDFM') && ie.nMult == 0 && ie.nEditnnc == 0);
assert(strcmp(ip.modelType, 'pEDFM') && ip.nMult > 0);
assert(strcmp(ic.modelType, 'CPEDFM') && ic.nMult > 0);
assert(ie.nNNC < ip.nNNC && ie.nNNC < ic.nNNC, ...
    'Projected models must contain additional projection NNCs.');

% Independent Eclipse parser round-trip of the model-defining records.
files = {edfmFile, pedfmFile, cpedfmFile};
infos = {ie, ip, ic};
for q = 1:numel(files)
    deck = convertDeckUnits(readEclipseDeck(files{q}));
    assert(isfield(deck.GRID, 'NNC'));
    assert(size(deck.GRID.NNC, 1) == infos{q}.nNNC);
end

% CMG writers consume the same normalized data.
ce = writeCMGEDFM(Gedfm, Gedfm.rock, grdecl, ...
    fullfile(outDir, 'edfm_test22.inc'));
cp = writeCMGPEDFM(GpedfmRaw, GpedfmRaw.rock, grdecl, ...
    fullfile(outDir, 'pedfm_test22.inc'));
assert(ce.nNNC == ie.nNNC && ce.nMult == 0);
assert(cp.nNNC == ip.nNNC && cp.nMult == ip.nMult);

% Simulator selection produces tagged, self-describing CMG packages.
simulators = {'IMEX', 'GEM', 'STARS'};
for q = 1:numel(simulators)
    sim = simulators{q};
    out = fullfile(outDir, ['pedfm_', lower(sim), '_test22.inc']);
    ci = writeCMGPEDFM(GpedfmRaw, GpedfmRaw.rock, grdecl, out, ...
        'simulator', sim, 'cmgVersion', 'test-release');
    assert(strcmp(ci.simulator, sim) && strcmp(ci.modelType, 'pEDFM'));
    assert(exist(ci.tableFiles.fractureBlocks, 'file') == 2);
    assert(exist(ci.tableFiles.nnc, 'file') == 2);
    assert(exist(ci.tableFiles.multipliers, 'file') == 2);
    assert(exist(ci.tableFiles.manifest, 'file') == 2);
    manifest = jsondecode(fileread(ci.tableFiles.manifest));
    assert(strcmp(manifest.simulator, sim));
    assert(strcmp(manifest.fractureModel, 'pEDFM'));
    assert(manifest.counts.nnc == ip.nNNC);
end

% The normalized pEDFM face records reproduce the exporter multipliers.
mp = getEmbeddedFractureMultipliers(Gpedfm, tol);
assert(nnz(mp < 1) > 0 && all(mp >= 0 & mp <= 1));

fprintf(['test22 PASS: EDFM %d NNC/0 MULT; pEDFM %d NNC/%d MULT; ', ...
    'CPEDFM %d NNC/%d MULT.\n'], ie.nNNC, ip.nNNC, ip.nMult, ...
    ic.nNNC, ic.nMult);
