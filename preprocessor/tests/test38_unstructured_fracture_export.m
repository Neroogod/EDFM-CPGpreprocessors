%% test38_unstructured_fracture_export
% Unstructured grid carrying fractures, exported as a connection-list model.
%
% CMG has no unstructured grid type, so the cells become an Nc x 1 x 1
% pseudo-grid whose entire connectivity is explicit. That representation
% turns out to suit fractures better than the structured one: there are no
% padding layers, because every cell already has its own block address, and
% projection seals fold into the connection transmissibilities instead of
% needing a separate multiplier keyword.
%
% What must hold:
%   * every cell of the assembled grid becomes a block
%   * pore volume transfers exactly, matrix and fractures together
%   * the connection list holds every internal matrix face plus every NNC
%   * sealing removes connections rather than over-assigning area

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
% CPEDFM-U supplies EDFMGridU, the unstructured gridder. EDFMshalegrid and
% EDFMGridCPG both assume 8-corner cells and fail on a PEBI grid inside
% checkplaneAABBintersect; EDFMGridU uses half-space clipping, which is exact
% for any convex cell.
uDir = fullfile(thisDir, '..', '..', 'CPEDFM-U');
if isfolder(uDir)
    addpath(uDir);
    if exist('startupCPEDFMU', 'file') == 2, evalc('startupCPEDFMU()'); end
end
% After the startup, which adds hfm and would otherwise shadow hfm_v2.
mrstModule add hfm_v2 deckformat

outDir = fullfile(thisDir, 'output_cmg');
if ~isfolder(outDir), mkdir(outDir); end

failures = {};
report = @(ok, name, detail) fprintf('%s %-34s %s\n', ...
    subsref({'FAIL', 'PASS'}, substruct('{}', {ok + 1})), name, detail);

%% ---- an unstructured grid -------------------------------------------
% A PEBI footprint extruded into layers: unstructured in map view, which is
% what makes cmgDetectGrid classify it as such (no cartDims).
tol = 1e-5;
rng(7);
pts = 100*rand(60, 2);
G2 = pebi(triangleGrid(pts));
G2.nodes.coords = G2.nodes.coords;             % already in metres
Gm = makeLayeredGrid(G2, [8 8 8]);             % three layers, 24 m thick
Gm = computeGeometry(Gm);
Gm.nodes.coords(:,3) = Gm.nodes.coords(:,3) + 2000;
Gm = computeGeometry(Gm);
Gm.rock = makeRock(Gm, 100*milli*darcy, 0.2);

gi = cmgDetectGrid(Gm);
fprintf('grid: %d cells, family %s\n', Gm.cells.num, gi.family);
assert(strcmp(gi.family, 'unstructured'), 'expected an unstructured grid');

%% ---- put a fracture through it ---------------------------------------
frac = struct('points', [30 10 2000; 70 90 2000; 70 90 2024; 30 10 2024], ...
    'aperture', 5*milli*meter, 'poro', 0.4, 'perm', 1*darcy);

[Gf, frac] = EDFMGridU(Gm, frac, 'Tolerance', tol, ...
    'fracturelist', 1, 'Verbose', false); %#ok<ASGLU>
Gf = fracturematrixShaleNNC3D(Gf, tol);

fi = cmgDetectFractureModel(Gf);
nMat = fi.nMatrix; nFrac = fi.nFrac;
fprintf('after EDFM: %d matrix + %d fracture cells, %d NNCs\n\n', ...
    nMat, nFrac, fi.nNNC);

okDetect = strcmp(fi.model, 'EDFM') && nFrac > 0;
report(okDetect, 'A. detection', sprintf('%s, %d fracture cells', fi.model, nFrac));
if ~okDetect, failures{end+1} = 'A'; end

%% ---- export -----------------------------------------------------------
info = writeCMG(Gf, fullfile(outDir, 'unstr_frac.dat'), ...
    'rock', Gf.rock, 'simulator', 'IMEX', 'cmgVersion', '2024.30', ...
    'tol', tol, 'hostDeck', 'runnable', 'title', 'UNSTR_FRAC');

model = buildCMGModel(Gf, Gf.rock, 'tol', tol);
c = model.meta.counts;

% Every cell becomes a block; no padding.
okBlocks = c.blocks == Gf.cells.num && c.paddingLayers == 0 && c.inactive == 0;
report(okBlocks, 'B. one block per cell', sprintf( ...
    '%d blocks = %d matrix + %d fracture, %d padding layers', ...
    c.blocks, c.matrixCells, c.fractureCells, c.paddingLayers));
if ~okBlocks, failures{end+1} = 'B'; end

% Pore volume, matrix and fractures together.
pvMRST = sum(poreVolume(Gf, Gf.rock));
pvCMG = sum(model.props.porvol);
relPV = abs(pvCMG - pvMRST)/pvMRST;
okPV = relPV < 1e-12;
report(okPV, 'C. pore volume', sprintf('%.4f m3, rel diff %.2e', pvMRST, relPV));
if ~okPV, failures{end+1} = 'C'; end

% Connectivity: internal matrix faces (less any fully sealed) plus the NNCs.
nInternal = nnz(all(Gf.Matrix.faces.neighbors > 0, 2));
nNNC = size(Gf.nnc.cells, 1);
okConn = c.connections == (nInternal - c.sealedFaces) + nNNC;
report(okConn, 'D. connection list complete', sprintf( ...
    '%d = %d matrix faces - %d sealed + %d NNCs', ...
    c.connections, nInternal, c.sealedFaces, nNNC));
if ~okConn, failures{end+1} = 'D'; end

% No connection may carry a non-finite or negative transmissibility.
T = model.conn(:,7);
okT = all(isfinite(T)) && all(T > 0);
report(okT, 'E. transmissibilities sane', sprintf( ...
    '%.3e .. %.3e md.m', min(T), max(T)));
if ~okT, failures{end+1} = 'E'; end

fprintf('\nwrote %s\n', info.dataset);

%% ---- the same model in Eclipse ----------------------------------------
outEcl = fullfile(thisDir, 'output_eclipse');
if ~isfolder(outEcl), mkdir(outEcl); end

iE = writeEclipse(Gf, fullfile(outEcl, 'unstr_frac.DATA'), ...
    'rock', Gf.rock, 'tol', tol, 'title', 'UNSTR_FRAC');
dE = readEclipseDeck(iE.dataset);

mdm = milli*darcy*meter;
ecl = centi*poise*meter^3/(day*barsa);

% Tolerance here is the deck's own precision, not the model's. PORV is
% written as %.8e - nine significant figures - so a value that has been
% through the file cannot agree more closely than about 1e-9 relative, and
% the sum over 216 blocks lands near 1e-10. Check C above compares the model
% in memory and is exact; this one cannot be.
okEclNNC = isfield(dE.GRID, 'NNC') && size(dE.GRID.NNC, 1) == c.connections;
okEclPV  = isfield(dE.GRID, 'PORV') && ...
    abs(sum(dE.GRID.PORV) - pvMRST)/pvMRST < 1e-8;
okEclT = false;
if okEclNNC
    want = model.conn(:,7)*(mdm/ecl);
    okEclT = max(abs(dE.GRID.NNC(:,7) - want)./max(abs(want), eps)) < 1e-8;
end
okEcl = okEclNNC && okEclPV && okEclT;
report(okEcl, 'F. Eclipse connection list', sprintf( ...
    'NNC %d/%d, PORV rel %.1e', size(dE.GRID.NNC, 1), c.connections, ...
    abs(sum(dE.GRID.PORV) - pvMRST)/pvMRST));
if ~okEcl
    failures{end+1} = 'F';
    fprintf('     nnc=%d pv=%d trans=%d\n', okEclNNC, okEclPV, okEclT);
end

%% ---- summary -----------------------------------------------------------
fprintf('\n');
if isempty(failures)
    fprintf('=== test38: ALL CHECKS PASSED ===\n');
else
    fprintf('=== test38: %d FAILURES: %s ===\n', numel(failures), ...
        strjoin(failures, ', '));
end
