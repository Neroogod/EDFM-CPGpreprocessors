%% test34_cornerpoint_fracture_cmg
% Fracture export on a corner-point grid: the padding layers are appended to
% ZCORN with corners bracketing each fracture cell's true depth, so depth
% comes from the geometry and the matrix blocks' own corner-derived depths
% are left untouched.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
cpgDir = fullfile(thisDir, '..', '..', 'CPEDFM-CPG');
if isfolder(cpgDir)
    addpath(cpgDir);
    if exist('startupCPEDFMCPG', 'file') == 2, startupCPEDFMCPG(); end
end
% hfm_v2, not hfm: EDFMGridCPG calls convexpolygonarea with two outputs, and
% only the hfm_v2 copy declares them. Loading hfm here fails inside the
% gridder with "Too many output arguments" at pebiAABBintersectCPGv5:648.
mrstModule add deckformat hfm_v2

outDir = fullfile(thisDir, 'output_cmg');
if ~isfolder(outDir), mkdir(outDir); end

failures = {};
tol = 1e-5;

%% ---- a faulted corner-point grid --------------------------------------
% Same configuration test20 exercises on the Eclipse side, translated to a
% realistic depth. Other fracture-plane placements hit a concave-intersection
% branch in the upstream gridder where duplicate copies of
% concavePolygonArea disagree on their output count.
datum = 2000;
dims = [20, 10, 8];
physdim = [200, 100, 40];
grdecl = simpleGrdecl(dims, 0.12);
grdecl.COORD(1:3:end) = grdecl.COORD(1:3:end)*physdim(1);
grdecl.COORD(2:3:end) = grdecl.COORD(2:3:end)*physdim(2);
grdecl.COORD(3:3:end) = grdecl.COORD(3:3:end)*physdim(3) + datum;
grdecl.ZCORN = grdecl.ZCORN*physdim(3) + datum;

G = computeGeometry(processGRDECL(grdecl));
G.rock = makeRock(G, 100*milli*darcy, 0.25);

fprintf('corner-point grid: %d cells of %d\n', G.cells.num, prod(dims));

%% ---- one fracture plane through it ------------------------------------
fracplanes = struct('points', ...
    [70 -5 -5; 110 105 -5; 150 105 45; 110 -5 45] + [0 0 datum], ...
    'aperture', 10*milli*meter, 'poro', 0.4, 'perm', 1*nano*darcy);

[Gf, fracplanes] = EDFMGridCPG(G, fracplanes, 'Tolerance', tol, ...
    'plotgrid', false, 'fracturelist', 1); %#ok<ASGLU>
Gf = fracturematrixShaleNNC3D(Gf, tol);

nMat = Gf.Matrix.cells.num;
nFrac = Gf.cells.num - nMat;
fprintf('after EDFM: %d matrix + %d fracture cells, %d NNCs\n', ...
    nMat, nFrac, size(Gf.nnc.cells, 1));

%% ---- detection ---------------------------------------------------------
gi = cmgDetectGrid(Gf);
fi = cmgDetectFractureModel(Gf);
okDetect = strcmp(gi.family, 'corner') && strcmp(fi.model, 'EDFM');
fprintf('%s detection: family=%s model=%s\n', ...
    subsref({'FAIL','PASS'}, substruct('{}', {okDetect+1})), gi.family, fi.model);
if ~okDetect, failures{end+1} = 'detect'; end

%% ---- export ------------------------------------------------------------
datFile = fullfile(outDir, 'frac_corner.dat');
info = writeCMG(Gf, datFile, 'rock', Gf.rock, 'grdecl', grdecl, ...
    'simulator', 'IMEX', 'cmgVersion', '2024.30', 'tol', tol, ...
    'hostDeck', 'runnable', 'title', 'FRAC_CORNER');

model = buildCMGModel(Gf, Gf.rock, 'grdecl', grdecl, 'tol', tol);
c = model.meta.counts;

%% ---- checks ------------------------------------------------------------
% ZCORN must have grown by exactly the padding layers.
nzTot = model.grid.dims(3);
expectZ = 8*model.grid.dims(1)*model.grid.dims(2)*nzTot;
okZ = numel(model.grid.zcorn) == expectZ;
fprintf('%s ZCORN extended: %d values for %dx%dx%d (expected %d)\n', ...
    subsref({'FAIL','PASS'}, substruct('{}', {okZ+1})), ...
    numel(model.grid.zcorn), model.grid.dims, expectZ);
if ~okZ, failures{end+1} = 'zcorn'; end

% Pore volume must transfer exactly through *PORVOL.
pvCMG = sum(model.props.porvol(model.props.active));
pvMRST = sum(poreVolume(Gf, Gf.rock));
relPV = abs(pvCMG - pvMRST)/pvMRST;
okPV = relPV < 1e-12;
fprintf('%s pore volume: MRST %.6f vs CMG %.6f m3 (rel %.2e)\n', ...
    subsref({'FAIL','PASS'}, substruct('{}', {okPV+1})), pvMRST, pvCMG, relPV);
if ~okPV, failures{end+1} = 'porevolume'; end

% Padding blocks must sit at the fracture depths, not stacked below the base.
fracZ = Gf.cells.centroids(nMat+1:end, 3);
zc = reshape(model.grid.zcorn, 2*model.grid.dims(1), 2*model.grid.dims(2), 2*nzTot);
depthErr = 0;
for q = 1:size(model.grid.fracIJK, 1)
    ijk = model.grid.fracIJK(q,:);
    ii = (2*ijk(1)-1):(2*ijk(1));
    jj = (2*ijk(2)-1):(2*ijk(2));
    zTop = mean(mean(zc(ii, jj, 2*ijk(3)-1)));
    zBot = mean(mean(zc(ii, jj, 2*ijk(3))));
    depthErr = max(depthErr, abs(0.5*(zTop+zBot) - fracZ(q)));
end
okDepth = depthErr < 1e-9;
fprintf('%s fracture block depths match centroids (max err %.3e m)\n', ...
    subsref({'FAIL','PASS'}, substruct('{}', {okDepth+1})), depthErr);
if ~okDepth, failures{end+1} = 'depth'; end

% One connection per NNC.
okConn = c.connections == size(Gf.nnc.cells, 1);
fprintf('%s connections: %d for %d NNCs\n', ...
    subsref({'FAIL','PASS'}, substruct('{}', {okConn+1})), ...
    c.connections, size(Gf.nnc.cells, 1));
if ~okConn, failures{end+1} = 'connections'; end

fprintf('\nblocks %d (%d inactive), padding layers %d\n', ...
    c.blocks, c.inactive, c.paddingLayers);
fprintf('wrote %s\n', info.dataset);

%% ---- summary -----------------------------------------------------------
fprintf('\n');
if isempty(failures)
    fprintf('=== test34: ALL CHECKS PASSED ===\n');
else
    fprintf('=== test34: %d FAILURES: %s ===\n', numel(failures), ...
        strjoin(failures, ', '));
end
