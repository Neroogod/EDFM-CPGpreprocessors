%% test31_eagleford_matrix_cmg
% Phase 1 acceptance test: export the Eagle Ford matrix grid (the locally
% refined tensorGrid that meshHFsystem builds for eagleFordEDFMpre.m) to a
% CMG dataset. Matrix only - fractures are added in phases 3-5.
%
% The point of this case is the non-uniform spacing: nxRefine=9 near the
% hydraulic fractures and nxRefine_small=4 between clusters, which is what
% forces *GRID *VARI rather than *CART.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
mrstModule add hfm ad-core ad-props upr shale
% createMultiStageHFs/meshHFsystem live in the shale module's utils folder,
% which mrstModule does not put on the path by itself.
shaleUtils = fullfile(mrstPath('shale'), 'utils');
if isfolder(shaleUtils), addpath(shaleUtils); end

outDir = fullfile(thisDir, 'output_cmg');
if ~isfolder(outDir), mkdir(outDir); end

%% ---- Eagle Ford geometry, matching eagleFordEDFMpre.m -----------------
tol = 1e-5;
physdim = [150 300 38.1];
fracHalfLength = 180;
fracHeight = 75*ft;
fractureSpacing = 150;
fracPerStage = 3;
NumStages = 1;
clusterSpacing = 10;
fracAperture = 3.0*milli*meter;
MatrixPoro = 0.12;
MatrixPerm = 900*nano*darcy;
fracCentroid1 = [65, fracHalfLength/2, physdim(3)-fracHeight/2];

[fracplanes, frac_centroid_s] = createMultiStageHFs( ...
    'numStages', NumStages, 'fracSpacing', fractureSpacing, ...
    'numFracsPerStage', fracPerStage, 'fracHalfLength', fracHalfLength, ...
    'fracHeight', fracHeight, 'clusterSpacing', clusterSpacing, ...
    'fracCentroid1', fracCentroid1, 'isStencil', 1); %#ok<ASGLU>

G = meshHFsystem(physdim, frac_centroid_s, 'numStages', NumStages, ...
    'numFracsPerStage', fracPerStage, 'ny', 20, 'nz', 5, ...
    'nxRefine_small', 4, 'nxRefine', 9, 'fracSpacing', fractureSpacing, ...
    'aperture', fracAperture);
G = computeGeometry(G);
G.rock = makeRock(G, MatrixPerm, MatrixPoro);

fprintf('Eagle Ford matrix grid: %d cells, cartDims %s\n', ...
    G.cells.num, mat2str(G.cartDims));

%% ---- detection --------------------------------------------------------
gi = cmgDetectGrid(G);
fi = cmgDetectFractureModel(G);
fprintf('detected grid family : %s (uniform=%d)\n', gi.family, gi.uniform);
fprintf('detected fracture    : %s (%s)\n', fi.model, fi.reason);
fprintf('dx range             : %.4g .. %.4g m over %d columns\n', ...
    min(gi.di), max(gi.di), numel(gi.di));

assert(strcmp(gi.family, 'vari'), ...
    'Eagle Ford must export as *VARI: its spacing is deliberately uneven.');
assert(strcmp(fi.model, 'matrix'), 'No fractures expected at this stage.');

%% ---- export -----------------------------------------------------------
datFile = fullfile(outDir, 'eagleford_matrix.dat');
info = writeCMG(G, datFile, 'simulator', 'IMEX', 'cmgVersion', '2024.30', ...
    'datum', 3000, 'title', 'EAGLEFORD_MATRIX', 'hostDeck', 'runnable');

fprintf('\nwrote %s\n', info.dataset);
fprintf('      %s\n', info.include);
fprintf('blocks %d (%d inactive), connections %d\n', ...
    info.counts.blocks, info.counts.inactive, info.counts.connections);

%% ---- conservation check ----------------------------------------------
% The *VARI block sizes must reproduce the grid's bulk volume exactly, and
% pore volume follows from the porosity array.
model = buildCMGModel(G, G.rock, 'datum', 3000);
bulkCMG = sum(model.grid.di) * sum(model.grid.dj) * sum(model.grid.dk);
bulkMRST = sum(G.cells.volumes);
relBulk = abs(bulkCMG - bulkMRST)/bulkMRST;

pvMRST = sum(poreVolume(G, G.rock));
pvCMG = bulkCMG * MatrixPoro;
relPV = abs(pvCMG - pvMRST)/pvMRST;

fprintf('bulk volume  MRST %.6f m3 vs CMG %.6f m3  (rel %.3e)\n', ...
    bulkMRST, bulkCMG, relBulk);
fprintf('pore volume  MRST %.6f m3 vs CMG %.6f m3  (rel %.3e)\n', ...
    pvMRST, pvCMG, relPV);

ok = relBulk < 1e-12 && relPV < 1e-12 && exist(info.dataset, 'file') == 2;
if ok
    fprintf('\n=== test31: PASS (export written, volume conserved) ===\n');
else
    fprintf('\n=== test31: FAILURES (relBulk %.3e relPV %.3e) ===\n', ...
        relBulk, relPV);
end
