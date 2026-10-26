%% 3D, two-phase simulation of the Johansen formation
% This example shows the creation of a 3D grid for the Johansen formation using
% the inbuilt mrst utility function makeJohansenVEgrid. This utility function
% creates both the 3D and the VE grid for Johansen but in this case we only use
% the 3D grid. We then run a two-phase water-gas simulation with 100 years of CO2 injection
% into a brine filled reservoir followed by 1000 years post injection modelling.
% This example is similar to basic_example_3D.m but with a more realistic grid.

%% Load modules
clc;clear;close all;
mrstModule add ad-core;
mrstModule add ad-props;
mrstModule add ad-blackoil;
mrstModule add co2lab-common co2lab-ve co2lab-spillpoint;
mrstModule add shale;
mrstModule add hfm_v2;
mrstModule add mrst-gui;
mrstModule add coarsegrid;
mrstModule add vemmech deckformat

% EDFMGridCPG and the corner-point EDFM utilities are not part of any
% mrstModule; CPEDFM-CPG's startup puts them on the path. Without this the
% script only works in a session where some earlier run happened to add them.
cpgDir = fullfile(mrstPath('preprocessor'), '..', 'CPEDFM-CPG');
if isfolder(cpgDir)
    addpath(cpgDir);
    if exist('startupCPEDFMCPG', 'file') == 2, startupCPEDFMCPG(); end
end

% ORDER MATTERS. startupCPEDFMCPG does "mrstModule add hfm", and hfm's
% convexpolygonarea returns a single output while EDFMGridCPG calls it with
% two. Re-adding hfm_v2 afterwards puts the two-output copy back in front;
% without this the gridder dies with "Too many output arguments" at
% pebiAABBintersectCPGv5 line 648.
mrstModule add hfm_v2

%% setup default option values

% setup default option values
opt = struct('nonlinearTolerance' , 1e-6            , ...
             'splittingTolerance' , 1e-3            , ...
             'verbose'            , false           , ...
             'splittingVerbose'   , false);  

%% Grid and rock

% Load Johansen model
[G, rock, bcIx] = makeJohansenVEgridFull();

topLayer_global = false(G.cartDims);
topLayer_global(1:end, 1:end, 1) = true;
topLayer = find(topLayer_global(G.cells.indexMap));

topLayer5_global = false(G.cartDims);
topLayer5_global(1:end, 1:end, 5) = true;
layer_5 = find(topLayer5_global(G.cells.indexMap));

topLayer6_global = false(G.cartDims);
topLayer6_global(1:end, 1:end, 6) = true;
layer_6 = find(topLayer6_global(G.cells.indexMap));

% modify the permeability in the x and y direction of the layer 1 and in
% the z direction of layers 1 - 5
rock.perm(topLayer, 1:2) = rock.perm(topLayer, 1:2)*1e4; 

%Grabbed the faceIDs of layer5 and layer6
layer5faceIDs = [];   
layer6faceIDs = [];

for i = 1:numel(layer_5)
    
    layer5cellFaceData = G.cells.faces(G.cells.facePos(layer_5(i)):G.cells.facePos(layer_5(i)+1)-1, :);
    layer6cellFaceData = G.cells.faces(G.cells.facePos(layer_6(i)):G.cells.facePos(layer_6(i)+1)-1, :);
    layer5faceIDs = [layer5faceIDs; layer5cellFaceData(:, 1)];
    layer6faceIDs = [layer6faceIDs; layer6cellFaceData(:, 1)];
end


%Get the common faces between layers 5 and 6
commonfaces = intersect(layer5faceIDs, layer6faceIDs);

%% Grab all the boundary faces in the reservoir domain

nx = G.cartDims(1); ny=G.cartDims(2); nz=G.cartDims(3);

backdir = searchForBoundaryFaces(G, 'BACK');
leftdir = searchForBoundaryFaces(G, 'LEFT');
rightdir = searchForBoundaryFaces(G, 'RIGHT');
frontdir = searchForBoundaryFaces(G, 'FRONT');

figure,
plotGrid(G, 'faceColor', 'none', 'EdgeAlpha', 0.1)
% plotFaces(G, backdir, 'b');
% plotFaces(G, frontdir, 'b');
plotFaces(G, leftdir, 'y');
plotFaces(G, rightdir, 'y');
% plotFaces(G, topfaces, 'r');
view(30,60), axis tight off;

% Grab the faces of some bottom cell layers to fix the missing cells 

cellsx1_global = false(G.cartDims);
cellsx1_global(48, 7:13, 11) = true;
cellsx1 = find(cellsx1_global(G.cells.indexMap));

cellsx2_global = false(G.cartDims);
cellsx2_global(46, 10:13, 11) = true;
cellsx2 = find(cellsx2_global(G.cells.indexMap));

cellsx3_global = false(G.cartDims);
cellsx3_global(45, 10, 11) = true;
cellsx3 = find(cellsx3_global(G.cells.indexMap));

cellsy4_global = false(G.cartDims);
cellsy4_global(44:47, 6, 11) = true;
cellsy4 = find(cellsy4_global(G.cells.indexMap));

cellsy5_global = false(G.cartDims);
cellsy5_global(47, 14, 11) = true;
cellsy5 = find(cellsy5_global(G.cells.indexMap));

cellsy6_global = false(G.cartDims);
cellsy6_global(45:46, 10, 11) = true;
cellsy6 = find(cellsy6_global(G.cells.indexMap));

cellsy7_global = false(G.cartDims);
cellsy7_global(44, 11, 11) = true;
cellsy7 = find(cellsy7_global(G.cells.indexMap));

% concatenate cells that have the same faces of interest
x1cells = [cellsx1; cellsx3];
x2cells = cellsx2;
y3cells = [cellsy5; cellsy6; cellsy7];
y4cells = cellsy4;


% Grab the faces of the cells of interest above and plot
x1cellsfacesAll = G.cells.faces(mcolon(G.cells.facePos(x1cells), G.cells.facePos(x1cells + 1) - 1), :);
x1cellsfaces = x1cellsfacesAll(x1cellsfacesAll(:, 2) == 1  , 1);

figure,
plotGrid(G, 'faceColor', 'none', 'EdgeAlpha', 0.1);
plotFaces(G, x1cellsfaces, 'b');

x2cellsfacesAll = G.cells.faces(mcolon(G.cells.facePos(x2cells), G.cells.facePos(x2cells + 1) - 1), :);
x2cellsfaces = x2cellsfacesAll(x2cellsfacesAll(:, 2) == 2  , 1);

figure,
plotGrid(G, 'faceColor', 'none', 'EdgeAlpha', 0.1);
plotFaces(G, x2cellsfaces, 'b');

y3cellsfacesAll = G.cells.faces(mcolon(G.cells.facePos(y3cells), G.cells.facePos(y3cells + 1) - 1), :);
y3cellsfaces = y3cellsfacesAll(y3cellsfacesAll(:, 2) == 3  , 1);

figure,
plotGrid(G, 'faceColor', 'none', 'EdgeAlpha', 0.1);
plotFaces(G, y3cellsfaces, 'b');

y4cellsfacesAll = G.cells.faces(mcolon(G.cells.facePos(y4cells), G.cells.facePos(y4cells + 1) - 1), :);
y4cellsfaces = y4cellsfacesAll(y4cellsfacesAll(:, 2) == 4  , 1);

figure,
plotGrid(G, 'faceColor', 'none', 'EdgeAlpha', 0.1);
plotFaces(G, y4cellsfaces, 'b');


% Grab the cells in layer 10 directly above the missing cells in layer 11 

cellsLayer10int1_global = false(G.cartDims);
cellsLayer10int1_global(44:47, 7:9, 10) = true;
cellsLayer10int1 = find(cellsLayer10int1_global(G.cells.indexMap));

cellsLayer10int2_global = false(G.cartDims);
cellsLayer10int2_global(44, 10, 10) = true;
cellsLayer10int2 = find(cellsLayer10int2_global(G.cells.indexMap));

cellsLayer10int3_global = false(G.cartDims);
cellsLayer10int3_global(47, 10:13, 10) = true;
cellsLayer10int3 = find(cellsLayer10int3_global(G.cells.indexMap));

cellsLayer10int = [cellsLayer10int1; cellsLayer10int2; cellsLayer10int3];

figure,
plotGrid(G, 'faceColor', 'none', 'EdgeAlpha', 0.1);
plotGrid(G, cellsLayer10int, 'FaceColor','r');


%% FIND FAULT - set inner boundary
% The main fault is assumed to be sealing and must therefore be represented
% as an inner boundary in the 2D grid.  Start by locating cells on each
% side of fault, which is defined to be between cells with index i=43 and
% i=44 for j<44.

[ijk{1:3}] = ind2sub(G.cartDims, G.cells.indexMap(:));
ijk        = [ijk{:}];

% idInIndexMap = sub2ind(G.cartDims,i3,j3,k3);

faultcells_right = find(ijk(:,1) == 43 & ijk(:,2) <= 50 & ijk(:,3) <= 5);
faultcells_left = find(ijk(:,1) == 44 & ijk(:,2) <= 50 & ijk(:,3) <= 5);
faultcells = [faultcells_right;faultcells_left];

Allfaultcellsright = find(ijk(:, 1) == 43 & ijk(:,2) <= 50 & ijk(:,3) >= 6);
Allfaultcellsleft = find(ijk(:, 1) == 44 & ijk(:,2) <= 50 & ijk(:,3) >= 6);

Allfaultcells = [Allfaultcellsright; Allfaultcellsleft];

G_mult = G;

%% Export the matrix grid to CMG
% Done here, before the EDFM stage below replaces G with a grid that also
% carries fracture cells - writeCMG does not export those yet.
%
% makeJohansenVEgridFull builds its grdecl internally and then clears it
% (see its line 105), so it has to be rebuilt here. Corner-point grids need
% it: MRST discards COORD/ZCORN inside processGRDECL and they cannot be
% recovered from G. Reconstructed exactly as that function does.

mrstModule add preprocessor

jdir   = fullfile(mrstPath('co2lab-common'), 'data', 'johansen');
sector = fullfile(jdir, 'NPD5');
grdecl = readGRDECL([sector '.grdecl']);
K = reshape(load([sector, '_Permeability.txt'])', [], 1);
p = reshape(load([sector, '_Porosity.txt'])',     [], 1);
grdecl.PERMX = K; grdecl.PERMY = K; grdecl.PERMZ = 0.1*K;
grdecl.PORO  = p;

info = writeCMG(G, 'Johansen_co2.dat', ...
    'rock',   rock, ...
    'grdecl', grdecl, ...
    'title',  'JOHANSEN_CO2');

fprintf('CMG dataset : %s\n', info.dataset);
fprintf('CMG include : %s\n', info.include);
fprintf('blocks      : %d (%d inactive)\n', info.counts.blocks, info.counts.inactive);

%% Setting up EDFM
G.rock=rock;
try
    Globals,          % ADFNE initialisation; not needed for the manually
catch                 % constructed fracture planes below
    warning('adfne:missing', 'ADFNE not on the path; continuing without it.');
end
rng(123467)
tol=1e-5;

%Fault 1 part 1
[i j k] = ind2sub(G.cartDims, G.cells.indexMap);
subGrid=extractSubgrid(G, i==43 & j>0 & j < 22 & k==6);

% figure,
% plotGrid(G,'FaceColor','none');
% plotGrid(G,subGrid.cells.global,'FaceColor','r');

fcellCntrd=G.cells.centroids(subGrid.cells.global,:);

pt1=fcellCntrd(1,:);
pt2=fcellCntrd(end,:);

fracplanes(1).points = [pt1-[-500,0,2000];
  pt2-[-500,0,2000];
  pt2+[500,0,2000];
  pt1+[500,0,2000]];

%Fault 1 part 2
[i j k] = ind2sub(G.cartDims, G.cells.indexMap);
subGrid=extractSubgrid(G, i==43 & j>19 & j < 48 & k==6);

% figure,
% plotGrid(G,'FaceColor','none');
% plotGrid(G,subGrid.cells.global,'FaceColor','r');

fcellCntrd=G.cells.centroids(subGrid.cells.global,:);

pt1=fcellCntrd(1,:);
pt2=fcellCntrd(end,:);

fracplanes(2).points = [pt1-[-500,0,2000];
  pt2-[400,0,2000];
  pt2+[-400,0,2000];
  pt1+[500,0,2000]];

%Fault 2 part 1
[i j k] = ind2sub(G.cartDims, G.cells.indexMap);
subGrid=extractSubgrid(G, i==45 & j>0 & j < 22 & k==6);

% figure,
% plotGrid(G,'FaceColor','none');
% plotGrid(G,subGrid.cells.global,'FaceColor','r');

fcellCntrd=G.cells.centroids(subGrid.cells.global,:);

pt1=fcellCntrd(1,:);
pt2=fcellCntrd(end,:);

fracplanes(3).points = [pt1-[-500,0,2000];
  pt2-[-500,0,2000];
  pt2+[500,0,2000];
  pt1+[500,0,2000]];

%Fault 2 part 2
[i j k] = ind2sub(G.cartDims, G.cells.indexMap);
subGrid=extractSubgrid(G, i==45 & j>19 & j < 48 & k==6);

% figure,
% plotGrid(G,'FaceColor','none');
% plotGrid(G,subGrid.cells.global,'FaceColor','r');

fcellCntrd=G.cells.centroids(subGrid.cells.global,:);

pt1=fcellCntrd(1,:);
pt2=fcellCntrd(end,:);

fracplanes(4).points = [pt1-[-500,0,2000];
  pt2-[-600,0,2000];
  pt2+[600,0,2000];
  pt1+[500,0,2000]];

numNFplanes=numel(fracplanes);
for i=1:numNFplanes
    fracplanes(i).aperture = 1;
    fracplanes(i).poro=0.5;
    fracplanes(i).perm=[100*darcy, 100*darcy, 10*darcy];
end

[G,fracplanes]=EDFMGridCPG(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
G = fracturematrixShaleNNC3D(G,tol);
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);

%set up TPFA operators
TPFAoperators = setupShaleEDFMOpsTPFA(G, G.rock, tol);


%% Export the fractured model to CMG
% G now carries the four fracture planes as extra cells with
% matrix-fracture and fracture-fracture NNCs. They are exported as padding
% k-layers appended to ZCORN, with the connections written as *SCONNECT.
%
% Note the rock argument: after EDFM, 'rock' describes only the matrix
% cells while G has matrix + fracture cells, so G.rock is the matching one.

fi = cmgDetectFractureModel(G);
fprintf('\nfracture model detected: %s (%s)\n', fi.model, fi.reason);
fprintf('matrix cells %d, fracture cells %d, NNCs %d\n', ...
    fi.nMatrix, fi.nFrac, fi.nNNC);

infoF = writeCMG(G, 'Johansen_co2_edfm.dat', ...
    'rock',   G.rock, ...
    'grdecl', grdecl, ...
    'title',  'JOHANSEN_CO2_EDFM');

fprintf('CMG dataset : %s\n', infoF.dataset);
fprintf('blocks      : %d (%d inactive), padding layers %d\n', ...
    infoF.counts.blocks, infoF.counts.inactive, infoF.counts.paddingLayers);
fprintf('connections : %d\n', infoF.counts.connections);

pvMRST = sum(poreVolume(G, G.rock));
fprintf('pore volume : %.4f m3 (matrix + fractures)\n', pvMRST);


