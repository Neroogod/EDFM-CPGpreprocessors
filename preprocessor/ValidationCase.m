close all
clear
clc
%% LOAD NECESSARY MODULES
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
mrstModule add ad-blackoil;     % AD blackoil solver
mrstModule add ad-core;         % AD core module
mrstModule add ad-props;        % AD properties
mrstModule add co2lab;
%% SET UP A STRUCTURED MATRIX GRID
tol=1e-5;
mrstModule add shale;

physdim=[1000,600,270];
fracPerm=1000*darcy;
fracPoro=0.5;
FracAperture=0.01*ft;
fracHalfLength=300;
MatPerm= 100*milli*darcy;

wellRadius=0.328084*ft;
pRef = 5000*psia;
pwf=1000*psia;
s0 = [0, 1, 0];
qinj=0.01;

fluid = initSimpleADIFluid('phases' , 'WOG', ...
    'mu' , [   1,  5, 0.2] * centi*poise     , ...
    'rho', [1000, 700, 250] * kilogram/meter^3, ...
    'c',   [1e-8, 1e-5, 1e-3] / barsa, ...
    'n'  , [   2,   2, 2], ...
    'pRef' , pRef);


totTime = 50*year;
nSteps = 15;
dt = rampupTimesteps(totTime, 90*day, nSteps);


%% tEDFM
% [grdecl,~,~]=makeModel3ExplicitFault(physdim);
celldims=[33,14,8];
[grdecl,~,~] = genericFaultCPGTrial(physdim, 'aperture', FracAperture,'fracHalfLength',fracHalfLength,...
    'fracHeight', 60,'tipNY',5,'nxRefine',60,...
    'ny',20,'nz',11,'FaultLoc','mid','gridType','cartNoref','dims',celldims);


G = processGRDECL(grdecl);
G = computeGeometry(G);

permMat =MatPerm*ones(G.cells.num,3);
permMat(:,3) = 0.01* permMat(:,3); 
G.rock=makeRock(G,permMat,0.36);
G.rock=makeRock(G,permMat,0.36);

G_flowac=G;

offset2=5;
offset=10;
offx=-0;
fracplanes(1).points = [physdim(1)/2-offx 0 0-offset2;
physdim(1)/2-offx 0 physdim(3)-30+offset;
physdim(1)/2-offx 600 physdim(3)-30+offset;
physdim(1)/2-offx 600 0-offset2];
X=fracplanes.points(:,1);
Y=fracplanes.points(:,2);
Z=fracplanes.points(:,3);
dims=G.cartDims;

rng(123)
%HR EDIT:  Make perturbation in z-direction
a   = rand ([numel(Z), 1]);
s   = randn([numel(Z), 2]);
fun = @(x, y, a, s) -5 *sin(pi * (x - 1/2)) - ...
    3 *sin(pi * (y + 2*x)) + ...
    a/2*sin(2*pi*x + s(1)) .* sin(4*pi*y + s(2));

for k = 1 : numel(Z)
    xi       = X(k) ./ physdim(1);
    eta      = Y(k) ./ physdim(2);
    Z(k) = Z(k) -  fun(xi, eta, a(k), s(k,:));
end
fracplanes(1).points=[X Y Z];
fracplanes(1).aperture = 0.01*ft;
fracplanes(1).poro=fracPoro;
fracplanes(1).perm=[fracPerm fracPerm fracPerm];

% Plane normals
points = [fracplanes(1).points;fracplanes(1).points(1,:)];
diffp = diff(points,1);
normal = cross(diffp(1,:), diffp(2,:));
fracplanes(1).normal = normal/norm(normal);

[G,fracplanes]=EDFMGridCPG(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
  %-Frac-Matrix NNCs
G = fracturematrixShaleNNC3D(G,tol);
  %-Frac-Frac NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);


fi = cmgDetectFractureModel(G);
fprintf('\nfracture model detected: %s (%s)\n', fi.model, fi.reason);
fprintf('matrix cells %d, fracture cells %d, NNCs %d\n', ...
    fi.nMatrix, fi.nFrac, fi.nNNC);

infoF = writeCMG(G, 'val_co2_edfm.dat', ...
    'rock',   G.rock, ...
    'grdecl', grdecl, ...
    'title',  'VAL_CO2_EDFM');

fprintf('CMG dataset : %s\n', infoF.dataset);
fprintf('blocks      : %d (%d inactive), padding layers %d\n', ...
    infoF.counts.blocks, infoF.counts.inactive, infoF.counts.paddingLayers);
fprintf('connections : %d\n', infoF.counts.connections);

pvMRST = sum(poreVolume(G, G.rock));
fprintf('pore volume : %.4f m3 (matrix + fractures)\n', pvMRST);


