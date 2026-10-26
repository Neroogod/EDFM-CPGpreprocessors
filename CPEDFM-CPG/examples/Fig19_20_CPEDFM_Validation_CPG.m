cpedfmCPGDir = fileparts(fileparts(mfilename('fullpath'))); addpath(cpedfmCPGDir); startupCPEDFMCPG();

%% It shows the limitation of pEDFM in modeling slanted sealing fractures

clc;clear;close all;

%% Load necessary modules
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
mrstModule add compositional;   % replace ad-blackoil with compositional module
mrstModule add ad-blackoil;     % AD blackoil solver
mrstModule add ad-core;         % AD core module
mrstModule add ad-props;        % AD properties
mrstModule add upr;             % AD upr
mrstModule add mpfa;

%% Set up a matrix grid
FracPerm=1e-6*nano*darcy;

FracPoro=0.001;
FracAperture=5*milli*meter;
MatrixPerm=1*darcy;
MatrixPoro=0.1;
Injrate=1/500;
pwf=10e6;


WellXloc=350;
WellRadius=0.01;

celldim = [50,50, 3];
physdim = [1000,1000, 150];
[G, grdecl] = makeExampleCornerPointGrid(celldim, physdim); %#ok<NASGU>
G.rock=makeRock(G,MatrixPerm,MatrixPoro);

% Plot Matrix Grid
% figure;
% plotGrid(G,'facealpha',0,'edgealpha',0.5);
% view(30,45);

%% Set up fracture planes

% % Fracture plane 1
fracplanes(1).points = [300 0 0;
    600 1000 0;
    600 1000 150;
    300 0 150];

fracplanes(1).aperture = FracAperture;
fracplanes(1).poro=FracPoro;
fracplanes(1).perm=FracPerm;

xmaxFS = physdim;
data.HF={[fracplanes(1).points(1,1),fracplanes(1).points(1,2);fracplanes(1).points(2,1),fracplanes(1).points(2,2)]};
%% 2D fractures and grid processing
physdimScaled= [1,1,1]; % Scaled dimension for proper gridding
xmax=physdimScaled(1:2);
sdfScaled2D= cellfun(@(x) (x./xmaxFS(1:2)), data.HF, 'UniformOutput', false); % Scaled SDF in 2D
Globals;
set1=[];
% Add NF to SDF
fractureScaled2D=sdfScaled2D;
for i=1:size(set1,1)
    fractureScaled2D=[fractureScaled2D {[set1(i,1),set1(i,2);set1(i,3),set1(i,4)]}];
end

%% GRID
dx   = 0.03; % Approximate cell size
protD = {@(p) 0.0005*ones(size(p,1),1)};
ApertureSD=FracAperture;
G2D = pebiGrid2D(dx, xmax(1:2), ...
    'cellConstraints', fractureScaled2D, ... % Fractures
    'protLayer',true,....
    'protD',protD,...
    'CCRefinement'   , true , ... % Refine fractures
    'CCEps',1/50,...
    'CCFactor' , 0.05  );    % Rel fracture cell size

G2DFD=G2D;
G2DFD.nodes.coords=G2DFD.nodes.coords.*xmaxFS(1:2);
% figure,
% plotGrid(G2DFD,G2DFD.cells.tag, 'faceColor', [1,1,1]*0.8,'LineWidth',1); axis tight
% plotGrid(G2DFD, 'faceColor', 'none', 'edgeAlpha', 0.1);

%% Make 3D layered grid
layers = diff(linspace(0, xmaxFS(3), celldim(3)+1));
G = makeLayeredGrid(G2DFD, layers);
G = computeGeometry(G);

% Convert the fracture to full dim 3D
FractureFullDim=cellfun(@(x) (x.*xmaxFS(1:2)), fractureScaled2D, 'UniformOutput', false);

FractureFulldim3D=[];
for i=1:numel(FractureFullDim)
    Curfracdata=FractureFullDim{1,i};
    FractureFulldim3D=[FractureFulldim3D {[Curfracdata(1,:) 0;Curfracdata(2,:) 0;Curfracdata(2,:) xmaxFS(3);Curfracdata(1,:) xmaxFS(3)]}];
end
G.cells.tag = repmat(G2DFD.cells.tag,celldim(3), 1);
G.cells.HF=G.cells.tag;
fracture_cells=G.cells.tag;
%% Well trajectories
CordInj=[physdim(1)/2-WellXloc, physdim(2)/2, physdim(3)/2];
inj=[];
for i=1:size(CordInj,1)
    distances1=sqrt(sum((G.cells.centroids(:,:)-CordInj(i,:)).^2,2));
    [~, distances1] = sort(distances1,'ascend');
    inj(i)=distances1(1);
end

%Production well
ProdInj=[physdim(1)/2+WellXloc physdim(2)/2 physdim(3)/2];
% To Do:   Find the CordInj and ProdInj automatically
prod=[];
for i=1:size(ProdInj,1)
    distances1=sqrt(sum((G.cells.centroids(:,:)-ProdInj(i,:)).^2,2));
    [~, distances1] = sort(distances1,'ascend');
    prod(i)=distances1(1);
end

% Make rock
rock = makeRock(G, MatrixPerm, MatrixPoro);
rock.perm(fracture_cells) = FracPerm; % Fracture permeability
FracPoreVolm=3.28;
rock.poro(fracture_cells) =FracPoreVolm./G.cells.volumes(fracture_cells);

fluid = initSimpleADIFluid('mu' , [   1,  2, 0.2] .* centi*poise     , ...
    'rho', [1000, 700, 250] .* kilogram/meter^3, ...
    'n'  , [   2,   2, 2]);

% Add compressibility to fluid
pRef = 20e6;
c_w = 4e-10;
c_o = 1e-5/barsa;
c_g = 1e-3/barsa;

fluid.bW = @(p) exp((p - pRef)*c_w);
fluid.bO = @(p) exp((p - pRef)*c_o);
fluid.bG = @(p) exp((p - pRef)*c_g);

%Surface Conditions
p_sc = 101325;  %atmospheric pressure
T_sc = 288.706; % 60 Farenheit


model = ThreePhaseBlackOilModel(G, rock, fluid, 'disgas', false, 'vapoil', false);
totTime = 15*year;

%% Make wells
rw=WellRadius; % Well radius
W = addWell([], G, rock, inj, ...
    'type', 'rate', 'val',Injrate,'rR',0.4238,'WI',3.17770105725592e-14,'InnerProduct', 'ip_tpf', 'Comp_i', [1, 0, 0], 'name', 'inj','Radius', rw,'Dir','z','Sign', 1);

W = addWell(W, G, rock, prod, ...
    'type', 'bhp', 'val', pwf,'rR',0.4238,'WI',3.29144440310462e-14, 'InnerProduct', 'ip_tpf', 'Comp_i', [1, 1, 0] ./ 2, 'name', 'prod_1','Radius', rw,'Dir','z','Sign', -1);


s0 = [0.2, 0.8, 0];   
state  = initResSol(G, pRef, s0);
dt = rampupTimesteps(totTime, 90*day, 8);
schedule = simpleSchedule(dt, 'W', W);

% Plot Grid
figure();
plotGrid(G, G.cells.HF , 'faceColor', 'y', 'edgeColor', 'none');
plotGrid(G,'faceColor', 'none', 'edgeAlpha', 0.1);
plotWell(G, schedule.control(1).W, 'color', 'r', 'fontSize', 18);
axis equal;
G_pebi=G;

%% Simulate problem
[ws_pebi, states_pebi, ~] = simulateScheduleAD(state, model, schedule);


%% Fracture-Fracture NNCs
% This calculates the transmissibilities between connected fracture and
% fracture grid blocks. Information is saved under G.nnc.
[G, grdecl] = makeExampleCornerPointGrid(celldim, physdim); %#ok<NASGU>
G.rock=makeRock(G,MatrixPerm,MatrixPoro);
poroMultiplyer=1.0001;
fracplanes(1).poro=FracPoro/poroMultiplyer;

tol=1e-5;
[G,fracplanes]=EDFMGridCPG(G,fracplanes,...
    'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));

%% Fracture-Matrix non-Neighbouring Connections (NNC)
% This calculates the transmissibilities between connected matrix and
% fracture grid blocks. Information is saved under G.nnc.

tol=1e-5;
%-Frac-Matrix NNCs
G = fracturematrixShaleNNC3D(G,tol);
%-Frac-Frac NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);
G_cpedfm=G;
G_cpedfm = pMatFracNNCsCPEDFM_CPG(G_cpedfm,tol);

gravity reset off
model_pedfm = ThreePhaseBlackOilModel(G_cpedfm, [], fluid, 'disgas', false, 'vapoil', false);
%% Setup TPFA Operators
% Generate operators that incorporates all the EDFM NNCs
G_cpedfm.rock.poro(G_cpedfm.Matrix.cells.num+1:end) =FracPoreVolm./G_cpedfm.cells.volumes(G_cpedfm.Matrix.cells.num+1:end);
TPFAoperators_pedfm = setupCPEDFMOpsTPFA_CPG(G_cpedfm, G_cpedfm.rock, tol);
model_pedfm.operators = TPFAoperators_pedfm;
%% Add wells
wellRadius = 0.01;
[nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));

inj=[];
for i=1:size(CordInj,1)
    distances1=sqrt(sum((G.cells.centroids(:,:)-CordInj(i,:)).^2,2));
    [~, distances1] = sort(distances1,'ascend');
    inj(i)=distances1(1);
end

%Production well
prod=[];
for i=1:size(ProdInj,1)
    distances1=sqrt(sum((G.cells.centroids(:,:)-ProdInj(i,:)).^2,2));
    [~, distances1] = sort(distances1,'ascend');
    prod(i)=distances1(1);
end




% Injector
cellinj = inj;
W_pedfm   = addWell([], G_cpedfm.Matrix, G_cpedfm.Matrix.rock, cellinj,'InnerProduct', 'ip_tpf','Type', ...
    'rate', 'Val', Injrate, 'Radius', wellRadius, 'Comp_i', [1, 0, 0], 'Name', 'inj', ...
    'rR',0.4238,'WI',3.17770105725592e-14,'Sign', 1);
% Producer
cellprod = prod;
W_pedfm   = addWell(W_pedfm, G_cpedfm.Matrix, G_cpedfm.Matrix.rock, cellprod, 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', pwf, 'Radius', wellRadius, 'Comp_i', [1, 1, 0] ./ 2, 'Name', 'prod_1', ...
    'rR',0.4238,'WI',3.29144440310462e-14,'Sign', -1);


% Plot setup
figure,
plotEDFMgrid(G);
hold on;
plotWell(G,W_pedfm);

%% Set up initial state and schedule

s0 = [0.2, 0.8, 0];  
state_pedfm  = initResSol(G_cpedfm, pRef, s0);
schedule_pedfm = simpleSchedule(dt, 'W', W_pedfm);

%% Simulate problem
[ws_pedfm, states_pedfm, ~] = simulateScheduleAD(state_pedfm, model_pedfm, schedule_pedfm);

ws = {ws_pebi, ws_pedfm};
states = {states_pebi, states_pedfm};

names = {'explicit','pEDFM'};
plotWellSols(ws, cumsum(schedule.step.val), 'datasetnames', names,'linestyles',{'r','--g'})


figure,
plotToolbar(G_pebi,states_pebi)
colormap jet;
axis tight equal;

figure,
plotToolbar(G_cpedfm,states_pedfm)
colormap jet;
axis tight equal;
