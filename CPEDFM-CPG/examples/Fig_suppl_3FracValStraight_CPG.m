cpedfmCPGDir = fileparts(fileparts(mfilename('fullpath'))); addpath(cpedfmCPGDir); startupCPEDFMCPG();

%% This will generate Figures 14, 15, 17 and 18 in SPE-201243-PA
%
clear 
clc
close all
opt = struct('nkr',        2, ...
    'shouldPlot', 1 );

%% Load necessary modules, etc
mrstModule add hfm;             % hybrid fracture module
mrstModule add ad-core;
mrstModule add ad-props ad-core
mrstModule add mrst-gui;        % plotting routines
mrstModule add compositional; %mrstModule add ad-blackoil;
mrstModule add upr;

%% Basic Parameters
tol=1e-3;
physdim = [600 600 45];

wellRadius = 0.25*ft;

matrixPoro = 0.16;
matrixPerm =1*milli*darcy;

w_f = 0.01*ft;
k_f = 50*darcy;
x_f = 990*ft;


%% Create Box Reservoir Grid 

celldim = [50,50, 3];
[G_matrix, grdecl] = makeExampleCornerPointGrid(celldim, physdim); %#ok<NASGU>
GG_matrix = G_matrix;
   
   
   

if (opt.shouldPlot)
    figure,
    plotGrid(G_matrix)%, view(5,45)
end

G_matrix.rock = makeShaleRock(G_matrix, matrixPerm, matrixPoro);

%% Create Fracture System

 %-Hydraulic Fracture Creation
fracplanes = struct;

tlt=0;
fracplanes(1).points = [270-tlt 100 0; 270+tlt 500 0; 270+tlt 500 45; 270-tlt 100 45];
fracplanes(2).points = [306-tlt 150 0; 306+tlt 450 0; 306+tlt 450 45; 306-tlt 150 45];
fracplanes(3).points = [342-tlt 100 0; 342+tlt 500 0; 342+tlt 500 45; 342-tlt 100 45];



barrierPerm = 1e-3*nano*darcy;
%barrier 1
fracplanes(1).aperture = w_f;
fracplanes(1).poro= matrixPoro;
fracplanes(1).perm= barrierPerm;

%fracture
fracplanes(2).aperture = w_f;
fracplanes(2).poro= matrixPoro;
fracplanes(2).perm= k_f;

%barrier 2
fracplanes(3).aperture = w_f;
fracplanes(3).poro= matrixPoro;
fracplanes(3).perm= barrierPerm;

if (opt.shouldPlot)
    plotfracongrid(G_matrix,fracplanes); % visualize to check before pre-process
end

%% Create Wells

wells = struct;
%We use a short horizontal well to allow us find fracture/well intersection
halfToe2HeelDist = 5;
midpoint = [306,300,22.5];
wells(1).points=[midpoint(1)-halfToe2HeelDist,midpoint(2),midpoint(3); midpoint(1)+halfToe2HeelDist,midpoint(2),midpoint(3)]*meter;
wells(1).radius=wellRadius;

%% visualize to check before pre-process
G=G_matrix;

%% EDFM PreProcessing
[G,fracplanes]=EDFMGridCPG(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
%-Frac-Matrix NNCs
G = fracturematrixShaleNNC3D(G,tol);
%-Frac-Frac NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);
%-Well-Frac NNCs
[G,wells] = wellfractureNNCs3D(G,fracplanes,wells,tol);
G_pedfm = G;
G_edfm=G;

%
G = pMatFracNNCsCPEDFM_CPG(G,tol);
% G = pMatFracNNCs3D(G,tol);
%%
% MRST includes both natural variables and overall composition. 
useNatural = true;

% Name of problem and pressure range
casename = 'c1c2c3';
% casename = 'onlydecane';

pwf = 1000*psia;

%% Set up pEDFM operators
TPFAoperators = setupCPEDFMOpsTPFA_CPG(G, G.rock, tol);

%% Define three-phase compressible flow model
nkr=2;
[fluid, info] = getShaleCompFluidCase(casename);
info.pressure = 5700*psia;  %this is redundant for the bakken case
info.temp = 352.594;

eosname = 'prcorr';  %'srk','rk','prcorr'
G1cell = cartGrid([1 1],[1 1]);
G1cell = computeGeometry(G1cell);
EOSModel = EquationOfStateModel(G1cell, fluid, eosname);
%Surface Conditions
p_sc = 101325; %atmospheric pressure
T_sc = 288.706;% 60 Farenheit
[~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);
flowfluid = initSimpleADIFluid('n', [nkr, nkr, nkr], 'rho', [1000, rhoO_S, rhoG_S]);  

gravity reset off
model_cpedfm = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', false);
model_cpedfm.operators = TPFAoperators;

s0 = [1 0];%1;
    
W = [];
for wi=1:numel(wells)
    W = addWellShale(W, G, G.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [1,0],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir', 'z');
end
for wi=1:numel(wells)
    W(wi).components = info.initial;
end

%% Well trajectories
plotWell(G,W);

M = csvread('CMG_timestep2.csv',1);    % Benchmark timesteps
dt_list=M(1:50,1)*day;
time_list=cumsum(convertTo(dt_list,day));

schedule_cpedfm = simpleSchedule(dt_list, 'W', W);

%% Impose initial pressure equilibrium
state_cpedfm = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model_cpedfm.EOSModel);

%% Run simulations
[ws_cpedfm , states_cpedfm, report_cpedfm] = simulateScheduleAD(state_cpedfm, model_cpedfm, schedule_cpedfm);
G_cpedfm = G;

%% ------------------------------------------------------------------------
%                        pEDFM Model Modifications
%  ------------------------------------------------------------------------
G=G_pedfm;  % replace G with the previously-stored EDFM grid
G = pMatFracNNCsPEDFM_CPG(G,tol);
%% Set up EDFM operators
TPFAoperators = setupCPEDFMOpsTPFA_CPG(G, G.rock, tol);

%% Define three-phase compressible flow model
model_pedfm = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', false);
model_pedfm.operators = TPFAoperators;
    
W = [];
for wi=1:numel(wells)
    W = addWellShale(W, G, G.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [1,0],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir', 'z');
end
for wi=1:numel(wells)
    W(wi).components = info.initial;
end

schedule_pedfm = simpleSchedule(dt_list, 'W', W);

%% Impose initial pressure equilibrium
state_pedfm = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model_pedfm.EOSModel);

%% Run simulations
[ws_pedfm , states_pedfm, report_pedfm] = simulateScheduleAD(state_pedfm, model_pedfm, schedule_pedfm);
G_pedfm = G;



%% ------------------------------------------------------------------------
%                        EDFM Model Modifications
%  ------------------------------------------------------------------------
G=G_edfm;  % replace G with the previously-stored EDFM grid

%% Set up EDFM operators
TPFAoperators = setupShaleEDFMOpsTPFA(G, G.rock, tol);

%% Define three-phase compressible flow model
model_edfm = NaturalVariablesCompositionalModel(G, [], flowfluid, fluid, 'water', false);
model_edfm.operators = TPFAoperators;
    
W = [];
for wi=1:numel(wells)
    W = addWellShale(W, G, G.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [1,0],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir', 'z');
end
for wi=1:numel(wells)
    W(wi).components = info.initial;
end

schedule_edfm = simpleSchedule(dt_list, 'W', W);

%% Impose initial pressure equilibrium
state_edfm = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model_edfm.EOSModel);

%% Run simulations
[ws_edfm , states_edfm, report_edfm] = simulateScheduleAD(state_edfm, model_edfm, schedule_edfm);
%% Refference solution MRST Pebi

xmaxFS = physdim;
data.HF={[fracplanes(1).points(1,1),fracplanes(1).points(1,2);fracplanes(1).points(2,1),fracplanes(1).points(2,2)],...
    [fracplanes(2).points(1,1),fracplanes(2).points(1,2);fracplanes(2).points(2,1),fracplanes(2).points(2,2)],...
    [fracplanes(3).points(1,1),fracplanes(3).points(1,2);fracplanes(3).points(2,1),fracplanes(3).points(2,2)]};

% 2D fractures and grid processing
physdimScaled= [1,1,1]; % Scaled dimension for proper gridding
xmax=physdimScaled(1:2);
sdfScaled2D= cellfun(@(x) (x./xmaxFS(1:2)), data.HF, 'UniformOutput', false); % Scaled SDF in 2D
fractureScaled2D=sdfScaled2D;


% GRID
dx   = 0.03; % Approximate cell size
protD = {@(p) (0.00005/3)*ones(size(p,1),1)};
ApertureSD=1/25;
G2D = pebiGrid2D(dx, xmax(1:2), ...
    'cellConstraints', fractureScaled2D, ... % Fractures
    'protLayer',true,....
    'protD',protD,...
    'CCRefinement'   , true , ... % Refine fractures
    'CCEps',1/50,...
    'CCFactor' , 0.1  );    % Rel fracture cell size

G2DFD=G2D;
G2DFD.nodes.coords=G2DFD.nodes.coords.*xmaxFS(1:2);
figure,
plotGrid(G2DFD,G2DFD.cells.tag, 'faceColor', [1,1,1]*0.8,'LineWidth',1); axis tight
plotGrid(G2DFD, 'faceColor', 'none', 'edgeAlpha', 0.1);

% Make 3D layered grid
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


% Well trajectories
%Production well
ProdInj=[physdim(1)/2 physdim(2)/2 physdim(3)/2];
prod=[];
for i=1:size(ProdInj,1)
    distances1=sqrt(sum((G.cells.centroids(fracture_cells,:)-ProdInj(i,:)).^2,2));
    [~, distances1] = sort(distances1,'ascend');
    prod(i)=distances1(1);
end
   map = find(fracture_cells);
prod = map(prod);

%Slot Drill Fractures
pts=G.cells.centroids(:,:);
o = Distance(pts,FractureFulldim3D(2));
nearcell=(o<ApertureSD);
SdfTagInit=zeros(size(nearcell(:,1)));
for i=1:size(nearcell,2)
    SdfTagInit=SdfTagInit | nearcell(:,i);
end


% Natural fractures
TempTag=G.cells.tag;
[rmvID,ID]=find(SdfTagInit);
TempTag(rmvID)=0;
NFtagInit=TempTag;


G.cells.NF=NFtagInit;
G.cells.SDF=SdfTagInit;




% Make rock
rock = makeRock(G, matrixPerm, matrixPoro);
rock.perm(SdfTagInit) = k_f; % Fracture permeability
rock.perm(NFtagInit) = barrierPerm; % Fracture permeability
% FracPoreVolm=3.28;
% rock.poro(fracture_cells) =FracPoreVolm./G.cells.volumes(fracture_cells);
fracPoroPebi=sum(G_cpedfm.cells.volumes(G_cpedfm.Matrix.cells.num+1:end))*matrixPoro/sum(G.cells.volumes(fracture_cells));
rock.poro(fracture_cells) =fracPoroPebi;
model = NaturalVariablesCompositionalModel(G, rock, flowfluid, fluid, 'water', false);



%% Make wells

W = addWell([], G, rock, prod, ...
    'type', 'bhp', 'val', pwf, 'Comp_i', [1, 0] , 'WI',2.649022003063879e-13,'Name', ['Prod',num2str(1)],'Radius', wellRadius,'Dir','z','Sign', -1);

W.components = info.initial;

state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, model.EOSModel);

schedule = simpleSchedule(dt_list, 'W', W);

% Plot Grid
figure();
plotGrid(G, G.cells.HF , 'faceColor', 'y', 'edgeColor', 'none');
plotGrid(G,'faceColor', 'none', 'edgeAlpha', 0.1);
plotWell(G, schedule.control(1).W, 'color', 'r', 'fontSize', 18);
axis equal;
G_pebi=G;

%% Simulate problem
[ws_pebi, states_pebi, ~] = simulateScheduleAD(state, model, schedule);



%% Plot the pEDFM, EDFM and Explicit Solution Profiles 

figure, 
plotToolbar(G_cpedfm, states_cpedfm)
view(40,30);
axis tight equal;
set(gca,'CLim',[8000000 12000000],'DataAspectRatio',[1 1 1],'ZDir',...
    'reverse');
figure, 
plotToolbar(G_pedfm, states_pedfm)
view(40,30);
axis tight equal;

figure, 
plotToolbar(G_edfm, states_edfm)
view(40,30);
axis tight equal;

figure, 
plotToolbar(G_pebi, states_pebi)
view(40,30);
axis tight equal;

% ws = {ws, ws_pedfm, ws_edfm};
% states = {states, states_pedfm, states_edfm};

ws = {ws_cpedfm, ws_pedfm,ws_edfm,ws_pebi};
states = {states_cpedfm, states_pedfm,states_edfm,states_pebi};

names = {'CPEDFM', 'pEDFM','EDFM','Reference'};
plotWellSols(ws, cumsum(schedule.step.val), 'datasetnames', names,'field','qOs','linestyles',{'g','b','r','--k'})
