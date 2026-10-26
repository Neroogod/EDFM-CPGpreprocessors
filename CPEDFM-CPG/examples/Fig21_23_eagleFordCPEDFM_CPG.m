cpedfmCPGDir = fileparts(fileparts(mfilename('fullpath'))); addpath(cpedfmCPGDir); startupCPEDFMCPG();

clc;
clear;
close all;

try
    %% This code requires ADFNE and the Statistics & Machine Learning Toolbox
    % The ADFNE folder should be put within the shale module folder
    Globals;
    warning('on');
catch
    %% Print warning when ADFNE is not available
    warning('adfne:missing', 'ADFNE or Statistics & Machine Learning Toolbox is not available!');
    warning('off', 'adfne:missing');
end


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
tol=1e-5;
physdim = [150 300 38.1];
fracHalfLength=180;
fracHeight=75*ft;
fractureSpacing=150;
fracPerStage=3;
NumStages=1;
clusterSpacing=10;
fracAperture=3.0*milli*meter;
fracPoro=0.5;
fracPerm=1.0*darcy;
MatrixPoro=0.12;
MatrixPerm=900*nano*darcy;
% MatrixPerm=10*micro*darcy;
wellLength=150;
wellRadius=0.328084*ft;
heelCooord=[0, 0, physdim(3)-tol];
ToeCoord=[heelCooord(1)+wellLength,0, physdim(3)-tol];
BHP=1000*psia;
fracCentroid1= [65, fracHalfLength/2, physdim(3)-fracHeight/2]; % Centroid of the first fracture
%% Create Fracture System
[fracplanes,frac_centroid_s] = createMultiStageHFs('numStages',NumStages,'fracSpacing', fractureSpacing,...
      'numFracsPerStage', fracPerStage,'fracHalfLength', fracHalfLength,'fracHeight',fracHeight, ...
      'clusterSpacing', clusterSpacing,'fracCentroid1',fracCentroid1,'isStencil',1);
  
for i=1:numel(fracplanes)
     fracplanes(i).aperture = fracAperture; 
     fracplanes(i).poro = fracPoro;
     fracplanes(i).perm = fracPerm;
end
fracplanes(4).points = [45 0 0; 51 50 0; 51 50 38.1; 45 0 38.1];
fracplanes(5).points = [50 20 0; 42 80 0; 42 80 38.1; 50 20 38.1];
% fracplanes(6).points = [35 80 0; 60 140 0; 60 140 38.1; 35 80 38.1];
% fracplanes(7).points = [60 110 0; 40 180 0; 40 180 38.1; 60 110 38.1];

AddAmnt=56;
fracplanes(6).points = [45+AddAmnt 0 0; 51+AddAmnt 50 0; 51+AddAmnt 50 38.1; 45+AddAmnt 0 38.1];
fracplanes(7).points = [50+AddAmnt 20 0; 42+AddAmnt 80 0; 42+AddAmnt 80 38.1; 50+AddAmnt 20 38.1];
% fracplanes(10).points = [35+AddAmnt 80 0; 60+AddAmnt 140 0; 60+AddAmnt 140 38.1; 35+AddAmnt 80 38.1];
% fracplanes(11).points = [60+AddAmnt 110 0; 40+AddAmnt 180 0; 40+AddAmnt 180 38.1; 60+AddAmnt 110 38.1];



% fracplanes(5).points = [100 0 0; 131 300 0; 131 300 38.1; 100 0 38.1];
% fracplanes(4).points = [57.78 140+47.37 0; 99.44 140+69 0; 99.44 140+69 38.1; 57.78 140+47.37 38.1];
%  fracplanes(4).points = [30 55 0; 120 55 0; 120 55 38.1; 30 55 38.1];
barrierPerm = 1e-3*nano*darcy;
fracplanes(4).aperture = fracAperture;
fracplanes(4).poro= fracPoro;
fracplanes(4).perm= barrierPerm; 
fracplanes(5).aperture = fracAperture;
fracplanes(5).poro= fracPoro;
fracplanes(5).perm= barrierPerm; 
fracplanes(6).aperture = fracAperture;
fracplanes(6).poro= fracPoro;
fracplanes(6).perm= barrierPerm; 
fracplanes(7).aperture = fracAperture;
fracplanes(7).poro= fracPoro;
fracplanes(7).perm= barrierPerm; 


% fracplanes(8).aperture = fracAperture;
% fracplanes(8).poro= fracPoro;
% fracplanes(8).perm= barrierPerm; 
% fracplanes(9).aperture = fracAperture;
% fracplanes(9).poro= fracPoro;
% fracplanes(9).perm= barrierPerm; 
% fracplanes(10).aperture = fracAperture;
% fracplanes(10).poro= fracPoro;
% fracplanes(10).perm= barrierPerm; 
% fracplanes(11).aperture = fracAperture;
% fracplanes(11).poro= fracPoro;
% fracplanes(11).perm= barrierPerm; 


G_layout = meshHFsystem(physdim,frac_centroid_s,'numStages',NumStages,'numFracsPerStage',fracPerStage,...
    'ny',20,'nz',5,'nxRefine_small',4,'nxRefine',9,'fracSpacing', fractureSpacing,'aperture',fracAperture);
[G_matrix, grdecl] = makeExampleCornerPointGrid(G_layout.cartDims, physdim); %#ok<NASGU>
clear G_layout

G_matrix.rock=makeRock(G_matrix,MatrixPerm,MatrixPoro);
if (opt.shouldPlot)
    plotfracongrid(G_matrix,fracplanes,'label',false); % visualize to check before pre-process
end
view(40,30);


%% Stochastic Natural Fracture Generation using ADFNE
rng(12345); 
set1 = []; set2 = []; set3 = [];

try
    %% This code requires ADFNE and the Statistics & Machine Learning Toolbox
    set1 = Field(DFN('dim',3,'n',50,'dir',60,'ddir',1e3,'minl',10,'mu',20, ...
        'maxl',20,'bbx',[tol,tol,tol,40-tol, ...
        physdim(2)-tol,physdim(3)-tol],'dip',45,'ddip',-100, ...
        'shape','l','q',4),'Poly');
    set2 = Field(DFN('dim',3,'n',50,'dir',-60,'ddir',1e3,'minl',10,'mu',20, ...
        'maxl',20,'bbx',[108,tol,tol,physdim(1)-tol, ...
        physdim(2)-tol,physdim(3)-tol],'dip',45,'ddip',-100, ...
        'shape','l','q',4),'Poly');
    set3 = Field(DFN('dim',3,'n',50,'dir',-60,'ddir',1e3,'minl',10,'mu',20, ...
        'maxl',20,'bbx',[tol,80,tol,physdim(1)-tol, ...
        physdim(2)-tol,physdim(3)-tol],'dip',45,'ddip',-100, ...
        'shape','l','q',4),'Poly');
catch ME
    %% Load generated fractures when ADFNE is not available
    warning('EagleFord:ADFNEFallback', ...
        'ADFNE generation stopped early (%s). Using available/cached sets.', ME.message);
    if isfile('EagleFordFracs.mat'), load('EagleFordFracs.mat'); end
end
    
[set1_,nonPlanarSets1,fracArea1] = processOptionalStochFracSet(set1);
[set2_,nonPlanarSets2,fracArea2] = processOptionalStochFracSet(set2);
[set3_,nonPlanarSets3,fracArea3] = processOptionalStochFracSet(set3);
fracArea = [fracArea1;fracArea2;fracArea3];

fprintf('%d of %d set1 fracs were OK while %d of %d set2 fracs were OK \n',...
    numel(set1_),numel(set1),numel(set2_),numel(set2));

fprintf('Number of nonplanar sets in sets 1 and 2 are : %d and %d respectively\n',...
    numel(nonPlanarSets1),numel(nonPlanarSets2));

% figure; 
if (opt.shouldPlot) %ADFNE Plot funcs     
    try      
    % Plot stereonet if ADFNE is available
        figure;          
        Draw('ply',set1_);    
        Draw('ply',set2_);view(45,30)             
        o = Orientation([set1_;set2_]);          
        figure,          
        Stereonet([o.Dip],[o.Dir],'density',true,'marker','*','ndip',6,...              
            'ndir',24,'cmap',@jet,'color','y');   
    catch 
    % Don’t plot if ADFNE is not available             
    end    
end 


fracSet = [set1_ ;set2_;set3_]; 

numHFplanes = numel(fracplanes);
numNFplanes = numel(fracSet); 
totalNumFracPlanes = numHFplanes + numNFplanes;
startIdx = numHFplanes + 1;
for i=1:numNFplanes
    idxGlobal = numHFplanes + i;
    fracplanes(idxGlobal).points = fracSet{i}(1:end-1,:);
    fracplanes(idxGlobal).aperture = 100*nano*meter; %1*micro*meter;
    fracplanes(idxGlobal).poro=0.5;
    if(mod(i,2)==false)%even NF idx
       fracplanes(idxGlobal).perm=1e-3*nano*darcy;%100*micro*darcy;%0.01*darcy;
    else%odd NF idx
       fracplanes(idxGlobal).perm=100*micro*darcy;%100*micro*darcy;%0.01*darcy;
    end
end
G_matrix.numHFplanes=numHFplanes;
G_matrix.numNFplanes=numNFplanes;

%% Create Wells
wells = struct;
wells(1).points=[heelCooord; ToeCoord];
wells(1).radius=wellRadius;

%% visualize to check before pre-process
G=G_matrix;

%% pEDFM PreProcessing
% The preprocessing takes a while. The user can  remove line 147 to 159 and
% instead load the preprocessed grid at line 162

[G,fracplanes]=EDFMGridCPG(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
  %-Frac-Matrix NNCs
G = fracturematrixShaleNNC3D(G,tol);
  %-Frac-Frac NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);
  %-Well-Fracs NNCs
[G,wells] = wellfractureNNCs3D_CPG(G,fracplanes,wells,tol);
  % HR: Continuous-Projection-based NNCs
G = pMatFracNNCsCPEDFM_CPG(G,tol);

% Set up pEDFM operators
TPFAoperators = setupCPEDFMOpsTPFA_CPG(G, G.rock, tol);

% Load preprocessed grid
% load('EagleFordGrid.mat');

if (opt.shouldPlot)
    figure,
    plotfracSystem_NF(G,fracplanes,3,wells,'label',false)
end
view(40,30);

% Name of problem and pressure range
casename = 'eagleford';
pwf = BHP;




%% Define three-phase compositional flow model
% We define a three-phase compositional model with or without water.
% This is done by first instantiating the compositional model, and then
% manually passing in the internal transmissibilities and the topological
% neighborship from the embedded fracture grid.

[fluid, info] = getShaleCompFluidCase(casename);

eosname = 'prcorr';  %'srk','rk','prcorr'
G1cell = cartGrid([1 1],[1 1]);
G1cell = computeGeometry(G1cell);
EOSModel = EquationOfStateModel(G1cell, fluid, eosname);

%Surface Conditions
p_sc = 101325; %atmospheric pressure
T_sc = 288.706;% 60 Farenheit
[~, ~, ~, ~, ~, rhoO_S, rhoG_S] = standaloneFlash(p_sc, T_sc, info.initial, EOSModel);

flowfluid = initSimpleADIFluid('phases', 'WOG', 'n', [opt.nkr, opt.nkr, opt.nkr], 'rho', [1000, rhoO_S, rhoG_S]); 
gravity reset on

% model = NatVarsShaleModel(G, [], flowfluid, fluid, 'water', true);
% model.operators = TPFAoperators;
arg = {G, [], flowfluid, fluid, 'water', true};

sparse_backend = SparseAutoDiffBackend();
useNatural=true;
if useNatural
    constructor = @NatVarsShaleModel;
else
    constructor = @GenericOverallCompositionModel;
end
modelSparseAD = constructor(arg{:}, 'AutoDiffBackend', sparse_backend);

modelSparseAD.operators = TPFAoperators;

%% Pick linear solver
% Plain sparse backslash deliberately bypasses selectLinearSolverAD's
% cell-reduction and ordering metadata, which are incompatible with the
% custom NatVarsShaleModel linearized-problem representation.
linsolve = BackslashSolverAD();
disp(linsolve)
nls = NonLinearSolver('LinearSolver', linsolve);
%% Set up initial state and schedule

ncomp = fluid.getNumberOfComponents();
s0 = [0.17, 0.74, 0.09];   

state = initCompositionalState(G, info.pressure, info.temp, s0, info.initial, modelSparseAD.EOSModel);

W = [];
for wi=1:numel(wells)
    W = addWellShale(W, G, G.Matrix.rock, wells(wi).XFracCellIDs, ...
        'comp_i', [0.17, 0.74, 0.09],'Name', ['Prod',num2str(wi)], 'Val',...
        pwf, 'sign', -1, 'Type', 'bhp','Radius', wells(wi).radius,'Dir','x');
end
for wi=1:numel(wells)
    W(wi).components = info.initial;
end

totTime = 30*year;
nSteps = 15;
dt = rampupTimesteps(totTime, 30*day, nSteps);

schedule = simpleSchedule(dt, 'W', W);
%% Simulate problem
% [ws, states, reports] = simulateScheduleAD(state, model, schedule, 'Verbose', true);
[ws, states, reports] = simulateScheduleAD(state, modelSparseAD, schedule, 'nonlinearsolver', nls, 'Verbose', true);
wsflashed = flashWS(ws, 'recoveryMethod', 'CGEOR','casename', 'eagleford'); 
%% plotting

figure, 
plotToolbar(G, states)
view(40,30);
axis tight equal;

figure,
plotWellSols(wsflashed,cumsum(schedule.step.val))

figure,
plotWellSols(ws,cumsum(schedule.step.val))
