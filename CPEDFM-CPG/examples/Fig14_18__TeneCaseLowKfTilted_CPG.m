cpedfmCPGDir = fileparts(fileparts(mfilename('fullpath'))); addpath(cpedfmCPGDir); startupCPEDFMCPG();

%% This will generate Figures 10 to 13 in SPE-201243-PA when the fractures 
%  are sealing
%% PEDFM Validation case of + fracture of single phase flow
% This code is used to generate the plots in Figures 10 to 13 of 
% Olorode et al., 2020, which are based on Tene et al., 2017 Figures 3 & 4
clear 
clc
close all
opt = struct('nkr',       2, ...
    'shouldPlot', 1 );
%     opt = merge_options(opt, varargin{:});
%% Load necessary modules
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
mrstModule add ad-blackoil;     % AD blackoil solver
mrstModule add ad-core;         % AD core module
mrstModule add ad-props;        % AD properties
mrstModule add upr;        % AD upr
%% Basic Parameters
tol=1e-5;
physdim = [9 9 3]*meter;
griddim = [63 63 3];
cellsize = physdim./griddim;

%% Create Box Reservoir Grid
[G_matrix, grdecl] = makeExampleCornerPointGrid(griddim, physdim); %#ok<NASGU>
G_matrix.rock=makeShaleRock(G_matrix,1*darcy,0.3);

%% Create Fracture System
fracplanes = struct;
% Tilt=cellsize(1)*2;
Tilt=0;
fracplanes(1).points=...
    [2,4.5-Tilt,0;
     7,4.5+Tilt,0;
     7,4.5+Tilt,3;
     2,4.5-Tilt,3];
% Tilt=0;
Tilt=cellsize(1)*2;
fracplanes(2).points=...
   [4.5+Tilt,2,0;
    4.5-Tilt,7,0;
    4.5-Tilt,7,3;
    4.5+Tilt,2,3];

% fracplanes(2).points=...
%    [4.3,6.8,0;
%     4.8,2.15,0;
%     4.8,2.15,3;
%     4.5-Tilt,2,3];


  fracPerm= 1e-9*darcy;
for i=1:numel(fracplanes)
     fracplanes(i).aperture = 0.003; %1*micro*meter;
     fracplanes(i).poro = 0.3;
     fracplanes(i).perm =fracPerm;
end


 %% Create Wells
wells = struct;

%% visualize to check before pre-process
checkIfCoplanar(fracplanes)

G=G_matrix;

%% EDFM PreProcessing
[G,fracplanes]=EDFMGridCPG(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
  %-Frac-Matrix NNCs
G = fracturematrixShaleNNC3D(G,tol);
  %-Frac-Frac NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);
  %-Well-Fracs NNCs
[G,wells] = wellfractureNNCs3D(G,fracplanes,wells,tol);

  % OMO: Projection-based NNCs % Switch from EDFM to pEDFM
G_edfm = G; %EDFM
% G_pedfm=G;
G_pedfm = pMatFracNNCsPEDFM_CPG(G,tol); % pEDFM on CPG
G_cpedfm = pMatFracNNCsCPEDFM_CPG(G,tol); %pEDFM
figure,
plotFracSystem(G_matrix,fracplanes,wells);

%% Define resolved explict grid
xmaxFS = physdim;
data.HF={[fracplanes(1).points(1,1),fracplanes(1).points(1,2);fracplanes(1).points(2,1),fracplanes(1).points(2,2)],...
    [fracplanes(2).points(1,1),fracplanes(2).points(1,2);fracplanes(2).points(2,1),fracplanes(2).points(2,2)]};
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
protD = {@(p) 0.0003*ones(size(p,1),1)};
% ApertureSD=FracAperture;
% G2D = pebiGrid2D(dx, xmax(1:2), ...
%     'cellConstraints', fractureScaled2D, ... % Fractures
%     'protLayer',true,....
%     'protD',protD,...
%     'CCRefinement'   , true , ... % Refine fractures
%     'CCEps',1/50,...
%     'CCFactor' , 0.05  );    % Rel fracture cell size
gS = [0.01,0.01];
pdims=[1,1];
G2D = compositePebiGrid2D(gS, pdims,'cellConstraints',fractureScaled2D,'CCFactor',1/8,'protLayer',true,'protD',protD);

G2DFD=G2D;
G2DFD.nodes.coords=G2DFD.nodes.coords.*xmaxFS(1:2);
figure,
plotGrid(G2DFD,G2DFD.cells.tag, 'faceColor', [1,1,1]*0.8,'LineWidth',1); axis tight
plotGrid(G2DFD, 'faceColor', 'none', 'edgeAlpha', 0.1);

%% Make 3D layered grid
layers = diff(linspace(0, xmaxFS(3), griddim(3)+1));
G = makeLayeredGrid(G2DFD, layers);
G = computeGeometry(G);

% Convert the fracture to full dim 3D
FractureFullDim=cellfun(@(x) (x.*xmaxFS(1:2)), fractureScaled2D, 'UniformOutput', false);

FractureFulldim3D=[];
for i=1:numel(FractureFullDim)
    Curfracdata=FractureFullDim{1,i};
    FractureFulldim3D=[FractureFulldim3D {[Curfracdata(1,:) 0;Curfracdata(2,:) 0;Curfracdata(2,:) xmaxFS(3);Curfracdata(1,:) xmaxFS(3)]}];
end
G.cells.tag = repmat(G2DFD.cells.tag,griddim(3), 1);
G.cells.HF=G.cells.tag;
fracture_cells=G.cells.tag;

leftFace=find(G.faces.centroids(:,1)==0);
rightFace=find(G.faces.centroids(:,1)==physdim(1));
% Make rock
MatrixPerm=1*darcy;
MatrixPoro=0.3;
rock = makeRock(G, MatrixPerm, MatrixPoro);
rock.perm(fracture_cells) = fracPerm; % Fracture permeability
% FracPoreVolm=3.28;
% rock.poro(fracture_cells) =FracPoreVolm./G.cells.volumes(fracture_cells);
rock.poro(fracture_cells)=0.3;

G_pebi=G;
% Plot Grid
figure,
plotGrid(G_pebi, G_pebi.cells.HF , 'faceColor', 'y', 'edgeColor', 'none');
plotGrid(G_pebi,'faceColor', 'none', 'edgeAlpha', 0.1);
axis tight;
%% Simulate problem


    
%% Define fluid properties
% Define a single-phase fluid model without capillarity.

fluid = initSimpleADIFluid('phases','W',       ... % Fluid phase: water
                           'mu',  1*centi*poise, ... % Viscosity
                           'rho', 1000,          ... % Surface density [kg/m^3]
                           'c',   1e-4/barsa,    ... % Fluid compressibility
                           'cR',  1e-5/barsa     ... % Rock compressibility
                           );

%% Define single-phase compressible flow model with pEDFM operators
% We define a three-phase black-oil model without dissolved gas or vaporized
% oil. This is done by first instantiating the blackoil model, and then
% manually passing in the internal transmissibilities and the topological
% neighborship from the embedded fracture grid.

gravity reset off
model_pedfm = WaterModel(G_pedfm,[], fluid);
model_pedfm.operators = setupCPEDFMOpsTPFA_CPG(G_pedfm, G_pedfm.rock, tol);

model_cpedfm = WaterModel(G_cpedfm,[], fluid);
model_cpedfm.operators = setupCPEDFMOpsTPFA_CPG(G_cpedfm, G_cpedfm.rock, tol);

model_edfm = WaterModel(G_edfm,[], fluid);
model_edfm.operators = setupShaleEDFMOpsTPFA(G_edfm, G_edfm.rock, tol);

model_pebi = WaterModel(G_pebi,rock, fluid);



s0 = (1);
bc_pedfm = [];
bc_pedfm  = pside(bc_pedfm, G_pedfm, 'LEFT', 10*barsa,'sat', s0);
bc_pedfm  = pside(bc_pedfm, G_pedfm, 'RIGHT', 0*barsa,'sat', s0);

bc_cpedfm = [];
bc_cpedfm  = pside(bc_cpedfm, G_cpedfm, 'LEFT', 10*barsa,'sat', s0);
bc_cpedfm  = pside(bc_cpedfm, G_cpedfm, 'RIGHT', 0*barsa,'sat', s0);

bc_edfm = [];
bc_edfm  = pside(bc_edfm, G_edfm, 'LEFT', 10*barsa,'sat', s0);
bc_edfm  = pside(bc_edfm, G_edfm, 'RIGHT', 0*barsa,'sat', s0);

bc_pebi = [];
% bc  = pside(bc, G2DFD, 'left', 1e6, 'sat', 1); 
% bc  = pside(bc, G_pebi, 'right', 1e5, 'sat', 1);
bc_pebi = addBC(bc_pebi, leftFace, 'pressure', 10*barsa, 'sat', s0);
bc_pebi = addBC(bc_pebi, rightFace, 'pressure', 0*barsa*barsa, 'sat', s0);

clf, figure
plotGrid(G_matrix,'FaceColor', 'none'); view(3);
plotFaces(G_matrix, bc_pedfm.face(bc_pedfm.value==0), 'b','facealpha',0.7);
plotFaces(G_matrix, bc_pedfm.face(bc_pedfm.value==1000000), 'r','facealpha',0.7);
hold on;
% 
% clf, figure
% plotGrid(G_matrix,'FaceColor', 'none'); view(3);
% plotFaces(G_matrix, bc_pedfm.value==0, 'r','facealpha',0.7);
% plotFaces(G_matrix, bc_pedfm.value==1000000, 'b','facealpha',0.7);
% hold on;

hf={};
for i = 1:length(fracplanes)
    X=fracplanes(i).points(:,1);
    Y=fracplanes(i).points(:,2);
    Z=fracplanes(i).points(:,3);
    hf{end+1}=fill3(X,Y,Z,'y');
    set(hf{end},'facealpha',0.7);
end
set(gca, 'CameraPosition', [-24.1107  -33.7559  -19.0289]);


%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water and oil initially. We also set up a simple-time step strategy that
% ramps up gradually towards 30 day time-steps.

s0 = (1);
pRef=5*barsa;
state_pedfm  = initResSol(G_pedfm, pRef, s0);
state_cpedfm = initResSol(G_cpedfm, pRef, s0);
state_edfm = initResSol(G_edfm, pRef, s0);
state_pebi = initResSol(G_pebi, pRef, s0);

totTime = 60*day;
nSteps = 10;
dt = rampupTimesteps(totTime, 30*day, nSteps);

schedule_pedfm = simpleSchedule(dt, 'bc', bc_pedfm);
schedule_cpedfm = simpleSchedule(dt, 'bc', bc_cpedfm);
schedule_edfm = simpleSchedule(dt, 'bc', bc_edfm);
schedule_pebi = simpleSchedule(dt, 'bc', bc_pebi);






%% Simulate problem
[ws_pedfm, states_pedfm, reports_pedfm] = simulateScheduleAD(state_pedfm, model_pedfm, schedule_pedfm);

[ws_cpedfm, states_cpedfm, reports_cpedfm] = simulateScheduleAD(state_cpedfm, model_cpedfm, schedule_cpedfm);

[ws_edfm, states_edfm, reports_edfm] = simulateScheduleAD(state_edfm, model_edfm, schedule_edfm);

[ws_pebi, states_pebi, reports_pebi] = simulateScheduleAD(state_pebi, model_pebi, schedule_pebi);

%% plotting
figure, 
plotToolbar(G_pedfm, states_pedfm)
view(40,30);
axis tight equal;

figure, 
plotToolbar(G_cpedfm, states_cpedfm)
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


%% Pressrue plot
% figure,
% plotCellData(G, convertTo(states{end,1}.pressure, barsa()), 'EdgeColor', 'k');
% xlabel('x'), ylabel('y'), zlabel('Depth');
% view(3); axis tight;colormap('jet');
% h=colorbar; caxis([0 10]); 
% h.Label.String = 'Pressure (bar)';
% set(gca, 'CameraPosition', [-24.1107  -33.7559  -19.0289]);


%% Pressure contour comparsion

% Figure and figure settings
figure('Position',[200 460 1200 400]);

[X,Y]=meshgrid(cellsize(1)/2:cellsize(1):physdim(1),cellsize(2)/2:cellsize(2):physdim(2));
nLayer=1;
NumGridLayer=griddim(1)*griddim(2);
startId=(nLayer-1)*NumGridLayer+1;
endId=(nLayer)*NumGridLayer;

pres_pedfm=reshape(states_pedfm{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);

subplot(1,4,1);
surf(X,Y,convertTo(pres_pedfm, barsa()),'EdgeAlpha',0.1), view(110,20)
colormap('jet');axis tight;
title('pEDFM solution pressure')

pres_cpedfm=reshape(states_cpedfm{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);

subplot(1,4,2);
surf(X,Y,convertTo(pres_cpedfm, barsa()),'EdgeAlpha',0.1), view(110,20)
colormap('jet');axis tight;
title('CPEDFM solution pressure')

pres_edfm=reshape(states_edfm{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);

subplot(1,4,3);
surf(X,Y,convertTo(pres_edfm, barsa()),'EdgeAlpha',0.1), view(110,20)
colormap('jet');axis tight;
title('EDFM solution pressure')

startId=1;
endId=G_pebi.cells.num/3;
% pres_pebi=reshape(states_pebi{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);
% pres_pebi2=reshape(states_pebi{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);

subplot(1,4,4);
%surf(X,Y,convertTo(pres_pebi, barsa()),'EdgeAlpha',0.1), view(130,20)
T = delaunay(G_pebi.cells.centroids(1:endId,1),G_pebi.cells.centroids(1:endId,2));
% figure;
trisurf(T,G_pebi.cells.centroids(1:endId,1),G_pebi.cells.centroids(1:endId,2),convertTo(states_pebi{end,1}.pressure(1:endId), barsa()),'EdgeAlpha',0.1), view(20,20)
colormap('jet');axis tight;
title('PEBI solution pressure')

% Topface=find(G_pebi.faces.centroids(:,3)==0);
% TopCell=1:G2DFD.cells.num;
% TopPres=states_pebi{12, 1}.pressure (TopCell);
% X=G_pebi.faces.centroids(Topface,1);
% Y=G_pebi.faces.centroids(Topface,2);
% surf(X,Y,convertTo(TopPres, barsa()),'EdgeAlpha',0.1), view(130,20)
%% Pressure overline comparsion
IJ_hFrac=ones(griddim(1),2)*((griddim(1)-1)/2+1); IJ_hFrac(:,1)=1:griddim(1);
% CenterCellIDs = sub2ind(G_pedfm.cartDims(1:2), IJ_hFrac(:,1), IJ_hFrac(:,2));
CenterCellIDs=(griddim(1)*(griddim(2)-1)/2+1:griddim(1)*(griddim(2)-1)/2+griddim(1))';
CenterCellIDs=CenterCellIDs-2*griddim(2);
% IJ_hFrac_r=ones(225,2)*113; IJ_hFrac_r(:,1)=1:225;
% CenterCellIDs_r = sub2ind(Gr.cartDims(1:2), IJ_hFrac_r(:,1), IJ_hFrac_r(:,2));

x_line=IJ_hFrac(:,1)*cellsize(1)-cellsize(1)/2;

pres_line_pedfm=states_pedfm{end,1}.pressure(CenterCellIDs);
pres_line_cpedfm=states_cpedfm{end,1}.pressure(CenterCellIDs);
pres_line_edfm=states_edfm{end,1}.pressure(CenterCellIDs);

% x_line_r=IJ_hFrac_r(:,1)*cellsize_r(1)-cellsize_r(1)/2;
% pres_line_r=states_r{end,1}.pressure(CenterCellIDs_r);
CenterCellIDsPEBI=sort(find(G_pebi.cells.centroids(:,2)==4.5 & G_pebi.cells.centroids(:,3)==1.5));
pres_line_pebi=states_pebi{end,1}.pressure(CenterCellIDsPEBI);
x_line_pebi=(G_pebi.cells.centroids(CenterCellIDsPEBI,1));
figure('rend','painters','pos',[10 10 640 500]);
title('Pressure solution over the horizontal center line');
% plot(x_line_r, convertTo(pres_line_r, barsa()),'k-', 'LineWidth', 2,'DisplayName','Explicit');
% hold on;
plot(x_line, convertTo(pres_line_pedfm, barsa()),'-','Color','r', ...
    'LineWidth', 2,'DisplayName','pEDFM');
hold on
plot(x_line, convertTo(pres_line_cpedfm, barsa()),'--','Color','g', ...
    'LineWidth', 2,'DisplayName','CPEDFM');
hold on
plot(x_line, convertTo(pres_line_edfm, barsa()),'-*','Color','b', ...
    'LineWidth', 2,'DisplayName','EDFM');
hold on
plot(x_line_pebi, convertTo(pres_line_pebi, barsa()),'o','Color','k', ...
    'LineWidth', 2,'DisplayName','Reference');
hold off;

xlim([-1e-3 9]);
ylim([-1e-3 10]);

set(gca,'FontSize',15);
xlabel('X, m');
ylabel('Pressure, Bar');
legend; axis tight;
