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
Tilt=0;
% Tilt=cellsize(1)*2;
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

%% Define resolved explicit reference grid
% This fine tensor grid is the fracture-conforming benchmark, not the
% embedded reservoir grid. The CPEDFM reservoir grid above is GRDECL CPG.
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
%% Define resolved explict grid
griddim_r = [225 225 2];
cellsize_r = physdim./griddim_r;

Gr = tensorGrid(0:cellsize_r(1):physdim(1),...
                      0:cellsize_r(2):physdim(2),...
                      0:cellsize_r(3):physdim(3));
Gr = computeGeometry(Gr);

Gr.rock=makeShaleRock(Gr,1*darcy,0.3);


%Set fracture perm
IJ_hFrac=ones(175-51+1,2)*113; IJ_hFrac(:,1)=51:175;
IJ_vFrac=ones(175-51+1,2)*113; IJ_vFrac(:,2)=51:175;
IJ_Fracs=[IJ_hFrac;IJ_vFrac];
FracCellIDs_2D = sub2ind(Gr.cartDims(1:2), IJ_Fracs(:,1), IJ_Fracs(:,2));

fractureCells=[];
NumGridLayer=griddim_r(1)*griddim_r(2);
for layerID = 1:griddim_r(3)
    fractureCells=[fractureCells; FracCellIDs_2D+(layerID-1)*NumGridLayer];
end
Gr.rock.perm(fractureCells) = fracplanes(1).perm; %Frac perm

figure,
show = false([Gr.cells.num, 1]);
show(fractureCells) = true;% Hide well cell
clf; 
plotCellData(Gr,convertTo(Gr.rock.perm,darcy),show,'EdgeColor', 'None','facealpha',0.5);
plotCellData(Gr,convertTo(Gr.rock.perm,darcy),~show,'EdgeColor', 'None','facealpha',0.0);
view(3); colorbar; axis equal tight

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


model_r = WaterModel(Gr,[], fluid);
model_r.operators = setupOperatorsTPFA(Gr, Gr.rock);



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

bcr = [];
bcr  = pside(bcr, Gr, 'LEFT', 10*barsa,'sat', s0);
bcr  = pside(bcr, Gr, 'RIGHT', 0*barsa,'sat', s0);

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
state_r = initResSol(Gr, pRef, s0);

totTime = 60*day;
nSteps = 10;
dt = rampupTimesteps(totTime, 30*day, nSteps);

schedule_pedfm = simpleSchedule(dt, 'bc', bc_pedfm);
schedule_cpedfm = simpleSchedule(dt, 'bc', bc_cpedfm);
schedule_edfm = simpleSchedule(dt, 'bc', bc_edfm);
schedule_r = simpleSchedule(dt, 'bc', bcr);






%% Simulate problem
[ws_pedfm, states_pedfm, reports_pedfm] = simulateScheduleAD(state_pedfm, model_pedfm, schedule_pedfm);

[ws_cpedfm, states_cpedfm, reports_cpedfm] = simulateScheduleAD(state_cpedfm, model_cpedfm, schedule_cpedfm);

[ws_edfm, states_edfm, reports_edfm] = simulateScheduleAD(state_edfm, model_edfm, schedule_edfm);

[ws_r, states_r, reports_r] = simulateScheduleAD(state_r, model_r, schedule_r, ...
   'afterStepFn', getPlotAfterStep(state_r, model_r, schedule_r));

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
plotToolbar(Gr, states_r)
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
% 
[Xr,Yr]=meshgrid(cellsize_r(1)/2:cellsize_r(1):physdim(1),cellsize_r(2)/2:cellsize_r(2):physdim(2));
nLayer=1;
NumGridLayer=griddim_r(1)*griddim_r(2);
startId=(nLayer-1)*NumGridLayer+1;
endId=(nLayer)*NumGridLayer;


pres_r=reshape(states_r{end,1}.pressure(startId:endId),[griddim_r(1),griddim_r(2)]);

[X,Y]=meshgrid(cellsize(1)/2:cellsize(1):physdim(1),cellsize(2)/2:cellsize(2):physdim(2));
nLayer=1;
NumGridLayer=griddim(1)*griddim(2);
startId=(nLayer-1)*NumGridLayer+1;
endId=(nLayer)*NumGridLayer;

pres_edfm=reshape(states_edfm{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);

subplot(1,4,1);
surf(X,Y,convertTo(pres_edfm, barsa()),'EdgeAlpha',0.1), view(130,20)
colormap('jet');axis tight;
title('EDFM')



pres_pedfm=reshape(states_pedfm{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);

subplot(1,4,2);
surf(X,Y,convertTo(pres_pedfm, barsa()),'EdgeAlpha',0.1), view(130,20)
colormap('jet');axis tight;
title('pEDFM')

pres_cpedfm=reshape(states_cpedfm{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);

subplot(1,4,3);
surf(X,Y,convertTo(pres_cpedfm, barsa()),'EdgeAlpha',0.1), view(130,20)
colormap('jet');axis tight;
title('CPEDFM')



subplot(1,4,4);
surf(Xr,Yr,convertTo(pres_r, barsa()),'EdgeAlpha',0.1), view(130,20)
colormap('jet');axis tight;
title('Reference')

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
IJ_hFrac_r=ones(225,2)*113; IJ_hFrac_r(:,1)=1:225;
CenterCellIDs_r = sub2ind(Gr.cartDims(1:2), IJ_hFrac_r(:,1), IJ_hFrac_r(:,2));

x_line=IJ_hFrac(:,1)*cellsize(1)-cellsize(1)/2;

pres_line_pedfm=states_pedfm{end,1}.pressure(CenterCellIDs);
pres_line_cpedfm=states_cpedfm{end,1}.pressure(CenterCellIDs);
pres_line_edfm=states_edfm{end,1}.pressure(CenterCellIDs);

x_line_r=IJ_hFrac_r(:,1)*cellsize_r(1)-cellsize_r(1)/2;
pres_line_r=states_r{end,1}.pressure(CenterCellIDs_r);
% CenterCellIDsPEBI=sort(find(G_pebi.cells.centroids(:,2)==4.5 & G_pebi.cells.centroids(:,3)==1.5));
% pres_line_pebi=states_pebi{end,1}.pressure(CenterCellIDsPEBI);
% x_line_pebi=(G_pebi.cells.centroids(CenterCellIDsPEBI,1));
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
plot(x_line_r, convertTo(pres_line_r, barsa()),'o','Color','k', ...
    'LineWidth', 2,'DisplayName','Reference');
hold off;

xlim([-1e-3 9]);
ylim([-1e-3 10]);

set(gca,'FontSize',15);
xlabel('X, m');
ylabel('Pressure, Bar');
legend; axis tight;
