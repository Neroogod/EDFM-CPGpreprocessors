cpedfmCPGDir = fileparts(fileparts(mfilename('fullpath'))); addpath(cpedfmCPGDir); startupCPEDFMCPG();

%% 
clc;clear;close all;

%% Load necessary modules
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
mrstModule add compositional;   % replace ad-blackoil with compositional module
mrstModule add ad-blackoil;     % AD blackoil solver
mrstModule add ad-core;         % AD core module
mrstModule add ad-props;        % AD properties

%% Set up a matrix grid
tol=1e-5;
Globals;

celldim = [30, 30, 5];
physdim = [60, 60, 10];
[G, grdecl] = makeExampleCornerPointGrid(celldim, physdim); %#ok<NASGU>
G.rock=makeRock(G,1*milli*darcy,0.1);


%% Set up fracture planes

% Case 7
% fracplanes(1).points = [20 0 0;
%                         38 60 10;
%                         36 60 20;
%                         18 0 10];


% Case 8
fracplanes(1).points = [10 0 0;
                        50 60 0;
                        50 60 10;
                        10 0 10];


                    
fracplanes(1).aperture = 10*milli*meter;
fracplanes(1).poro=0.4;
fracplanes(1).perm=1*nano*darcy;  



 
%% Test with random fracture

% rng(56471);  %(simple NF)
% tol4domain = tol*1e3;
% 
% 
% 
% fracplanes=[];
% 
% set1_1 = Field(DFN('dim',3,'n',10,'dip',45,'ddip',0,'dir',90,'ddir',0,'minl',6,...
%            'mu',7,'maxl',9,'bbx',[tol4domain,tol4domain,tol4domain,physdim(1)-tol4domain,physdim(2)-tol4domain,physdim(3)-tol4domain],...
%            'shape','q'),'Poly');
%  [set1,nonPlanarSets1,fracArea1] = processStochFracs(set1_1); 
% 
% numHFplanes = numel(fracplanes);
% numNFplanes = numel(set1); 
% totalNumFracPlanes = numHFplanes + numNFplanes;
% startIdx = numHFplanes + 1;
% for i=1:numNFplanes
%     idxGlobal = numHFplanes + i;
%     fracplanes(idxGlobal).points = set1{i}(1:end-1,:);
%     fracplanes(idxGlobal).aperture = 100*nano*meter; %1*micro*meter;
%     fracplanes(idxGlobal).poro=0.5;
%     if(mod(i,2)==false)%even NF idx
%        fracplanes(idxGlobal).perm=100*micro*darcy;%100*micro*darcy;%0.01*darcy;
%     else%odd NF idx
%        fracplanes(idxGlobal).perm=1*nano*darcy;%100*micro*darcy;%0.01*darcy;
%     end
% end
%  
 %%
 


% Plot fracture planes
figure,
plotfracongrid(G,fracplanes); % visualize to check before pre-process
view(30,45)

%% Construct fracture grid


[G,fracplanes]=EDFMGridCPG(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));

% Plot Fracture grid
figure;
plotGrid(cartGrid([1 1 1],physdim),'facealpha',0);
hold on;
plotGrid(G,G.Matrix.cells.num+1:G.cells.num);
axis tight equal
title('Fracture Grid')
view(30,45);

%% Fracture-Matrix non-Neighbouring Connections (NNC)

  %-Frac-Matrix NNCs
G = fracturematrixShaleNNC3D(G,tol);
  %-Frac-Frac NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);
G = pMatFracNNCsCPEDFM_CPG(G,tol);
% Plot NNCs

LR1=find(G.nnc.nncFaces(:,2)==1 | G.nnc.nncFaces(:,2)==2);
LRF1=G.nnc.nncFaces(LR1,1);
FB1=find(G.nnc.nncFaces(:,2)==3 | G.nnc.nncFaces(:,2)==4);
FBF1=G.nnc.nncFaces(FB1,1);
TB1=find(G.nnc.nncFaces(:,2)==5 | G.nnc.nncFaces(:,2)==6);
TBF1=G.nnc.nncFaces(TB1,1);
figure,
plotGrid(G,'facealpha',0,'edgealpha',0.1);
plotFaces(G,LRF1,'b','facealpha',0.5,'edgealpha',1);
plotFaces(G,FBF1,'r','facealpha',0.5,'edgealpha',1);
% plotFaces(G,TBF1,'k','facealpha',0.2,'edgealpha',1);
plotGrid(G,G.Matrix.cells.num+1:G.cells.num,'facealpha',0.1,'edgealpha',0);
axis tight equal;
title('Matrix-Fracture NNCs Processed')
view(30,45);


% LR2=find(G.nnc.nncFacesOld(:,2)==1 | G.nnc.nncFacesOld(:,2)==2);
% LRF2=G.nnc.nncFacesOld(LR2,1);
% FB2=find(G.nnc.nncFacesOld(:,2)==3 | G.nnc.nncFacesOld(:,2)==4);
% FBF2=G.nnc.nncFacesOld(FB2,1);
% TB2=find(G.nnc.nncFacesOld(:,2)==5 | G.nnc.nncFacesOld(:,2)==6);
% TBF2=G.nnc.nncFacesOld(TB2,1);
% 
% figure,
% plotfracongrid(cartGrid([1 1 1],physdim),fracplanes);
% plotGrid(G,G.nnc.originalcells(:,1),'facealpha',0,'edgealpha',1);
% plotFaces(G,LRF2,'b','facealpha',1,'edgealpha',1);
% plotFaces(G,FBF2,'r','facealpha',1,'edgealpha',1);
% plotFaces(G,TBF2,'g','facealpha',1,'edgealpha',1);
% plotGrid(G,G.Matrix.cells.num+1:G.cells.num,'facealpha',0.2);
% axis tight equal
% title('Matrix-Fracture NNCs Original')
% view(30,45);





