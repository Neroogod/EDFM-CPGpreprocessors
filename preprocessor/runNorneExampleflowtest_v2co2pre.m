%% Example: Poroelasticity simulation applied to the Norne case.
% 
% The simulation options are gathered in the opt structure. If opt=[] the
% simulation is run with the default options defined below
%
% **  Summary of the options ** 
%
% option 'norne_case' :
%
%     * 'full'       : 7392 cells
%     * 'mini Norne' :  605 cells
%
% option 'bc_case' :
%
%     * 'no displacement' : All nodes belonging to external faces have displacement
%                           equal to zero
%     * 'bottom fixed'    : The nodes that belong to the bottom have zero
%                           displacement, while a given pressure is imposed on
%                           the external faces that are not bottom faces.
%
% option 'method' :
%
%     * 'fully coupled'          : The mechanical and flow equations are solved fully coupled.
%     * 'fixed stress splitting' : The mechanical and flow equations are solved
%                                  sequentially using a fixed stress splitting
%
% option 'fluid_model' :
%
%     * 'blackoil'  : blackoil model is used for the fluid (gas is injected, see
%                     schedule below)
%     * 'oil water' : Two phase oil-water
%     * 'water'     : water model is used for the fluid

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

clc;clear;close all;

mrstModule add ad-core;
mrstModule add ad-props;
mrstModule add ad-blackoil;
mrstModule add co2lab-common co2lab-ve co2lab-spillpoint;
mrstModule add shale;
mrstModule add mrst-gui;
mrstModule add deckformat;
mrstModule add coarsegrid;

%% Load Norne grid

grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
grdecl = readGRDECL(grdecl);
fullCartDims = grdecl.cartDims;
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);
grdecl = cutGrdecl(grdecl, [10 25; 35 55; 1 22]); 

grdecl.ACTNUM = ones(size(grdecl.ACTNUM));
G = processGRDECL(grdecl);
G = G(1);
G = computeGeometry(G);

figure,
plotGrid(G, 'facecolor', 'none','edgeAlpha',0.1);  
axis tight
view(-64,70)

% Perm and poro values

perm = [grdecl.PERMX, grdecl.PERMY, grdecl.PERMZ];
rock.perm = perm(G.cells.indexMap, :);
rock.poro = max(grdecl.PORO(G.cells.indexMap), 0.1);

%% Setup rock parameters (for flow)

perm = [grdecl.PERMX, grdecl.PERMY, grdecl.PERMZ];
rock.perm = perm(G.cells.indexMap, :);
rock.poro = max(grdecl.PORO(G.cells.indexMap), 0.1);


figure,
plotGrid(G, 'facecolor', 'none','edgeAlpha',0.1);
axis tight
view(-64,70)


%% Export the grid to CMG
% writeCMG detects the grid family from G (processGRDECL -> *GRID *CORNER
% here) and writes two files:
%   norne_co2.inc  the reservoir description: *GRID/*COORD/*ZCORN/
%                  *PERMI/*PERMJ/*PERMK/*POR, straight from G and rock
%   norne_co2.dat  a minimal dataset that *INCLUDEs it - grid only, no
%                  fluid model, no tables, no schedule
%
% Open the .dat in Builder and populate PVT, rock-fluid, initialization and
% the RUN section there. Nothing in either file is invented: every number
% comes from the MRST model.
%
% Corner-point grids need grdecl passed as well - MRST discards COORD/ZCORN
% inside processGRDECL, so they cannot be recovered from G alone.

mrstModule add preprocessor

info = writeCMG(G, 'norne_co2.dat', ...
    'rock',   rock, ...
    'grdecl', grdecl, ...
    'title',  'NORNE_CO2_PRE');

fprintf('CMG dataset : %s\n', info.dataset);
fprintf('CMG include : %s\n', info.include);
fprintf('blocks      : %d (%d inactive)\n', info.counts.blocks, info.counts.inactive);

% To drop the grid into the existing GEM CO2 deck instead of opening the stub
% in Builder, emit the include on its own and replace that deck's
% COORD/ZCORN/ACTNUM/PERMI/PERMJ/PERMK/POR includes with one *INCLUDE line:
%
%   info = writeCMG(G, 'norne_co2.dat', 'rock', rock, 'grdecl', grdecl, ...
%                   'hostDeck', 'none');


