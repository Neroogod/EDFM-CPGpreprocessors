%% test33_sconnect_is_live
% Proves the exported *SCONNECT records actually carry flow in CMG.
%
% Deck acceptance does not prove this: IMEX runs to Normal Termination even
% when a *SCONNECT record addresses a null block, so a silently ignored
% connection list would look identical to a working one.
%
% The test builds an EDFM model and completes a producer in a FRACTURE cell
% only. Every geometric connection of the padding region is zeroed by the
% export, so that cell's sole path to the reservoir is the matrix-fracture
% *SCONNECT list. Two decks are written:
%
%   live     the export as generated
%   severed  identical, with the *SCONNECT block commented out
%
% If the connections are live, the well drains the matrix and produces far
% more than the fracture's own pore volume. If they are being ignored, the
% severed and live runs are indistinguishable.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
cpgDir = fullfile(thisDir, '..', '..', 'CPEDFM-CPG');
if isfolder(cpgDir)
    addpath(cpgDir);
    if exist('startupCPEDFMCPG', 'file') == 2, startupCPEDFMCPG(); end
end
% hfm_v2 is the canonical EDFM module here; see the note in test32.
mrstModule add deckformat hfm_v2

outDir = fullfile(thisDir, 'output_cmg', 'sconnect');
if ~isfolder(outDir), mkdir(outDir); end

%% ---- EDFM model --------------------------------------------------------
tol = 1e-5;
dims = [12, 10, 2];
physdim = [120, 100, 20];

G0 = computeGeometry(cartGrid(dims, physdim));
G0.rock = makeRock(G0, 100*milli*darcy, 0.2);

frac = struct('points', [35 0 0; 85 100 0; 85 100 20; 35 0 20], ...
    'aperture', 5*milli*meter, 'poro', 0.4, 'perm', 1*darcy);

[G, frac] = EDFMshalegrid(G0, frac, 'Tolerance', tol, ...
    'plotgrid', false, 'fracturelist', 1); %#ok<ASGLU>
G = fracturematrixShaleNNC3D(G, tol);
G = prepareEDFMExport(G);

nMat = G.Matrix.cells.num;
nFrac = G.cells.num - nMat;
pvAll = poreVolume(G, G.rock);
pvMatrix = sum(pvAll(1:nMat));
pvFrac = sum(pvAll(nMat+1:end));

fprintf('matrix cells %d (pv %.2f m3), fracture cells %d (pv %.4f m3)\n', ...
    nMat, pvMatrix, nFrac, pvFrac);
fprintf('fracture pore volume is %.5f%% of the matrix\n', 100*pvFrac/pvMatrix);

%% ---- a producer completed in one fracture cell only -------------------
fracCell = nMat + 1;                       % first fracture cell
W = struct('name', 'FRACPROD', 'cells', fracCell, ...
    'type', 'bhp', 'val', 50*barsa, 'WI', 1e-11, ...
    'compi', [0 1], 'sign', -1, 'dir', 'z');

info = writeCMG(G, fullfile(outDir, 'live.dat'), 'rock', G.rock, ...
    'simulator', 'IMEX', 'cmgVersion', '2024.30', 'datum', 2000, ...
    'wells', W, 'timesteps', repmat(5*day, 6, 1), ...
    'hostDeck', 'runnable', 'title', 'SCONNECT_LIVE');

fprintf('exported: %d blocks, %d connections, well at UBA %s\n', ...
    info.counts.blocks, info.counts.connections, ...
    mat2str(reshape(info.counts.fractureCells > 0, 1, [])));

%% ---- severed twin: same deck, *SCONNECT commented out -----------------
liveInc = fullfile(outDir, 'live.inc');
sevInc  = fullfile(outDir, 'severed.inc');
sevDat  = fullfile(outDir, 'severed.dat');

txt = fileread(liveInc);
lines = strsplit(txt, newline);
inBlock = false;
for q = 1:numel(lines)
    L = strtrim(lines{q});
    if strcmp(L, '*SCONNECT')
        inBlock = true;
        lines{q} = ['** ' lines{q}];
        continue
    end
    if inBlock
        if isempty(L) || startsWith(L, '*') || startsWith(L, '**')
            inBlock = false;
        else
            lines{q} = ['** ' lines{q}];
        end
    end
end
fid = fopen(sevInc, 'w'); fprintf(fid, '%s', strjoin(lines, newline)); fclose(fid);

datTxt = fileread(fullfile(outDir, 'live.dat'));
datTxt = strrep(datTxt, 'live.inc', 'severed.inc');
datTxt = strrep(datTxt, 'SCONNECT_LIVE', 'SCONNECT_SEVERED');
fid = fopen(sevDat, 'w'); fprintf(fid, '%s', datTxt); fclose(fid);

fprintf('\nwrote %s\n     %s\n', fullfile(outDir, 'live.dat'), sevDat);
fprintf(['\nRun both in IMEX. If the connections are live, cumulative oil\n' ...
         'from live.dat must exceed the fracture pore volume (%.4f m3) by\n' ...
         'a wide margin, while severed.dat cannot.\n'], pvFrac);

save(fullfile(outDir, 'expected.mat'), 'pvFrac', 'pvMatrix');
fprintf('\n=== test33 export complete ===\n');
