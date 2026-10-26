%% test20_eclipse_roundtrip
% Verify the Eclipse deck exporter by reading the deck back with MRST's
% deckformat parser and reproducing every CPEDFM data structure from the
% deck alone.
%
% Case: the skewed, faulted corner-point grid with a slanted sealing
% fracture from CPEDFM-CPG test01.
%
% Checks:
%   A. GRID round-trip: the extended corner-point grid parses, and its
%      active cell count equals matrix cells + fracture cells.
%   B. NNC records: count and transmissibilities match G.nnc exactly
%      (after METRIC -> SI unit conversion by the parser).
%   C. MULT seals: the parsed MULTX/MULTY/MULTZ arrays reproduce the
%      face multipliers of transmultpEDFM_CPG on every logical-neighbor
%      sealed face.
%   D. EDITNNC seals: fault-juxtaposed sealed pairs are exported with
%      the correct multipliers.
%   E. Fracture-cell overrides: PORV and DEPTH match the MRST pore
%      volumes and centroid depths.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', '..', 'CPEDFM-CPG', 'startupCPEDFMCPG.m'));
addpath(fullfile(thisDir, '..'));
mrstModule add deckformat

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
failures = {};

%% Build the CPEDFM model (test01 configuration)
dims = [20, 10, 8];
grdecl = simpleGrdecl(dims, 0.12);
physdim = [200, 100, 40];
grdecl.COORD(1:3:end) = grdecl.COORD(1:3:end) * physdim(1);
grdecl.COORD(2:3:end) = grdecl.COORD(2:3:end) * physdim(2);
grdecl.COORD(3:3:end) = grdecl.COORD(3:3:end) * physdim(3);
grdecl.ZCORN = grdecl.ZCORN * physdim(3);

G = processGRDECL(grdecl);
G = computeGeometry(G);
G.rock = makeRock(G, 100*milli*darcy, 0.25);

fp = struct;
fp(1).points = [ 70  -5  -5; 110 105  -5; 150 105  45; 110  -5  45];
fp(1).aperture = 10*milli*meter;
fp(1).poro = 0.4;
fp(1).perm = 1*nano*darcy;

[G, fp] = EDFMGridCPG(G, fp, 'Tolerance', tol, 'plotgrid', false, ...
    'fracturelist', 1);
G = fracturematrixShaleNNC3D(G, tol);
Gp = pMatFracNNCsCPEDFM_CPG(G, tol);

nMat = G.Matrix.cells.num;
nFrac = Gp.cells.num - nMat;
multRef = transmultpEDFM_CPG(Gp, tol);

%% Export
[iW, ~, ~] = ind2sub(dims, G.Matrix.cells.indexMap(1:nMat));
wells = struct('name', {'INJ', 'PROD'}, ...
    'type', {'rate', 'bhp'}, ...
    'val',  {50*meter^3/day, 100*barsa}, ...
    'comp', {'WATER', 'OIL'}, ...
    'cells', {find(iW == 1, 4), find(iW == dims(1), 4)});

deckFile = fullfile(outDir, 'CPEDFM_TEST20.DATA');
info = writeEclipseCPEDFM(Gp, Gp.rock, grdecl, wells, deckFile, 'tol', tol);
fprintf('exported: %d NNCs, %d MULT seals, %d EDITNNC seals, dims %s\n', ...
    info.nNNC, info.nMult, info.nEditnnc, mat2str(info.dims));

%% Read back
deck = readEclipseDeck(deckFile);
deck = convertDeckUnits(deck);

%% Check A: grid round-trip
gr = deck.GRID;
gr.cartDims = deck.RUNSPEC.DIMENS;
G2 = processGRDECL(gr);
nActive = G2.cells.num;
fprintf('A. parsed grid: %d active cells (expected %d)\n', nActive, nMat + nFrac);
if nActive == nMat + nFrac
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('A: active cells %d vs %d', nActive, nMat+nFrac);
    fprintf('   FAIL\n');
end

%% Check B: NNC records
okB = isfield(deck.GRID, 'NNC');
if okB
    nnc = deck.GRID.NNC;
    okB = size(nnc, 1) == info.nNNC;
    if okB
        % match by (ijk1, ijk2) rows; parser preserves order
        dIJK = max(max(abs(nnc(:,1:6) - info.nncIJK)));
        relT = max(abs(nnc(:,7) - Gp.nnc.T) ./ max(abs(Gp.nnc.T), eps));
        fprintf('B. NNC: %d records, ijk diff %d, trans rel diff %.3e\n', ...
            size(nnc,1), dIJK, relT);
        okB = dIJK == 0 && relT < 1e-6;
    end
end
if okB
    fprintf('   PASS\n');
else
    failures{end+1} = 'B: NNC records mismatch';
    fprintf('B. NNC: FAIL\n');
end

%% Check C: MULT seals on logical-neighbor faces
names = {'MULTX', 'MULTY', 'MULTZ'};
okC = true; nChecked = 0; maxDiff = 0;
for q = 1:size(info.multRec, 1)
    d = info.multRec(q, 1);
    ijk = info.multRec(q, 2:4);
    val = info.multRec(q, 5);
    if ~isfield(deck.GRID, names{d}), okC = false; break; end
    arr = deck.GRID.(names{d});
    lin = sub2ind(info.dims, ijk(1), ijk(2), ijk(3));
    maxDiff = max(maxDiff, abs(arr(lin) - val));
    nChecked = nChecked + 1;
end
fprintf('C. MULT seals: %d/%d records verified, max diff %.3e\n', ...
    nChecked, size(info.multRec,1), maxDiff);
if okC && nChecked == size(info.multRec,1) && maxDiff < 1e-10
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('C: MULT mismatch (max diff %.3e)', maxDiff);
    fprintf('   FAIL\n');
end

%% Check D: EDITNNC seals
if info.nEditnnc > 0
    % MRST's deck parser does not implement EDITNNC (it is reported in
    % deck.UnhandledKeywords.EDIT); Eclipse and OPM Flow support it.
    % Verify the written records by parsing the deck text directly.
    txt = fileread(deckFile);
    blk = regexp(txt, 'EDITNNC\s*\n(.*?)\n/', 'tokens', 'once');
    okD = ~isempty(blk);
    if okD
        rows = sscanf(strrep(blk{1}, '/', ''), ' %d %d %d %d %d %d %f', [7, inf])';
        okD = size(rows,1) == info.nEditnnc && ...
            max(max(abs(rows(:,1:6) - info.editnncRec(:,1:6)))) == 0 && ...
            max(abs(rows(:,7) - info.editnncRec(:,7))) < 1e-10;
    end
    unh = isfield(deck, 'UnhandledKeywords') && ...
        isfield(deck.UnhandledKeywords, 'EDIT') && ...
        any(strcmp(deck.UnhandledKeywords.EDIT, 'EDITNNC'));
    fprintf('D. EDITNNC: %d records (text-verified; MRST parser skips keyword: %d)\n', ...
        info.nEditnnc, unh);
    if okD
        fprintf('   PASS\n');
    else
        failures{end+1} = 'D: EDITNNC mismatch';
        fprintf('   FAIL\n');
    end
else
    fprintf('D. EDITNNC: none required for this case (no fault-juxtaposed seals)\n');
end

%% Check E: fracture PORV and DEPTH overrides
pvRef = poreVolume(Gp, Gp.rock);
pvRef = pvRef(nMat+1:end);
zRef  = Gp.cells.centroids(nMat+1:end, 3);
% MRST merges EDIT-section keywords into deck.GRID
okE = isfield(deck.GRID, 'PORV') && isfield(deck.GRID, 'DEPTH');
if okE
    fracLin = sub2ind(info.dims, info.fracIJK(:,1), info.fracIJK(:,2), info.fracIJK(:,3));
    pvDeck = deck.GRID.PORV(fracLin);
    zDeck  = deck.GRID.DEPTH(fracLin);
    relPV = max(abs(pvDeck - pvRef) ./ max(pvRef, eps));
    dZ    = max(abs(zDeck - zRef));
    fprintf('E. fracture overrides: PORV rel diff %.3e, DEPTH max diff %.3e m\n', ...
        relPV, dZ);
    okE = relPV < 1e-6 && dZ < 1e-4;
end
if okE
    fprintf('   PASS\n');
else
    failures{end+1} = 'E: PORV/DEPTH mismatch';
    fprintf('E. fracture overrides: FAIL\n');
end

%% Summary
if isempty(failures)
    fprintf('\n=== test20: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test20: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test20 failed');
end
