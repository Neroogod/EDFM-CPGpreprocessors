%% test35_eclipse_matrix_export
% Eclipse matrix-grid export, verified by reading the deck back with MRST's
% own reader and reproducing the grid from the file alone.
%
% Writing a deck that merely looks right is easy; this checks that an
% independent parser recovers the same geometry and properties. Reading the
% .DATA rather than the include exercises the RUNSPEC/GRID skeleton too, so a
% deck that a simulator would reject cannot pass.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
mrstModule add deckformat

outDir = fullfile(thisDir, 'output_eclipse');
if ~isfolder(outDir), mkdir(outDir); end

failures = {};
report = @(ok, name, detail) fprintf('%s %-34s %s\n', ...
    subsref({'FAIL', 'PASS'}, substruct('{}', {ok + 1})), name, detail);

%% ---- A: corner-point ---------------------------------------------------
dims = [8 6 4]; physdim = [400 300 40]; datum = 2500;
grdecl = simpleGrdecl(dims, 0.12);
grdecl.COORD(1:3:end) = grdecl.COORD(1:3:end)*physdim(1);
grdecl.COORD(2:3:end) = grdecl.COORD(2:3:end)*physdim(2);
grdecl.COORD(3:3:end) = grdecl.COORD(3:3:end)*physdim(3) + datum;
grdecl.ZCORN          = grdecl.ZCORN*physdim(3) + datum;
grdecl.ACTNUM         = ones(prod(dims), 1);

Ga = computeGeometry(processGRDECL(grdecl));
[ai, aj, ak] = ind2sub(dims, (1:Ga.cells.num)');
permA = [100 + 5*ai + 3*aj + 20*ak, 100 + 5*ai + 3*aj + 20*ak, ...
         10 + 0.5*ai + 0.3*aj + 2*ak];
poroA = 0.18 + 0.001*ai + 0.002*ak;
rockA = makeRock(Ga, permA*milli*darcy, poroA);

ia = writeEclipse(Ga, fullfile(outDir, 'ecl_corner.DATA'), ...
    'rock', rockA, 'grdecl', grdecl, 'title', 'ECL_CORNER');

dA = readEclipseDeck(ia.dataset);
okA = isequal(reshape(dA.RUNSPEC.DIMENS, 1, []), dims) && ...
    max(abs(dA.GRID.COORD(:) - grdecl.COORD(:))) < 1e-6 && ...
    max(abs(dA.GRID.ZCORN(:) - grdecl.ZCORN(:))) < 1e-6 && ...
    max(abs(dA.GRID.PERMX(:) - permA(:,1))) < 1e-6 && ...
    max(abs(dA.GRID.PERMZ(:) - permA(:,3))) < 1e-6 && ...
    max(abs(dA.GRID.PORO(:)  - poroA(:)))   < 1e-9;
report(okA, 'A. corner-point round-trip', ...
    sprintf('%d blocks; COORD/ZCORN/PERMX/PERMZ/PORO recovered', prod(dims)));
if ~okA, failures{end+1} = 'A'; end

%% ---- B: uniform Cartesian ---------------------------------------------
cdims = [10 7 5]; cphys = [500 350 50]; cdatum = 1800;
Gb = computeGeometry(cartGrid(cdims, cphys));
[bi, ~, bk] = ind2sub(cdims, (1:Gb.cells.num)');
permB = [80 + 4*bi + 15*bk, 80 + 4*bi + 15*bk, 8 + 0.4*bi];
rockB = makeRock(Gb, permB*milli*darcy, 0.15 + 0.002*bi);

ib = writeEclipse(Gb, fullfile(outDir, 'ecl_cart.DATA'), ...
    'rock', rockB, 'datum', cdatum, 'title', 'ECL_CART');

dz = cphys(3)/cdims(3);
dB = readEclipseDeck(ib.dataset);
okB = isequal(reshape(dB.RUNSPEC.DIMENS, 1, []), cdims) && ...
    max(abs(dB.GRID.DXV(:) - cphys(1)/cdims(1))) < 1e-9 && ...
    max(abs(dB.GRID.DYV(:) - cphys(2)/cdims(2))) < 1e-9 && ...
    max(abs(dB.GRID.DZV(:) - dz)) < 1e-9 && ...
    numel(dB.GRID.TOPS) == Gb.cells.num && ...
    max(abs(dB.GRID.TOPS(:) - (cdatum + (bk-1)*dz))) < 1e-9 && ...
    max(abs(dB.GRID.PERMX(:) - permB(:,1))) < 1e-6;
report(okB, 'B. uniform Cartesian round-trip', ...
    sprintf('DXV/DYV/DZV %g/%g/%g; TOPS %g..%g over %d blocks', ...
    dB.GRID.DXV(1), dB.GRID.DYV(1), dB.GRID.DZV(1), ...
    min(dB.GRID.TOPS), max(dB.GRID.TOPS), numel(dB.GRID.TOPS)));
if ~okB, failures{end+1} = 'B'; end

%% ---- C: variable Cartesian --------------------------------------------
xs = [0 12 30 40 70 85 130]; ys = [0 20 45 100]; zs = [0 4 12 18];
Gc = computeGeometry(tensorGrid(xs, ys, zs));
rockC = makeRock(Gc, 120*milli*darcy, 0.21);

ic = writeEclipse(Gc, fullfile(outDir, 'ecl_vari.DATA'), ...
    'rock', rockC, 'datum', 2200, 'title', 'ECL_VARI');

dC = readEclipseDeck(ic.dataset);
okC = max(abs(dC.GRID.DXV(:) - diff(xs)')) < 1e-9 && ...
    max(abs(dC.GRID.DYV(:) - diff(ys)')) < 1e-9 && ...
    max(abs(dC.GRID.DZV(:) - diff(zs)')) < 1e-9 && ...
    numel(dC.GRID.TOPS) == Gc.cells.num;
report(okC, 'C. variable Cartesian round-trip', ...
    sprintf('DXV = [%s]', strtrim(num2str(dC.GRID.DXV(:)', '%g '))));
if ~okC, failures{end+1} = 'C'; end

%% ---- D: the grid rebuilds from the deck alone -------------------------
gBack = readGRDECL(ia.include);
Grebuilt = computeGeometry(processGRDECL(gBack));
rockBack = grdecl2Rock(gBack, Grebuilt.cells.indexMap);
rockBack.perm = rockBack.perm*milli*darcy;      % the deck is METRIC, so mD
relPV = abs(sum(poreVolume(Grebuilt, rockBack)) - ...
            sum(poreVolume(Ga, rockA))) / sum(poreVolume(Ga, rockA));
okD = Grebuilt.cells.num == Ga.cells.num && relPV < 1e-9;
report(okD, 'D. grid rebuilt from the deck', ...
    sprintf('%d cells, pore volume rel diff %.2e', Grebuilt.cells.num, relPV));
if ~okD, failures{end+1} = 'D'; end

%% ---- E: an unknown grid kind is refused -------------------------------
% All four kinds now have an Eclipse form (test36 covers fractures on
% structured and corner-point, test38 the connection list). What must still
% raise is a kind the writer does not recognise, so a future grid family
% cannot slip through as a silently wrong deck.
okE = false;
try
    writeEclipseModel(struct( ...
        'meta', struct('fractureModel', 'matrix'), ...
        'grid', struct('kind', 'radial')), ...
        fullfile(outDir, 'should_not_exist.grdecl'));
catch ME
    okE = strcmp(ME.identifier, 'writeEclipseModel:fracGridKind');
end
report(okE, 'E. unknown grid kind refused', ...
    'raises rather than writing a wrong deck');
if ~okE, failures{end+1} = 'E'; end

%% ---- summary -----------------------------------------------------------
fprintf('\n');
if isempty(failures)
    fprintf('=== test35: ALL CHECKS PASSED ===\n');
else
    fprintf('=== test35: %d FAILURES: %s ===\n', numel(failures), ...
        strjoin(failures, ', '));
end
