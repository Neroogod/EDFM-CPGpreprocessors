%% test30_cmg_grid_export
% Phase 0-2 acceptance test for the MRST -> CMG preprocessor:
% grid-family detection, fracture-method detection, and matrix-only export
% for uniform Cartesian, variable Cartesian, corner-point and unstructured
% grids.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
cpgDir = fullfile(thisDir, '..', '..', 'CPEDFM-CPG');
if isfolder(cpgDir)
    addpath(cpgDir);
    if exist('startupCPEDFMCPG', 'file') == 2, startupCPEDFMCPG(); end
end
mrstModule add deckformat upr hfm

outDir = fullfile(thisDir, 'output_cmg');
if ~isfolder(outDir), mkdir(outDir); end

failures = {};
report = @(ok, name, detail) fprintf('%s %-34s %s\n', ...
    subsref({'FAIL', 'PASS'}, substruct('{}', {ok + 1})), name, detail);

%% ---------------------------------------------------------------- A: cart
Ga = computeGeometry(cartGrid([8 6 4], [80 60 20]));
Ga.rock = makeRock(Ga, 100*milli*darcy, 0.2);
ia = cmgDetectGrid(Ga);
okA = strcmp(ia.family, 'cart') && isequal(ia.dims, [8 6 4]) && ia.uniform;
report(okA, 'A. cartGrid -> *CART', sprintf('family=%s uniform=%d', ia.family, ia.uniform));
if ~okA, failures{end+1} = 'A'; end

%% ---------------------------------------------------------------- B: vari
% Deliberately uneven spacing, as meshHFsystem produces for Eagle Ford.
xs = [0 5 15 20 40 45 70];
ys = [0 10 25 60];
zs = [0 3 9 12];
Gb = computeGeometry(tensorGrid(xs, ys, zs));
Gb.rock = makeRock(Gb, 50*milli*darcy, 0.15);
ib = cmgDetectGrid(Gb);
okB = strcmp(ib.family, 'vari') && ~ib.uniform && ...
    numel(ib.di) == numel(xs) - 1 && ...
    max(abs(ib.di - diff(xs))) < 1e-9;
report(okB, 'B. tensorGrid -> *VARI', ...
    sprintf('family=%s di=[%s]', ib.family, num2str(ib.di, '%g ')));
if ~okB, failures{end+1} = 'B'; end

%% -------------------------------------------------------------- C: corner
grdecl = simpleGrdecl([6 5 3], 0.12);
Gc = computeGeometry(processGRDECL(grdecl));
Gc.rock = makeRock(Gc, 200*milli*darcy, 0.25);
ic = cmgDetectGrid(Gc);
okC = strcmp(ic.family, 'corner');
report(okC, 'C. processGRDECL -> *CORNER', sprintf('family=%s', ic.family));
if ~okC, failures{end+1} = 'C'; end

%% -------------------------------------------------------- D: unstructured
rng(42);
Gd = computeGeometry(pebi(triangleGrid(rand(40, 2)*100)));
Gd = makeLayeredGrid(Gd, 2);
Gd = computeGeometry(Gd);
Gd.rock = makeRock(Gd, 80*milli*darcy, 0.18);
id = cmgDetectGrid(Gd);
okD = strcmp(id.family, 'unstructured') && id.dims(1) == Gd.cells.num;
report(okD, 'D. pebi -> connection list', ...
    sprintf('family=%s cells=%d', id.family, Gd.cells.num));
if ~okD, failures{end+1} = 'D'; end

%% ------------------------------------------- E: fracture-method detection
% The empty-pMMfaces trap: prepareEDFMExport writes empty projection arrays
% onto a normalized EDFM grid, so an isfield test would call it CPEDFM.
tol = 1e-5;
dimsF = [12 10 2];
physdim = [120 100 20];
G0 = computeGeometry(cartGrid(dimsF, physdim));
G0.rock = makeRock(G0, 100*milli*darcy, 0.2);
frac = struct('points', [35 0 0; 85 100 0; 85 100 20; 35 0 20], ...
    'aperture', 5*milli*meter, 'poro', 0.4, 'perm', 1*nano*darcy);

haveHFM = exist('EDFMshalegrid', 'file') == 2;
if haveHFM
    [Gf, frac] = EDFMshalegrid(G0, frac, 'Tolerance', tol, ...
        'plotgrid', false, 'fracturelist', 1); %#ok<ASGLU>
    Gf = fracturematrixShaleNNC3D(Gf, tol);

    dPlain  = cmgDetectFractureModel(G0);
    dEdfm   = cmgDetectFractureModel(prepareEDFMExport(Gf));
    dPedfm  = cmgDetectFractureModel(preparePEDFMExport(pMatFracNNCs3D(Gf, tol)));
    dCpedfm = cmgDetectFractureModel(pMatFracNNCsCPEDFM_CPG(Gf, tol));

    okE = strcmp(dPlain.model, 'matrix') && strcmp(dEdfm.model, 'EDFM') && ...
        strcmp(dPedfm.model, 'pEDFM') && strcmp(dCpedfm.model, 'CPEDFM');
    report(okE, 'E. fracture-method detection', sprintf('%s / %s / %s / %s', ...
        dPlain.model, dEdfm.model, dPedfm.model, dCpedfm.model));
    if ~okE, failures{end+1} = 'E'; end
else
    fprintf('SKIP E. fracture-method detection      (hfm module unavailable)\n');
end

%% ------------------------------------------------------- F: matrix export
cases = {'cart', Ga, []; 'vari', Gb, []; 'corner', Gc, grdecl; 'unstr', Gd, []};
for q = 1:size(cases, 1)
    name = cases{q, 1};
    Gq = cases{q, 2};
    args = {};
    if ~isempty(cases{q, 3}), args = {'grdecl', cases{q, 3}}; end

    % 'runnable' adds a placeholder fluid model so the deck can be executed
    % by the IMEX acceptance harness. The default 'stub' writes grid only.
    dat = fullfile(outDir, sprintf('mx_%s.dat', name));
    iq = writeCMG(Gq, dat, 'simulator', 'IMEX', 'cmgVersion', '2024.30', ...
        'datum', 2000, 'hostDeck', 'runnable', args{:});

    pvMRST = sum(poreVolume(Gq, Gq.rock));
    okF = exist(iq.include, 'file') == 2 && exist(iq.dataset, 'file') == 2 && ...
        strcmp(iq.fractureModel, 'matrix');

    % For the connection-list route the deck must reproduce the fine pore
    % volume through VOLMOD, and carry one connection per internal face.
    detail = sprintf('blocks=%d conn=%d pv=%.6g m3', iq.counts.blocks, ...
        iq.counts.connections, pvMRST);
    if strcmp(iq.gridFamily, 'unstructured')
        mq = buildCMGModel(Gq, Gq.rock);
        pvCMG = sum(mq.props.porvol);       % *PORVOL carries it directly
        relPV = abs(pvCMG - pvMRST)/pvMRST;
        nInternal = nnz(all(Gq.faces.neighbors > 0, 2));
        okF = okF && relPV < 1e-12 && iq.counts.connections == nInternal;
        detail = sprintf('%s pvRel=%.2e conn/faces=%d/%d', detail, relPV, ...
            iq.counts.connections, nInternal);
    end

    report(okF, sprintf('F.%d %s matrix export', q, name), detail);
    if ~okF, failures{end+1} = sprintf('F.%d', q); end
end

%% ---------------------------------------------------------------- summary
fprintf('\n');
if isempty(failures)
    fprintf('=== test30: ALL CHECKS PASSED ===\n');
else
    fprintf('=== test30: %d FAILURES: %s ===\n', numel(failures), ...
        strjoin(failures, ', '));
end
