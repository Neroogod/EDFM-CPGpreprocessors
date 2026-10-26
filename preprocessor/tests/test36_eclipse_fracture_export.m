%% test36_eclipse_fracture_export
% Eclipse export of embedded-fracture models, verified by reading the deck
% back with MRST's parser.
%
% The CMG and Eclipse writers consume the same cmgModel struct, so this also
% checks that one model really does feed both: the connection counts, pore
% volume and block counts must agree between the two decks.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
cpgDir = fullfile(thisDir, '..', '..', 'CPEDFM-CPG');
if isfolder(cpgDir)
    addpath(cpgDir);
    if exist('startupCPEDFMCPG', 'file') == 2, startupCPEDFMCPG(); end
end
% hfm_v2 after the CPG startup; see the note in test32.
mrstModule add deckformat hfm_v2

outDir = fullfile(thisDir, 'output_eclipse');
if ~isfolder(outDir), mkdir(outDir); end

failures = {};
report = @(ok, name, detail) fprintf('%s %-32s %s\n', ...
    subsref({'FAIL', 'PASS'}, substruct('{}', {ok + 1})), name, detail);

%% ---- one shared embedded-fracture grid --------------------------------
tol = 1e-5;
G0 = computeGeometry(cartGrid([12 10 2], [120 100 20]));
G0.rock = makeRock(G0, 100*milli*darcy, 0.2);
frac = struct('points', [35 0 0; 85 100 0; 85 100 20; 35 0 20], ...
    'aperture', 5*milli*meter, 'poro', 0.4, 'perm', 1*nano*darcy);

[Gb, frac] = EDFMshalegrid(G0, frac, 'Tolerance', tol, ...
    'plotgrid', false, 'fracturelist', 1); %#ok<ASGLU>
Gb = fracturematrixShaleNNC3D(Gb, tol);

cases = struct( ...
    'name', {'EDFM', 'pEDFM', 'CPEDFM'}, ...
    'G', {prepareEDFMExport(Gb), ...
          preparePEDFMExport(pMatFracNNCs3D(Gb, tol)), ...
          pMatFracNNCsCPEDFM_CPG(Gb, tol)});

mdm = milli*darcy*meter;
ecl = centi*poise*meter^3/(day*barsa);

for q = 1:numel(cases)
    name = cases(q).name;
    Gq = cases(q).G;

    info = writeEclipse(Gq, fullfile(outDir, sprintf('frac_%s.DATA', lower(name))), ...
        'rock', Gq.rock, 'datum', 2000, 'tol', tol, 'title', ['FRAC_' name]);

    model = buildCMGModel(Gq, Gq.rock, 'datum', 2000, 'tol', tol);
    deck  = readEclipseDeck(info.dataset);

    nNNC = size(model.conn, 1);
    okCount = isfield(deck.GRID, 'NNC') && size(deck.GRID.NNC, 1) == nNNC;

    % Transmissibilities must survive the md.m -> cP.m3/day/bar conversion.
    okTrans = false;
    if okCount && nNNC > 0
        want = model.conn(:,7)*(mdm/ecl);
        got  = deck.GRID.NNC(:, 7);
        okTrans = max(abs(got - want)./max(abs(want), eps)) < 1e-8;
    end

    % PORV must reproduce the MRST pore volume exactly.
    pvMRST = sum(poreVolume(Gq, Gq.rock));
    okPV = isfield(deck.GRID, 'PORV') && ...
        abs(sum(deck.GRID.PORV) - pvMRST)/pvMRST < 1e-12;

    ok = okCount && okTrans && okPV;
    report(ok, sprintf('%s Eclipse export', name), ...
        sprintf('NNC %d/%d, PORV rel %.1e, blocks %d', ...
        size(deck.GRID.NNC, 1), nNNC, ...
        abs(sum(deck.GRID.PORV) - pvMRST)/pvMRST, model.meta.counts.blocks));
    if ~ok
        failures{end+1} = name; %#ok<AGROW>
        fprintf('     count=%d trans=%d pv=%d\n', okCount, okTrans, okPV);
    end
end

%% ---- the two writers agree on the same model --------------------------
Gc = cases(3).G;
mc = buildCMGModel(Gc, Gc.rock, 'datum', 2000, 'tol', tol);

iCMG = writeCMG(Gc, fullfile(outDir, 'agree_cmg.dat'), ...
    'rock', Gc.rock, 'datum', 2000, 'tol', tol, 'title', 'AGREE');
iECL = writeEclipse(Gc, fullfile(outDir, 'agree_ecl.DATA'), ...
    'rock', Gc.rock, 'datum', 2000, 'tol', tol, 'title', 'AGREE');

okAgree = isequal(iCMG.counts.blocks, iECL.counts.blocks) && ...
    isequal(iCMG.counts.connections, size(mc.conn, 1)) && ...
    iECL.nNNC == size(mc.conn, 1) && ...
    iECL.nMult == size(mc.mult, 1);
report(okAgree, 'CMG and Eclipse agree', ...
    sprintf('%d blocks, %d connections, %d seals in both', ...
    iCMG.counts.blocks, iECL.nNNC, iECL.nMult));
if ~okAgree, failures{end+1} = 'agreement'; end

%% ---- summary -----------------------------------------------------------
fprintf('\n');
if isempty(failures)
    fprintf('=== test36: ALL CHECKS PASSED ===\n');
else
    fprintf('=== test36: %d FAILURES: %s ===\n', numel(failures), ...
        strjoin(failures, ', '));
end
