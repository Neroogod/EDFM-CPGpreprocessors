%% test50_cmg_runnable_ow
% The point of preprocessorAll: an MRST case exported and run directly.
%
% A two-phase oil-water waterflood is built in MRST, simulated there, then
% exported as a complete IMEX dataset. The dataset is run by IMEX 2024.30 and
% the two answers are compared. This is the acceptance test the original
% preprocessor's README lists as missing ("an actual Eclipse / OPM Flow / CMG
% run of the exported deck").
%
% Checks:
%   A. buildFullDeck produces tables with no fidelity warnings for a fluid
%      the deck formats can represent exactly.
%   B. Sampled tables reproduce the MRST fluid functions to round-off.
%   C. The dataset and its include are written.
%   D. IMEX runs the dataset to Normal Termination with zero errors.
%   E. IMEX's material balance is unity and its cumulative water injection
%      matches the schedule exactly.
%   F. IMEX's oil production and average pressure agree with the MRST run.
%
% Running IMEX needs CMG on the machine; when it is absent, D-F report SKIP
% rather than failing, and A-C still run.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
mrstModule add preprocessor ad-core ad-props ad-blackoil

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

failures = {};

%% ---- The MRST case -----------------------------------------------------
gravity reset on
dims = [10, 10, 3];
physdim = [200, 200, 30];
G = cartGrid(dims, physdim);
G.nodes.coords(:,3) = G.nodes.coords(:,3) + 2000;   % genuine depth
G = computeGeometry(G);
rock = makeRock(G, 200*milli*darcy, 0.22);

fluid = initSimpleADIFluid('phases', 'WO', ...
    'mu',  [1, 5]*centi*poise, ...
    'rho', [1000, 800], ...
    'c',   [4.4e-5, 1e-4]/barsa, ...
    'n',   [2, 2], ...
    'pRef', 250*barsa);

model = GenericBlackOilModel(G, rock, fluid, ...
    'water', true, 'oil', true, 'gas', false);

state0 = initResSol(G, 250*barsa, [0.15, 0.85]);

rateInj = 100*meter^3/day;
bhpProd = 200*barsa;
cInj  = find(G.cells.centroids(:,1) < 25);
cProd = find(G.cells.centroids(:,1) > 175);
W = addWell([], G, rock, cInj,  'Type','rate','Val', rateInj, ...
    'Comp_i',[1 0], 'Name','INJ',  'refDepth', 2000);
W = addWell(W, G, rock, cProd, 'Type','bhp', 'Val', bhpProd, ...
    'Comp_i',[0 1], 'Name','PROD', 'refDepth', 2000);

totTime = 360*day;
schedule = simpleSchedule(rampupTimesteps(totTime, 90*day, 3), 'W', W);

fprintf('MRST case: %d cells, %d steps, injection %.4g m3/day\n', ...
    G.cells.num, numel(schedule.step.val), convertTo(rateInj, meter^3/day));

%% ---- Run it in MRST ----------------------------------------------------
[wsMRST, statesMRST] = simulateScheduleAD(state0, model, schedule);

pvAll = poreVolume(G, rock);
pMeanMRST = sum(statesMRST{end}.pressure .* pvAll) / sum(pvAll);
qoMRST = zeros(numel(wsMRST),1);
for k = 1:numel(wsMRST)
    qoMRST(k) = abs(wsMRST{k}(2).qOs);
end
cumOilMRST = sum(qoMRST .* schedule.step.val);
fprintf('MRST: pore-volume-weighted mean p = %.4f bar, cum oil = %.4f m3\n', ...
    pMeanMRST/barsa, cumOilMRST);

%% ---- Check A: build the deck ------------------------------------------
deck = buildFullDeck(model, state0, schedule);
fprintf('A. deck built: %s, %d PVT samples, %d sat samples, %d warnings\n', ...
    deck.meta.gridFamily, numel(deck.tab.p), numel(deck.tab.sw), numel(deck.warnings));
for q = 1:numel(deck.warnings), fprintf('     %s\n', deck.warnings{q}); end
okA = isempty(deck.warnings);
if okA
    fprintf('   PASS (no fidelity loss for this fluid)\n');
else
    failures{end+1} = sprintf('A: %d fidelity warnings', numel(deck.warnings));
    fprintf('   FAIL\n');
end

%% ---- Check B: tables reproduce the MRST fluid --------------------------
pT = deck.tab.p;
errBw  = max(abs(deck.tab.water.B  - 1./fluid.bW(pT))  ./ (1./fluid.bW(pT)));
errBo  = max(abs(deck.tab.oil.B    - 1./fluid.bO(pT))  ./ (1./fluid.bO(pT)));
errMuo = max(abs(deck.tab.oil.mu   - fluid.muO(pT))    ./ fluid.muO(pT));
sT = deck.tab.sw;
% A two-phase fluid exposes krO; krOW appears only in three-phase fluids.
if isfield(fluid, 'krOW'), krowFn = fluid.krOW; else, krowFn = fluid.krO; end
errKrw  = max(abs(deck.tab.krw  - fluid.krW(sT)));
errKrow = max(abs(deck.tab.krow - krowFn(1-sT)));
fprintf('B. table fidelity: Bw %.2e Bo %.2e muO %.2e krw %.2e krow %.2e\n', ...
    errBw, errBo, errMuo, errKrw, errKrow);
okB = max([errBw errBo errMuo errKrw errKrow]) < 1e-12;
if okB
    fprintf('   PASS\n');
else
    failures{end+1} = 'B: sampled tables do not reproduce the fluid functions';
    fprintf('   FAIL\n');
end

%% ---- Check C: write the dataset ---------------------------------------
datFile = fullfile(outDir, 'ow_waterflood.dat');
info = writeCMGDeck(deck, datFile, 'title', 'MRST_OW_WATERFLOOD');
okC = isfile(info.dataset) && isfile(info.include);
fprintf('C. dataset written: %s (%.1f kB) + include\n', ...
    info.dataset, localBytes(info.dataset)/1024);
if okC
    fprintf('   PASS\n');
else
    failures{end+1} = 'C: dataset or include missing';
    fprintf('   FAIL\n');
end

%% ---- Checks D-F: run IMEX and compare ----------------------------------
imex = localFindIMEX();
if isempty(imex)
    fprintf('D-F. SKIP: IMEX not found on this machine.\n');
else
    fprintf('D. running IMEX: %s\n', imex);
    outFile = strrep(info.dataset, '.dat', '.out');
    if isfile(outFile), delete(outFile); end
    cmd = sprintf('cd /d "%s" && "%s" -f "%s" -log', outDir, imex, info.dataset);
    [st, ~] = system(cmd);
    logTxt = '';
    if isfile(outFile), logTxt = fileread(outFile); end

    nErr = numel(regexp(logTxt, 'ERROR \(from subroutine', 'start'));
    normalEnd = contains(logTxt, 'Normal Termination');
    fprintf('   exit=%d, IMEX errors=%d, Normal Termination=%d\n', st, nErr, normalEnd);
    okD = normalEnd && nErr == 0;
    if okD
        fprintf('   PASS (the exported deck runs as it stands)\n');
    else
        failures{end+1} = sprintf('D: IMEX errors=%d normalEnd=%d', nErr, normalEnd);
        fprintf('   FAIL\n');
        err = regexp(logTxt, 'ERROR \(from subroutine[^\n]*\n[^\n]*\n[^\n]*', 'match');
        for q = 1:min(3, numel(err)), fprintf('     %s\n', err{q}); end
    end

    % E. material balance and injected volume
    mb = regexp(logTxt, 'Material Balances \(owg\):\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)', 'tokens', 'once');
    cumInj = localGrabCum(logTxt, 'Cumulative Injection');
    expInj = convertTo(rateInj, meter^3/day) * convertTo(totTime, day) / 1000; % MM3 = 1e3 m3
    okE = false;
    if ~isempty(mb)
        mbv = cellfun(@str2double, mb);
        relInj = abs(cumInj - expInj)/max(expInj, eps);
        fprintf('E. material balance %s | cum water inj %.4f vs schedule %.4f (rel %.2e)\n', ...
            mat2str(mbv, 5), cumInj, expInj, relInj);
        okE = all(abs(mbv - 1) < 1e-3) && relInj < 1e-6;
    else
        fprintf('E. material balance not found in the log\n');
    end
    if okE
        fprintf('   PASS\n');
    else
        failures{end+1} = 'E: material balance or injected volume off';
        fprintf('   FAIL\n');
    end

    % F. agreement with MRST
    cumOilIMEX = localGrabCum(logTxt, 'Cumulative Production') * 1000;  % -> m3
    pIMEX = localGrabPressure(logTxt);
    relOil = abs(cumOilIMEX - cumOilMRST)/max(cumOilMRST, eps);
    relP   = abs(pIMEX - pMeanMRST/barsa*100)/max(pMeanMRST/barsa*100, eps);
    fprintf('F. cum oil  IMEX %.3f m3 vs MRST %.3f m3   (rel %.3f)\n', ...
        cumOilIMEX, cumOilMRST, relOil);
    fprintf('   mean p    IMEX %.2f kPa vs MRST %.2f kPa (rel %.3f)\n', ...
        pIMEX, pMeanMRST/barsa*100, relP);
    okF = relOil < 0.05 && relP < 0.02;
    if okF
        fprintf('   PASS (exported deck reproduces the MRST answer)\n');
    else
        failures{end+1} = sprintf('F: oil rel %.3f, pressure rel %.3f', relOil, relP);
        fprintf('   FAIL\n');
    end
end

%% ---- Summary ------------------------------------------------------------
save(fullfile(outDir, 'test50_results.mat'), 'pMeanMRST', 'cumOilMRST');
if isempty(failures)
    fprintf('\n=== test50: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test50: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test50 failed');
end

%% ---- local helpers ------------------------------------------------------
function n = localBytes(f)
    d = dir(f); n = d.bytes;
end

function exe = localFindIMEX()
    exe = '';
    root = getenv('CMG_ROOT');
    if isempty(root), root = 'C:\Program Files (x86)\CMG'; end
    c = dir(fullfile(root, 'IMEX', '*', 'Win_x64', 'EXE', 'mx*.exe'));
    if ~isempty(c)
        [~, ix] = max([c.datenum]);
        exe = fullfile(c(ix).folder, c(ix).name);
    end
end

function v = localGrabCum(txt, label)
% Pull the water column from IMEX's "Field Total" block. Columns are
% Oil, Gas, Water; production takes water, injection takes water too.
    v = NaN;
    pat = [label '\s+([\d.eE+-]+|NA)\s+([\d.eE+-]+|NA)\s+([\d.eE+-]+|NA)'];
    t = regexp(txt, pat, 'tokens', 'once');
    if isempty(t), return; end
    if contains(label, 'Injection')
        v = str2double(t{3});
    else
        v = str2double(t{1});     % oil column
    end
end

function p = localGrabPressure(txt)
    p = NaN;
    t = regexp(txt, 'Average reservoir pressure[^:]*:\s*([\d.]+)', 'tokens', 'once');
    if ~isempty(t), p = str2double(t{1}); end
end
