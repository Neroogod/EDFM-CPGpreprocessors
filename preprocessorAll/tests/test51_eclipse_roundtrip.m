%% test51_eclipse_roundtrip
% The Eclipse side of the full-deck export: write a complete .DATA deck from
% an MRST case, read it back with MRST's own deckformat parser, and check
% that every section survived with the values it went in with.
%
% No Eclipse or OPM Flow binary is available on this machine, so the parser
% is the strongest available check on the Eclipse path. The CMG path is
% verified by an actual simulator run in test50; that asymmetry is stated in
% the README rather than glossed over.
%
% Checks:
%   A. Deck and grid file are written.
%   B. The deck parses, and RUNSPEC carries the right phases and dimensions.
%   C. PROPS: PVTW, PVDO, SWOF and DENSITY round-trip to the sampled values.
%   D. SOLUTION: PRESSURE and SWAT round-trip to state0.
%   E. SCHEDULE: one COMPDAT row per perforation with MRST's own connection
%      factors, and TSTEP summing to the schedule's total time.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
mrstModule add preprocessor deckformat ad-core ad-props ad-blackoil

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

failures = {};

%% ---- The MRST case -----------------------------------------------------
gravity reset on
dims = [8, 8, 3];
G = cartGrid(dims, [160, 160, 30]);
G.nodes.coords(:,3) = G.nodes.coords(:,3) + 1500;
G = computeGeometry(G);
rock = makeRock(G, 150*milli*darcy, 0.25);

fluid = initSimpleADIFluid('phases', 'WO', ...
    'mu',  [0.8, 3]*centi*poise, ...
    'rho', [1020, 780], ...
    'c',   [4.0e-5, 8e-5]/barsa, ...
    'n',   [2, 3], ...
    'pRef', 180*barsa);

model = GenericBlackOilModel(G, rock, fluid, ...
    'water', true, 'oil', true, 'gas', false);
state0 = initResSol(G, 180*barsa, [0.2, 0.8]);

cInj  = find(G.cells.centroids(:,1) < 25);
cProd = find(G.cells.centroids(:,1) > 135);
W = addWell([], G, rock, cInj,  'Type','rate','Val', 60*meter^3/day, ...
    'Comp_i',[1 0], 'Name','I1', 'refDepth', 1500);
W = addWell(W, G, rock, cProd, 'Type','bhp','Val', 150*barsa, ...
    'Comp_i',[0 1], 'Name','P1', 'refDepth', 1500);
totTime = 300*day;
schedule = simpleSchedule(rampupTimesteps(totTime, 75*day, 2), 'W', W);

deck = buildFullDeck(model, state0, schedule);
fprintf('deck built: %d cells, %d steps, %d warnings\n', ...
    G.cells.num, deck.meta.nSteps, numel(deck.warnings));

%% ---- Check A: write ----------------------------------------------------
dataFile = fullfile(outDir, 'OW_ECL.DATA');
info = writeEclipseDeck(deck, dataFile, 'title', 'MRST_OW_ECL');
okA = isfile(info.deck) && isfile(info.grid);
fprintf('A. written: %s (%.1f kB) + %s\n', info.deck, localBytes(info.deck)/1024, info.grid);
if okA, fprintf('   PASS\n'); else
    failures{end+1} = 'A: deck or grid file missing'; fprintf('   FAIL\n');
end

%% ---- Check B: parse ----------------------------------------------------
okB = false; dk = [];
try
    dk = readEclipseDeck(info.deck);
    dk = convertDeckUnits(dk);
    nx = dk.RUNSPEC.cartDims;
    hasW = isfield(dk.RUNSPEC,'WATER') && dk.RUNSPEC.WATER;
    hasO = isfield(dk.RUNSPEC,'OIL')   && dk.RUNSPEC.OIL;
    hasG = isfield(dk.RUNSPEC,'GAS')   && dk.RUNSPEC.GAS;
    okB = isequal(nx(:)', dims) && hasW && hasO && ~hasG;
    fprintf('B. parsed: dims %s, phases W=%d O=%d G=%d\n', mat2str(nx(:)'), hasW, hasO, hasG);
catch me
    fprintf('B. parse FAILED: %s\n', me.message);
end
if okB, fprintf('   PASS\n'); else
    failures{end+1} = 'B: deck did not parse or RUNSPEC wrong'; fprintf('   FAIL\n');
end

%% ---- Check C: PROPS ----------------------------------------------------
okC = false;
if ~isempty(dk) && isfield(dk, 'PROPS')
    P = dk.PROPS;
    % PVTW: [pRef Bw Cw muW Cvw] in SI after convertDeckUnits
    pvtw = P.PVTW(1,:);
    dBw = abs(pvtw(2) - deck.tab.water.Bref)/deck.tab.water.Bref;
    dMuw = abs(pvtw(4) - deck.tab.water.muRef)/deck.tab.water.muRef;
    % PVDO: columns [p Bo muO]
    pvdo = P.PVDO{1};
    dBo = max(abs(pvdo(:,2) - deck.tab.oil.B)./deck.tab.oil.B);
    dMuo = max(abs(pvdo(:,3) - deck.tab.oil.mu)./deck.tab.oil.mu);
    % SWOF: [sw krw krow pcow]
    swof = P.SWOF{1};
    dKrw  = max(abs(swof(:,2) - deck.tab.krw));
    dKrow = max(abs(swof(:,3) - deck.tab.krow));
    % DENSITY: [oil water gas]
    dens = P.DENSITY(1,:);
    dRho = max(abs(dens(1:2) - [deck.tab.rhoS(2), deck.tab.rhoS(1)]));
    fprintf('C. PROPS: Bw %.2e muW %.2e Bo %.2e muO %.2e krw %.2e krow %.2e rho %.2e\n', ...
        dBw, dMuw, dBo, dMuo, dKrw, dKrow, dRho);
    okC = max([dBw dMuw dBo dMuo]) < 1e-6 && max([dKrw dKrow]) < 1e-6 && dRho < 1e-6;
end
if okC, fprintf('   PASS\n'); else
    failures{end+1} = 'C: PROPS did not round-trip'; fprintf('   FAIL\n');
end

%% ---- Check D: SOLUTION -------------------------------------------------
okD = false;
if ~isempty(dk) && isfield(dk, 'SOLUTION')
    S = dk.SOLUTION;
    nTot = prod(dims);
    dP  = max(abs(S.PRESSURE(1:nTot) - state0.pressure))/mean(state0.pressure);
    dSw = max(abs(S.SWAT(1:nTot) - state0.s(:,1)));
    fprintf('D. SOLUTION: pressure %.2e, SWAT %.2e\n', dP, dSw);
    okD = dP < 1e-6 && dSw < 1e-9;
end
if okD, fprintf('   PASS\n'); else
    failures{end+1} = 'D: SOLUTION did not round-trip'; fprintf('   FAIL\n');
end

%% ---- Check E: SCHEDULE -------------------------------------------------
okE = false;
if ~isempty(dk) && isfield(dk, 'SCHEDULE')
    txt = fileread(info.deck);
    nPerfExpected = numel(cInj) + numel(cProd);
    compBlock = regexp(txt, 'COMPDAT(.*?)\n/', 'tokens', 'once');
    nRows = 0; cfMaxErr = NaN;
    if ~isempty(compBlock)
        rows = regexp(compBlock{1}, '''(\w+)''\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+''OPEN''\s+1\*\s+([\d.eE+-]+)', 'tokens');
        nRows = numel(rows);
        ecl = centi*poise*meter^3/(day*barsa);
        wiAll = [deck.wells(1).W.wi];
        cf = cellfun(@(r) str2double(r{6}), rows);
        if numel(cf) == numel(wiAll)
            cfMaxErr = max(abs(cf(:) - wiAll(:)/ecl) ./ max(abs(wiAll(:)/ecl), eps));
        end
    end
    tstepTotal = sum(deck.steps.dtDays);
    dt = 0;
    if isfield(dk.SCHEDULE, 'step') && isfield(dk.SCHEDULE.step, 'val')
        dt = sum(convertTo(dk.SCHEDULE.step.val, day));
    end
    fprintf('E. SCHEDULE: COMPDAT rows %d (expected %d), CF rel.err %.2e, TSTEP %.4f vs %.4f days\n', ...
        nRows, nPerfExpected, cfMaxErr, dt, tstepTotal);
    okE = nRows == nPerfExpected && (isnan(cfMaxErr) || cfMaxErr < 1e-6) && ...
          abs(dt - tstepTotal) < 1e-6*max(tstepTotal,1);
end
if okE, fprintf('   PASS\n'); else
    failures{end+1} = 'E: SCHEDULE did not round-trip'; fprintf('   FAIL\n');
end

%% ---- Summary ------------------------------------------------------------
if isempty(failures)
    fprintf('\n=== test51: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test51: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test51 failed');
end

function n = localBytes(f)
    d = dir(f); n = d.bytes;
end
