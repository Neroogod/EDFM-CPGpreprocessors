%% demoWaterfloodExport
% End-to-end demonstration of modules/preprocessorAll:
%
%   1. build an ordinary MRST case  (heterogeneous quarter-five-spot waterflood)
%   2. run it in MRST
%   3. export it as a COMPLETE CMG dataset
%   4. run that dataset in IMEX
%   5. read IMEX's results back out of the SR3
%   6. plot MRST and CMG side by side
%
% Nothing about the case is special to the exporter: it is written the way
% any MRST study is written, and the export happens in two lines at the end.
%
% Figures are written to examples/output/ as PNGs and the paths are printed.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
mrstModule add preprocessor ad-core ad-props ad-blackoil

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

%% ======================================================================
%  1. An ordinary MRST case
%  ======================================================================
gravity reset on
dims    = [20, 20, 3];
physdim = [400, 400, 30];

G = cartGrid(dims, physdim);
G.nodes.coords(:,3) = G.nodes.coords(:,3) + 2000;    % reservoir at 2000 m
G = computeGeometry(G);

% Heterogeneous permeability, so the flood front has something to react to.
rng(7);
% Kernel dimensions must be odd: gaussianField uses padSz = (sz-1)/2.
poro = gaussianField(dims, [0.12, 0.32], [5 5 3]);
poro = poro(:);
% Kozeny-Carman form, scaled to a realistic reservoir range (tens to a few
% hundred mD). Without the 1e4 factor this returns ~0.01 mD and the injector
% cannot deliver its rate, which shows up as thousands of bar.
perm = 1e4 * (poro.^3) ./ ((1 - poro).^2) * milli*darcy;
rock = makeRock(G, perm, poro);

fluid = initSimpleADIFluid('phases', 'WO', ...
    'mu',   [0.5, 4]*centi*poise, ...
    'rho',  [1010, 830], ...
    'c',    [4.4e-5, 1.0e-4]/barsa, ...
    'n',    [2, 2], ...
    'pRef', 250*barsa);

model  = GenericBlackOilModel(G, rock, fluid, ...
    'water', true, 'oil', true, 'gas', false);
state0 = initResSol(G, 250*barsa, [0.15, 0.85]);

% Quarter five-spot: injector at one corner, producer at the opposite one.
[ii, jj] = deal(G.cells.centroids(:,1), G.cells.centroids(:,2));
cInj  = find(ii < 30  & jj < 30);
cProd = find(ii > 370 & jj > 370);
rateInj = 150*meter^3/day;
bhpProd = 200*barsa;

W = addWell([], G, rock, cInj,  'Type','rate','Val', rateInj, ...
    'Comp_i',[1 0], 'Name','INJ',  'refDepth', 2000, 'Radius', 0.1);
W = addWell(W, G, rock, cProd, 'Type','bhp', 'Val', bhpProd, ...
    'Comp_i',[0 1], 'Name','PROD', 'refDepth', 2000, 'Radius', 0.1);

totTime  = 3*year;
schedule = simpleSchedule(rampupTimesteps(totTime, 90*day, 4), 'W', W);

fprintf('=== MRST case ===\n');
fprintf('  %d x %d x %d = %d cells, %d m x %d m x %d m at 2000 m depth\n', ...
    dims, G.cells.num, physdim);
fprintf('  perm %.1f - %.1f mD, poro %.3f - %.3f\n', ...
    min(perm)/(milli*darcy), max(perm)/(milli*darcy), min(poro), max(poro));
fprintf('  INJ %.0f m3/day, PROD %.0f bar, %d report steps over %.1f years\n\n', ...
    convertTo(rateInj, meter^3/day), convertTo(bhpProd, barsa), ...
    numel(schedule.step.val), convertTo(totTime, year));

%% ======================================================================
%  2. Run it in MRST
%  ======================================================================
fprintf('=== running MRST ===\n');
tM = tic;
[wsMRST, statesMRST] = simulateScheduleAD(state0, model, schedule);
fprintf('  done in %.1f s\n\n', toc(tM));

tDaysM  = cumsum(convertTo(schedule.step.val, day));
qoMRST  = cellfun(@(x) abs(x(2).qOs), wsMRST) * day;    % m3/day
qwMRST  = cellfun(@(x) abs(x(2).qWs), wsMRST) * day;
bhpMRST = cellfun(@(x) x(1).bhp, wsMRST) / barsa;
swMRST  = statesMRST{end}.s(:,1);
pMRST   = statesMRST{end}.pressure / barsa;

%% ======================================================================
%  3. Export as a complete CMG dataset  -- the two lines that matter
%  ======================================================================
fprintf('=== exporting to CMG ===\n');
deck = buildFullDeck(model, state0, schedule);
info = writeCMGDeck(deck, fullfile(outDir, 'waterflood.dat'), ...
    'title', 'MRST_WATERFLOOD_DEMO');

fprintf('  dataset : %s\n', info.dataset);
fprintf('  include : %s\n', info.include);
fprintf('  fidelity warnings: %d\n', numel(deck.warnings));
for q = 1:numel(deck.warnings), fprintf('    %s\n', deck.warnings{q}); end
fprintf('\n');

%% ======================================================================
%  4. Run the exported dataset in IMEX
%  ======================================================================
imex = localFindIMEX();
haveCMG = ~isempty(imex);
sr3 = strrep(info.dataset, '.dat', '.sr3');
if haveCMG
    fprintf('=== running IMEX ===\n  %s\n', imex);
    outLog = strrep(info.dataset, '.dat', '.out');
    for f = {outLog, sr3}
        if isfile(f{1}), delete(f{1}); end
    end
    tC = tic;
    system(sprintf('cd /d "%s" && "%s" -f "%s" -log', outDir, imex, info.dataset));
    logTxt = ''; if isfile(outLog), logTxt = fileread(outLog); end
    nErr = numel(regexp(logTxt, 'ERROR \(from subroutine', 'start'));
    ok = contains(logTxt, 'Normal Termination') && nErr == 0;
    fprintf('  %s in %.1f s (%d errors)\n\n', ...
        localTern(ok, 'Normal Termination', 'ABNORMAL termination'), toc(tC), nErr);
else
    fprintf('=== IMEX not found: showing the MRST side only ===\n\n');
end

%% ======================================================================
%  5. Read IMEX results back
%  ======================================================================
haveSr3 = haveCMG && isfile(sr3);
if haveSr3
    cmg = readCMGSr3(sr3);
    swCMG = cmg.props.SW(:, end);
    pCMG  = cmg.props.PRES(:, end) / 100;               % kPa -> bar

    pick = @(list, name) find(strcmpi(strtrim(list), name), 1);
    iProd = pick(cmg.wells.names, 'PROD');
    iInj  = pick(cmg.wells.names, 'INJ');
    iOil  = pick(cmg.wells.vars,  'OILRATSC');
    iWat  = pick(cmg.wells.vars,  'WATRATSC');
    iBhp  = pick(cmg.wells.vars,  'BHP');

    haveWells = ~isempty(iProd) && ~isempty(iInj) && ...
                ~isempty(iOil) && ~isempty(iWat) && ~isempty(iBhp);
    if ~haveWells
        fprintf(['  NOTE: a well curve could not be located in the SR3.\n' ...
                 '        wells: %s\n        vars : %s\n'], ...
            strjoin(strtrim(cmg.wells.names), ', '), ...
            strjoin(strtrim(cmg.wells.vars),  ', '));
        tDaysC = []; qoCMG = []; qwCMG = []; bhpCMG = [];
    else
        tDaysC = cmg.wells.timeDays;
        qoCMG  = squeeze(cmg.wells.data(iProd, iOil, :));
        qwCMG  = squeeze(cmg.wells.data(iProd, iWat, :));
        bhpCMG = squeeze(cmg.wells.data(iInj,  iBhp, :)) / 100;   % kPa -> bar
    end

    fprintf('=== MRST vs CMG at %.0f days ===\n', tDaysM(end));
    fprintf('  mean Sw        MRST %.4f    CMG %.4f    (diff %.4f)\n', ...
        mean(swMRST), mean(swCMG), abs(mean(swMRST)-mean(swCMG)));
    fprintf('  mean pressure  MRST %.2f bar CMG %.2f bar (diff %.2f)\n', ...
        mean(pMRST), mean(pCMG), abs(mean(pMRST)-mean(pCMG)));
    if haveWells
        fprintf('  final oil rate MRST %.2f m3/d CMG %.2f m3/d\n', ...
            qoMRST(end), qoCMG(end));
    end
    fprintf('  cell-by-cell Sw: max |diff| %.4f, RMS %.4f\n\n', ...
        max(abs(swMRST - swCMG)), sqrt(mean((swMRST - swCMG).^2)));
end

%% ======================================================================
%  6. Figures
%  ======================================================================
figFiles = {};

% ---- Figure 1: the MRST case ------------------------------------------
f1 = figure('Name','1. MRST case','Position',[60 60 1250 420],'Color','w');
subplot(1,3,1);
plotCellData(G, poro, 'EdgeAlpha', 0.05); view(-35, 55); axis tight equal off
colorbar('southoutside'); title('porosity (MRST input)');
subplot(1,3,2);
plotCellData(G, swMRST, 'EdgeAlpha', 0.05); view(-35, 55); axis tight equal off
caxis([0.15 0.85]); colorbar('southoutside'); title('MRST water saturation, 3 years');
subplot(1,3,3);
plotCellData(G, pMRST, 'EdgeAlpha', 0.05); view(-35, 55); axis tight equal off
colorbar('southoutside'); title('MRST pressure [bar], 3 years');
figFiles{end+1} = fullfile(outDir, 'demo_1_mrst_case.png');
exportgraphics(f1, figFiles{end}, 'Resolution', 140);

% ---- Figure 2: MRST vs CMG saturation ---------------------------------
if haveSr3
    f2 = figure('Name','2. MRST vs CMG','Position',[60 60 1250 430],'Color','w');
    cl = [0.15 0.85];
    subplot(1,3,1);
    plotCellData(G, swMRST, 'EdgeAlpha', 0.05); view(-35,55); axis tight equal off
    caxis(cl); colorbar('southoutside'); title('MRST  S_w');
    subplot(1,3,2);
    plotCellData(G, swCMG, 'EdgeAlpha', 0.05); view(-35,55); axis tight equal off
    caxis(cl); colorbar('southoutside'); title('CMG IMEX  S_w  (from the exported deck)');
    subplot(1,3,3);
    plotCellData(G, swCMG - swMRST, 'EdgeAlpha', 0.05); view(-35,55); axis tight equal off
    colorbar('southoutside'); title('difference  (CMG - MRST)');
    figFiles{end+1} = fullfile(outDir, 'demo_2_mrst_vs_cmg_sw.png');
    exportgraphics(f2, figFiles{end}, 'Resolution', 140);
end

% ---- Figure 3: well curves --------------------------------------------
f3 = figure('Name','3. Well response','Position',[60 60 1150 400],'Color','w');
subplot(1,3,1); hold on; grid on
plot(tDaysM, qoMRST, '-o', 'LineWidth', 1.6, 'DisplayName', 'MRST');
if haveSr3 && haveWells, plot(tDaysC, qoCMG, '--s', 'LineWidth', 1.4, 'DisplayName', 'CMG IMEX'); end
xlabel('time [days]'); ylabel('oil rate [m^3/day]'); title('PROD oil rate'); legend('Location','best');

subplot(1,3,2); hold on; grid on
plot(tDaysM, qwMRST, '-o', 'LineWidth', 1.6, 'DisplayName', 'MRST');
if haveSr3 && haveWells, plot(tDaysC, qwCMG, '--s', 'LineWidth', 1.4, 'DisplayName', 'CMG IMEX'); end
xlabel('time [days]'); ylabel('water rate [m^3/day]'); title('PROD water rate'); legend('Location','best');

subplot(1,3,3); hold on; grid on
plot(tDaysM, bhpMRST, '-o', 'LineWidth', 1.6, 'DisplayName', 'MRST');
if haveSr3 && haveWells, plot(tDaysC, bhpCMG, '--s', 'LineWidth', 1.4, 'DisplayName', 'CMG IMEX'); end
xlabel('time [days]'); ylabel('BHP [bar]'); title('INJ bottom-hole pressure'); legend('Location','best');
figFiles{end+1} = fullfile(outDir, 'demo_3_well_curves.png');
exportgraphics(f3, figFiles{end}, 'Resolution', 140);

%% ======================================================================
fprintf('=== figures ===\n');
for k = 1:numel(figFiles), fprintf('  %s\n', figFiles{k}); end
fprintf('\n=== CMG files (open these in Results / Builder) ===\n');
fprintf('  dataset : %s\n', info.dataset);
if haveSr3, fprintf('  results : %s\n', sr3); end
fprintf('\nDONE\n');

%% ---- local helpers ------------------------------------------------------
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

function s = localTern(cond, a, b)
    if cond, s = a; else, s = b; end
end
