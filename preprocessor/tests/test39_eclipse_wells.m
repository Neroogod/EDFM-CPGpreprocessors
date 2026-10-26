%% test39_eclipse_wells
% Wells in the Eclipse writer.
%
% They are emitted commented out, because they belong in SCHEDULE and the
% stub deck stops at the end of GRID. The test therefore reads the deck as
% text: a commented block is invisible to a parser, which is exactly the
% property that makes it safe to write, and exactly what makes it easy to get
% silently wrong.
%
% The value that matters is the connection factor. On a padding-cell or
% connection-list grid the block geometry is a placeholder, so Eclipse's
% radial well model would produce a meaningless index; MRST's WI is written
% explicitly instead, converted from md.m to cP.m3/(day.bar).

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
mrstModule add deckformat incomp

outDir = fullfile(thisDir, 'output_eclipse');
if ~isfolder(outDir), mkdir(outDir); end

failures = {};
report = @(ok, name, detail) fprintf('%s %-32s %s\n', ...
    subsref({'FAIL', 'PASS'}, substruct('{}', {ok + 1})), name, detail);

%% ---- a grid with two wells --------------------------------------------
dims = [10 8 4];
G = computeGeometry(cartGrid(dims, [500 400 40]));
rock = makeRock(G, 100*milli*darcy, 0.2);

[ii, jj, ~] = ind2sub(dims, (1:G.cells.num)');
injCol  = find(ii == 2 & jj == 2);
prodCol = find(ii == 9 & jj == 7);

W = addWell([], G, rock, injCol, 'Type', 'rate', 'Val', 300*meter^3/day, ...
    'Radius', 0.1, 'Name', 'INJ', 'Comp_i', [1 0], 'sign', 1);
W = addWell(W, G, rock, prodCol, 'Type', 'bhp', 'Val', 150*barsa, ...
    'Radius', 0.1, 'Name', 'PROD', 'Comp_i', [0 1], 'sign', -1);

info = writeEclipse(G, fullfile(outDir, 'ecl_wells.DATA'), ...
    'rock', rock, 'wells', W, 'datum', 2000, 'title', 'ECL_WELLS');

txt = fileread(info.dataset);
lines = strsplit(txt, newline);

%% ---- A: both wells appear ---------------------------------------------
okA = info.nWells == 2 && contains(txt, 'WELSPECS') && contains(txt, 'COMPDAT') ...
    && contains(txt, '''INJ''') && contains(txt, '''PROD''');
report(okA, 'A. wells written', sprintf('%d wells, WELSPECS + COMPDAT present', info.nWells));
if ~okA, failures{end+1} = 'A'; end

%% ---- B: one COMPDAT row per completion --------------------------------
compRows = sum(cellfun(@(L) ~isempty(regexp(L, '^--\s+''(INJ|PROD)''\s+\d+\s+\d+\s+\d+\s+\d+\s+''OPEN''', 'once')), lines));
nPerf = numel(W(1).cells) + numel(W(2).cells);
okB = compRows == nPerf;
report(okB, 'B. one row per completion', sprintf('%d rows for %d perforations', compRows, nPerf));
if ~okB, failures{end+1} = 'B'; end

%% ---- C: connection factors are MRST's WI, converted -------------------
mdm = milli*darcy*meter;
ecl = centi*poise*meter^3/(day*barsa);
model = buildCMGModel(G, rock, 'datum', 2000, 'wells', W);

got = [];
for k = 1:numel(lines)
    t = regexp(lines{k}, '^--\s+''(?:INJ|PROD)''\s+\d+\s+\d+\s+\d+\s+\d+\s+''OPEN''\s+1\*\s+([\d.eE+-]+)', 'tokens', 'once');
    if ~isempty(t), got(end+1,1) = str2double(t{1}); end %#ok<AGROW>
end
want = [model.wells(1).wi; model.wells(2).wi]*(mdm/ecl);
okC = numel(got) == numel(want) && max(abs(got - want)./max(abs(want), eps)) < 1e-7;
report(okC, 'C. connection factors', sprintf( ...
    '%d values, max rel diff %.2e, range %.4g..%.4g', numel(got), ...
    max(abs(got - want)./max(abs(want), eps)), min(got), max(got)));
if ~okC, failures{end+1} = 'C'; end

%% ---- D: the well block stays invisible to the parser ------------------
% Everything well-related is commented, so the deck must still parse and must
% not gain a SCHEDULE section.
% readEclipseDeck always creates a SCHEDULE struct; what shows whether any
% well data was read is SCHEDULE.control, which stays empty while the block
% is commented.
deck = readEclipseDeck(info.dataset);
okD = isfield(deck, 'SCHEDULE') && isempty(deck.SCHEDULE.control) && ...
    isequal(reshape(deck.RUNSPEC.DIMENS, 1, []), dims);
report(okD, 'D. deck still parses', sprintf( ...
    'DIMENS %s recovered, %d well controls parsed (block is inert)', ...
    mat2str(reshape(deck.RUNSPEC.DIMENS, 1, [])), numel(deck.SCHEDULE.control)));
if ~okD, failures{end+1} = 'D'; end

%% ---- E: controls follow the MRST well type ----------------------------
okE = contains(txt, 'WCONINJE') && contains(txt, 'RATE') && ...
      contains(txt, 'WCONPROD') && contains(txt, 'BHP');
report(okE, 'E. controls match', 'rate injector, BHP producer');
if ~okE, failures{end+1} = 'E'; end

fprintf('\nwrote %s\n', info.dataset);

%% ---- summary -----------------------------------------------------------
fprintf('\n');
if isempty(failures)
    fprintf('=== test39: ALL CHECKS PASSED ===\n');
else
    fprintf('=== test39: %d FAILURES: %s ===\n', numel(failures), ...
        strjoin(failures, ', '));
end
