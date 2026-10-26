function ok = webSelfTest(workDir)
% Run the seven pipeline stages exactly the way the backend will: through
% webStage, with JSON option files and JSON results. Nothing here calls a
% stage function directly, so a dispatcher that forgets a case fails here.
%
% Two chains are exercised:
%   A  Cartesian -> fractures -> EDFM gridding -> NNCs -> pEDFM -> CMG
%   B  the same grid exported straight from stage 1 (matrix only, Eclipse)
%
% Chain B is the "spit out the matrix grid first" path, and it is a separate
% chain rather than an assertion because that is how the web app uses it: an
% export can branch off any completed stage.

if nargin < 1 || isempty(workDir)
    workDir = fullfile(tempdir, 'webSelfTest');
end
if isfolder(workDir), rmdir(workDir, 's'); end
mkdir(workDir);

fprintf('web self-test\n  work: %s\n\n', workDir);
ok = true;

%% ---- chain A ---------------------------------------------------------
opts1 = struct('source', 'cartesian', 'dims', [12 10 6], ...
    'physdim', [120 100 30], 'perm', 100, 'poro', 0.2, 'datum', 2000);
r1 = localRun(1, '', fullfile(workDir, 's1.mat'), opts1, workDir);
ok = ok && localCheck('1  grid', r1, @(s) s.cells == 720 && strcmp(s.family, 'cart'));

% Two crossing planes: one gives fracture-matrix NNCs, the pair gives
% fracture-fracture NNCs, so stage 5 has something real to do.
%
% The z values are 5..25, not 2005..2025. cartGrid builds at z=0 and the
% datum is applied by the writer at export time, so a plane placed at
% reservoir depth would miss the grid entirely.
%
% The offsets (35, 55) matter too. Blocks are 10 m, so a plane at x=30 lies
% exactly on a cell face: every connection comes out as "fracmat boundary"
% and pMatFracNNCs3D, which only projects "fracmat interior" pairs, finds
% nothing to do. Cutting through cell interiors is what makes stage 6 real.
p1 = struct('points', [35 10 5; 35 90 5; 35 90 25; 35 10 25], ...
    'aperture', 1e-3, 'poro', 0.5, 'perm', 1e4);
p2 = struct('points', [10 55 5; 110 55 5; 110 55 25; 10 55 25], ...
    'aperture', 1e-3, 'poro', 0.5, 'perm', 1e4);
opts2 = struct('mode', 'explicit', 'planes', {{p1, p2}});
r2 = localRun(2, fullfile(workDir, 's1.mat'), fullfile(workDir, 's2.mat'), ...
    opts2, workDir);
ok = ok && localCheck('2  fractures', r2, @(s) s.planes == 2 && s.outsideGrid == 0);

% No tol: stage 3 derives it from the grid extent, and that path is the
% one the web app uses, so it is the one worth testing.
r3 = localRun(3, fullfile(workDir, 's2.mat'), fullfile(workDir, 's3.mat'), ...
    struct(), workDir);
ok = ok && localCheck('3  EDFM grid', r3, ...
    @(s) s.fractureCells > 0 && strcmp(s.gridder, 'EDFMshalegrid'));

r4 = localRun(4, fullfile(workDir, 's3.mat'), fullfile(workDir, 's4.mat'), ...
    struct('tol', 1e-5), workDir);
ok = ok && localCheck('4  frac-matrix', r4, @(s) s.added > 0);

r5 = localRun(5, fullfile(workDir, 's4.mat'), fullfile(workDir, 's5.mat'), ...
    struct('tol', 1e-5), workDir);
ok = ok && localCheck('5  frac-frac', r5, @(s) s.total >= r4.stats.total);

r6 = localRun(6, fullfile(workDir, 's5.mat'), fullfile(workDir, 's6.mat'), ...
    struct('method', 'pEDFM', 'tol', 1e-5), workDir);
ok = ok && localCheck('6  projection', r6, ...
    @(s) strcmpi(s.detected, 'pEDFM') && s.projectionRecords > 0);

opts7 = struct('format', 'cmg', 'outDir', workDir, 'caseName', 'chainA', ...
    'title', 'WEB SELF-TEST A', 'hostDeck', 'stub');
r7 = localRun(7, fullfile(workDir, 's6.mat'), '', opts7, workDir);
ok = ok && localCheck('7  CMG export', r7, ...
    @(s) strcmpi(s.fractureModel, 'pEDFM') && s.connections > 0);

%% ---- chain C: CPEDFM off the same stage-5 checkpoint ------------------
% Branching from s5 rather than rerunning stages 1-5 keeps the third method
% covered for the cost of two stages.
rC6 = localRun(6, fullfile(workDir, 's5.mat'), fullfile(workDir, 's6c.mat'), ...
    struct('method', 'CPEDFM'), workDir);
ok = ok && localCheck('6C CPEDFM', rC6, ...
    @(s) strcmpi(s.detected, 'CPEDFM') && s.projectionRecords > 0);

optsC = struct('format', 'eclipse', 'outDir', workDir, 'caseName', 'chainC', ...
    'title', 'WEB SELF-TEST C', 'hostDeck', 'stub');
rC7 = localRun(7, fullfile(workDir, 's6c.mat'), '', optsC, workDir);
ok = ok && localCheck('7C CPEDFM Eclipse', rC7, ...
    @(s) strcmpi(s.fractureModel, 'CPEDFM') && s.connections > 0);

%% ---- chain B: matrix only, straight off stage 1 -----------------------
optsB = struct('format', 'eclipse', 'outDir', workDir, 'caseName', 'chainB', ...
    'title', 'WEB SELF-TEST B', 'hostDeck', 'stub');
rB = localRun(7, fullfile(workDir, 's1.mat'), '', optsB, workDir);
ok = ok && localCheck('7B matrix Eclipse', rB, ...
    @(s) strcmpi(s.fractureModel, 'matrix') && s.blocks == 720);

%% ---- the exports must actually exist and be non-trivial ---------------
for f = {r7.files.dataset, r7.files.include, rB.files.dataset, rB.files.include}
    d = dir(f{1});
    good = ~isempty(d) && d.bytes > 100;
    ok = ok && localReport(sprintf('   file %s', localBase(f{1})), good, ...
        localSize(d));
end

fprintf('\n%s\n', localVerdict(ok));
end

%--------------------------------------------------------------------------
function r = localRun(stage, inFile, outFile, opts, workDir)
optsFile   = fullfile(workDir, sprintf('opts%d.json', stage));
resultFile = fullfile(workDir, sprintf('result%d_%s.json', stage, ...
    datestr(now, 'HHMMSSFFF'))); %#ok<TNOW1,DATST>
fid = fopen(optsFile, 'w');
fprintf(fid, '%s', jsonencode(opts));
fclose(fid);
webStage(stage, inFile, outFile, optsFile, resultFile);
r = jsondecode(fileread(resultFile));
end

%--------------------------------------------------------------------------
function ok = localCheck(label, r, pred)
ok = r.ok && pred(r.stats);
localReport(label, ok, localSummary(r));
end

%--------------------------------------------------------------------------
function ok = localReport(label, ok, note)
if ok, mark = 'PASS'; else, mark = 'FAIL'; end
fprintf('  [%s] %-20s %s\n', mark, label, note);
end

%--------------------------------------------------------------------------
function s = localSummary(r)
f = fieldnames(r.stats);
parts = cell(1, 0);
for k = 1:min(numel(f), 5)
    v = r.stats.(f{k});
    if ischar(v)
        parts{end+1} = sprintf('%s=%s', f{k}, v); %#ok<AGROW>
    elseif isscalar(v)
        parts{end+1} = sprintf('%s=%.6g', f{k}, v); %#ok<AGROW>
    end
end
s = strjoin(parts, ' ');
if ~isempty(r.warnings)
    s = sprintf('%s  [%d warning(s)]', s, numel(r.warnings));
end
end

%--------------------------------------------------------------------------
function s = localSize(d)
if isempty(d), s = '(missing)'; else, s = sprintf('%d bytes', d.bytes); end
end

%--------------------------------------------------------------------------
function b = localBase(p)
[~, n, e] = fileparts(p); b = [n e];
end

%--------------------------------------------------------------------------
function v = localVerdict(ok)
if ok, v = 'ALL STAGES PASSED'; else, v = 'SELF-TEST FAILED'; end
end
