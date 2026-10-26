function out = webStage7Export(inFile, opts)
% Stage 7: write the deck. Runs after any of stages 1, 3, 4, 5 or 6, so a
% matrix-only export is not a special case - it is simply this stage reached
% early.
%
% opts.format    'cmg' (default) or 'eclipse'
% opts.simulator 'IMEX' | 'GEM' | 'STARS'   (CMG only)
% opts.hostDeck  'stub' | 'none' | 'runnable'
% opts.outDir, opts.caseName, opts.title

S = load(inFile);
G = S.G;
rock = G.rock;                       % after EDFM this covers fracture cells too
if ~isfield(G, 'rock') && isfield(S, 'rock'), rock = S.rock; end

format   = lower(localField(opts, 'format', 'cmg'));
outDir   = localField(opts, 'outDir', pwd);
caseName = localField(opts, 'caseName', 'case');
titleStr = localField(opts, 'title', upper(caseName));
hostDeck = localField(opts, 'hostDeck', 'stub');
tol      = localField(opts, 'tol', 1e-5);
datum    = [];
if isfield(S, 'datum'), datum = S.datum; end
datum = localField(opts, 'datum', datum);

grdecl = [];
if isfield(S, 'grdecl'), grdecl = S.grdecl; end

if ~isfolder(outDir), mkdir(outDir); end

args = {'rock', rock, 'grdecl', grdecl, 'datum', datum, ...
        'tol', tol, 'title', titleStr, 'hostDeck', hostDeck};

switch format
    case 'cmg'
        info = writeCMG(G, fullfile(outDir, [caseName '.dat']), ...
            'simulator', localField(opts, 'simulator', 'IMEX'), ...
            'cmgVersion', localField(opts, 'cmgVersion', '2024.30'), args{:});
        files = struct('dataset', info.dataset, 'include', info.include, ...
            'manifest', info.manifest);

    case 'eclipse'
        assert(~strcmp(hostDeck, 'runnable'), 'webStage7Export:hostDeck', ...
            ['The Eclipse writer has no runnable host deck: it stops at the ' ...
             'end of GRID. Use ''stub'' or ''none''.']);
        info = writeEclipse(G, fullfile(outDir, [caseName '.DATA']), args{:});
        files = struct('dataset', info.dataset, 'include', info.include);

    otherwise
        error('webStage7Export:format', ...
            'Unknown output format "%s"; expected cmg or eclipse.', format);
end

fi = cmgDetectFractureModel(G);
out.warnings = {};
out.stats = struct( ...
    'format', format, ...
    'gridFamily', info.gridFamily, ...
    'fractureModel', info.fractureModel, ...
    'detection', fi.reason, ...
    'blocks', info.counts.blocks, ...
    'inactive', info.counts.inactive, ...
    'connections', info.counts.connections, ...
    'poreVolume', sum(poreVolume(G, rock)));
if isfield(info.counts, 'paddingLayers')
    out.stats.paddingLayers = info.counts.paddingLayers;
end
out.files = files;
end

%--------------------------------------------------------------------------
function v = localField(s, name, dflt)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = dflt;
end
end
