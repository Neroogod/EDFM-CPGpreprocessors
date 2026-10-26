function webStage(stageNum, inFile, outFile, optsFile, resultFile)
% Single entry point the web backend calls for every pipeline stage.
%
%   webStage(1, '',            'stage1.mat', 'opts.json', 'result.json')
%   webStage(7, 'stage1.mat',  '',           'opts.json', 'result.json')
%
% Options arrive as a JSON *file* rather than a string: embedding JSON in a
% MATLAB command line means quoting it through a shell, and one apostrophe
% in a title would break the call.
%
% The result JSON is always written - including on failure, with ok=false
% and the MATLAB message - so the backend never has to guess why a stage
% produced nothing.

result = struct('ok', false, 'stage', stageNum, 'errors', {{}}, ...
    'warnings', {{}}, 'stats', struct(), 'files', struct());

try
    paths = webSetupPaths();
    result.environment = paths;

    if isempty(optsFile) || exist(optsFile, 'file') ~= 2
        opts = struct();
    else
        opts = jsondecode(fileread(optsFile));
    end

    t0 = tic;
    switch stageNum
        case 1
            out = webStage1Grid(outFile, opts);
        case 2
            out = webStage2Fractures(inFile, outFile, opts);
        case 3
            out = webStage3EDFM(inFile, outFile, opts);
        case {4, 5, 6}
            out = webStage456(stageNum, inFile, outFile, opts);
        case 7
            out = webStage7Export(inFile, opts);
        otherwise
            error('webStage:unknownStage', ...
                'Stage %d is not a pipeline stage; expected 1 to 7.', stageNum);
    end
    result.elapsed_s = toc(t0);

    if isfield(out, 'stats'), result.stats = out.stats; end
    if isfield(out, 'files'), result.files = out.files; end
    if isfield(out, 'warnings'), result.warnings = out.warnings; end
    result.ok = true;

catch ME
    result.errors = {ME.message};
    for k = 1:min(numel(ME.stack), 5)
        result.errors{end+1} = sprintf('  in %s (line %d)', ...
            ME.stack(k).name, ME.stack(k).line); %#ok<AGROW>
    end
end

fid = fopen(resultFile, 'w');
if fid > 0
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s\n', jsonencode(result, 'PrettyPrint', true));
end

if ~result.ok
    error('webStage:failed', '%s', strjoin(result.errors, ' | '));
end
end
