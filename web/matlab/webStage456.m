function out = webStage456(stage, inFile, outFile, opts)
% Stages 4, 5 and 6 of the pipeline. They share this file because each is a
% single upstream call plus bookkeeping, and keeping them together makes the
% sequence readable:
%
%   4  fracturematrixShaleNNC3D      matrix-fracture connections
%   5  fracturefractureShaleNNCs3D   fracture-fracture connections
%   6  projection                    pEDFM or CPEDFM, or nothing for EDFM
%
% Stage 5 is a no-op with a single plane, and stage 6 is a no-op for plain
% EDFM. Both say so rather than failing, because "nothing to do" is a normal
% outcome the user should see reported.

S = load(inFile);
tol = localField(opts, 'tol', 1e-5);
warnings = {};

switch stage
    case 4
        before = size(localNNC(S.G), 1);
        t0 = tic;
        S.G = fracturematrixShaleNNC3D(S.G, tol);
        elapsed = toc(t0);
        added = size(localNNC(S.G), 1) - before;
        stats = struct('step', 'fracture-matrix NNCs', ...
            'added', added, 'total', size(localNNC(S.G), 1), 'seconds', elapsed);
        if added == 0
            warnings{end+1} = ...
                'No matrix-fracture connections were created.';
        end

    case 5
        if numel(S.fracplanes) < 2
            stats = struct('step', 'fracture-fracture NNCs', ...
                'added', 0, 'total', size(localNNC(S.G), 1), 'seconds', 0);
            warnings{end+1} = ...
                'Only one fracture plane: there are no plane pairs to intersect.';
        else
            before = size(localNNC(S.G), 1);
            t0 = tic;
            [S.G, S.fracplanes] = fracturefractureShaleNNCs3D(S.G, S.fracplanes, tol);
            elapsed = toc(t0);
            stats = struct('step', 'fracture-fracture NNCs', ...
                'added', size(localNNC(S.G), 1) - before, ...
                'total', size(localNNC(S.G), 1), 'seconds', elapsed);
        end

    case 6
        method = upper(localField(opts, 'method', 'EDFM'));
        before = size(localNNC(S.G), 1);
        t0 = tic;
        switch method
            case 'EDFM'
                S.G = prepareEDFMExport(S.G);
            case 'PEDFM'
                S.G = preparePEDFMExport(pMatFracNNCs3D(S.G, tol));
            case 'CPEDFM'
                gi = cmgDetectGrid(S.G);
                if strcmp(gi.family, 'unstructured')
                    assert(exist('pMatFracNNCsCPEDFM_U', 'file') == 2, ...
                        'webStage456:noCPEDFMU', ...
                        'pMatFracNNCsCPEDFM_U is unavailable (CPEDFM-U not on the path).');
                    S.G = pMatFracNNCsCPEDFM_U(S.G, tol);
                else
                    S.G = pMatFracNNCsCPEDFM_CPG(S.G, tol);
                end
            otherwise
                error('webStage456:method', ...
                    'Unknown method "%s"; expected EDFM, pEDFM or CPEDFM.', method);
        end
        elapsed = toc(t0);
        fi = cmgDetectFractureModel(S.G);
        stats = struct('step', 'projection', 'method', method, ...
            'detected', fi.model, ...
            'added', size(localNNC(S.G), 1) - before, ...
            'total', size(localNNC(S.G), 1), ...
            'projectionRecords', max(fi.nProjFace, fi.nProjPair), ...
            'seconds', elapsed);
        if ~strcmpi(fi.model, method)
            warnings{end+1} = sprintf( ...
                'Requested %s but the grid now reports %s.', method, fi.model);
        end
        if ~strcmp(method, 'EDFM') && stats.projectionRecords == 0
            % The projection routines only act on "fracmat interior" pairs.
            % A fracture lying exactly on cell faces produces none, so the
            % export is EDFM in all but name - and silently, unless we say so.
            warnings{end+1} = [ ...
                'No projection records were produced. Every fracture-matrix ' ...
                'connection is a boundary one - the planes may lie on cell ' ...
                'faces - so this export is equivalent to EDFM.'];
        end

    otherwise
        error('webStage456:stage', 'Stage %d is not handled here.', stage);
end

save(outFile, '-struct', 'S', '-v7.3');

out.warnings = warnings;
out.stats = stats;
out.files = struct('checkpoint', outFile);
end

%--------------------------------------------------------------------------
function c = localNNC(G)
c = zeros(0, 2);
if isfield(G, 'nnc') && isfield(G.nnc, 'cells') && ~isempty(G.nnc.cells)
    c = G.nnc.cells;
end
end

%--------------------------------------------------------------------------
function v = localField(s, name, dflt)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = dflt;
end
end
