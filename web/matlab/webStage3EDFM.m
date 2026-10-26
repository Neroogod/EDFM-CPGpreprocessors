function out = webStage3EDFM(inFile, outFile, opts)
% Stage 3: cut the fracture planes into the grid, creating fracture cells.
%
% The gridder is chosen from the grid family, because they are not
% interchangeable:
%
%   cart / vari   EDFMshalegrid   hexahedral cells
%   corner        EDFMGridCPG     corner-point hexahedra
%   unstructured  EDFMGridU       half-space clipping, any convex cell
%
% EDFMshalegrid and EDFMGridCPG both assume eight corners per cell and fail
% inside checkplaneAABBintersect on a PEBI grid, so picking by family is not
% a convenience - it is the difference between working and a confusing error.

S = load(inFile);
G = S.G;
assert(isfield(S, 'fracplanes') && ~isempty(S.fracplanes), ...
    'webStage3EDFM:noFractures', 'Stage 2 has not defined any fractures.');

gridInfo = cmgDetectGrid(G);
warnings = {};

% The gridders' Tolerance is an absolute length in metres, and their default
% of 1e-5 is tuned for shale-scale models a few hundred metres across. On a
% field-scale grid the clipping runs out of room and aborts with "Tolerance
% too tight" - which says nothing about the actual cause. Scale the starting
% value with the model instead, and if it still fails, loosen and retry.
extent = max(max(G.nodes.coords, [], 1) - min(G.nodes.coords, [], 1));
tol = localField(opts, 'tol', max(1e-5, 1e-7*extent));

attempt = 0;
while true
    attempt = attempt + 1;
    t0 = tic;
    try
        [Gf, fracplanes, gridder] = localGrid(G, S.fracplanes, gridInfo, tol);
        elapsed = toc(t0);
        break
    catch ME
        if ~contains(ME.message, 'Tolerance too tight') || attempt >= 4
            rethrow(ME);
        end
        tol = tol*10;
        warnings{end+1} = sprintf( ...
            'The gridder could not resolve the geometry; retrying with tolerance %g m.', ...
            tol); %#ok<AGROW>
    end
end

nMat = Gf.Matrix.cells.num;
nFrac = Gf.cells.num - nMat;

if nFrac == 0
    warnings{end+1} = ...
        'No fracture cells were created: the planes may not intersect the grid.';
end

S.G = Gf;
S.fracplanes = fracplanes; %#ok<STRNU>
save(outFile, '-struct', 'S', '-v7.3');

out.warnings = warnings;
out.stats = struct( ...
    'gridder', gridder, ...
    'family', gridInfo.family, ...
    'matrixCells', nMat, ...
    'fractureCells', nFrac, ...
    'planes', numel(fracplanes), ...
    'tolerance', tol, ...
    'seconds', elapsed);
out.files = struct('checkpoint', outFile);
end

%--------------------------------------------------------------------------
function [Gf, fracplanes, gridder] = localGrid(G, fracplanes, gridInfo, tol)
list = 1:numel(fracplanes);
switch gridInfo.family
    case {'cart', 'vari'}
        gridder = 'EDFMshalegrid';
        [Gf, fracplanes] = EDFMshalegrid(G, fracplanes, ...
            'Tolerance', tol, 'plotgrid', false, 'fracturelist', list);

    case 'corner'
        gridder = 'EDFMGridCPG';
        [Gf, fracplanes] = EDFMGridCPG(G, fracplanes, ...
            'Tolerance', tol, 'plotgrid', false, 'fracturelist', list);

    case 'unstructured'
        gridder = 'EDFMGridU';
        assert(exist('EDFMGridU', 'file') == 2, 'webStage3EDFM:noGridderU', ...
            ['EDFMGridU is unavailable; the CPEDFM-U module is not on the ' ...
             'path. Unstructured grids cannot use the hexahedral gridders.']);
        [Gf, fracplanes] = EDFMGridU(G, fracplanes, ...
            'Tolerance', tol, 'fracturelist', list, 'Verbose', false);

    otherwise
        error('webStage3EDFM:family', ...
            'No EDFM gridder for grid family "%s".', gridInfo.family);
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
