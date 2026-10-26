function out = webStage2Fractures(inFile, outFile, opts)
% Stage 2: define the fracture planes. No gridding happens here - this stage
% only builds and validates the fracplanes struct, so a bad plane is caught
% in milliseconds rather than after a five-minute intersection run.
%
% opts.mode is one of
%   'explicit'   planes: array of {points [4x3 or Nx3], aperture, poro, perm}
%   'multistage' numStages, fracsPerStage, spacing, clusterSpacing,
%                halfLength, height, centroid1  (createMultiStageHFs)
%
% Common: aperture [m], poro, perm [mD] applied to any plane that does not
% carry its own.

S = load(inFile);
G = S.G;

mode = lower(localField(opts, 'mode', 'explicit'));
apDef = localField(opts, 'aperture', 1e-3);
poDef = localField(opts, 'poro', 0.5);
pmDef = localField(opts, 'perm', 1e4);        % mD

warnings = {};

switch mode
    case 'explicit'
        planes = localField(opts, 'planes', []);
        assert(~isempty(planes), 'webStage2Fractures:noPlanes', ...
            'No fracture planes were supplied.');
        if ~iscell(planes), planes = num2cell(planes); end
        fracplanes = struct('points', {}, 'aperture', {}, 'poro', {}, 'perm', {});
        for k = 1:numel(planes)
            p = planes{k};
            pts = reshape(localField(p, 'points', []), [], 3);
            assert(size(pts, 1) >= 3, 'webStage2Fractures:points', ...
                'Plane %d needs at least three points.', k);
            fracplanes(k).points   = pts; %#ok<AGROW>
            fracplanes(k).aperture = localField(p, 'aperture', apDef); %#ok<AGROW>
            fracplanes(k).poro     = localField(p, 'poro', poDef); %#ok<AGROW>
            fracplanes(k).perm     = localField(p, 'perm', pmDef)*milli*darcy; %#ok<AGROW>
        end

    case 'multistage'
        shaleUtils = fullfile(mrstPath('shale'), 'utils');
        if isfolder(shaleUtils), addpath(shaleUtils); end
        assert(exist('createMultiStageHFs', 'file') == 2, ...
            'webStage2Fractures:noShale', ...
            'createMultiStageHFs is unavailable; the shale module is missing.');
        [fracplanes, ~] = createMultiStageHFs( ...
            'numStages',        localField(opts, 'numStages', 1), ...
            'fracSpacing',      localField(opts, 'spacing', 100), ...
            'numFracsPerStage', localField(opts, 'fracsPerStage', 3), ...
            'fracHalfLength',   localField(opts, 'halfLength', 100), ...
            'fracHeight',       localField(opts, 'height', 40), ...
            'clusterSpacing',   localField(opts, 'clusterSpacing', 10), ...
            ... jsondecode turns a JSON array into a column; the generator
            ... wants a row triple, and silently misreads a column.
            'fracCentroid1',    reshape(localField(opts, 'centroid1', ...
                                    [50 50 20]), 1, []), ...
            'isStencil', 1);
        for k = 1:numel(fracplanes)
            fracplanes(k).aperture = apDef;
            fracplanes(k).poro     = poDef;
            fracplanes(k).perm     = pmDef*milli*darcy;
        end

    otherwise
        error('webStage2Fractures:mode', 'Unknown fracture mode "%s".', mode);
end

% ---- validate against the grid ----------------------------------------
lo = min(G.nodes.coords, [], 1);
hi = max(G.nodes.coords, [], 1);
areas = zeros(numel(fracplanes), 1);
outside = 0;
for k = 1:numel(fracplanes)
    pts = fracplanes(k).points;
    areas(k) = localPolyArea(pts);
    if any(max(pts, [], 1) < lo - 1e-9) || any(min(pts, [], 1) > hi + 1e-9)
        outside = outside + 1;
    end
    assert(fracplanes(k).aperture > 0, 'webStage2Fractures:aperture', ...
        'Plane %d has a non-positive aperture.', k);
    assert(all(fracplanes(k).perm > 0), 'webStage2Fractures:perm', ...
        'Plane %d has a non-positive permeability.', k);
end
if outside > 0
    warnings{end+1} = sprintf( ...
        ['%d plane(s) lie entirely outside the grid bounding box and will ' ...
         'produce no fracture cells.'], outside);
end

S.fracplanes = fracplanes; %#ok<STRNU>
save(outFile, '-struct', 'S', '-v7.3');

out.warnings = warnings;
out.stats = struct( ...
    'mode', mode, ...
    'planes', numel(fracplanes), ...
    'totalArea', sum(areas), ...
    'minArea', min(areas), ...
    'maxArea', max(areas), ...
    'apertureMM', apDef*1000, ...
    'outsideGrid', outside);
out.files = struct('checkpoint', outFile);
end

%--------------------------------------------------------------------------
function A = localPolyArea(pts)
% Fan triangulation from the first vertex; exact for a planar polygon.
A = 0;
for t = 2:size(pts, 1)-1
    A = A + 0.5*norm(cross(pts(t,:) - pts(1,:), pts(t+1,:) - pts(1,:)));
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
