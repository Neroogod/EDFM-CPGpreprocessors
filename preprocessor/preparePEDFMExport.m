function Gout = preparePEDFMExport(G, varargin)
% Convert standard pEDFM projection records to the face-indexed contract.
%
% OPTIONS
%   ambiguousFacePolicy - 'error' (default) or 'areaWeighted'. A pEDFM
%       cell pair may correspond to several fault-split physical faces.
%       The default refuses to guess. 'areaWeighted' distributes the
%       projected area over the candidate faces in proportion to their
%       physical areas.

opt = struct('ambiguousFacePolicy', 'error');
opt = merge_options(opt, varargin{:});
policy = validatestring(opt.ambiguousFacePolicy, {'error', 'areaWeighted'});

assert(isfield(G, 'Matrix') && isfield(G, 'FracGrid'), ...
    'G must be an assembled pEDFM grid.');
assert(isfield(G, 'nnc') && isfield(G.nnc, 'cells') && isfield(G.nnc, 'T'), ...
    'G must contain G.nnc.cells and G.nnc.T.');
assert(isfield(G.nnc, 'pMMneighs'), ...
    'Standard pEDFM data G.nnc.pMMneighs is missing.');

Gout = G;
pairs = G.nnc.pMMneighs;
nProj = size(pairs, 1);

if nProj == 0
    Gout.nnc.pMMfaces = zeros(0, 1);
    Gout.nnc.pMMarea = zeros(0, 1);
    Gout.nnc.exportModelType = 'pEDFM';
    return
end

assert(size(pairs, 2) == 2 && all(pairs(:) >= 1) && ...
    all(pairs(:) <= G.Matrix.cells.num), ...
    'pMMneighs must contain valid matrix-cell pairs.');

projArea = localProjectedAreas(G, nProj);
faces = zeros(0, 1);
areas = zeros(0, 1);

for q = 1:nProj
    pair = pairs(q, :);
    nb = G.Matrix.faces.neighbors;
    candidate = find((nb(:,1) == pair(1) & nb(:,2) == pair(2)) | ...
                     (nb(:,1) == pair(2) & nb(:,2) == pair(1)));
    candidate = candidate(all(nb(candidate, :) > 0, 2));

    if isempty(candidate)
        error('preparePEDFMExport:noFace', ...
            'No physical matrix face was found for projected pair [%d %d].', pair);
    elseif numel(candidate) == 1
        faces(end+1, 1) = candidate; %#ok<AGROW>
        areas(end+1, 1) = projArea(q); %#ok<AGROW>
    elseif strcmp(policy, 'error')
        error('preparePEDFMExport:ambiguousFaultFace', ...
            ['Projected pair [%d %d] maps to %d fault-split faces. ', ...
             'Record face IDs during pEDFM projection or use ', ...
             '''ambiguousFacePolicy'',''areaWeighted''.'], ...
            pair(1), pair(2), numel(candidate));
    else
        faceArea = G.Matrix.faces.areas(candidate);
        weights = faceArea ./ sum(faceArea);
        faces = [faces; candidate]; %#ok<AGROW>
        areas = [areas; projArea(q).*weights]; %#ok<AGROW>
    end
end

Gout.nnc.pMMfaces = faces;
Gout.nnc.pMMarea = areas;
Gout.nnc.exportModelType = 'pEDFM';

% Legacy pMatFracNNCs3D appends transmissibilities before collapsing
% duplicate projected pairs, but collapses cells/areas without applying
% the same unique indices to G.nnc.T. Rebuild the final projected NNC
% transmissibilities from the retained cells and retained union areas.
nCells = size(Gout.nnc.cells, 1);
nBase = nCells - nProj;
assert(nBase >= 0 && numel(Gout.nnc.T) >= nBase, ...
    'The pEDFM NNC arrays cannot be reconciled.');
projectionCells = Gout.nnc.cells(nBase+1:end, :);
assert(size(projectionCells, 1) == nProj, ...
    'The final pEDFM NNC rows do not match pMMneighs.');

k1 = Gout.rock.perm(projectionCells(:,1), 1);
k2 = Gout.rock.perm(projectionCells(:,2), 1);
kHarm = zeros(nProj, 1);
positive = (k1 + k2) > 0;
kHarm(positive) = 2*k1(positive).*k2(positive)./(k1(positive) + k2(positive));
distance = sqrt(sum((Gout.cells.centroids(projectionCells(:,1), :) - ...
                     Gout.cells.centroids(projectionCells(:,2), :)).^2, 2));
distance = max(distance, eps);
projectionT = kHarm.*projArea./distance;
Gout.nnc.T = [Gout.nnc.T(1:nBase); projectionT];

assert(numel(Gout.nnc.pMMfaces) == numel(Gout.nnc.pMMarea));
assert(all(isfinite(Gout.nnc.pMMarea)) && all(Gout.nnc.pMMarea >= 0), ...
    'pEDFM projected areas must be finite and nonnegative.');
assert(size(Gout.nnc.cells, 1) == numel(Gout.nnc.T), ...
    'pEDFM NNC cells and transmissibilities have inconsistent lengths.');
end

function projectedArea = localProjectedAreas(G, nProj)
if isfield(G.nnc, 'pMMarea') && numel(G.nnc.pMMarea) == nProj
    projectedArea = G.nnc.pMMarea(:);
    return
end

assert(isfield(G.nnc, 'area') && numel(G.nnc.area) >= nProj, ...
    'Cannot recover pEDFM projected areas from G.nnc.area.');

% pMatFracNNCs3D appends one projected matrix-fracture NNC, and hence one
% projected area, for every pMMneighs record. These are the final nProj
% entries even when the original NNC normal array has a different length.
projectedArea = G.nnc.area(end-nProj+1:end);
projectedArea = projectedArea(:);
end
