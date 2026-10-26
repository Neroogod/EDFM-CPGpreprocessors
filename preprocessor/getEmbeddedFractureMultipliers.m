function multiplier = getEmbeddedFractureMultipliers(G, tol)
% Assemble matrix-face multipliers for normalized EDFM-family exports.

if nargin < 2, tol = 1e-5; end
multiplier = ones(G.faces.num, 1);

if ~isfield(G.nnc, 'pMMfaces')
    if isfield(G.nnc, 'pMMneighs') && ~isempty(G.nnc.pMMneighs)
        error('getEmbeddedFractureMultipliers:notNormalized', ...
            ['Pair-indexed pEDFM projections must be normalized with ', ...
             'preparePEDFMExport before writing.']);
    end
    return
end

if isempty(G.nnc.pMMfaces)
    return
end

assert(isfield(G.nnc, 'pMMarea') && ...
    numel(G.nnc.pMMfaces) == numel(G.nnc.pMMarea), ...
    'Face-indexed projection arrays pMMfaces and pMMarea are inconsistent.');

face = G.nnc.pMMfaces(:);
area = G.nnc.pMMarea(:);
assert(all(face >= 1 & face <= G.faces.num & face == floor(face)), ...
    'Projection face IDs are outside the global face range.');
assert(all(isfinite(area) & area >= 0), ...
    'Projected face areas must be finite and nonnegative.');

[uniqueFace, ~, group] = unique(face);
projectedArea = accumarray(group, area);
faceArea = G.faces.areas(uniqueFace);
openArea = max(faceArea - projectedArea, 0);
openArea(openArea < tol) = 0;
multiplier(uniqueFace) = openArea ./ faceArea;
end
