function Gout = prepareEDFMExport(G)
% Normalize an EDFM grid for the common embedded-fracture deck writers.
%
% Projection NNCs are removed when a projected grid is supplied, making
% the requested model type explicit instead of accidentally exporting a
% pEDFM/CPEDFM model under an EDFM filename.

assert(isfield(G, 'Matrix') && isfield(G, 'FracGrid'), ...
    'G must be an assembled EDFM grid with Matrix and FracGrid fields.');
assert(isfield(G, 'nnc') && isfield(G.nnc, 'cells') && isfield(G.nnc, 'T'), ...
    'G must contain G.nnc.cells and G.nnc.T.');

Gout = G;
n = size(Gout.nnc.cells, 1);
keep = true(n, 1);
originalT = Gout.nnc.T(:);

if isfield(Gout.nnc, 'type') && numel(Gout.nnc.type) == n
    nncType = string(Gout.nnc.type(:));
    keep = ~contains(lower(nncType), 'fracpmat') & ...
           ~contains(lower(nncType), 'projection');
end

fields = fieldnames(Gout.nnc);
for i = 1:numel(fields)
    name = fields{i};
    value = Gout.nnc.(name);
    if ~isempty(value) && size(value, 1) == n
        Gout.nnc.(name) = value(keep, :);
    end
end

% Legacy pEDFM grids can have more T entries than NNC cell rows because
% duplicate projections were collapsed incompletely. Original EDFM NNCs
% precede all projection entries, so retain that leading block explicitly.
nKeep = nnz(keep);
assert(numel(originalT) >= nKeep, 'Not enough transmissibilities for EDFM NNCs.');
Gout.nnc.T = originalT(1:nKeep);

Gout.nnc.pMMfaces = zeros(0, 1);
Gout.nnc.pMMarea = zeros(0, 1);
Gout.nnc.pMMneighs = zeros(0, 2);
Gout.nnc.exportModelType = 'EDFM';

assert(size(Gout.nnc.cells, 1) == numel(Gout.nnc.T), ...
    'EDFM NNC cells and transmissibilities have inconsistent lengths.');
assert(all(isfinite(Gout.nnc.T)), ...
    'EDFM NNC transmissibilities must be finite before export.');
end
