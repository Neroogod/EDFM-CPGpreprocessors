function info = cmgDetectGrid(G)
% Classify an MRST grid into the CMG grid family that can represent it.
%
% SYNOPSIS:
%   info = cmgDetectGrid(G)
%
% RETURNS:
%   info - struct with
%     .family    'cart' | 'vari' | 'corner' | 'unstructured'
%     .gridType  the MRST G.type stamp, joined
%     .dims      [ni nj nk] for structured families; [nCells 1 1] otherwise
%     .di,.dj,.dk  block sizes [m] per index direction (cart/vari only)
%     .uniform   true when every direction has constant spacing
%     .ztop      depth to the top of the grid [m] (cart/vari only)
%     .nCells    number of active cells
%
% METHOD:
%   Structured grids are the ones that carry cartDims. MRST has many
%   unstructured constructors (pebi, clippedPebi2D, triangleGrid,
%   tetrahedralGrid, ...) whose type strings vary, but none of them set
%   cartDims, so that field - not the type string - is the discriminator.
%   Within the structured families G.type separates corner-point from
%   tensor-product grids, and uniform spacing selects *CART over *VARI.
%
% SEE ALSO:
%   cmgDetectFractureModel, buildCMGModel, writeCMG

assert(isstruct(G) && isfield(G, 'cells') && isfield(G, 'faces'), ...
    'cmgDetectGrid:notAGrid', 'G does not look like an MRST grid.');

% On an assembled embedded-fracture grid the fracture cells contribute nodes
% of their own, so the global grid is no longer an axis-aligned tensor
% product even when the reservoir grid is. The family, dimensions and block
% spacings all come from the matrix grid; the fracture cells are addressed
% separately as padding layers.
nCellsTotal = G.cells.num;
if isfield(G, 'Matrix') && isfield(G, 'FracGrid') && ...
        isfield(G.Matrix, 'cells') && isfield(G.Matrix, 'faces')
    G = G.Matrix;
end

if isfield(G, 'type')
    gridType = strjoin(cellstr(G.type(:))', ',');
else
    gridType = '';
end

info = struct('family', '', 'gridType', gridType, 'dims', [], ...
    'di', [], 'dj', [], 'dk', [], 'uniform', false, ...
    'ztop', 0, 'zbot', 0, 'nCells', nCellsTotal, ...
    'nMatrixCells', G.cells.num);

% Depth range is needed by every family - it sets the initialization datum
% and the contact depths - so compute it before any family-specific return.
if size(G.nodes.coords, 2) > 2
    info.ztop = min(G.nodes.coords(:, 3));
    info.zbot = max(G.nodes.coords(:, 3));
end

% ---- unstructured: no cartDims, so no I-J-K addressing exists ----------
if ~isfield(G, 'cartDims') || isempty(G.cartDims)
    info.family = 'unstructured';
    info.dims = [G.cells.num, 1, 1];
    return
end

dims = reshape(double(G.cartDims), 1, []);
if numel(dims) == 2, dims = [dims, 1]; end
info.dims = dims;

% ---- corner-point ------------------------------------------------------
if ~isempty(gridType) && contains(gridType, 'processGRDECL')
    info.family = 'corner';
    return
end

% ---- tensor product: recover the per-direction block sizes -------------
% A tensorGrid (and hence cartGrid) is axis aligned, so the distinct node
% coordinates in each direction are exactly the block boundaries.
[info.di, okI] = localSpacing(G.nodes.coords(:,1), dims(1));
[info.dj, okJ] = localSpacing(G.nodes.coords(:,2), dims(2));
if size(G.nodes.coords, 2) > 2
    [info.dk, okK] = localSpacing(G.nodes.coords(:,3), dims(3));
else
    info.dk = 1; okK = true;
end

assert(okI && okJ && okK, 'cmgDetectGrid:notTensor', ...
    ['Grid carries cartDims but its nodes are not an axis-aligned tensor ', ...
     'product. Export it as corner-point or unstructured instead.']);

info.uniform = localIsConstant(info.di) && localIsConstant(info.dj) && ...
    localIsConstant(info.dk);

if info.uniform
    info.family = 'cart';
else
    info.family = 'vari';
end
end

%--------------------------------------------------------------------------
function [d, ok] = localSpacing(coord, n)
% Distinct coordinates along one axis -> n block sizes.
c = uniquetol(coord(:), 1e-9, 'DataScale', max(1, max(abs(coord(:)))));
ok = numel(c) == n + 1;
if ok
    d = diff(sort(c(:)))';
else
    d = [];
end
end

%--------------------------------------------------------------------------
function tf = localIsConstant(d)
tf = isscalar(d) || (max(d) - min(d)) <= 1e-9*max(1, max(abs(d)));
end
