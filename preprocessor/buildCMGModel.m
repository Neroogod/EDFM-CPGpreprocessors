function model = buildCMGModel(G, rock, varargin)
% Assemble the simulator-neutral CMG model description from an MRST grid.
%
% SYNOPSIS:
%   model = buildCMGModel(G, rock)
%   model = buildCMGModel(G, rock, 'simulator', 'GEM')
%
% PARAMETERS:
%   G    - any MRST grid: cartGrid, tensorGrid, processGRDECL, or an
%          unstructured grid, with or without an assembled fracture grid.
%   rock - rock of the grid (defaults to G.rock when omitted/empty).
%
% OPTIONAL ('key', value):
%   'simulator' - 'IMEX' (default), 'GEM' or 'STARS'
%   'cmgVersion'- release label recorded in the model metadata
%   'model'     - 'auto' (default) or a forced fracture method
%   'grdecl'    - corner-point description, for the *CORNER fast path
%   'datum'     - depth [m] of the top of the grid; overrides the grid's own
%                 z when the model was built at or near zero depth
%   'tol'       - tolerance handed to the multiplier assembly (default 1e-5)
%
% RETURNS:
%   model - struct with .meta, .grid, .props, .conn, .mult (see below).
%           .conn is nConn x 7 [i1 j1 k1 i2 j2 k2 trans], trans in md.m,
%           and is empty for a matrix-only grid.
%
% This function performs no file I/O. writeCMGModel turns the result into
% CMG cards; keeping the two apart is what lets every grid family and every
% fracture method share one writer.
%
% SEE ALSO:
%   cmgDetectGrid, cmgDetectFractureModel, writeCMGModel, writeCMG

opt = struct('simulator', 'IMEX', 'cmgVersion', 'unspecified', ...
             'model', 'auto', 'grdecl', [], 'datum', [], 'tol', 1e-5, ...
             'wells', []);
opt = merge_options(opt, varargin{:});

if nargin < 2 || isempty(rock)
    assert(isfield(G, 'rock'), 'buildCMGModel:noRock', ...
        'No rock supplied and G.rock is absent.');
    rock = G.rock;
end

gridInfo = cmgDetectGrid(G);
fracInfo = cmgDetectFractureModel(G);

if ~strcmpi(opt.model, 'auto')
    forced = validatestring(opt.model, {'matrix', 'EDFM', 'pEDFM', 'CPEDFM'});
    if ~strcmpi(forced, fracInfo.model)
        fracInfo.reason = sprintf('%s; forced to %s by caller', ...
            fracInfo.reason, forced);
    end
    fracInfo.model = forced;
end

profile = getCMGSimulatorProfile(opt.simulator, opt.cmgVersion);

% ---- units -------------------------------------------------------------
% CMG SI datasets take permeability in mD and lengths in m. Transmissibility
% for *SCONNECT is md.m, which is exactly what the existing CMG writers use.
permUnit = milli*darcy;

model = struct();
model.meta = struct( ...
    'simulator', profile.simulator, ...
    'cmgVersion', profile.version, ...
    'unitSystem', 'SI', ...
    'transUnit', 'md.m', ...
    'gridFamily', gridInfo.family, ...
    'gridType', gridInfo.gridType, ...
    'fractureModel', fracInfo.model, ...
    'detection', fracInfo.reason, ...
    'profile', profile);

model.grid = gridInfo;
model.detect = fracInfo;

% =======================================================================
% Matrix-only grids: the cells of the CMG grid are the cells of G.
% =======================================================================
if ~fracInfo.hasFrac
    model = localBuildMatrixOnly(model, G, rock, gridInfo, opt, permUnit);
    model.wells = localBuildWells(G, opt.wells, model);
    model.meta.counts.wells = numel(model.wells);
    return
end

% =======================================================================
% Fracture grids: matrix cells keep their addresses, fracture cells are
% appended as padding k-layers below the model.
% =======================================================================
if strcmp(gridInfo.family, 'unstructured')
    % No padding layers are needed here: a connection-list model already
    % addresses every cell individually, so a fracture cell is simply
    % another block and its NNCs are simply more connections.
    model = localBuildFracConnList(model, G, rock, fracInfo, opt, permUnit);
else
    assert(any(strcmp(gridInfo.family, {'cart', 'vari', 'corner'})), ...
        'buildCMGModel:fracGridFamily', ...
        'Fracture export has no form for grid family "%s".', gridInfo.family);
    if strcmp(gridInfo.family, 'corner')
        assert(~isempty(opt.grdecl), 'buildCMGModel:needGrdecl', ...
            ['Corner-point fracture export needs the grdecl (COORD/ZCORN) ', ...
             'so the padding layers can be appended. Pass ''grdecl'', grdecl.']);
    end
    model = localBuildWithFractures(model, G, rock, gridInfo, fracInfo, opt, permUnit);
end
model.wells = localBuildWells(G, opt.wells, model);
model.meta.counts.wells = numel(model.wells);
end

%--------------------------------------------------------------------------
function model = localBuildWithFractures(model, G, rock, gridInfo, fracInfo, ...
    opt, permUnit)
% The non-intrusive EDFM representation, in CMG terms.
%
%   * fracture cells occupy padding k-layers below the matrix grid; their
%     geometry is a placeholder and unused slots are *NULL
%   * true fracture pore volume enters through *VOLMOD, true depth through
%     *DEPTH-TOP
%   * every geometric connection of the padding, and the matrix-padding
%     interface, is zeroed, so fracture cells talk only through *SCONNECT
%   * every MRST NNC becomes one *SCONNECT record
%   * pEDFM/CPEDFM seals become *TRANS?*MOD on logical neighbours and
%     negative *SCONNECT (a multiplier) on irregular pairs
%
% The grid is always emitted as *VARI, even when the matrix grid is uniform:
% padding blocks sit at depths unrelated to the layer they occupy, and only
% *VARI lets blocks of one layer lie at different depths.

nMat = fracInfo.nMatrix;
nFrac = fracInfo.nFrac;
NX = gridInfo.dims(1); NY = gridInfo.dims(2); NZ = gridInfo.dims(3);

dzPad = 1.0;                                  % placeholder padding thickness
% Provisional layer count; the real one falls out of the column-based
% placement below, where it is the deepest stack over any single column.
nzExtra = max(ceil(nFrac/(NX*NY)), 1);
dims = [NX, NY, NZ + nzExtra];
nTot = prod(dims);
isCorner = strcmp(gridInfo.family, 'corner');

model.grid.dims = dims;
if ~isempty(opt.datum) && ~isCorner
    shift = opt.datum - model.grid.ztop;
    model.grid.ztop = opt.datum;
    model.grid.zbot = model.grid.zbot + shift;
    zShift = shift;
else
    zShift = 0;
end
model.grid.nzMatrix = NZ;
model.grid.nzPad = nzExtra;

if isCorner
    model.grid.kind = 'corner';
else
    % Padding blocks sit at depths unrelated to the layer they occupy, and
    % only *VARI lets blocks of one layer lie at different depths.
    model.grid.kind = 'vari';
    model.grid.di = gridInfo.di;
    model.grid.dj = gridInfo.dj;
    model.grid.dk = [gridInfo.dk(:)', repmat(dzPad, 1, nzExtra)];
end

% ---- addressing --------------------------------------------------------
if isfield(G.Matrix.cells, 'indexMap')
    matSlot = G.Matrix.cells.indexMap(:);
else
    matSlot = (1:nMat)';
end
[mi, mj, mk] = ind2sub([NX NY NZ], matSlot);
matIJK = [mi(:), mj(:), mk(:)];

% Each fracture cell is placed in the padding stack directly beneath the
% matrix column it intersects, not in an arbitrary slot. With sequential
% slots the padding blocks take their (i,j) from a running counter and their
% depth from the fracture, so they appear scattered across the map and
% float above or below the reservoir - the model is numerically identical
% but impossible to inspect, and geometrically it interleaves padding blocks
% with unrelated columns.
fracIJK = localPlaceFractureCells(G, nMat, nFrac, matIJK, NX, NY, NZ);
nzExtra = max(fracIJK(:,3)) - NZ;
dims(3) = NZ + nzExtra;
nTot = prod(dims);
model.grid.dims = dims;
model.grid.nzPad = nzExtra;
if ~isCorner
    model.grid.dk = [gridInfo.dk(:)', repmat(dzPad, 1, nzExtra)];
end

matLin  = sub2ind(dims, matIJK(:,1),  matIJK(:,2),  matIJK(:,3));
fracLin = sub2ind(dims, fracIJK(:,1), fracIJK(:,2), fracIJK(:,3));

model.grid.matIJK = matIJK;
model.grid.fracIJK = fracIJK;

% ---- properties --------------------------------------------------------
perm = rock.perm;
if size(perm, 2) == 1, perm = repmat(perm, 1, 3); end

permi = zeros(nTot,1); permj = permi; permk = permi; poro = permi;
permi(matLin) = convertTo(perm(1:nMat, 1), permUnit);
permj(matLin) = convertTo(perm(1:nMat, min(2,end)), permUnit);
permk(matLin) = convertTo(perm(1:nMat, min(3,end)), permUnit);
poro(matLin)  = rock.poro(1:nMat);

pf = perm(nMat+1:end, :);
permi(fracLin) = convertTo(pf(:,1), permUnit);
permj(fracLin) = convertTo(pf(:,min(2,end)), permUnit);
permk(fracLin) = convertTo(pf(:,min(3,end)), permUnit);
% Pore volume comes from *PORVOL, but CMG still wants a porosity for every
% active block, and zero would be rejected.
poro(fracLin) = max(rock.poro(nMat+1:end), 1e-6);

active = false(nTot,1);
active(matLin) = true;
active(fracLin) = true;

% ---- pore volume: *PORVOL sets it directly ------------------------------
% *PORVOL overrides the porosity-times-bulk-volume calculation, so the pore
% volume MRST computed transfers exactly for every block - no dependence on
% how CMG evaluates the bulk volume of a placeholder padding cell, and none
% of the accumulation hazard *VOLMOD carries when a block is repeated.
pvAll = poreVolume(G, rock);
porvol = zeros(nTot,1);
porvol(matLin)  = pvAll(1:nMat);
porvol(fracLin) = pvAll(nMat+1:end);

zFrac = G.cells.centroids(nMat+1:end, 3) + zShift;

model.props = struct('permi', permi, 'permj', permj, 'permk', permk, ...
    'poro', poro, 'active', active, 'nTot', nTot, ...
    'volmod', [], 'porvol', porvol, 'depth', []);

if isCorner
    % Corner-point: give each padding block corners that bracket the true
    % fracture depth. Depth then comes from the geometry, so *DEPTH-TOP is
    % not needed and the matrix blocks' own corner-derived depths are left
    % untouched. Padding blocks need not be conforming - every connection
    % among them is removed below.
    [model.grid.coord, model.grid.zcorn] = localExtendCornerPoint( ...
        opt.grdecl, NX, NY, NZ, nzExtra, fracIJK, zFrac, dzPad);
else
    % Structured: per-block depth through *DEPTH-TOP.
    layerTop = model.grid.ztop + [0, cumsum(model.grid.dk(1:end-1))];
    depthTop = zeros(nTot,1);
    depthTop(matLin) = layerTop(matIJK(:,3));
    depthTop(fracLin) = zFrac - 0.5*dzPad;
    [~, ~, kAll] = ind2sub(dims, (1:nTot)');
    fill = ~active;
    depthTop(fill) = layerTop(kAll(fill));
    model.props.depth = depthTop;
end

% ---- connections: every NNC becomes one *SCONNECT ----------------------
transUnit = milli*darcy*meter;
nncCells = G.nnc.cells;
uba = @(c) localUBA(c, nMat, matIJK, fracIJK);
nNNC = size(nncCells, 1);
conn = zeros(nNNC, 7);
for q = 1:nNNC
    conn(q, 1:3) = uba(nncCells(q,1));
    conn(q, 4:6) = uba(nncCells(q,2));
end
conn(:,7) = convertTo(G.nnc.T(:), transUnit);
model.conn = conn;

% ---- seals from the projection stage -----------------------------------
[model.mult, sealConn] = localSeals(G, matIJK, opt.tol);
if ~isempty(sealConn)
    model.conn = [model.conn; sealConn];
end

% ---- isolate the padding region ----------------------------------------
% Fracture cells must communicate only through the explicit connections, so
% every automatic connection inside the padding, and across the interface
% between the last matrix layer and the first padding layer, is zeroed.
model.isolate = struct('keyword', {}, 'range', {}, 'value', {});
model.isolate(end+1) = struct('keyword', '*TRANSK', ...
    'range', [1 NX 1 NY NZ NZ], 'value', 0);
for kw = {'*TRANSI', '*TRANSJ', '*TRANSK'}
    model.isolate(end+1) = struct('keyword', kw{1}, ...
        'range', [1 NX 1 NY NZ+1 dims(3)], 'value', 0);
end

model.meta.counts = struct( ...
    'cells', G.cells.num, ...
    'blocks', nTot, ...
    'inactive', nTot - nnz(active), ...
    'matrixCells', nMat, ...
    'fractureCells', nFrac, ...
    'paddingLayers', nzExtra, ...
    'connections', size(model.conn,1), ...
    'multipliers', size(model.mult,1));
end

%--------------------------------------------------------------------------
function model = localBuildFracConnList(model, G, rock, fracInfo, opt, permUnit)
% Unstructured grid carrying fractures, as a connection-list pseudo-grid.
%
% This is the simplest of the fracture paths. The padding-layer machinery
% exists only because a structured grid has no spare addresses; a
% connection-list model gives every cell its own block already, so matrix
% and fracture cells are numbered straight through and the whole
% connectivity - internal matrix faces plus every NNC - is one *SCONNECT
% list.
%
% Projection seals fold in rather than being emitted separately. In the
% structured path a sealed matrix face becomes a *TRANS?*MOD multiplier
% because the face is implied by the grid; here the face IS an explicit
% connection, so the multiplier is applied to its transmissibility directly.
% Same physics, one mechanism fewer.

nMat = fracInfo.nMatrix;
nTot = G.cells.num;
transUnit = milli*darcy*meter;

model.grid.kind = 'connlist';
model.grid.dims = [nTot, 1, 1];
model.grid.di = 1; model.grid.dj = 1; model.grid.dk = 1;

perm = rock.perm;
if size(perm, 2) == 1, perm = repmat(perm, 1, 3); end

model.props = struct( ...
    'permi', convertTo(perm(:,1), permUnit), ...
    'permj', convertTo(perm(:,min(2,end)), permUnit), ...
    'permk', convertTo(perm(:,min(3,end)), permUnit), ...
    'poro',  max(rock.poro(:), 1e-6), ...
    'active', true(nTot,1), ...
    'nTot', nTot, ...
    'volmod', [], ...
    'porvol', poreVolume(G, rock), ...
    'depth', G.cells.centroids(:,3) - 0.5);

model.grid.ztop = min(model.props.depth);
model.grid.zbot = max(model.props.depth) + 1;

% ---- matrix-matrix connections, with any projection seals applied -------
Nb = G.Matrix.faces.neighbors;
internal = all(Nb > 0, 2);
Tface = localFaceTrans(G.Matrix, struct('perm', rock.perm(1:nMat, :)));

eta = getEmbeddedFractureMultipliers(G, opt.tol);
if numel(eta) >= G.Matrix.faces.num
    Tface = Tface .* eta(1:G.Matrix.faces.num);
end

pairs = Nb(internal, :);
Tm = Tface(internal);
keep = Tm > 0;                      % a fully sealed face carries no flow
mmConn = [pairs(keep,1), ones(nnz(keep),2), ...
          pairs(keep,2), ones(nnz(keep),2), convertTo(Tm(keep), transUnit)];

% ---- fracture connections, straight from the NNC list -------------------
nncCells = G.nnc.cells;
fracConn = [nncCells(:,1), ones(size(nncCells,1),2), ...
            nncCells(:,2), ones(size(nncCells,1),2), ...
            convertTo(G.nnc.T(:), transUnit)];

model.conn = [mmConn; fracConn];
model.mult = zeros(0,5);            % folded into the transmissibilities above
model.isolate = struct('keyword', {'*TRANSI'}, ...
    'range', {[1 nTot 1 1 1 1]}, 'value', {0});

model.meta.counts = struct( ...
    'cells', nTot, 'blocks', nTot, 'inactive', 0, ...
    'matrixCells', nMat, 'fractureCells', nTot - nMat, ...
    'paddingLayers', 0, ...
    'connections', size(model.conn,1), ...
    'multipliers', 0, ...
    'sealedFaces', nnz(internal) - nnz(keep));
end

%--------------------------------------------------------------------------
function Tface = localFaceTrans(Gm, rockM)
% Two-point face transmissibility: the harmonic average of the two
% half-transmissibilities either side of the face.
T = computeTrans(Gm, rockM);
faceOfHalf = Gm.cells.faces(:, 1);
Tface = 1 ./ accumarray(faceOfHalf, 1 ./ max(T, eps), [Gm.faces.num, 1]);
end

%--------------------------------------------------------------------------
function fracIJK = localPlaceFractureCells(G, nMat, nFrac, matIJK, NX, NY, NZ)
% Give every fracture cell an (i,j) under its host matrix column.
%
% An EDFM fracture cell is created where a fracture plane cuts one matrix
% cell, and the matrix-fracture NNC records that pairing. Placing the
% padding block in the host's column keeps it where the fracture actually
% is, so the exported grid can be inspected in Builder. Cells sharing a
% column stack downwards through successive padding layers.

host = zeros(nFrac, 1);
if isfield(G, 'nnc') && isfield(G.nnc, 'cells') && ~isempty(G.nnc.cells)
    c1 = G.nnc.cells(:,1);
    c2 = G.nnc.cells(:,2);
    % matrix-fracture pairs, either orientation
    isMF = (c1 <= nMat) & (c2 > nMat);
    isFM = (c2 <= nMat) & (c1 > nMat);
    host(c2(isMF) - nMat) = c1(isMF);
    host(c1(isFM) - nMat) = c2(isFM);
end

% A fracture cell ALWAYS takes its host column. Cells sharing a column stack
% downwards through successive padding layers, so the padding-layer count is
% the deepest such stack.
%
% Spilling a displaced cell sideways into a neighbouring free column looks
% attractive - it keeps the padding to one layer - but it is wrong. The block
% still carries the fracture's true depth, while the column it landed in sits
% at a different depth, so the block floats above or below the reservoir.
% Keeping the host column guarantees the depth lies inside that column's own
% vertical extent. The extra layers are the cost of a grid you can look at.

used = zeros(NX, NY);
fracIJK = zeros(nFrac, 3);
fallback = 0;
for q = 1:nFrac
    if host(q) > 0
        i = matIJK(host(q), 1);
        j = matIJK(host(q), 2);
    else
        % Fracture cells with no matrix contact (fracture-fracture only)
        % still need a valid unique address.
        fallback = fallback + 1;
        lin = mod(fallback - 1, NX*NY) + 1;
        j = ceil(lin/NX);
        i = lin - (j-1)*NX;
    end
    used(i,j) = used(i,j) + 1;
    fracIJK(q, :) = [i, j, NZ + used(i,j)];
end
end

%--------------------------------------------------------------------------
function [coord, zcorn] = localExtendCornerPoint(grdecl, NX, NY, NZ, ...
    nzExtra, fracIJK, zFrac, dzPad)
% Append nzExtra padding layers to a corner-point description.
%
% Each used padding block gets corners at zFrac +/- dzPad/2, so its centre
% is the true depth of the fracture cell it represents. Unused slots are
% stacked below the reservoir base and are *NULL in the deck.
%
% ZCORN is indexed (2*NX, 2*NY, 2*NZ): for block (i,j,k) the four top
% corners are (2i-1:2i, 2j-1:2j, 2k-1) and the four bottom ones 2k.

zcorn0 = reshape(grdecl.ZCORN(:), 2*NX, 2*NY, 2*NZ);
zbot = zcorn0(:,:,end);

nzTot = NZ + nzExtra;
zc = zeros(2*NX, 2*NY, 2*nzTot);
zc(:,:,1:2*NZ) = zcorn0;

% Default filler for the padding: a flat stack below the reservoir base.
for L = 1:nzExtra
    zc(:,:,2*NZ + 2*L - 1) = zbot + (L-1)*dzPad;
    zc(:,:,2*NZ + 2*L    ) = zbot +  L   *dzPad;
end

% Used slots override that filler with the fracture's own depth.
half = 0.5*dzPad;
for q = 1:size(fracIJK, 1)
    i = fracIJK(q,1); j = fracIJK(q,2); k = fracIJK(q,3);
    ii = (2*i-1):(2*i);
    jj = (2*j-1):(2*j);
    zc(ii, jj, 2*k-1) = zFrac(q) - half;
    zc(ii, jj, 2*k  ) = zFrac(q) + half;
end
zcorn = zc(:);

% Pillars must span the deepest corner now in use, otherwise the padding
% corners fall outside the pillar segment that defines them.
coord = grdecl.COORD(:);
zTopPillar = coord(3:6:end);
zBotPillar = coord(6:6:end);
deepest = max(zcorn) + dzPad;
zBotPillar = max(zBotPillar, deepest);
% Keep pillars non-degenerate.
zBotPillar = max(zBotPillar, zTopPillar + dzPad);
coord(6:6:end) = zBotPillar;
end

%--------------------------------------------------------------------------
function ijk = localUBA(c, nMat, matIJK, fracIJK)
if c <= nMat
    ijk = matIJK(c, :);
else
    ijk = fracIJK(c - nMat, :);
end
end

%--------------------------------------------------------------------------
function [mult, sealConn] = localSeals(G, matIJK, tol)
% pEDFM/CPEDFM matrix-face seals.
%
% Logical neighbours become a *TRANS?*MOD on the lower-indexed block.
% Anything else - a fault-juxtaposed pair - has no direction to modify, so
% it becomes a *SCONNECT carrying a NEGATIVE value, which CMG reads as a
% transmissibility multiplier on the existing connection.
mult = zeros(0,5);
sealConn = zeros(0,7);

eta = getEmbeddedFractureMultipliers(G, tol);
sealed = find(eta < 1 - 1e-12);
sealed = sealed(sealed <= G.Matrix.faces.num);

for q = 1:numel(sealed)
    f = sealed(q);
    nb = G.Matrix.faces.neighbors(f, :);
    if any(nb == 0), continue; end
    a = matIJK(nb(1), :);
    b = matIJK(nb(2), :);
    d = b - a;
    if nnz(d) == 1 && max(abs(d)) == 1
        dir = find(d ~= 0);
        if d(dir) > 0, low = a; else, low = b; end
        mult(end+1, :) = [dir, low, eta(f)]; %#ok<AGROW>
    else
        sealConn(end+1, :) = [a, b, -eta(f)]; %#ok<AGROW>
    end
end
end

%--------------------------------------------------------------------------
function model = localBuildMatrixOnly(model, G, rock, gridInfo, opt, permUnit)

switch gridInfo.family
    case {'cart', 'vari'}
        nTot = prod(gridInfo.dims);
        model.grid.kind = gridInfo.family;
        if ~isempty(opt.datum)
            % Shift the whole section, not just its top, so the contact
            % depths derived from zbot stay below the reservoir.
            shift = opt.datum - model.grid.ztop;
            model.grid.ztop = opt.datum;
            model.grid.zbot = model.grid.zbot + shift;
        end

    case 'corner'
        nTot = prod(gridInfo.dims);
        model.grid.kind = 'corner';
        assert(~isempty(opt.grdecl), 'buildCMGModel:needGrdecl', ...
            ['Corner-point export needs the grdecl (COORD/ZCORN); MRST ', ...
             'discards them in processGRDECL. Pass ''grdecl'', grdecl.']);
        model.grid.coord = opt.grdecl.COORD(:);
        model.grid.zcorn = opt.grdecl.ZCORN(:);

    case 'unstructured'
        % Connection-list pseudo-grid: CMG has no unstructured grid type, so
        % the cells become an Nc x 1 x 1 placeholder and every internal face
        % is emitted as an explicit connection.
        nTot = G.cells.num;
        model.grid.kind = 'connlist';
        model.grid.dims = [nTot, 1, 1];

    otherwise
        error('buildCMGModel:family', 'Unhandled grid family %s', gridInfo.family);
end

% ---- active-cell map ---------------------------------------------------
if strcmp(model.grid.kind, 'connlist')
    cellSlot = (1:G.cells.num)';
    active = true(nTot, 1);
else
    if isfield(G.cells, 'indexMap')
        cellSlot = G.cells.indexMap(:);
    else
        cellSlot = (1:G.cells.num)';
    end
    active = false(nTot, 1);
    active(cellSlot) = true;
end

% ---- property arrays over the full index space -------------------------
perm = rock.perm;
if size(perm, 2) == 1, perm = repmat(perm, 1, 3); end

permi = zeros(nTot, 1); permj = permi; permk = permi; poro = permi;
permi(cellSlot) = convertTo(perm(:, 1), permUnit);
permj(cellSlot) = convertTo(perm(:, min(2, end)), permUnit);
permk(cellSlot) = convertTo(perm(:, min(3, end)), permUnit);
poro(cellSlot)  = rock.poro(:);

model.props = struct('permi', permi, 'permj', permj, 'permk', permk, ...
    'poro', poro, 'active', active, 'nTot', nTot);

% ---- pore-volume and depth handling for the placeholder pseudo-grid ----
model.props.volmod = [];
model.props.porvol = [];
model.props.depth = [];
if strcmp(model.grid.kind, 'connlist')
    % Placeholder geometry carries no physics: *PORVOL restores the true
    % pore volume directly (overriding porosity x bulk volume) and an
    % explicit depth restores the gravity head.
    dz = 1.0;
    model.props.porvol = poreVolume(G, rock);
    % *DEPTH-TOP wants the depth of each block's top face, not its centre.
    model.props.depth = G.cells.centroids(:, 3) - 0.5*dz;
    model.grid.di = 1; model.grid.dj = 1; model.grid.dk = dz;
    model.grid.ztop = min(model.props.depth);
    model.grid.zbot = max(model.props.depth) + dz;
end

% ---- connections -------------------------------------------------------
if strcmp(model.grid.kind, 'connlist')
    model.conn = localInternalFaceConnections(G, rock);
else
    model.conn = zeros(0, 7);
end
model.mult = zeros(0, 5);
model.isolate = struct('keyword', {}, 'range', {}, 'value', {});

model.meta.counts = struct( ...
    'cells', G.cells.num, ...
    'blocks', nTot, ...
    'inactive', nTot - nnz(active), ...
    'connections', size(model.conn, 1), ...
    'multipliers', size(model.mult, 1));
end

%--------------------------------------------------------------------------
function wells = localBuildWells(G, W, model)
% Translate an MRST well structure into CMG completion records.
%
% The well index is carried across explicitly. CMG's *PERF *WI takes the
% geometric part of the well index in md.m - the same unit and the same
% definition (q = wi * mobility * dp) that MRST's W.WI uses - so the value
% transfers directly. Letting CMG recompute it from *GEOMETRY would apply a
% radial-flow model to whatever block shape the export produced, which is
% wrong for corner-point cells and meaningless for a placeholder pseudo-grid.

wells = struct('name', {}, 'isInjector', {}, 'type', {}, 'val', {}, ...
    'comp', {}, 'uba', {}, 'wi', {});
if isempty(W)
    return
end

transUnit = milli*darcy*meter;
dims = model.grid.dims;
isConnList = strcmp(model.grid.kind, 'connlist');

% On a fracture grid the address of a cell depends on which side of the
% matrix/fracture split it falls, so reuse the same map the connections use.
hasFrac = isfield(model.grid, 'matIJK') && isfield(model.grid, 'fracIJK');
if hasFrac
    matIJK = model.grid.matIJK;
    fracIJK = model.grid.fracIJK;
    nMat = size(matIJK, 1);
end

if isfield(G.cells, 'indexMap')
    indexMap = G.cells.indexMap(:);
else
    indexMap = (1:G.cells.num)';
end

for k = 1:numel(W)
    w = W(k);
    cells = w.cells(:);

    if isConnList
        uba = [cells, ones(numel(cells), 2)];
    elseif hasFrac
        uba = zeros(numel(cells), 3);
        for c = 1:numel(cells)
            uba(c, :) = localUBA(cells(c), nMat, matIJK, fracIJK);
        end
    else
        [ii, jj, kk] = ind2sub(dims, indexMap(cells));
        uba = [ii(:), jj(:), kk(:)];
    end

    if isfield(w, 'sign') && ~isempty(w.sign)
        isInj = w.sign > 0;
    else
        isInj = false;
    end

    comp = 'WATER';
    if isfield(w, 'compi') && numel(w.compi) >= 2 && w.compi(2) > 0.5
        comp = 'OIL';
    end

    if isfield(w, 'WI') && ~isempty(w.WI)
        wi = convertTo(abs(w.WI(:)), transUnit);
    else
        wi = -ones(numel(cells), 1);   % '-' lets CMG compute it
    end

    entry = struct('name', localWellName(w, k), 'isInjector', isInj, ...
        'type', lower(char(string(w.type))), 'val', w.val, ...
        'comp', comp, 'uba', uba, 'wi', wi);
    wells(end+1) = entry; %#ok<AGROW>
end
end

%--------------------------------------------------------------------------
function name = localWellName(w, k)
if isfield(w, 'name') && ~isempty(w.name)
    name = char(string(w.name));
else
    name = sprintf('WELL%d', k);
end
% CMG well names are quoted strings; keep them free of the quote character.
name = strrep(name, '''', '_');
end

%--------------------------------------------------------------------------
function conn = localInternalFaceConnections(G, rock)
% Every internal face of an unstructured grid becomes one explicit
% connection, carrying the same two-point transmissibility MRST would use.
N = G.faces.neighbors;
internal = all(N > 0, 2);
pairs = N(internal, :);

T = computeTrans(G, rock);              % one half-trans per cell-face
faceOfHalf = G.cells.faces(:, 1);
% Harmonic average of the two half-transmissibilities = face transmissibility.
Tface = 1 ./ accumarray(faceOfHalf, 1 ./ max(T, eps), [G.faces.num, 1]);
Tface = Tface(internal);

transUnit = milli*darcy*meter;
conn = [pairs(:, 1), ones(size(pairs, 1), 2), ...
        pairs(:, 2), ones(size(pairs, 1), 2), ...
        convertTo(Tface, transUnit)];
conn = conn(:, [1 2 3 4 5 6 7]);
end
