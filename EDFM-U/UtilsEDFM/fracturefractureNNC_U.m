function [G, fracplanes] = fracturefractureNNC_U(G, fracplanes, tol, varargin)
% Fracture-fracture NNCs for an unstructured EDFM grid.
%
% SYNOPSIS:
%   [G, fracplanes] = fracturefractureNNC_U(G, fracplanes, tol)
%
% PARAMETERS:
%   G          - global EDFM grid after EDFMgrid_U and fracturematrixNNC_U.
%   fracplanes - the struct array returned by EDFMgrid_U (carries .normal).
%   tol        - geometric tolerance.
%
% RETURNS:
%   G          - with 'fracfrac' rows appended to G.nnc.
%   fracplanes - with .intersects listing each plane's intersecting planes.
%
% WHY THIS IS NOT A WRAPPER
%   The Cartesian routine fracturefractureShaleNNCs3D cannot run on an
%   unstructured grid. It builds its candidate filter with
%   hfm/findneighbours, which opens with assert(isfield(G,'cartDims')) and
%   constructs a 26-cell 3x3x3 cloud from logical indexing - neither is
%   available on a PEBI or tetrahedral grid. Its inner routine fgridNNCs3D
%   additionally assumes that row i of matrix_connection.cells belongs to
%   fracture cell mod(i,Nf), a convention that only holds when every
%   intersection is of the same kind; EDFMgrid_U records the mapping
%   explicitly in matrix_connection.fraccells so that face-coplanar and
%   interior intersections can coexist.
%
%   So the two grid-dependent pieces are replaced here - the cloud by the
%   node-sharing adjacency of cellNeighbourhood_U, the implicit ordering by
%   the explicit fraccells map - while the geometry is the same as the
%   Cartesian release: PEBIPEBIintersection3D for the intersection segment
%   and its average distances, and the Moinfar et al. (2013) harmonic
%   transmissibility
%
%       T_i = k_i * w_i * L / df_i,   T = T_i*T_j / (T_i + T_j).
%
% SEE ALSO:
%   EDFMgrid_U, fracturematrixNNC_U, cellNeighbourhood_U

opt = struct('Verbose', true);
opt = merge_options(opt, varargin{:});
verbose = opt.Verbose;

assert(isfield(G, 'FracGrid') && isfield(G, 'Matrix'), ...
    'G must be a global EDFM grid (run EDFMgrid_U first).');
assert(isfield(G.nnc, 'cells'), ...
    'Run fracturematrixNNC_U before fracturefractureNNC_U.');

t1 = tic;
fracNames = fieldnames(G.FracGrid);
nFrac = numel(fracNames);
if nFrac < 2
    dispif(verbose, 'fracturefractureNNC_U: fewer than two fracture grids; nothing to do.\n');
    return
end

dispif(verbose, 'Processing fracture-fracture NNCs (%d planes)...\n', nFrac);

% node-sharing adjacency of the matrix grid, built once
A = cellNeighbourhood_U(G.Matrix);

newCells = zeros(0,2);
newT     = zeros(0,1);

for a = 1:(nFrac-1)
    Fa = G.FracGrid.(fracNames{a});
    ka = Fa.fracplanenumber;

    for b = (a+1):nFrac
        Fb = G.FracGrid.(fracNames{b});
        kb = Fb.fracplanenumber;

        if ismember(b, G.FracGrid.(fracNames{a}).fracgridnnc)
            dispif(verbose, '  %s vs %s checked previously.\n', ...
                fracNames{a}, fracNames{b});
            continue
        end
        G.FracGrid.(fracNames{a}).fracgridnnc(end+1) = b;
        G.FracGrid.(fracNames{b}).fracgridnnc(end+1) = a;

        % same-set fractures are taken not to intersect (as in the
        % Cartesian release)
        if isfield(fracplanes, 'SetID') && ...
                isequal(fracplanes(ka).SetID, fracplanes(kb).SetID)
            continue
        end

        % parallel planes cannot intersect, unless they belong to the same
        % piecewise-planar fracture
        sameFracture = isfield(fracplanes, 'fractureID') && ...
            isequal(fracplanes(ka).fractureID, fracplanes(kb).fractureID);
        if ~sameFracture && ...
                norm(cross(fracplanes(ka).normal, fracplanes(kb).normal)) < tol
            continue
        end

        % do the two planes intersect at all?
        if ~PEBIPEBIintersection3D(fracplanes(ka).points, fracplanes(kb).points, tol)
            continue
        end

        [cells_ab, T_ab] = planePairNNCs(Fa, Fb, ...
            fracplanes(ka).aperture, fracplanes(kb).aperture, A, tol);

        if ~isempty(T_ab)
            newCells = [newCells; cells_ab]; %#ok<AGROW>
            newT     = [newT;     T_ab];     %#ok<AGROW>
            fracplanes(ka).intersects = appendUnique(fracplanes, ka, kb);
            fracplanes(kb).intersects = appendUnique(fracplanes, kb, ka);
            dispif(verbose, '  %s vs %s: %d NNCs\n', ...
                fracNames{a}, fracNames{b}, numel(T_ab));
        end
    end
end

if isempty(newT)
    dispif(verbose, 'fracturefractureNNC_U: no fracture-fracture NNCs found.\n');
    return
end

% de-duplicate (a pair can be reached from either side)
[newCells, ia] = unique(newCells, 'rows');
newT = newT(ia);
n = numel(newT);

G.nnc.cells      = [G.nnc.cells;      newCells];
G.nnc.T          = [G.nnc.T;          newT];
G.nnc.type       = [G.nnc.type;       repmat({'fracfrac'}, n, 1)];
G.nnc.area       = [G.nnc.area;       zeros(n,1)];
G.nnc.normal     = [G.nnc.normal;     zeros(n,3)];
G.nnc.planepoint = [G.nnc.planepoint; zeros(n,3)];
if isfield(G.nnc, 'CI'),   G.nnc.CI   = [G.nnc.CI;   nan(n,1)]; end
if isfield(G.nnc, 'davg'), G.nnc.davg = [G.nnc.davg; nan(n,1)]; end

dispif(verbose, 'Fracture-fracture NNCs done: %d connections in %.3g s.\n', ...
    n, toc(t1));
end

% =========================================================================
function [cells, T] = planePairNNCs(Fa, Fb, apA, apB, A, tol)
% NNCs between the cells of two intersecting fracture grids.

% host matrix cell of each fracture cell (any one will do for the filter:
% a face-coplanar fracture cell has two, and both are adjacent to it)
hostA = hostCells(Fa);
hostB = hostCells(Fb);

% candidate fracture-cell pairs: their host matrix cells must touch
nA = Fa.cells.num;
pairs = zeros(0,2);
for ca = 1:nA
    if isnan(hostA(ca)), continue; end
    cloud = find(A(:, hostA(ca)));
    cand  = find(ismember(hostB, cloud));
    if ~isempty(cand)
        pairs = [pairs; repmat(ca, numel(cand), 1), cand(:)]; %#ok<AGROW>
    end
end
if isempty(pairs)
    cells = zeros(0,2); T = zeros(0,1);
    return
end
pairs = unique(pairs, 'rows');

% rigorous step: intersect the fracture-cell polygons
polyA = cellPolygons(Fa);
polyB = cellPolygons(Fb);

np = size(pairs,1);
cells = zeros(np, 2);
T     = zeros(np, 1);
keep  = false(np, 1);

for m = 1:np
    ca = pairs(m,1); cb = pairs(m,2);
    [xsect, ~, xlength, df] = PEBIPEBIintersection3D(polyA{ca}, polyB{cb}, tol);
    if ~xsect || xlength <= tol || any(df <= tol)
        continue
    end
    Ta = Fa.rock.perm(ca) * apA * xlength / df(1);
    Tb = Fb.rock.perm(cb) * apB * xlength / df(2);

    cells(m,:) = sort([ca + Fa.cells.start - 1, cb + Fb.cells.start - 1]);
    T(m)       = (Ta*Tb)/(Ta + Tb);
    keep(m)    = true;
end

cells = cells(keep,:);
T     = T(keep);
end

% -------------------------------------------------------------------------
function h = hostCells(F)
% One representative host matrix cell per fracture cell.
h = nan(F.cells.num, 1);
mc = F.matrix_connection;
h(mc.fraccells) = mc.cells;      % later rows overwrite, either is valid
end

% -------------------------------------------------------------------------
function P = cellPolygons(F)
% Mid-plane polygon of each fracture cell, as PEBIPEBIintersection3D wants.
% Fracture grids are one layer thick, so a cell's node list is the bottom
% ring followed by the top ring; nodes2D holds the un-extruded mid-plane.
[cn, cpos] = gridCellNodes(F, 1:F.cells.num);
P = cell(F.cells.num, 1);
for c = 1:F.cells.num
    idx = cn(cpos(c):(cpos(c+1)-1));
    P{c} = F.nodes2D.coords(idx(1:numel(idx)/2), :);
end
end

% -------------------------------------------------------------------------
function v = appendUnique(fracplanes, k, other)
if isfield(fracplanes, 'intersects') && ~isempty(fracplanes(k).intersects)
    v = unique([fracplanes(k).intersects(:); other]);
else
    v = other;
end
end
