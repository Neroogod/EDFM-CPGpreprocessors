function [G, fracplanes] = EDFMgrid_U(G, fracplanes, varargin)
% EDFM fracture gridding on general unstructured grids.
%
% SYNOPSIS:
%   [G, fracplanes] = EDFMgrid_U(G, fracplanes)
%   [G, fracplanes] = EDFMgrid_U(G, fracplanes, 'Tolerance', 1e-5)
%
% PARAMETERS:
%   G          - matrix grid through computeGeometry, carrying G.rock.
%                Any grid whose cells are convex polyhedra: Cartesian,
%                PEBI/Voronoi (2D extruded or 3D), tetrahedral, prismatic,
%                corner-point.
%   fracplanes - struct array with fields
%                  points   - Nx3 ordered vertices of a planar polygon
%                  aperture - fracture aperture [m]
%                  poro     - fracture porosity
%                  perm     - fracture permeability
%
% OPTIONAL PARAMETERS:
%   'Tolerance'   - geometric tolerance, default 1e-5.
%   'fracturelist'- indices of fracplanes to process, default all.
%   'Verbose'     - progress printing, default true.
%
% RETURNS:
%   G          - global EDFM grid with G.Matrix, G.FracGrid.FracN and G.nnc,
%                assembled by assembleShaleGlobalGrid.
%   fracplanes - annotated with .normal, .fracgrid and .processed.
%
% METHOD:
%   Each fracture polygon is clipped against every candidate matrix cell
%   with clipPolyByCell_U (half-space clipping, exact for convex cells).
%   The surviving polygons become the cells of a 2D fracture grid built by
%   fractureplanegeneralgrid and extruded to the aperture.
%
% FACE-COPLANAR FRACTURES:
%   A fracture lying in the plane of a matrix face borders TWO matrix
%   cells. This routine creates ONE fracture cell for such a face and
%   registers it against BOTH neighbours (type 1, 'fracmat boundary'), so
%   the fracture surface area is split correctly between them. A fracture
%   strictly inside a cell borders that one cell on both of its sides and
%   is registered once (type 2, 'fracmat interior'); its surface area is
%   twice the intersection area. fracMatNNC_U applies that factor from the
%   recorded type rather than by counting rows, so mixed boundary/interior
%   fractures - which the Cartesian implementation cannot express - are
%   handled correctly.
%
% RETURNS (G.FracGrid.FracN.matrix_connection):
%   cells          - matrix cell of each fracture-matrix connection
%   fraccells      - local fracture-cell index of each connection
%   area           - intersection area of each connection
%   type           - 1 boundary (face-coplanar), 2 interior
%   fraccellpoints - clipped polygon of each FRACTURE CELL (not connection)
%   normal         - fracture plane normal, one row per connection
%   planepoint     - a point on the fracture plane, one row per connection
%
% SEE ALSO:
%   clipPolyByCell_U, fracturematrixNNC_U, fracMatNNC_U

opt = struct('Tolerance', 1e-5, ...
             'fracturelist', 1:numel(fracplanes), ...
             'Verbose', true);
opt = merge_options(opt, varargin{:});
tol     = opt.Tolerance;
verbose = opt.Verbose;

t1 = tic;

%% Restructure G on first call
if ~all(isfield(G, {'Matrix', 'FracGrid'}))
    Gtemp          = G;
    Gtemp.Matrix   = G;
    Gtemp.FracGrid = struct;
    Gtemp.nnc      = [];
    G              = Gtemp;
end

if ~isfield(fracplanes, 'processed')
    for i = 1:numel(fracplanes)
        fracplanes(i).processed = false;
    end
end

fracturelist = opt.fracturelist(~[fracplanes(opt.fracturelist).processed]);
if isempty(fracturelist)
    dispif(verbose, 'EDFMgrid_U: nothing to process.\n');
    return
end

%% Plane normals
for i = fracturelist
    P      = [fracplanes(i).points; fracplanes(i).points(1,:)];
    dP     = diff(P, 1);
    normal = cross(dP(1,:), dP(2,:));
    assert(norm(normal) > tol, ...
        'fracplane %d is degenerate: its first three vertices are collinear.', i);
    fracplanes(i).normal = normal / norm(normal);
end

Gm   = G.Matrix;
nMat = Gm.cells.num;

%% Per-cell node bounding boxes (computed once, reused by every fracture)
[cn, cpos] = gridCellNodes(Gm, 1:nMat);
cellNo     = rldecode((1:nMat)', diff(cpos));
X          = Gm.nodes.coords(cn, :);
cLo = zeros(nMat,3); cHi = zeros(nMat,3);
for k = 1:3
    cLo(:,k) = accumarray(cellNo, X(:,k), [nMat 1], @min);
    cHi(:,k) = accumarray(cellNo, X(:,k), [nMat 1], @max);
end

%% Intersect each fracture with the grid
Fstore = cell(numel(fracturelist), 1);
for k = 1:numel(fracturelist)
    i           = fracturelist(k);
    points      = fracplanes(i).points;
    planenormal = fracplanes(i).normal;
    aperture    = fracplanes(i).aperture;

    dispif(verbose, 'EDFMgrid_U: fracture %d of %d ... ', k, numel(fracturelist));

    % coarse filter on overlapping bounding boxes
    bbLo = min(points, [], 1) - tol;
    bbHi = max(points, [], 1) + tol;
    cand = find(all(cHi >= bbLo, 2) & all(cLo <= bbHi, 2));

    fcellpoly = cell(numel(cand), 1);   % polygon of each fracture cell
    fcellarea = zeros(numel(cand), 1);
    connCell  = zeros(2*numel(cand), 1);
    connFrac  = zeros(2*numel(cand), 1);
    connType  = zeros(2*numel(cand), 1);
    faceCell  = zeros(Gm.faces.num, 1); % face -> fracture cell already made
    nF = 0; nC = 0;

    for q = 1:numel(cand)
        c = cand(q);
        [xarea, poly, onFace] = clipPolyByCell_U(points, Gm, c, tol);
        if xarea <= 0
            continue
        end

        if onFace > 0
            if faceCell(onFace) > 0
                continue                       % other neighbour did it
            end
            nF = nF + 1;
            fcellpoly{nF} = poly;
            fcellarea(nF) = xarea;
            faceCell(onFace) = nF;

            nb = Gm.faces.neighbors(onFace,:);
            nb = nb(nb > 0);
            for b = 1:numel(nb)
                nC = nC + 1;
                connCell(nC) = nb(b);
                connFrac(nC) = nF;
                connType(nC) = 1;              % fracmat boundary
            end
        else
            nF = nF + 1;
            fcellpoly{nF} = poly;
            fcellarea(nF) = xarea;

            nC = nC + 1;
            connCell(nC) = c;
            connFrac(nC) = nF;
            connType(nC) = 2;                  % fracmat interior
        end
    end

    if nF == 0
        warning('EDFMgrid_U:noIntersection', ...
            'Fracture %d does not intersect the grid; skipped.', i);
        dispif(verbose, 'no intersection.\n');
        Fstore{k} = [];
        continue
    end

    fcellpoly = fcellpoly(1:nF);
    fcellarea = fcellarea(1:nF);
    connCell  = connCell(1:nC);
    connFrac  = connFrac(1:nC);
    connType  = connType(1:nC);

    % Build the fracture-plane grid from the clipped polygons
    V = vertcat(fcellpoly{:});
    C = cellfun(@(p) 1:size(p,1), fcellpoly, 'UniformOutput', false);
    for j = 2:numel(C)
        C{j} = C{j} + C{j-1}(end);
    end
    Fgrid = fractureplanegeneralgrid(V, C, points, planenormal, aperture, tol);

    assert(Fgrid.cells.num == nF, ...
        ['EDFMgrid_U: fracture %d produced %d grid cells for %d clipped ', ...
         'polygons. Try a tighter ''Tolerance''.'], i, Fgrid.cells.num, nF);

    % volume diagnostic: the extruded fracture grid should reproduce the
    % total clipped area times the aperture
    volErr = abs(sum(Fgrid.cells.volumes) - sum(fcellarea)*aperture);
    if volErr > tol * max(sum(fcellarea)*aperture, 1)
        warning('EDFMgrid_U:volumeMismatch', ...
            ['Fracture %d: gridded volume differs from area*aperture by ', ...
             '%.3g (%.2g%%). Consider a tighter ''Tolerance''.'], ...
            i, volErr, 100*volErr/(sum(fcellarea)*aperture));
    end

    Fgrid.matrix_connection.cells          = connCell;
    Fgrid.matrix_connection.fraccells      = connFrac;
    Fgrid.matrix_connection.area           = fcellarea(connFrac);
    Fgrid.matrix_connection.type           = connType;
    Fgrid.matrix_connection.fraccellpoints = fcellpoly;
    Fgrid.matrix_connection.normal         = repmat(planenormal, nC, 1);
    Fgrid.matrix_connection.planepoint     = repmat(points(1,:), nC, 1);

    Fstore{k} = Fgrid;
    dispif(verbose, '%d fracture cells, %d connections (%d boundary).\n', ...
        nF, nC, sum(connType == 1));
end

kept         = ~cellfun('isempty', Fstore);
Fstore       = Fstore(kept);
fracturelist = fracturelist(kept);

if isempty(Fstore)
    error('EDFMgrid_U:noFractures', 'No fracture intersected the grid.');
end

%% Append fracture grids to G.FracGrid
FracGrid_new = G.FracGrid;
cstart   = G.cells.num + 1;
fstart   = G.faces.num + 1;
nstart   = G.nodes.num + 1;
lastfrac = numel(fieldnames(G.FracGrid));

for j = 1:numel(fracturelist)
    i  = fracturelist(j);
    fn = ['Frac', num2str(j + lastfrac)];

    FracGrid_new.(fn)             = Fstore{j};
    FracGrid_new.(fn).cells.start = cstart;
    FracGrid_new.(fn).faces.start = fstart;
    FracGrid_new.(fn).nodes.start = nstart;
    cstart = cstart + FracGrid_new.(fn).cells.num;
    fstart = fstart + FracGrid_new.(fn).faces.num;
    nstart = nstart + FracGrid_new.(fn).nodes.num;

    nfc = FracGrid_new.(fn).cells.num;
    FracGrid_new.(fn).rock.perm     = ones(nfc,1) * fracplanes(i).perm;
    FracGrid_new.(fn).rock.poro     = ones(nfc,1) * fracplanes(i).poro;
    FracGrid_new.(fn).rock.isMatrix = zeros(nfc,1);
    if isfield(G.rock, 'ntg')
        ntg = 1;
        if isfield(fracplanes, 'ntg') && ~isempty(fracplanes(i).ntg)
            ntg = fracplanes(i).ntg;
        end
        FracGrid_new.(fn).rock.ntg = ones(nfc,1) * ntg;
    end

    FracGrid_new.(fn).points          = fracplanes(i).points;
    FracGrid_new.(fn).fracplanenumber = i;
    FracGrid_new.(fn).matrixnnc       = false;
    FracGrid_new.(fn).fracgridnnc     = [];
    FracGrid_new.(fn).gridtype        = FracGrid_new.(fn).type;
    FracGrid_new.(fn)                 = rmfield(FracGrid_new.(fn), 'type');

    fracplanes(i).fracgrid  = j + lastfrac;
    fracplanes(i).processed = true;
end
G.FracGrid = FracGrid_new;

%% Assemble the global grid
Gtemp          = G.Matrix;
Gtemp.FracGrid = G.FracGrid;
Gtemp.nnc      = G.nnc;
G              = assembleShaleGlobalGrid(Gtemp);

dispif(verbose, 'EDFMgrid_U: %d fracture grids assembled in %.3g s.\n', ...
    numel(fracturelist), toc(t1));
end
