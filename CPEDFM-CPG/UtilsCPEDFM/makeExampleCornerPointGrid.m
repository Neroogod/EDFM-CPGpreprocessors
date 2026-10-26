function [G, grdecl] = makeExampleCornerPointGrid(dims, physDims, varargin)
%MAKEEXAMPLECORNERPOINTGRID Construct a distorted GRDECL corner-point grid.
%
%   [G, GRDECL] = makeExampleCornerPointGrid(DIMS, PHYSDIMS)
%   constructs an NX-by-NY-by-NZ grid through COORD/ZCORN and
%   processGRDECL. Unlike cartGrid/tensorGrid, the default grid has skewed
%   pillars and non-planar layer surfaces, and therefore exercises the
%   corner-point geometry path used by EDFMGridCPG.
%
% OPTIONS
%   pillarSkew - lateral displacement at the bottom, as a fraction of the
%                mean lateral cell size (default 0.20).
%   layerWarp  - layer-surface relief as a fraction of the mean cell
%                thickness (default 0.12).
%   faultDrop  - vertical offset applied to pillars/corners on the second
%                half of the I direction (default zero).
%   actnum     - optional DIMS-sized logical/numeric active-cell array.
%   preserveCpNodes - pass PreserveCpNodes to processGRDECL (default true).

% The deformation vanishes on the lateral boundary and keeps layer order.

opt = struct('pillarSkew', 0.20, 'layerWarp', 0.12, ...
    'faultDrop', 0, 'actnum', [], 'preserveCpNodes', true);
opt = merge_options(opt, varargin{:});

dims = reshape(dims, 1, []);
physDims = reshape(physDims, 1, []);
assert(numel(dims) == 3 && all(dims >= 1) && all(dims == round(dims)), ...
    'dims must contain three positive integers.');
assert(numel(physDims) == 3 && all(physDims > 0), ...
    'physDims must contain three positive lengths.');

nx = dims(1); ny = dims(2); nz = dims(3);
xv = linspace(0, physDims(1), nx + 1);
yv = linspace(0, physDims(2), ny + 1);
zv = linspace(0, physDims(3), nz + 1);
[X, Y] = ndgrid(xv, yv);

% Pillar endpoints. The top remains at its logical location; bottom
% endpoints are shifted smoothly to make non-vertical, non-parallel
% pillars without moving the outer boundary.
xi = X/physDims(1);
yi = Y/physDims(2);
window = sin(pi*xi).*sin(pi*yi);
dx = opt.pillarSkew*(physDims(1)/nx)*window.*sin(2*pi*yi);
dy = opt.pillarSkew*(physDims(2)/ny)*window.*cos(2*pi*xi);

topX = X;             topY = Y;
botX = X + dx;        botY = Y + dy;
topZ = zeros(size(X)); botZ = physDims(3)*ones(size(X));

if opt.faultDrop ~= 0
    faultSide = X > physDims(1)/2;
    topZ(faultSide) = topZ(faultSide) + opt.faultDrop;
    botZ(faultSide) = botZ(faultSide) + opt.faultDrop;
end

nPillars = (nx + 1)*(ny + 1);
coord = zeros(nPillars, 6);
coord(:,1) = topX(:); coord(:,2) = topY(:); coord(:,3) = topZ(:);
coord(:,4) = botX(:); coord(:,5) = botY(:); coord(:,6) = botZ(:);

% ZCORN uses repeated logical corner indices. Apply a smooth relief that
% vanishes at top/bottom, and a consistent fault offset when requested.
[XI, YI, ZI] = ndgrid(1 + fix((1:2*nx)/2), ...
                      1 + fix((1:2*ny)/2), ...
                      1 + fix((1:2*nz)/2));
Xc = xv(XI); Yc = yv(YI); Zc = zv(ZI);
relief = opt.layerWarp*(physDims(3)/nz).* ...
    sin(pi*Xc/physDims(1)).*sin(pi*Yc/physDims(2)).* ...
    sin(pi*Zc/physDims(3));
Zc = Zc + relief;
if opt.faultDrop ~= 0
    Zc(Xc > physDims(1)/2) = Zc(Xc > physDims(1)/2) + opt.faultDrop;
end

grdecl = struct('cartDims', dims, 'COORD', reshape(coord.', [], 1), ...
    'ZCORN', Zc(:));
if isempty(opt.actnum)
    grdecl.ACTNUM = ones(prod(dims), 1, 'int32');
else
    assert(numel(opt.actnum) == prod(dims), 'actnum must match prod(dims).');
    grdecl.ACTNUM = int32(opt.actnum(:) ~= 0);
end

G = processGRDECL(grdecl, 'PreserveCpNodes', opt.preserveCpNodes);
G = computeGeometry(G);
G.type = [G.type, {'synthetic-corner-point'}];
end
