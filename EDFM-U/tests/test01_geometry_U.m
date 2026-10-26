%% test01_geometry_U
% Unit tests for the EDFM-U geometry kernel: the polygon-vs-cell clip
% (clipPolyByCell_U), the cell-vs-halfspace clip (clipCellByPlane_U) and
% the average normal distance <d> (calcDavg_U).
%
% Checks
%   A. <d> on a Cartesian cell, several plane positions and orientations,
%      against the closed-form integral of |n.(x-x0)| over the cube.
%   B. Volume/centroid of a clipped cell: the two sides reproduce the
%      grid's own volume and centroid, on Cartesian, PEBI and tetrahedral
%      cells.
%   C. Area conservation: the clipped areas of a fracture polygon over all
%      matrix cells sum to the polygon's own area.
%   D. Face-coplanar detection: a fracture placed exactly on an interior
%      face is detected as a boundary intersection by BOTH neighbours and
%      yields a single shared fracture cell.
%   E. Fracture grid volume: sum of fracture cell volumes = area*aperture.
%
% No flow simulation; runs in seconds.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupEDFMU.m'));
addpath(thisDir);                       % boxMeanAbsDist_U (test reference)

tol = 1e-9;
failures = {};

%% A. <d> on a Cartesian cell against the analytic integral
Gc = computeGeometry(cartGrid([1 1 1], [2 3 5]));   % deliberately anisotropic
L  = [2 3 5];

maxErr = 0;
cases = { [0 0 1], [0 0 2.5];
          [0 0 1], [0 0 1.0];
          [0 0 1], [0 0 0.0];
          [1 0 0], [0.7 0 0];
          [1 1 0], [1.0 1.5 0];
          [1 2 3], [1.0 1.5 2.5];
          [3 -1 2], [0.5 2.0 1.0] };
for i = 1:size(cases,1)
    n  = cases{i,1}; p0 = cases{i,2};
    d  = calcDavg_U(Gc, 1, n, p0, tol);
    ex = boxMeanAbsDist_U([0 0 0], L, n, p0);
    maxErr = max(maxErr, abs(d - ex)/max(ex, eps));
end
fprintf('A. Cartesian <d> vs analytic: max rel error %.3e\n', maxErr);
if maxErr < 1e-9
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('A: max rel err %.3e', maxErr); fprintf('   FAIL\n');
end

%% B. clipped volume + centroid reproduce the grid's own values
rng(4);
grids = struct('name', {}, 'G', {});

grids(1).name = 'Cartesian';
grids(1).G    = computeGeometry(cartGrid([4 4 3], [40 40 30]));

G2 = makeLayeredGrid(pebiGrid2D(9, [40 40]), 3);
G2.nodes.coords(:,3) = G2.nodes.coords(:,3) * 10;
grids(2).name = 'extruded PEBI';
grids(2).G    = computeGeometry(G2);

pts = [40*rand(60,2), 30*rand(60,1)];
pts = [pts; combineCorners([0 40], [0 40], [0 30])];
DT  = delaunayTriangulation(pts);
grids(3).name = 'tetrahedral';
grids(3).G    = computeGeometry(tetrahedralGrid(DT.Points, DT.ConnectivityList));

nrm = [0.5 0.7 0.5]; nrm = nrm/norm(nrm);
p0  = [20 20 15];
okB = true;
for g = 1:numel(grids)
    Gg = grids(g).G;
    eV = 0; eC = 0; nCut = 0;
    for c = 1:Gg.cells.num
        Vc = Gg.cells.volumes(c); Cc = Gg.cells.centroids(c,:);
        [Vm, Cm] = clipCellByPlane_U(Gg, c, nrm, p0, 1e-10);
        [Vp, Cp] = clipCellByPlane_U(Gg, c, -nrm, p0, 1e-10);
        eV = max(eV, abs(Vm + Vp - Vc)/Vc);
        if Vm > 0 && Vp > 0
            nCut = nCut + 1;
            Cmix = (Vm*Cm + Vp*Cp)/(Vm + Vp);
            eC = max(eC, norm(Cmix - Cc)/nthroot(Vc,3));
        end
    end
    fprintf('B. %-14s (%4d cells, %d cut): dV %.2e, dC %.2e\n', ...
        grids(g).name, Gg.cells.num, nCut, eV, eC);
    okB = okB && eV < 1e-9 && eC < 1e-9 && nCut > 0;
end
if okB
    fprintf('   PASS\n');
else
    failures{end+1} = 'B: clipped volume/centroid do not close'; fprintf('   FAIL\n');
end

%% C. area conservation of the polygon-vs-cell clip
poly = [ 8  8  4;
        32 12  4;
        32 12 26;
         8  8 26];
Aexact = polyArea3D_U(poly);
okC = true;
for g = 1:numel(grids)
    Gg = grids(g).G;
    Asum = 0; nHit = 0;
    for c = 1:Gg.cells.num
        a = clipPolyByCell_U(poly, Gg, c, 1e-9);
        Asum = Asum + a;
        nHit = nHit + (a > 0);
    end
    rel = abs(Asum - Aexact)/Aexact;
    fprintf('C. %-14s clipped area %.6f vs %.6f (%d cells), rel %.2e\n', ...
        grids(g).name, Asum, Aexact, nHit, rel);
    okC = okC && rel < 1e-9;
end
if okC
    fprintf('   PASS\n');
else
    failures{end+1} = 'C: clipped areas do not sum to the polygon area'; fprintf('   FAIL\n');
end

%% D. face-coplanar fracture -> boundary type on both neighbours
Gf = computeGeometry(cartGrid([4 4 3], [40 40 30]));
Gf.rock = makeRock(Gf, 100*milli*darcy, 0.2);

fp = struct('points', [20  2  2; 20 38  2; 20 38 28; 20  2 28], ...
            'aperture', 0.01, 'poro', 0.4, 'perm', 1*darcy);
Gf = EDFMgrid_U(Gf, fp, 'Tolerance', 1e-7, 'Verbose', false);

mc  = Gf.FracGrid.Frac1.matrix_connection;
nfc = Gf.FracGrid.Frac1.cells.num;
xc  = Gf.Matrix.cells.centroids(mc.cells, 1);

fprintf('D. face-coplanar fracture: %d fracture cells, %d connections\n', nfc, numel(mc.cells));
fprintf('   types: %d boundary, %d interior; host x-centroids on both sides: %d / %d\n', ...
    sum(mc.type == 1), sum(mc.type == 2), sum(xc < 20), sum(xc > 20));
% every connection is a boundary one, each fracture cell is shared by two
% matrix cells, one on each side of x = 20
shareCount = accumarray(mc.fraccells, 1, [nfc 1]);
if all(mc.type == 1) && all(shareCount == 2) && ...
        sum(xc < 20) == numel(mc.cells)/2 && sum(xc > 20) == numel(mc.cells)/2
    fprintf('   PASS\n');
else
    failures{end+1} = 'D: face-coplanar registration wrong'; fprintf('   FAIL\n');
end

%% E. fracture grid volume = clipped area * aperture
okE = true;
for g = 1:numel(grids)
    Gg = grids(g).G;
    Gg.rock = makeRock(Gg, 100*milli*darcy, 0.2);
    fpe = struct('points', poly, 'aperture', 0.01, 'poro', 0.4, 'perm', 1*darcy);
    Gg = EDFMgrid_U(Gg, fpe, 'Tolerance', 1e-7, 'Verbose', false);
    F = Gg.FracGrid.Frac1;
    rel = abs(sum(F.cells.volumes) - Aexact*0.01) / (Aexact*0.01);
    fprintf('E. %-14s fracture volume rel error %.2e (%d fracture cells)\n', ...
        grids(g).name, rel, F.cells.num);
    okE = okE && rel < 1e-6;
end
if okE
    fprintf('   PASS\n');
else
    failures{end+1} = 'E: fracture grid volume != area*aperture'; fprintf('   FAIL\n');
end

%% Summary
if isempty(failures)
    fprintf('\n=== test01: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test01: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test01 failed');
end

%% ------------------------------ helpers ---------------------------------
function P = combineCorners(x, y, z)
[X, Y, Z] = ndgrid(x, y, z);
P = [X(:), Y(:), Z(:)];
end
