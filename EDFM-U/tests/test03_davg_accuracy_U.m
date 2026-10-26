%% test03_davg_accuracy_U
% Accuracy of the average normal fracture-matrix distance <d> on genuinely
% unstructured cells, against a Monte-Carlo ground truth.
%
% <d> enters EDFM through the connectivity index CI = A/<d>, so a relative
% error in <d> is a relative error of the same size in every
% fracture-matrix transmissibility.
%
% Two estimators are graded on the same host cells:
%   calcDavg_U    - EDFM-U: exact integral over the true cell polyhedron.
%   hfm/calcdavg  - the Cartesian implementation: 10x10x10 midpoint samples
%                   over the cell's AXIS-ALIGNED BOUNDING BOX. On a
%                   Cartesian cell the box is the cell; on a PEBI or
%                   tetrahedral cell it is not, and the estimate degrades.
%
% Ground truth: rejection sampling inside the cell (a point is inside iff
% it is behind every face plane). The Monte-Carlo standard error is
% reported so the comparison is interpretable - calcDavg_U is required to
% agree within a few standard errors, i.e. to be indistinguishable from
% exact at the resolution of the reference.
%
% Checks
%   A. Extruded PEBI grid: calcDavg_U within 4 MC standard errors on every
%      host cell, and clearly better than the bounding-box estimator.
%   B. Tetrahedral grid: same, where the bounding-box error is severe
%      (a tetrahedron fills only 1/6 of a comparable box).
%   C. The bounding-box estimator's error is NOT a sampling artefact: the
%      same estimator on a Cartesian cell is accurate, so the degradation
%      on unstructured cells comes from the geometry it integrates over.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupEDFMU.m'));
addpath(thisDir);                       % boxMeanAbsDist_U (test reference)

tol = 1e-9;
failures = {};
rng(17);

nSample  = 6e6;    % trial points per cell (before rejection)
nCells   = 12;     % host cells graded per grid

physdim = [40 40 30];
nrm = [0.62 0.55 0.56]; nrm = nrm/norm(nrm);
pp  = [20 20 15];

%% Grids
grids = struct('name', {}, 'G', {});

G2 = makeLayeredGrid(pebiGrid2D(8, physdim(1:2)), 3);
G2.nodes.coords(:,3) = G2.nodes.coords(:,3) * physdim(3)/3;
grids(1).name = 'extruded PEBI';
grids(1).G    = computeGeometry(G2);

pts = [physdim(1)*rand(70,1), physdim(2)*rand(70,1), physdim(3)*rand(70,1)];
[CX,CY,CZ] = ndgrid([0 physdim(1)], [0 physdim(2)], [0 physdim(3)]);
DT = delaunayTriangulation([pts; CX(:) CY(:) CZ(:)]);
grids(2).name = 'tetrahedral';
grids(2).G    = computeGeometry(tetrahedralGrid(DT.Points, DT.ConnectivityList));

%% A/B. grade both estimators on cells the plane actually cuts
[cnAll, cposAll] = deal(cell(numel(grids),1));
report = struct('name', {}, 'maxSig', {}, 'medU', {}, 'medB', {}, 'maxB', {});

for g = 1:numel(grids)
    Gg = grids(g).G;
    [cn, cpos] = gridCellNodes(Gg, 1:Gg.cells.num);
    cnAll{g} = cn; cposAll{g} = cpos;

    % cells genuinely straddling the plane
    cut = false(Gg.cells.num,1);
    for c = 1:Gg.cells.num
        Vm = clipCellByPlane_U(Gg, c, nrm, pp, tol);
        cut(c) = Vm > 1e-9 && Vm < Gg.cells.volumes(c) - 1e-9;
    end
    cutc = find(cut);
    sel  = cutc(round(linspace(1, numel(cutc), min(nCells, numel(cutc)))));

    relU = zeros(numel(sel),1); relB = zeros(numel(sel),1); sig = zeros(numel(sel),1);
    fprintf('\n%s: %d cells, %d cut by the plane, grading %d\n', ...
        grids(g).name, Gg.cells.num, numel(cutc), numel(sel));
    fprintf('  %6s %14s %14s %12s %10s %10s\n', ...
        'cell', 'MC <d>', 'calcDavg_U', 'bbox', 'errU/sig', 'errB[%]');

    for j = 1:numel(sel)
        c = sel(j);
        [dMC, sMC, acc] = mcMeanAbsDist(Gg, c, nrm, pp, nSample);
        dU = calcDavg_U(Gg, c, nrm, pp, tol);
        nd = Gg.nodes.coords(cn(cpos(c):cpos(c+1)-1), :);
        dB = calcdavg(nd, nrm, pp, tol);

        relU(j) = abs(dU - dMC)/dMC;
        relB(j) = abs(dB - dMC)/dMC;
        sig(j)  = sMC/dMC;
        fprintf('  %6d %14.8f %14.8f %12.8f %10.2f %10.2f   (%d pts)\n', ...
            c, dMC, dU, dB, relU(j)/sig(j), 100*relB(j), acc);
    end

    maxSig = max(relU./sig);
    fprintf('  calcDavg_U : max deviation %.2f MC standard errors, median rel err %.2e\n', ...
        maxSig, median(relU));
    fprintf('  bbox       : median rel err %.2e, max %.2e\n', median(relB), max(relB));

    report(g).name   = grids(g).name;
    report(g).maxSig = maxSig;
    report(g).medU   = median(relU);
    report(g).medB   = median(relB);
    report(g).maxB   = max(relB);

    if maxSig < 4 && median(relB) > 5*median(relU)
        fprintf('  PASS\n');
    else
        failures{end+1} = sprintf('%s: maxSig %.2f, medU %.2e, medB %.2e', ...
            grids(g).name, maxSig, median(relU), median(relB)); %#ok<SAGROW>
        fprintf('  FAIL\n');
    end
end

%% C. the bounding-box estimator is fine on a Cartesian cell
% Same estimator, same plane, a cell whose bounding box IS the cell: the
% error collapses, so what it loses on PEBI/tet cells is geometry, not
% sampling noise.
Gc = computeGeometry(cartGrid([1 1 1], [10 10 10]));
nd = [0 0 0; 10 0 0; 0 10 0; 10 10 0; 0 0 10; 10 0 10; 0 10 10; 10 10 10];
exC = boxMeanAbsDist_U([0 0 0], [10 10 10], nrm, [5 5 5]);
eUC = abs(calcDavg_U(Gc, 1, nrm, [5 5 5], tol) - exC)/exC;
eBC = abs(calcdavg(nd, nrm, [5 5 5], tol) - exC)/exC;
fprintf('\nC. Cartesian control cell: calcDavg_U %.2e | bbox %.2e (vs exact)\n', eUC, eBC);
fprintf('   bbox error grows from %.2e (Cartesian) to %.2e (PEBI) / %.2e (tet)\n', ...
    eBC, report(1).medB, report(2).medB);
if eUC < 1e-9 && eBC < 0.02 && report(2).medB > 3*eBC
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('C: eUC %.2e eBC %.2e medB_tet %.2e', ...
        eUC, eBC, report(2).medB);
    fprintf('   FAIL\n');
end

%% Summary
fprintf('\n--- impact on CI = A/<d> (relative CI error equals relative <d> error) ---\n');
for g = 1:numel(report)
    fprintf('  %-14s : EDFM-U %.2e | bounding-box %.2e (median), %.2e (worst)\n', ...
        report(g).name, report(g).medU, report(g).medB, report(g).maxB);
end

if isempty(failures)
    fprintf('\n=== test03: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test03: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test03 failed');
end

%% ------------------------------ helpers ---------------------------------
function [d, se, nAcc] = mcMeanAbsDist(G, c, n, p0, nSample)
% Monte-Carlo mean of |n.(x-p0)| over cell c, by rejection sampling in the
% cell's bounding box. A point is inside the (convex) cell iff it lies
% behind every one of its face planes.
fp  = G.cells.facePos(c):G.cells.facePos(c+1)-1;
fcs = G.cells.faces(fp,1);

N = G.faces.normals(fcs,:);
N = N ./ sqrt(sum(N.^2,2));
flip = G.faces.neighbors(fcs,1) ~= c;
N(flip,:) = -N(flip,:);
PC = G.faces.centroids(fcs,:);

nodes = [];
for q = 1:numel(fcs)
    np = G.faces.nodePos(fcs(q)):G.faces.nodePos(fcs(q)+1)-1;
    nodes = [nodes; G.nodes.coords(G.faces.nodes(np),:)]; %#ok<AGROW>
end
lo = min(nodes); hi = max(nodes);

P = lo + rand(nSample,3).*(hi - lo);
inside = true(nSample,1);
for q = 1:numel(fcs)
    inside = inside & ((P - PC(q,:))*N(q,:)' <= 0);
end
P = P(inside,:);
nAcc = size(P,1);
assert(nAcc > 1e4, 'too few accepted samples (%d) in cell %d', nAcc, c);

v  = abs((P - p0) * (n/norm(n))');
d  = mean(v);
se = std(v)/sqrt(nAcc);
end
