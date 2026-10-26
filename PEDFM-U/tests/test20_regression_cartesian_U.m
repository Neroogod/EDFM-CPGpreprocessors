%% test20_regression_cartesian_U
% The unstructured pEDFM must reproduce the published Cartesian pEDFM when
% the grid IS Cartesian - a Cartesian cell is a convex polyhedron, so the
% general rule has to cover it.
%
% Oracle: pMatFracNNCs3D.m (modules/CPEDFM/UpdateMay27_2023/pEDFM_nnc), the
% reference implementation of Olorode, Wang & Rashid (2020), SPE-201243-PA.
%
% WHAT THIS PINS DOWN
%   The paper picks the nearer face of each of the three axis pairs. This
%   module has no axes to work with, so it pairs each face with the face
%   whose outward normal is most anti-parallel and keeps the nearer one.
%   On a hexahedron the antipodal pairs are exactly the axis pairs, so the
%   two rules must agree - this test is what makes that claim checkable
%   rather than asserted.
%
% Checks:
%   A. Identical set of projection matrix-matrix cell pairs (unordered).
%   B. Identical projected areas (Eq. 24 numerator).
%   C. Identical multiset of pM-F transmissibilities (Eqs. 24-25).
%   D. Identical resulting face multipliers (Eq. 27).
%   E. Every recorded face is incident on both cells of its pMMneighs pair.
%   F. G.nnc arrays stay equal length.
%
% Known, deliberate differences from the reference - the same two
% documented for PEDFM-CPG:
%   1. The reference emits zero-area NNCs where a projection direction is
%      orthogonal to the fracture normal; this module skips them. Check A
%      asserts every reference-only pair has reference area exactly 0.
%   2. The reference collapses colocated projections to one NNC per cell
%      pair, keeping only the first fracture cell; this module keeps one
%      per (projection cell, fracture cell). Check B compares raw areas on
%      pairs reached once and union-corrected totals on pairs reached twice.
%
% Pair ORDER is not compared: the reference encodes the grid side in the
% column order of pMMneighs, which has no meaning without logical
% directions, so this module always stores [host, proj]. Pairs are
% therefore compared unordered.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupPEDFMU.m'));
Globals; % ADFNE globals required by the EDFM gridding path

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
failures = {};

%% Cartesian grid
G = computeGeometry(cartGrid([9, 9, 3], [90, 90, 30]));
G.rock = makeRock(G, 100*milli*darcy, 0.3);

%% One planar fracture, inclined to all three axes
fracplanes = struct;
fracplanes(1).points = [12  0  0;
                        72 90  0;
                        90 90 30;
                        30  0 30];
fracplanes(1).aperture = 10*milli*meter;
fracplanes(1).poro = 0.4;
fracplanes(1).perm = 1*nano*darcy;

[G, fracplanes] = EDFMshalegrid(G, fracplanes, ...
    'Tolerance', tol, 'plotgrid', false, 'fracturelist', 1);
G = fracturematrixShaleNNC3D(G, tol);

isInt = strcmp(G.nnc.type, 'fracmat interior');
hostsI = G.nnc.cells(isInt, 1);
nrm = G.nnc.normal(find(isInt, 1), :);
fprintf('grid: %d matrix cells | %d interior NNCs | %d unique hosts\n', ...
    G.Matrix.cells.num, sum(isInt), numel(unique(hostsI)));
fprintf('fracture normal: [%.4f %.4f %.4f]\n', nrm);
assert(numel(unique(hostsI)) == sum(isInt), ...
    'regression regime needs one fracture cell per host');
assert(all(abs(nrm) > 0.1), 'fracture must be inclined to all three axes');

nnc0 = size(G.nnc.cells, 1);

%% Reference and this module, from the same starting grid
Gref = pMatFracNNCs3D(G, tol);
Gnew = pMatFracNNCsPEDFM_U(G, tol);

s = Gnew.nnc.pStats;
fprintf('PEDFM-U: %d hosts, %d projections, %d antipodal ties\n', ...
    s.nHosts, s.nProjections, s.nTies);

refPairs = sort(Gref.nnc.pMMneighs, 2);
refArea  = Gref.nnc.area(strcmp(Gref.nnc.type, 'fracpmat'));
newPairs = sort(Gnew.nnc.pMMneighs, 2);
newArea  = Gnew.nnc.area(nnc0+1:end);        % raw Eq. 24 areas
fprintf('projections: reference %d pairs | PEDFM-U %d rows\n', ...
    size(refPairs,1), size(newPairs,1));

%% Check A: same set of projection cell pairs
uNew = unique(newPairs, 'rows');
onlyRef = setdiff(refPairs, uNew, 'rows');
onlyNew = setdiff(uNew, refPairs, 'rows');
okA = isempty(onlyNew);
if ~isempty(onlyRef)
    [~, loc] = ismember(onlyRef, refPairs, 'rows');
    okA = okA && all(refArea(loc) < 1e-12);
end
if okA
    fprintf('A. pM-M pair sets identical (%d pairs; %d zero-area ref-only): PASS\n', ...
        size(uNew,1), size(onlyRef,1));
else
    failures{end+1} = sprintf('A: %d only-ref (nonzero area), %d only-new', ...
        size(onlyRef,1), size(onlyNew,1));
    fprintf('A. pM-M pairs: FAIL (%d only-ref, %d only-new)\n', ...
        size(onlyRef,1), size(onlyNew,1));
end

%% Check B: projected areas
[uP, ~, ic] = unique(newPairs, 'rows');
cnt = accumarray(ic, 1);
[tfR, locR] = ismember(uP, refPairs, 'rows');
[~, ~, gFace] = unique(Gnew.nnc.pMMfaces);
faceTotUnion = accumarray(gFace, Gnew.nnc.pMMarea);

maxRelSingle = 0; nSingle = 0; maxRelDup = 0; nDup = 0;
for k = 1:size(uP,1)
    if ~tfR(k), continue; end
    aRef = refArea(locR(k));
    rows = find(ic == k);
    if cnt(k) == 1
        nSingle = nSingle + 1;
        maxRelSingle = max(maxRelSingle, abs(newArea(rows) - aRef)/max(aRef, eps));
    else
        nDup = nDup + 1;
        maxRelDup = max(maxRelDup, ...
            abs(faceTotUnion(gFace(rows(1))) - aRef)/max(aRef, eps));
    end
end
okB = (maxRelSingle < 1e-10) && (maxRelDup < 1e-10);
if okB
    fprintf(['B. projected areas match: %d singleton pairs (raw, max rel. diff %.3e), ', ...
             '%d duplicated pairs (union, max rel. diff %.3e): PASS\n'], ...
        nSingle, maxRelSingle, nDup, maxRelDup);
else
    failures{end+1} = sprintf('B: singleton %.3e, duplicated %.3e', ...
        maxRelSingle, maxRelDup);
    fprintf('B. projected areas: FAIL (singleton %.3e, duplicated %.3e)\n', ...
        maxRelSingle, maxRelDup);
end

%% Check C: pM-F transmissibilities as a multiset
refT = Gref.nnc.T(nnc0+1:end);
newT = Gnew.nnc.T(nnc0+1:end);
if numel(refT) == numel(newT)
    relT = max(abs(sort(refT) - sort(newT))) / max(max(abs(refT)), eps);
    okC = relT < 1e-10;
    if okC
        fprintf('C. pM-F transmissibilities match as a multiset (%d values, max rel. diff %.3e): PASS\n', ...
            numel(refT), relT);
    else
        failures{end+1} = sprintf('C: transmissibility max rel. diff %.3e', relT);
        fprintf('C. pM-F transmissibilities: FAIL (max rel. diff %.3e)\n', relT);
    end
else
    failures{end+1} = sprintf('C: T lengths differ (ref %d, new %d)', ...
        numel(refT), numel(newT));
    fprintf('C. pM-F transmissibilities: FAIL (lengths %d vs %d)\n', ...
        numel(refT), numel(newT));
end

%% Check D: resulting face multipliers
multRef = transmultpEDFM(Gref, tol);
multNew = transmultPEDFM3D_CPG(Gnew, tol);
nMF = G.Matrix.faces.num;
dM = abs(multRef(1:nMF) - multNew(1:nMF));
okD = max(dM) < 1e-10;
if okD
    fprintf('D. face multipliers match (max abs. diff %.3e; %d faces reduced): PASS\n', ...
        max(dM), sum(multNew(1:nMF) < 1 - 1e-12));
else
    failures{end+1} = sprintf('D: multiplier max abs. diff %.3e', max(dM));
    fprintf('D. face multipliers: FAIL (max abs. diff %.3e on %d faces)\n', ...
        max(dM), sum(dM > 1e-10));
end

%% Check E: face incidence
inc = false(numel(Gnew.nnc.pMMfaces), 1);
for q = 1:numel(Gnew.nnc.pMMfaces)
    nbrs = Gnew.faces.neighbors(Gnew.nnc.pMMfaces(q), :);
    inc(q) = any(nbrs == Gnew.nnc.pMMneighs(q,1)) && ...
             any(nbrs == Gnew.nnc.pMMneighs(q,2));
end
okE = all(inc);
if okE
    fprintf('E. every pMMface incident on its cell pair: PASS\n');
else
    failures{end+1} = sprintf('E: %d faces not incident', sum(~inc));
    fprintf('E. face incidence: FAIL (%d of %d)\n', sum(~inc), numel(inc));
end

%% Check F: G.nnc internal consistency
okF = isequal(size(Gnew.nnc.cells,1), numel(Gnew.nnc.T), ...
              numel(Gnew.nnc.area), numel(Gnew.nnc.type));
fprintf('F. G.nnc consistency: cells=%d T=%d area=%d type=%d\n', ...
    size(Gnew.nnc.cells,1), numel(Gnew.nnc.T), numel(Gnew.nnc.area), ...
    numel(Gnew.nnc.type));
if okF
    fprintf('   PASS\n');
else
    failures{end+1} = 'F: nnc arrays inconsistent';
    fprintf('   FAIL\n');
end

%% Summary
save(fullfile(outDir, 'test20_results.mat'), 'refPairs', 'newPairs', ...
     'refArea', 'newArea', 'refT', 'newT', 's');
if isempty(failures)
    fprintf('\n=== test20: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test20: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test20 failed');
end
