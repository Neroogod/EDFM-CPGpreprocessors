%% test00_regression_cartesian
% The corner-point pEDFM must reproduce the published Cartesian pEDFM when
% the grid IS Cartesian.
%
% Oracle: pMatFracNNCs3D.m (modules/CPEDFM/UpdateMay27_2023/pEDFM_nnc), the
% reference implementation of Olorode, Wang & Rashid (2020), SPE-201243-PA.
%
% The fracture is a single plane inclined to all three axes, so that the
% distance test of Fig. 2 (not just the equidistance tie-break) is
% exercised, every direction carries a nonzero projected area, and each
% host matrix cell contains exactly one fracture cell.
%
% Checks:
%   A. Identical set of projection matrix-matrix cell pairs.
%   B. Identical projected areas (Eq. 24 numerator).
%   C. Identical multiset of pM-F transmissibilities (Eqs. 24-25).
%   D. Identical resulting face multipliers (Eq. 27) - the quantity that
%      actually reaches the flow solver.
%   E. Every recorded face is incident on both cells of its pMMneighs pair.
%   F. G.nnc stays internally consistent (see note below).
%
% TWO KNOWN, DELIBERATE DIFFERENCES FROM THE REFERENCE
%
%   1. Zero-area projections. For a fracture whose normal is orthogonal to
%      a grid direction, the projection onto that direction has exactly
%      zero area. The reference still emits an NNC with area 0 and
%      transmissibility 0; this module skips it. The check below asserts
%      that every pair present only in the reference has reference area 0,
%      so the difference can never hide a real missing projection.
%
%   2. Colocated projections. When two host cells project onto opposite
%      sides of the same interface, the reference collapses them to one NNC
%      per cell pair and keeps only the first fracture cell's connection;
%      this module keeps one pM-F NNC per (projection cell, fracture cell)
%      pair. Check B therefore requires exact agreement on pairs reached
%      once, and requires the reference value to be reproduced by one of
%      this module's rows on pairs reached twice.
%
%   Difference 2 exposes a genuine inconsistency in the reference:
%   pMatFracNNCs3D appends the PRE-collapse transmissibility vector to
%   G.nnc.T (line 519) but the POST-collapse cells/area/type arrays (lines
%   570-572), so G.nnc.T ends up longer than G.nnc.cells by the number of
%   collapsed duplicates. Check F asserts this module does not do that.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupPEDFMCPG.m'));
Globals; % ADFNE globals required by the EDFM gridding path

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
failures = {};

%% Cartesian grid
celldim = [9, 9, 3];
physdim = [90, 90, 30];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);
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
nInterior = sum(isInt);
hostsI = G.nnc.cells(isInt, 1);
nrm = G.nnc.normal(find(isInt, 1), :);
fprintf('grid: %d matrix cells | %d interior NNCs | %d unique hosts\n', ...
    G.Matrix.cells.num, nInterior, numel(unique(hostsI)));
fprintf('fracture normal: [%.4f %.4f %.4f] (inclined to all three axes)\n', nrm);
assert(nInterior > 10, 'test setup produced too few interior NNCs');
assert(numel(unique(hostsI)) == nInterior, ...
    'test setup has colocated fracture cells; regression needs one per host');
assert(all(abs(nrm) > 0.1), ...
    'fracture must be inclined to all three axes for this regression');

nnc0 = size(G.nnc.cells, 1);   % NNC count before any projection

%% Reference (Cartesian) and this module, from the same starting grid
Gref = pMatFracNNCs3D(G, tol);
Gnew = pMatFracNNCsPEDFM3D_CPG(G, tol);

s = Gnew.nnc.pStats;
fprintf('cases: C0=%d C1=%d C2=%d C3=%d | fallbacks=%d dropped=%d\n', ...
    s.nCase0, s.nCase1, s.nCase2, s.nCase3, s.nFallback, s.nDropped);
assert(s.nCase0 > 0, 'distance test never fired; geometry is degenerate');
assert(s.nFallback == 0, ...
    'node-containment fallback fired on a Cartesian grid; it should not');

refPairs = Gref.nnc.pMMneighs;
refArea  = Gref.nnc.area(strcmp(Gref.nnc.type, 'fracpmat'));
newPairs = Gnew.nnc.pMMneighs;
% Raw per-record projected area (Eq. 24), the same quantity the reference
% stores. G.nnc.pMMarea is this after the Eq. 29 union correction and is
% checked separately in B2.
newArea  = Gnew.nnc.area(nnc0+1:end);
fprintf('projections: reference %d pairs | PEDFM-CPG %d rows\n', ...
    size(refPairs,1), size(newPairs,1));

%% Check A: same set of projection cell pairs (modulo zero-area entries)
uNew = unique(newPairs, 'rows');
onlyRef = setdiff(refPairs, uNew, 'rows');
onlyNew = setdiff(uNew, refPairs, 'rows');

okA = isempty(onlyNew);
if ~isempty(onlyRef)
    % Permitted only if the reference itself assigned them zero area.
    [~, loc] = ismember(onlyRef, refPairs, 'rows');
    okA = okA && all(refArea(loc) < 1e-12);
end
if okA
    fprintf('A. pM-M pair sets identical (%d pairs; %d zero-area ref-only): PASS\n', ...
        size(uNew,1), size(onlyRef,1));
else
    failures{end+1} = sprintf('A: %d only-ref (nonzero area), %d only-new', ...
        size(onlyRef,1), size(onlyNew,1));
    fprintf('A. pM-M pairs: FAIL\n');
end

%% Check B: projected areas
[uP, ~, ic] = unique(newPairs, 'rows');
cnt = accumarray(ic, 1);
[tfR, locR] = ismember(uP, refPairs, 'rows');

% For a pair reached once, the raw Eq. 24 area must match the reference
% exactly. For a pair reached twice the reference applies its own Eq. 29
% union (pMatFracNNCs3D lines 545-561) before collapsing, so the
% like-for-like comparison there is against this module's union-corrected
% total on the carrying face.
faceTotUnion = accumarray(findgroups(Gnew.nnc.pMMfaces), Gnew.nnc.pMMarea);
[~, ~, gFace] = unique(Gnew.nnc.pMMfaces);

maxRelSingle = 0; nSingle = 0; maxRelDup = 0; nDup = 0;
for k = 1:size(uP,1)
    if ~tfR(k), continue; end
    aRef = refArea(locR(k));
    rows = find(ic == k);
    if cnt(k) == 1
        nSingle = nSingle + 1;
        maxRelSingle = max(maxRelSingle, ...
            abs(newArea(rows) - aRef)/max(aRef, eps));
    else
        nDup = nDup + 1;
        aUni = faceTotUnion(gFace(rows(1)));
        maxRelDup = max(maxRelDup, abs(aUni - aRef)/max(aRef, eps));
    end
end
okB = (maxRelSingle < 1e-10) && (maxRelDup < 1e-10);
if okB
    fprintf(['B. projected areas match: %d singleton pairs (raw, max rel. diff %.3e), ', ...
             '%d duplicated pairs (union, max rel. diff %.3e): PASS\n'], ...
        nSingle, maxRelSingle, nDup, maxRelDup);
else
    failures{end+1} = sprintf('B: singleton rel. diff %.3e, duplicated rel. diff %.3e', ...
        maxRelSingle, maxRelDup);
    fprintf('B. projected areas: FAIL (singleton %.3e, duplicated %.3e)\n', ...
        maxRelSingle, maxRelDup);
end

%% Check B2: the Eq. 29 union correction behaves
% On every face the union-corrected total must not exceed the face area,
% must not exceed the raw sum (a union is never larger than a sum), and
% must equal the raw area wherever a face carries a single projection.
pA = Gnew.nnc.pMMarea;
[uf, ~, fic] = unique(Gnew.nnc.pMMfaces);
sumRaw = accumarray(fic, newArea);
sumUni = accumarray(fic, pA);
nPerFace = accumarray(fic, 1);
isLone = nPerFace(fic) == 1;          % records sitting alone on their face
nSingleFace = sum(nPerFace == 1);
okB2 = all(sumUni <= G.faces.areas(uf) + 1e-9) && ...
       all(sumUni <= sumRaw + 1e-9) && ...
       all(abs(pA(isLone) - newArea(isLone)) < 1e-12);
if okB2
    fprintf(['B2. union correction sane (%d faces, %d single-projection; ', ...
             'union <= sum and <= face area): PASS\n'], numel(uf), nSingleFace);
else
    failures{end+1} = 'B2: union-corrected areas violate union <= sum <= face area';
    fprintf('B2. union correction: FAIL\n');
end

%% Check C: pM-F transmissibilities, compared as a multiset
% The reference's appended T block is PRE-collapse and therefore has the
% same length as this module's; it cannot be aligned to refPairs (see the
% header note), so compare the sorted values.
refT = Gref.nnc.T(nnc0+1:end);
newT = Gnew.nnc.T(nnc0+1:end);
if numel(refT) == numel(newT)
    dT = abs(sort(refT) - sort(newT));
    relT = max(dT) / max(max(abs(refT)), eps);
    okC = relT < 1e-10;
    if okC
        fprintf('C. pM-F transmissibilities match as a multiset (%d values, max rel. diff %.3e): PASS\n', ...
            numel(refT), relT);
    else
        failures{end+1} = sprintf('C: transmissibility max rel. diff %.3e', relT);
        fprintf('C. pM-F transmissibilities: FAIL (max rel. diff %.3e)\n', relT);
    end
else
    failures{end+1} = sprintf('C: T block lengths differ (ref %d, new %d)', ...
        numel(refT), numel(newT));
    fprintf('C. pM-F transmissibilities: FAIL (lengths %d vs %d)\n', ...
        numel(refT), numel(newT));
end

%% Check D: resulting face multipliers (what the flow solver sees)
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

%% Check E: recorded faces are incident on their cell pair
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
    failures{end+1} = sprintf('E: %d faces not incident on their pair', sum(~inc));
    fprintf('E. face incidence: FAIL (%d of %d)\n', sum(~inc), numel(inc));
end

%% Check F: G.nnc internal consistency
nC = size(Gnew.nnc.cells,1); nT = numel(Gnew.nnc.T);
nA = numel(Gnew.nnc.area);   nY = numel(Gnew.nnc.type);
okF = isequal(nC, nT, nA, nY);
fprintf('F. G.nnc consistency: this module cells=%d T=%d area=%d type=%d | reference cells=%d T=%d\n', ...
    nC, nT, nA, nY, size(Gref.nnc.cells,1), numel(Gref.nnc.T));
if okF
    fprintf('   PASS (reference is inconsistent by %d entries; see header)\n', ...
        numel(Gref.nnc.T) - size(Gref.nnc.cells,1));
else
    failures{end+1} = sprintf('F: nnc arrays inconsistent (cells=%d T=%d area=%d type=%d)', ...
        nC, nT, nA, nY);
    fprintf('   FAIL\n');
end

%% Summary
save(fullfile(outDir, 'test00_results.mat'), 'refPairs', 'newPairs', ...
     'refArea', 'newArea', 'refT', 'newT', 's');
if isempty(failures)
    fprintf('\n=== test00: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test00: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test00 failed');
end
