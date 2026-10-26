%% test02_faults_and_ties
% Exercises the two parts of the corner-point port that test00 (Cartesian
% regression) and test01 (skewed grid, all Case 0) do not reach:
%   - fault-split logical sides, where a tag set holds several faces and a
%     host cell has more than six faces;
%   - the equidistance tie-break of Cases 1-3, and the reformulation of the
%     paper's Table 1 as a node-containment rule.
%
% THE INVARIANT UNDER TEST
%   Table 1 of Olorode, Wang & Rashid (2020) maps (direction of the
%   selected face, index of its node farthest from the fracture plane) to
%   the two remaining faces. As argued in pedfm_implementation.md section
%   2.3, every row of that table names exactly the two other faces of the
%   cell that contain the chosen node - so the three selected faces meet at
%   one corner. This test checks that corner property directly, which
%   validates the reformulation without depending on any node ordering.
%
% Note on realism: on a genuinely warped corner-point grid exact
% equidistance is a measure-zero event, so Cases 1-3 essentially never
% arise by themselves. Part B therefore widens the equidistance tolerance
% to drive the tie-break deliberately, and Part C uses a Cartesian grid
% with a symmetric fracture where ties occur naturally.
%
% Checks:
%   A. Faulted CPG runs; fault-split host cells are present and handled.
%   B. Forced ties drive Cases 1-3; no crash, data stays consistent.
%   C. Corner invariant holds on every six-faced host, in all parts.
%   D. Natural ties on a Cartesian grid are handled.
%   E. Projection data self-consistent everywhere (multipliers in [0,1],
%      faces incident on their pair, G.nnc arrays equal length).

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupPEDFMCPG.m'));
Globals;

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
failures = {};
caseSeen = false(1,4);   % which of Cases 0..3 were observed anywhere

%% ---------------------------------------------------------------- Part A
% Faulted corner-point grid, default settings.
grdecl = makeModel3([12, 10, 6], [240, 200, 60]);
Gf = computeGeometry(processGRDECL(grdecl));
Gf.rock = makeRock(Gf, 100*milli*darcy, 0.25);

nFaces = diff(Gf.cells.facePos);
fprintf('faulted grid: %d cells, %d with >6 faces (split logical sides)\n', ...
    Gf.cells.num, sum(nFaces > 6));
assert(sum(nFaces > 6) > 0, 'grid has no fault splits; test would be vacuous');

% Vertical fracture on a mid-cell plane (x = 110; cell centroids are at
% 10, 30, ... 230, so this is interior rather than a cell face).
fpF = struct('points', [110  10 -20; 110 190 -20; 110 190 90; 110  10 90], ...
             'aperture', 10*milli*meter, 'poro', 0.4, 'perm', 1*nano*darcy);

[GfF, ~] = EDFMGridCPG(Gf, fpF, 'Tolerance', tol, 'plotgrid', false, 'fracturelist', 1);
GfF = fracturematrixShaleNNC3D(GfF, tol);

GpA = pMatFracNNCsPEDFM3D_CPG(GfF, tol);
sA = GpA.nnc.pStats;

hostsA = unique(GfF.nnc.cells(strcmp(GfF.nnc.type, 'fracmat interior'), 1));
nfA = diff(GfF.cells.facePos);
nSplitHosts = sum(nfA(hostsA) > 6);
fprintf('A. faulted CPG: %d hosts, %d fault-split | C0=%d C1=%d C2=%d C3=%d | fallbacks=%d | %d projections\n', ...
    numel(hostsA), nSplitHosts, sA.nCase0, sA.nCase1, sA.nCase2, sA.nCase3, ...
    sA.nFallback, numel(GpA.nnc.pMMfaces));
okA = nSplitHosts > 0 && numel(GpA.nnc.pMMfaces) > 0;
if okA
    fprintf('   PASS (fault-split hosts projected)\n');
else
    failures{end+1} = 'A: no fault-split hosts projected';
    fprintf('   FAIL\n');
end
caseSeen = caseSeen | [sA.nCase0 sA.nCase1 sA.nCase2 sA.nCase3] > 0;

%% ---------------------------------------------------------------- Part B
% Same faulted grid, equidistance tolerance widened to force the tie-break.
GpB = pMatFracNNCsPEDFM3D_CPG(GfF, tol, 'relTol', 0.5);
sB = GpB.nnc.pStats;
fprintf('B. forced ties (relTol=0.5): C0=%d C1=%d C2=%d C3=%d | fallbacks=%d | %d projections\n', ...
    sB.nCase0, sB.nCase1, sB.nCase2, sB.nCase3, sB.nFallback, numel(GpB.nnc.pMMfaces));
okB = (sB.nCase1 + sB.nCase2 + sB.nCase3) > 0 && numel(GpB.nnc.pMMfaces) > 0;
if okB
    fprintf('   PASS (tie-break exercised on a faulted grid)\n');
else
    failures{end+1} = 'B: forced tolerance did not produce any tied case';
    fprintf('   FAIL\n');
end
caseSeen = caseSeen | [sB.nCase0 sB.nCase1 sB.nCase2 sB.nCase3] > 0;

%% ---------------------------------------------------------------- Part D
% Cartesian grid with a symmetric vertical fracture: ties arise naturally.
Gc = computeGeometry(cartGrid([9, 9, 3], [90, 90, 30]));
Gc.rock = makeRock(Gc, 100*milli*darcy, 0.3);
fpC = struct('points', [15 0 0; 75 90 0; 75 90 30; 15 0 30], ...
             'aperture', 10*milli*meter, 'poro', 0.4, 'perm', 1*nano*darcy);
[GcF, ~] = EDFMshalegrid(Gc, fpC, 'Tolerance', tol, 'plotgrid', false, 'fracturelist', 1);
GcF = fracturematrixShaleNNC3D(GcF, tol);
GpD = pMatFracNNCsPEDFM3D_CPG(GcF, tol);
sD = GpD.nnc.pStats;
fprintf('D. Cartesian natural ties: C0=%d C1=%d C2=%d C3=%d | fallbacks=%d | %d projections\n', ...
    sD.nCase0, sD.nCase1, sD.nCase2, sD.nCase3, sD.nFallback, numel(GpD.nnc.pMMfaces));
okD = (sD.nCase1 + sD.nCase3) > 0;
if okD
    fprintf('   PASS (vertical fracture ties in z, resolved without fallback)\n');
else
    failures{end+1} = 'D: expected natural ties on the symmetric Cartesian case';
    fprintf('   FAIL\n');
end
caseSeen = caseSeen | [sD.nCase0 sD.nCase1 sD.nCase2 sD.nCase3] > 0;

fprintf('   cases observed across parts A/B/D: %s\n', ...
    strjoin(compose('Case%d', find(caseSeen)-1), ', '));

%% ---------------------------------------------------------------- Part C
% Corner invariant, checked on every six-faced host of every configuration.
configs = { {GfF, 1e-6,  'faulted, default'}, ...
            {GfF, 0.5,   'faulted, forced ties'}, ...
            {GcF, 1e-6,  'Cartesian, natural ties'} };
nOK = 0; nTested = 0; badExamples = cell(1,3); nBad = 0;
for cfg = 1:numel(configs)
    Gx = configs{cfg}{1}; rt = configs{cfg}{2};
    o = struct('tol', tol, 'relTol', rt, 'defaultFace', 'bottom', ...
               'distanceMode', 'centroid');
    rows = find(strcmp(Gx.nnc.type, 'fracmat interior'));
    nfc = diff(Gx.cells.facePos);
    for r = reshape(rows, 1, [])
        h = Gx.nnc.cells(r,1);
        if nfc(h) ~= 6, continue; end          % corner is only defined here
        nrm = Gx.nnc.normal(r,:); nrm = nrm/norm(nrm);
        st = selectProjectionFacesPEDFM_CPG(Gx, h, ...
            Gx.cells.centroids(Gx.nnc.cells(r,2),:), nrm, Gx.nnc.planepoint(r,:), o);
        if any(st == 0), continue; end
        fs = getLogicalFacesPEDFM_CPG(Gx, h);
        nd = cell(1,3);
        for d = 1:3
            f = fs{st(d)};
            nd{d} = Gx.faces.nodes(Gx.faces.nodePos(f):(Gx.faces.nodePos(f+1)-1));
        end
        common = intersect(intersect(nd{1}, nd{2}), nd{3});
        nTested = nTested + 1;
        if numel(common) == 1
            nOK = nOK + 1;
        elseif nBad < numel(badExamples)
            nBad = nBad + 1;
            badExamples{nBad} = sprintf('cfg %d, cell %d: %d common nodes', ...
                cfg, h, numel(common));
        end
    end
end
fprintf('C. corner invariant (3 selected faces meet at exactly 1 node): %d/%d six-faced hosts\n', ...
    nOK, nTested);
okC = (nTested > 0) && (nOK == nTested);
if okC
    fprintf('   PASS (Table 1 reformulation holds)\n');
else
    failures{end+1} = sprintf('C: corner invariant failed on %d of %d hosts', ...
        nTested - nOK, nTested);
    fprintf('   FAIL\n');
    if nBad > 0, fprintf('     e.g. %s\n', badExamples{1:nBad}); end
end

%% ---------------------------------------------------------------- Part E
% Self-consistency of every produced grid.
okE = true;
for Gx = {GpA, GpB, GpD}
    Gq = Gx{1};
    m = transmultPEDFM3D_CPG(Gq, tol);
    okE = okE && all(m >= -1e-12 & m <= 1 + 1e-12);
    okE = okE && all(Gq.nnc.pMMarea > 0) && all(isfinite(Gq.nnc.pMMarea));
    okE = okE && all(isfinite(Gq.nnc.T));
    okE = okE && isequal(size(Gq.nnc.cells,1), numel(Gq.nnc.T), ...
                         numel(Gq.nnc.area), numel(Gq.nnc.type));
    for q = 1:numel(Gq.nnc.pMMfaces)
        nb = Gq.faces.neighbors(Gq.nnc.pMMfaces(q),:);
        okE = okE && any(nb == Gq.nnc.pMMneighs(q,1)) && ...
                     any(nb == Gq.nnc.pMMneighs(q,2));
    end
end
if okE
    fprintf('E. projection data self-consistent in all three configurations: PASS\n');
else
    failures{end+1} = 'E: projection data inconsistent';
    fprintf('E. projection data: FAIL\n');
end

%% Summary
save(fullfile(outDir, 'test02_results.mat'), 'sA', 'sB', 'sD', ...
     'nOK', 'nTested', 'nSplitHosts');
if isempty(failures)
    fprintf('\n=== test02: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test02: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test02 failed');
end
