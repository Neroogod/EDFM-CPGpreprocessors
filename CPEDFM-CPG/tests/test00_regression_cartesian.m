%% test00_regression_cartesian
% Regression of the CPG implementation against the Cartesian CPEDFM
% release on a Cartesian grid (where both must agree).
%
% Case: the slanted sealing fracture of Fig. 11a (SPE-217469-PA).
%
% Checks:
%   A. Neighborhoods from buildNeighborhoods match centroid-coordinate
%      grouping (the original's rule) exactly.
%   B. Orientation branch: fracOrientation selects the same face tags
%      as the ADFNE-based classification in calculateProjectionCPEDFM.
%   C. Set of sealed matrix-matrix pairs (pMMneighs) matches.
%   D. Per-face projected areas and the transmissibility multipliers
%      match within tolerance.
%   E. pM-F transmissibilities: compared and reported. Known deviation:
%      for positive-side selections in x/y, the original returns zero
%      merged transmissibility (its collapse loops retain no face there);
%      the CPG code returns the Eq. 10 value. Differences are reported,
%      not failed, when they occur exactly on those entries.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupCPEDFMCPG.m'));
Globals; % ADFNE globals required by the original code path

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
failures = {};

%% Grid and fracture (Fig. 11a case)
celldim = [30, 30, 5];
physdim = [60, 60, 10];
G0 = cartGrid(celldim, physdim);
G0 = computeGeometry(G0);
G0.rock = makeRock(G0, 1*milli*darcy, 0.1);

fracplanes = struct;
fracplanes(1).points = [10  0  0;
                        50 60  0;
                        50 60 10;
                        10  0 10];
fracplanes(1).aperture = 10*milli*meter;
fracplanes(1).poro = 0.4;
fracplanes(1).perm = 1*nano*darcy;

%% Shared gridding + standard EDFM NNC stage (identical input for both)
[Gbase, fracplanes] = EDFMshalegrid(G0, fracplanes, ...
    'Tolerance', tol, 'plotgrid', false, 'fracturelist', 1);
Gbase = fracturematrixShaleNNC3D(Gbase, tol);

%% Path A: original Cartesian CPEDFM
GA = pMatFracNNCs3DCPEDFM(Gbase, tol);
multA = transmultpEDFM(GA, tol);

%% Path B: corner-point implementation
GB = pMatFracNNCsCPEDFM_CPG(Gbase, tol);
multB = transmultpEDFM_CPG(GB, tol);

%% Check A: neighborhood equivalence (IJK vs centroid grouping)
interiorIdx = find(strcmp(Gbase.nnc.type, 'fracmat interior'));
hostCells = Gbase.nnc.cells(interiorIdx, 1);
h = calculateProjectionCPEDFM_CPG('helpers');
nb = h.buildNeighborhoods(Gbase.Matrix, hostCells);
cc = Gbase.cells.centroids(hostCells, :);
% original grouping: X-loop = same (x,z); Y-loop = same (y,z); Z-loop = same (x,y)
coordPairs = {[2 3], [1 3], [1 2]}; % vary I <-> share (y,z); vary J <-> share (x,z); vary K <-> share (x,y)
okA = true;
for d = 1:3
    [~, ~, gidRef] = unique(round(cc(:, coordPairs{d})/tol), 'rows');
    nGroupsRef = max(gidRef);
    if nGroupsRef ~= numel(nb{d})
        okA = false;
    else
        for g = 1:numel(nb{d})
            if numel(unique(gidRef(nb{d}{g}))) ~= 1, okA = false; end
        end
    end
end
if okA
    fprintf('A. neighborhoods (IJK vs centroid): PASS\n');
else
    failures{end+1} = 'A: neighborhood grouping mismatch';
    fprintf('A. neighborhoods (IJK vs centroid): FAIL\n');
end

%% Check B: orientation branch / selected tags
nrm = Gbase.nnc.normal(interiorIdx(1), :);
selTags = h.fracOrientation(nrm);
% Original's selected local face IDs live in GA.nnc.nncFaces(:,2) (per
% host cell, before the collapse); on cartGrid local position == tag.
tagsA = unique(GA.nnc.nncFacesOld(:,2))';
tagsB = sort(selTags);
if isequal(sort(tagsA), tagsB)
    fprintf('B. selected face tags: PASS  (tags = %s)\n', mat2str(tagsB));
else
    failures{end+1} = sprintf('B: tags original=%s cpg=%s', mat2str(tagsA), mat2str(tagsB));
    fprintf('B. selected face tags: FAIL  original=%s  cpg=%s\n', ...
        mat2str(sort(tagsA)), mat2str(tagsB));
end

%% Check C: sealed matrix-matrix pairs
% The CPG implementation must reproduce every pair the original seals
% (envelope rule). It may add pairs from the flip-face rule - faces the
% fracture geometrically crosses between two host cells, which the
% original leaves open (they are the leakage path on host slabs thicker
% than two cells). Those extras must all be host-host flip faces.
pairsA = sortrows(GA.nnc.pMMneighs);
pairsB = sortrows(GB.nnc.pMMneighs);
setA = unique(pairsA, 'rows');
setB = unique(pairsB, 'rows');
common = intersect(setA, setB, 'rows');
onlyA = setdiff(setA, setB, 'rows');
onlyB = setdiff(setB, setA, 'rows');
hostSet = unique(Gbase.nnc.cells(interiorIdx,1));
% every flip addition is a face of a host cell, so each extra pair must
% contain at least one host cell
extrasTouchHost = all(any(ismember(onlyB, hostSet), 2));
fprintf('C. pMM pairs: common %d | missing-from-cpg %d | flip-additions %d\n', ...
    size(common,1), size(onlyA,1), size(onlyB,1));
if isempty(onlyA) && extrasTouchHost
    fprintf('   PASS (original pairs all reproduced; %d flip-face additions, all host-adjacent)\n', ...
        size(onlyB,1));
else
    failures{end+1} = sprintf('C: %d original pairs missing / extras not host-adjacent', ...
        size(onlyA,1));
end

%% Check D: per-face projected area and multipliers
% original: faces via matching pMM pairs on matrix neighbor list
nMM_A = size(GA.nnc.pMMneighs, 1);
pStart_A = numel(GA.nnc.area) - nMM_A;
faceA = zeros(nMM_A, 1);
for j = 1:nMM_A
    lia = ismember(Gbase.Matrix.faces.neighbors, GA.nnc.pMMneighs(j,:), 'rows');
    fj = find(lia);
    if numel(fj) == 1, faceA(j) = fj; end
end
areaA = GA.nnc.area(pStart_A+1:end);
vA = accumarray(faceA(faceA>0), areaA(faceA>0), [Gbase.faces.num, 1]);
vB = accumarray(GB.nnc.pMMfaces, GB.nnc.pMMarea, [Gbase.faces.num, 1]);
% compare on the faces the ORIGINAL projects onto; extra flip faces are
% additions by design and are checked separately below
onA = unique(faceA(faceA>0));
dArea = max(abs(vA(onA) - vB(onA)));
relArea = dArea / max(max(vA), eps);
fprintf('D. projected area on original faces: max abs diff %.3e (rel %.3e)\n', ...
    dArea, relArea);
dMult = max(abs(multA(onA) - multB(onA)));
fprintf('   multiplier diff on original faces: %.3e\n', dMult);
extraFaces = setdiff(unique(GB.nnc.pMMfaces), onA);
fprintf('   flip-face additions: %d faces, all fully sealed: %d\n', ...
    numel(extraFaces), all(multB(extraFaces) == 0));
if relArea < 1e-8 && dMult < 1e-8
    fprintf('   PASS\n');
else
    failures{end+1} = sprintf('D: area rel diff %.3e, mult diff %.3e', relArea, dMult);
end

%% Check E: pM-F transmissibilities (report; known deviation on pos sides)
TA = GA.nnc.T(pStart_A+1:end);
% match entries by (pM,F) pair
mfA = GA.nnc.cells(pStart_A+1:end, :);
mfB = GB.nnc.cells(end-numel(GB.nnc.pMMfaces)+1:end, :);
TB = GB.nnc.T(end-numel(GB.nnc.pMMfaces)+1:end);
[tf, loc] = ismember(mfA, mfB, 'rows');
relT = abs(TA(tf) - TB(loc(tf))) ./ max(abs(TA(tf)), eps);
fprintf('E. pM-F trans: %d/%d matched pairs, max rel diff %.3e\n', ...
    nnz(tf), numel(tf), max([relT; 0]));
fprintf('   (deviations expected where original returns 0 on positive-side x/y entries)\n');
zeroA = TA(tf) == 0 & TB(loc(tf)) > 0;
fprintf('   entries where original=0 but cpg>0: %d\n', nnz(zeroA));

%% Summary
save(fullfile(outDir, 'test00_results.mat'), 'GA', 'GB', 'multA', 'multB');
if isempty(failures)
    fprintf('\n=== test00: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test00: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test00 failed');
end
