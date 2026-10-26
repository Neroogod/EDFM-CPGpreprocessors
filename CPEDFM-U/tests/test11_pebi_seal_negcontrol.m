%% test11_pebi_seal_negcontrol
% CPEDFM on a genuinely unstructured grid: extruded PEBI (Voronoi)
% cells, where no structured concept (face tags, IJK, orientation
% branches) exists at all.
%
% Checks:
%   A. Pipeline runs end to end on the PEBI grid.
%   B. SEAL: a low-permeability fracture slanted with respect to all
%      axes and spanning the full domain cross-section splits the matrix
%      into 2 macro components; plain EDFM stays at 1 (leaks).
%   C. NEGATIVE CONTROL: a fracture spanning only the central part of
%      the domain does not compartmentalize (tip guard).
%   D. Sanity: multipliers in [0,1]; sealed faces consistent with their
%      cell pairs; projection figure written to tests/output.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupCPEDFMU.m'));
outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
failures = {};
physdim = [100, 100, 20];
nLayers = 6;

rng(42);
makeG = @() localPebiGrid(physdim, nLayers);

%% ---- Case B: full sealing fracture ------------------------------------
G = makeG();
fprintf('PEBI grid: %d cells, %d faces, faces/cell %d-%d\n', G.cells.num, ...
    G.faces.num, min(diff(G.cells.facePos)), max(diff(G.cells.facePos)));

fp = struct;
fp(1).points = [ 35  -5  -5;
                 55 105  -5;
                 75 105  25;
                 55  -5  25];       % slanted in x-y and x-z, full extent
fp(1).aperture = 10*milli*meter;
fp(1).poro = 0.4;
fp(1).perm = 1*nano*darcy;

[G, fp] = EDFMGridU(G, fp, 'Tolerance', tol, 'fracturelist', 1);
G = fracturematrixNNC_U(G, tol);
Gp = pMatFracNNCsCPEDFM_U(G, tol);
opsU = setupCPEDFMOpsTPFA_CPG(Gp, Gp.rock, tol);
opsE = setupShaleEDFMOpsTPFA(G, G.rock, tol);
fprintf('A. pipeline: PASS (%d host cells, %d sealed faces)\n', ...
    nnz(strcmp(G.nnc.type, 'fracmat interior')), numel(unique(Gp.nnc.pMMfaces)));

[nU, szU] = localComponents(G, opsU);
[nE, ~]   = localComponents(G, opsE);
nMat = G.Matrix.cells.num;
macroU = nnz(szU > 0.05*nMat);
fprintf('B. components: CPEDFM-U %d (%d macro) | EDFM %d\n', nU, macroU, nE);
if macroU == 2 && nE == 1 && sum(szU(szU > 0.05*nMat)) > 0.9*nMat
    fprintf('   PASS (2 macro compartments, EDFM leaks)\n');
else
    failures{end+1} = sprintf('B: macro U=%d, EDFM=%d', macroU, nE);
    fprintf('   FAIL\n');
end

% projection figure
try
    fig = figure('Visible', 'off', 'Position', [100 100 1100 480]);
    subplot(1,2,1); hold on
    plotGrid(G.Matrix, 'FaceColor', 'none', 'EdgeAlpha', 0.05);
    plotGrid(Gp, Gp.nnc.cells(strcmp(Gp.nnc.type,'fracmat interior'),2), ...
        'FaceColor', 'y', 'EdgeAlpha', 0.3);
    view(-35, 25); axis equal tight off; title('fracture cells (PEBI)');
    subplot(1,2,2); hold on
    plotGrid(G.Matrix, 'FaceColor', 'none', 'EdgeAlpha', 0.05);
    plotFaces(G.Matrix, unique(Gp.nnc.pMMfaces), 'FaceColor', 'b', 'EdgeAlpha', 0.3);
    view(-35, 25); axis equal tight off; title('continuous projection (sealed faces)');
    print(fig, fullfile(outDir, 'test11_projections.png'), '-dpng', '-r150');
    close(fig);
    fprintf('   projection figure saved to tests/output/test11_projections.png\n');
catch me
    fprintf('   figure generation skipped: %s\n', me.message);
end

%% Check D: sanity
mult = transmultpEDFM_CPG(Gp, tol);
okD = all(mult >= 0 & mult <= 1) && all(Gp.nnc.pMMarea > 0);
for q = 1:numel(Gp.nnc.pMMfaces)
    nb = Gp.faces.neighbors(Gp.nnc.pMMfaces(q),:);
    okD = okD && any(nb(:)' == Gp.nnc.pMMneighs(q,1)) ...
              && any(nb(:)' == Gp.nnc.pMMneighs(q,2));
end
if okD
    fprintf('D. sanity (multipliers, areas, face-pair consistency): PASS\n');
else
    failures{end+1} = 'D: sanity checks failed';
    fprintf('D. sanity: FAIL\n');
end

%% ---- Case C: partial fracture (negative control) -----------------------
G2 = makeG();
fp2 = struct;
fp2(1).points = [ 45  30   4;
                  55  70   4;
                  65  70  16;
                  55  30  16];      % central patch only
fp2(1).aperture = 10*milli*meter;
fp2(1).poro = 0.4;
fp2(1).perm = 1*nano*darcy;

[G2, fp2] = EDFMGridU(G2, fp2, 'Tolerance', tol, 'fracturelist', 1);
G2 = fracturematrixNNC_U(G2, tol);
G2p = pMatFracNNCsCPEDFM_U(G2, tol);
ops2 = setupCPEDFMOpsTPFA_CPG(G2p, G2p.rock, tol);
[n2, ~] = localComponents(G2, ops2);
fprintf('C. partial fracture: %d component(s), %d sealed faces\n', ...
    n2, numel(unique(G2p.nnc.pMMfaces)));
if n2 == 1
    fprintf('   PASS (no artificial compartmentalization)\n');
else
    failures{end+1} = sprintf('C: partial fracture gave %d components', n2);
    fprintf('   FAIL\n');
end

%% Summary
if isempty(failures)
    fprintf('\n=== test11: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test11: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test11 failed');
end

%% ---- local helpers ------------------------------------------------------
function G = localPebiGrid(physdim, nLayers)
% extruded PEBI grid scaled to physdim
G2D = pebiGrid2D(physdim(1)/12, physdim(1:2));
G = makeLayeredGrid(G2D, nLayers);
G.nodes.coords(:,3) = G.nodes.coords(:,3) * physdim(3)/nLayers;
G = computeGeometry(G);
G.rock = makeRock(G, 100*milli*darcy, 0.25);
end

function [nComp, compSizes] = localComponents(G, ops)
nMat = G.Matrix.cells.num;
inter = all(G.Matrix.faces.neighbors ~= 0, 2);
Nmm = G.Matrix.faces.neighbors(inter,:);
Tm = ops.T_all(1:G.Matrix.faces.num); Tm = Tm(inter);
keep = Tm > max(Tm)*1e-12;
g = graph(Nmm(keep,1), Nmm(keep,2), [], nMat);
comp = conncomp(g);
nComp = max(comp);
compSizes = accumarray(comp(:), 1);
end
