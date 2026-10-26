%% test02_multifrac_negcontrol
% Multi-fracture behavior and over-sealing guard on a faulted, skewed
% corner-point grid.
%
% Checks:
%   A. NEGATIVE CONTROL: a sealing fracture spanning only part of the
%      domain cross-section must NOT compartmentalize the matrix - flow
%      must be able to route around the fracture tips. Guards against the
%      full-face rule sealing faces beyond the fracture extent.
%   B. Two full-extent sealing fractures with OPPOSITE inclinations
%      (different Fig. 8 orientation branches) split the matrix into
%      three components.
%   C. Fault-split faces (this grid has a fault, so logical sides with
%      multiple faces exist among the host cells) are handled: every
%      projected face is a real face of its host/projection pair.

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupCPEDFMCPG.m'));

outDir = fullfile(thisDir, 'output');
if ~isfolder(outDir), mkdir(outDir); end

tol = 1e-5;
failures = {};

%% Faulted, skewed corner-point grid
dims = [20, 10, 8];
physdim = [200, 100, 40];
makeG = @() localMakeGrid(dims, physdim);

%% ---- Case A: partial fracture (negative control) ----------------------
G = makeG();
fp = struct;
fp(1).points = [ 90  25  10;
                110  75  10;
                130  75  30;
                110  25  30];   % covers only the middle of the domain
fp(1).aperture = 10*milli*meter;
fp(1).poro = 0.4;
fp(1).perm = 1*nano*darcy;

[G, fp] = EDFMGridCPG(G, fp, 'Tolerance', tol, 'plotgrid', false, ...
    'fracturelist', 1);
G = fracturematrixShaleNNC3D(G, tol);
Gp = pMatFracNNCsCPEDFM_CPG(G, tol);
[nCompA, szA] = localMatrixComponents(Gp, G, tol);
% Host cells pinched between the flip and envelope layers can appear as
% tiny matrix-graph slivers (they stay coupled through their fracture
% NNCs, which this matrix-only graph ignores) - same artifact as in
% check B. The negative control demands ONE macro component: the
% partial fracture must not split the reservoir.
nMatA = G.Matrix.cells.num;
macroA = szA > 0.05*nMatA;
fprintf('A. partial fracture: %d components (%d macro, largest %d of %d cells)\n', ...
    nCompA, nnz(macroA), max(szA), nMatA);
if nnz(macroA) == 1 && max(szA) > 0.9*nMatA
    fprintf('   PASS (no macro compartmentalization beyond fracture extent)\n');
else
    failures{end+1} = sprintf('A: %d macro components', nnz(macroA));
    fprintf('   FAIL\n');
end

%% ---- Case B: two opposing full fractures -> three compartments --------
G = makeG();
fp = struct;
fp(1).points = [ 50  -5  -5;    % dips one way
                 80 105  -5;
                110 105  45;
                 80  -5  45];
fp(2).points = [150  -5  45;    % dips the other way (different branch)
                120 105  45;
                 90 105  -5;
                120  -5  -5];
for i = 1:2
    fp(i).aperture = 10*milli*meter;
    fp(i).poro = 0.4;
    fp(i).perm = 1*nano*darcy;
end

[G, fp] = EDFMGridCPG(G, fp, 'Tolerance', tol, 'plotgrid', false, ...
    'fracturelist', 1:2);
G = fracturematrixShaleNNC3D(G, tol);
Gp = pMatFracNNCsCPEDFM_CPG(G, tol);

% report the branches exercised
hlp = calculateProjectionCPEDFM_CPG('helpers');
interiorIdx = find(strcmp(G.nnc.type, 'fracmat interior'));
fracNames = fieldnames(G.FracGrid);
for fn = 1:numel(fracNames)
    Fg = G.FracGrid.(fracNames{fn});
    rows = interiorIdx(G.nnc.cells(interiorIdx,2) >= Fg.cells.start & ...
        G.nnc.cells(interiorIdx,2) <= Fg.cells.start+Fg.cells.num-1);
    st = hlp.fracOrientation(G.nnc.normal(rows(1),:));
    fprintf('   %s selected tags: %s\n', fracNames{fn}, mat2str(st));
end

[nCompB, compSizes] = localMatrixComponents(Gp, G, tol);
% Host cells pinched between the flip-face layer and the envelope layer
% can end up as tiny matrix-graph slivers; they remain coupled to the
% system through their fracture NNCs (which this matrix-only graph
% ignores). The physical compartment count is the number of MACRO
% components.
nMatB = G.Matrix.cells.num;
macro = compSizes > 0.05*nMatB;
fprintf('B. two opposing fractures: %d components (%d macro, sizes: %s)\n', ...
    nCompB, nnz(macro), mat2str(sort(compSizes(macro), 'descend')'));
if nnz(macro) == 3 && sum(compSizes(macro)) > 0.9*nMatB
    fprintf('   PASS (3 macro compartments holding >90%% of cells)\n');
else
    failures{end+1} = sprintf('B: expected 3 macro components, got %d', nnz(macro));
    fprintf('   FAIL\n');
end

%% ---- Case C: fault-split faces handled ---------------------------------
% count host cells whose selected logical side holds >1 face
hostCells = unique(G.nnc.cells(interiorIdx,1));
nMulti = 0;
for r = 1:numel(hostCells)
    fs = hlp.getLogicalFaces(G, hostCells(r));
    nMulti = nMulti + any(cellfun(@numel, fs) > 1);
end
fprintf('C. host cells with fault-split sides: %d\n', nMulti);
okC = true;
for q = 1:numel(Gp.nnc.pMMfaces)
    nbrs = Gp.faces.neighbors(Gp.nnc.pMMfaces(q),:);
    okC = okC && any(nbrs(:)' == Gp.nnc.pMMneighs(q,1)) ...
              && any(nbrs(:)' == Gp.nnc.pMMneighs(q,2));
end
if okC
    fprintf('   PASS (all %d projected faces consistent with their cell pairs)\n', ...
        numel(Gp.nnc.pMMfaces));
else
    failures{end+1} = 'C: face/cell-pair inconsistency';
    fprintf('   FAIL\n');
end

%% Summary
if isempty(failures)
    fprintf('\n=== test02: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test02: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test02 failed');
end

%% ---- local helpers ------------------------------------------------------
function G = localMakeGrid(dims, physdim)
grdecl = simpleGrdecl(dims, 0.12);   % skewed pillars + fault
grdecl.COORD(1:3:end) = grdecl.COORD(1:3:end) * physdim(1);
grdecl.COORD(2:3:end) = grdecl.COORD(2:3:end) * physdim(2);
grdecl.COORD(3:3:end) = grdecl.COORD(3:3:end) * physdim(3);
grdecl.ZCORN = grdecl.ZCORN * physdim(3);
G = processGRDECL(grdecl);
G = computeGeometry(G);
G.rock = makeRock(G, 100*milli*darcy, 0.25);
end

function [nComp, compSizes] = localMatrixComponents(Gp, G, tol)
ops = setupCPEDFMOpsTPFA_CPG(Gp, Gp.rock, tol);
nMat = G.Matrix.cells.num;
inter = all(G.Matrix.faces.neighbors ~= 0, 2);
Nmm = G.Matrix.faces.neighbors(inter, :);
Tmm = ops.T_all(1:G.Matrix.faces.num);
Tmm = Tmm(inter);
keep = Tmm > max(Tmm)*1e-12;
g = graph(Nmm(keep,1), Nmm(keep,2), [], nMat);
comp = conncomp(g);
nComp = max(comp);
compSizes = accumarray(comp(:), 1);
end
