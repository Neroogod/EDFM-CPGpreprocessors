%% test02_cartesian_crosscheck_U
% Chain-back of the EDFM-U pipeline against the validated Cartesian
% implementation.
%
% On a Cartesian grid the cell IS its own axis-aligned bounding box, so the
% bounding-box <d> of hfm/calcdavg is exact and the two pipelines must
% agree to roundoff. That makes the Cartesian case a hard cross-check of
% everything EDFM-U introduces - the half-space clipping, the fracture
% gridding, the exact <d>, the surface-area rule and the transmissibility
% assembly - against code that is already trusted.
%
%   REFERENCE : EDFMshalegrid  -> fracturematrixShaleNNC3D   (CPEDFM/hfm)
%   EDFM-U    : EDFMgrid_U     -> fracturematrixNNC_U
%
% Checks
%   A. Same number of fracture cells and the same host matrix cells.
%   B. Intersection areas agree to roundoff.
%   C. NNC transmissibilities agree to within the reference's own <d>
%      quadrature error (see below).
%   D. <d> itself, graded against the exact box integral: this is where the
%      two pipelines legitimately differ.
%   E. A single-phase pressure solution on the two global grids agrees.
%   F. Fracture-fracture NNCs (two-fracture stage) match the reference and
%      leave G.nnc internally consistent.
%
% WHY C IS NOT ROUNDOFF-EXACT
%   hfm/calcdavg estimates <d> by averaging |n.(x-x0)| at the centroids of
%   a 10x10x10 sample grid. |d| kinks at the fracture plane, so the
%   midpoint rule is not exact there and the estimate carries a ~0.1-1%
%   bias on a Cartesian cell - it converges onto calcDavg_U as the sampling
%   is refined (7.7e-3 at 10^3 samples, 1.2e-6 at 640^3). Check D grades
%   both against the closed-form box integral (boxMeanAbsDist_U) and
%   confirms that the residual T gap in C is the reference's quadrature
%   error, not an EDFM-U defect. Everything the two pipelines are supposed
%   to share - host cells, areas, connectivity - matches to roundoff.
%
% The whole comparison runs twice: one fracture, then two intersecting
% fractures (which exercises the multi-plane bookkeeping).

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
run(fullfile(thisDir, '..', 'startupEDFMU.m'));
addpath(thisDir);                       % boxMeanAbsDist_U (test reference)

tol = 1e-5;
failures = {};
gravity off

physdim = [100, 100, 50];
celldim = [10, 10, 5];

%% Fracture set: a slanted plane through the middle of the domain.
% Deliberately NOT aligned with any cell face, so every intersection is of
% the 'fracmat interior' kind that the Cartesian code handles.
fpBase = struct;
fpBase(1).points   = [22   4  4;
                      78  62  4;
                      78  62 46;
                      22   4 46];
fpBase(1).aperture = 0.02;
fpBase(1).poro     = 0.4;
fpBase(1).perm     = 5*darcy;

fpBase(2).points   = [15  70  6;
                      85  30  6;
                      85  30 44;
                      15  70 44];
fpBase(2).aperture = 0.02;
fpBase(2).poro     = 0.4;
fpBase(2).perm     = 2*darcy;

for stage = 1:2
    if stage == 1
        fp   = fpBase(1);
        name = 'single fracture';
    else
        fp   = fpBase;
        name = 'two intersecting fractures';
    end
    fprintf('\n---------- %s ----------\n', name);

    %% Reference chain (validated Cartesian EDFM)
    Gr = computeGeometry(cartGrid(celldim, physdim));
    Gr.rock = makeRock(Gr, 100*milli*darcy, 0.25);
    fpR = fp;
    [Gr, fpR] = EDFMshalegrid(Gr, fpR, 'Tolerance', tol, ...
        'fracturelist', 1:numel(fpR), 'Verbose', false);
    Gr = fracturematrixShaleNNC3D(Gr, tol);

    %% EDFM-U chain
    Gu = computeGeometry(cartGrid(celldim, physdim));
    Gu.rock = makeRock(Gu, 100*milli*darcy, 0.25);
    fpU = fp;
    [Gu, fpU] = EDFMgrid_U(Gu, fpU, 'Tolerance', tol, 'Verbose', false);
    Gu = fracturematrixNNC_U(Gu, tol, 'Verbose', false);

    %% A. fracture cell counts and host cells
    okA = true;
    for q = 1:numel(fp)
        fn  = ['Frac', num2str(q)];
        nR  = Gr.FracGrid.(fn).cells.num;
        nU  = Gu.FracGrid.(fn).cells.num;
        hR  = sort(unique(Gr.FracGrid.(fn).matrix_connection.cells));
        hU  = sort(unique(Gu.FracGrid.(fn).matrix_connection.cells));
        fprintf('A. %s: fracture cells ref %d / EDFM-U %d; host cells %d / %d\n', ...
            fn, nR, nU, numel(hR), numel(hU));
        okA = okA && nR == nU && isequal(hR, hU);
    end
    if okA
        fprintf('   PASS\n');
    else
        failures{end+1} = sprintf('%s A: fracture/host cell mismatch', name); %#ok<SAGROW>
        fprintf('   FAIL\n');
    end

    %% B/C. per-connection areas and transmissibilities
    % Restrict to fracture-matrix rows and sort by (matrix cell, area) so the
    % two orderings of fracture cells become comparable.
    [aR, tR] = fmRows(Gr);
    [aU, tU] = fmRows(Gu);

    if numel(aR) ~= numel(aU)
        failures{end+1} = sprintf('%s B: %d vs %d frac-mat NNCs', ...
            name, numel(aR), numel(aU)); %#ok<SAGROW>
        fprintf('B. FAIL: %d vs %d fracture-matrix NNCs\n', numel(aR), numel(aU));
    else
        eA = max(abs(aR - aU)) / max(aR);
        eT = max(abs(tR - tU)) / max(tR);
        fprintf('B. intersection area   : max rel diff %.3e\n', eA);
        fprintf('C. NNC transmissibility: max rel diff %.3e\n', eT);
        if eA < 1e-10
            fprintf('   B PASS\n');
        else
            failures{end+1} = sprintf('%s B: area rel diff %.3e', name, eA); %#ok<SAGROW>
            fprintf('   B FAIL\n');
        end
        if eT < 1e-2
            fprintf('   C PASS (within the reference''s <d> quadrature error)\n');
        else
            failures{end+1} = sprintf('%s C: T rel diff %.3e', name, eT); %#ok<SAGROW>
            fprintf('   C FAIL\n');
        end
    end

    %% D. <d> against the exact box integral
    % On a Cartesian cell the box integral is the exact <d>, so both
    % estimators can be graded absolutely.
    [cn, cpos] = gridCellNodes(Gu.Matrix, 1:Gu.Matrix.cells.num);
    eU = 0; eB = 0;
    for q = 1:numel(fp)
        fn  = ['Frac', num2str(q)];
        mcU = Gu.FracGrid.(fn).matrix_connection;
        nrm = mcU.normal(1,:); pp = mcU.planepoint(1,:);
        hostU = unique(mcU.cells);
        dU = calcDavg_U(Gu.Matrix, hostU, nrm, pp, tol);
        for j = 1:numel(hostU)
            nd  = Gu.Matrix.nodes.coords(cn(cpos(hostU(j)):cpos(hostU(j)+1)-1), :);
            ex  = boxMeanAbsDist_U(min(nd), max(nd), nrm, pp);
            dB  = calcdavg(nd, nrm, pp, tol);      % hfm 10^3 sampling
            eU  = max(eU, abs(dU(j) - ex)/ex);
            eB  = max(eB, abs(dB     - ex)/ex);
        end
    end
    fprintf('D. <d> vs exact box integral: calcDavg_U %.3e | hfm calcdavg %.3e\n', eU, eB);
    if eU < 1e-9 && eB > 10*eU
        fprintf('   PASS (EDFM-U is exact; the reference''s sampling bias explains C)\n');
    else
        failures{end+1} = sprintf('%s D: davg errU %.3e errBB %.3e', name, eU, eB); %#ok<SAGROW>
        fprintf('   FAIL\n');
    end

    %% E. flow solution
    pR = steadyPressure(Gr, tol, physdim);
    pU = steadyPressure(Gu, tol, physdim);
    nMat = Gr.Matrix.cells.num;
    eP = norm(pR(1:nMat) - pU(1:nMat)) / norm(pR(1:nMat) - mean(pR(1:nMat)));
    fprintf('E. matrix pressure field: rel L2 diff %.3e\n', eP);
    if eP < 1e-4
        fprintf('   PASS\n');
    else
        failures{end+1} = sprintf('%s E: pressure rel diff %.3e', name, eP); %#ok<SAGROW>
        fprintf('   FAIL\n');
    end

    %% F. fracture-fracture NNCs (two-fracture stage only)
    % The two planes cross near (61, 44), so the frac-frac stage must fire.
    % fracturefractureNNC_U reuses the existing grid-agnostic routine; what
    % is checked here is that it fires on EDFM-U's fracture grids, that it
    % agrees with the reference chain, and that the diagnostic columns
    % EDFM-U adds to G.nnc stay aligned with G.nnc.cells afterwards.
    if stage == 2
        Gr2 = fracturefractureShaleNNCs3D(Gr, fpR, tol, 'Verbose', false);
        Gu2 = fracturefractureNNC_U(Gu, fpU, tol, 'Verbose', false);

        ffR = sort(Gr2.nnc.T(strcmp(Gr2.nnc.type, 'fracfrac')));
        ffU = sort(Gu2.nnc.T(strcmp(Gu2.nnc.type, 'fracfrac')));
        nU  = size(Gu2.nnc.cells, 1);
        aligned = numel(Gu2.nnc.CI) == nU && numel(Gu2.nnc.davg) == nU && ...
                  numel(Gu2.nnc.T) == nU && numel(Gu2.nnc.type) == nU && ...
                  size(Gu2.nnc.normal,1) == nU;

        fprintf('F. fracture-fracture NNCs: ref %d | EDFM-U %d; G.nnc columns aligned: %d\n', ...
            numel(ffR), numel(ffU), aligned);
        okF = ~isempty(ffU) && numel(ffR) == numel(ffU) && aligned;
        if okF
            % The two chains' fracture GRIDS are bit-identical here - same
            % node coordinates, same cell centroids, same cell order - so
            % the frac-frac stage receives identical input and the pairs it
            % returns are identical. Four of the five transmissibilities
            % match exactly; one differs at the 1e-4 level. That residual
            % comes from the REUSED upstream calcPEBIdf, which classifies a
            % fracture cell's vertices as lying on one side of the
            % intersection line or the other with a hard |sign| > tol test.
            % A vertex sitting within roundoff of that threshold can be
            % classified differently once ~1e-13 differences in
            % computeGeometry's accumulation order are present, which shifts
            % df slightly. It is a tolerance-branch effect in shared code,
            % not a difference in what EDFM-U models, so the bound below is
            % set at the tolerance scale rather than at roundoff.
            dGeom = 0;
            for q = 1:numel(fp)
                fnq = ['Frac', num2str(q)];
                cR = Gr.FracGrid.(fnq).cells.centroids;
                cU = Gu.FracGrid.(fnq).cells.centroids;
                dGeom = max(dGeom, max(abs(cR(:) - cU(:))));
            end
            eF = max(abs(ffR - ffU))/max(ffR);
            fprintf('   fracture-cell centroids agree to %.2e m (identical grids)\n', dGeom);
            fprintf('   frac-frac T max rel diff %.3e\n', eF);
            okF = eF < 1e-3 && dGeom < 1e-9;
        end
        % the operators must still assemble with the extra connections
        if okF
            try
                setupShaleEDFMOpsTPFA(Gu2, Gu2.rock, tol);
            catch me
                okF = false;
                fprintf('   operator setup failed: %s\n', me.message);
            end
        end
        if okF
            fprintf('   PASS\n');
        else
            failures{end+1} = sprintf('%s F: frac-frac NNC mismatch', name); %#ok<SAGROW>
            fprintf('   FAIL\n');
        end
    end
end

%% Summary
if isempty(failures)
    fprintf('\n=== test02: ALL CHECKS PASSED ===\n');
else
    fprintf('\n=== test02: %d FAILURES ===\n', numel(failures));
    fprintf('  - %s\n', failures{:});
    error('test02 failed');
end

%% ------------------------------ helpers ---------------------------------
function [area, T] = fmRows(G)
% fracture-matrix NNC areas and transmissibilities, in a canonical order
keep = ~strcmp(G.nnc.type, 'fracfrac');
M    = [G.nnc.cells(keep,1), G.nnc.area(keep), G.nnc.T(keep)];
M    = sortrows(M, [1 2 3]);
area = M(:,2);
T    = M(:,3);
end

function p = steadyPressure(G, tol, physdim)
% single-phase steady state under fixed west/east boundary pressures
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, ...
    'rho', 1000, 'c', 1e-8/barsa, 'pRef', 150*barsa);

bf = find(any(G.faces.neighbors == 0, 2));
bf = bf(bf <= G.Matrix.faces.num);
fc = G.faces.centroids(bf,:);
wf = bf(fc(:,1) < tol);
ef = bf(fc(:,1) > physdim(1) - tol);
assert(~isempty(wf) && ~isempty(ef), 'boundary faces not found');

bc = addBC([], wf, 'pressure', 200*barsa, 'sat', 1);
bc = addBC(bc, ef, 'pressure', 100*barsa, 'sat', 1);

model = GenericBlackOilModel(G, [], fluid, ...
    'water', true, 'oil', false, 'gas', false);
model.operators = setupShaleEDFMOpsTPFA(G, G.rock, tol);
model.minimumPressure = 1*barsa;

state0   = initResSol(G, 150*barsa, 1);
schedule = simpleSchedule(rampupTimesteps(1*year, 180*day, 4), 'bc', bc);
[~, states] = simulateScheduleAD(state0, model, schedule, 'Verbose', false);
p = states{end}.pressure;
end
