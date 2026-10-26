%% test37_projection_nonuniform
% Risk R1: are the pEDFM / CPEDFM projections correct on a grid with
% non-uniform spacing?
%
% Every other test exercises them on a UNIFORM cartGrid. Eagle Ford's
% meshHFsystem deliberately refines near the fractures (nxRefine=9 beside
% them, nxRefine_small=4 between clusters), so the projection runs on cells
% of very different sizes. A method that assumed uniform spacing would give
% plausible but wrong projected areas, and nothing downstream would notice:
% getEmbeddedFractureMultipliers clamps with max(faceArea - projArea, 0), so
% an over-assignment would silently become a fully sealed face.
%
% CHECK 1 - pEDFM against an analytic value. For a planar fracture of area A
% with unit normal n, the shadow gathered on faces normal to direction d is
% exactly A*|n_d|, whatever the grid. The fracture is held clear of every
% domain boundary: shadow reaching a boundary face has no cell pair, so it is
% never recorded, and the invariant would undercount through no fault of the
% projection.
%
% CHECK 2 - grid independence. The two grids distribute the same shadow over
% different numbers of faces; the totals must still agree.
%
% CHECK 3 - CPEDFM cannot be held to the same invariant. Its flip-face rule
% seals every host-adjacent face whose neighbour centroids straddle the
% plane, which is a cut of the connectivity graph rather than a shadow, so
% the sealed area legitimately depends on the grid. What must hold is that no
% face is over-assigned and every multiplier stays in [0,1].

clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
cpgDir = fullfile(thisDir, '..', '..', 'CPEDFM-CPG');
if isfolder(cpgDir)
    addpath(cpgDir);
    if exist('startupCPEDFMCPG', 'file') == 2, startupCPEDFMCPG(); end
end
mrstModule add hfm_v2

tol = 1e-5;
physdim = [120 100 20];

% Inset on every side, so no shadow can land on a domain boundary.
fracPts = [45 20 4; 75 80 4; 75 80 16; 45 20 16];
e1 = fracPts(2,:) - fracPts(1,:);
e2 = fracPts(4,:) - fracPts(1,:);
Afrac  = norm(cross(e1, e2));
nrm    = cross(e1, e2)/Afrac;
expect = Afrac*abs(nrm);

fprintf('fracture area   : %.6f m2, normal [%.4f %.4f %.4f]\n', Afrac, nrm);
fprintf('expected shadow : x %.6f  y %.6f  z %.6f\n\n', expect);

grids = struct('name', {}, 'G', {});
grids(1).name = 'uniform';
grids(1).G = computeGeometry(cartGrid([12 10 8], physdim));

xs = [0 20 25 30 40 50 60 70 80 90 100 110 120];
ys = [0 25 30 35 45 55 65 75 85 100];
zs = [0 2 4 6 10 12 14 16 20];
grids(2).name = 'non-uniform';
grids(2).G = computeGeometry(tensorGrid(xs, ys, zs));

fprintf('uniform     : %d cells, uniform spacing\n', grids(1).G.cells.num);
fprintf('non-uniform : %d cells, dx contrast %.1f:1, dz contrast %.1f:1\n\n', ...
    grids(2).G.cells.num, max(diff(xs))/min(diff(xs)), ...
    max(diff(zs))/min(diff(zs)));

failures = {};
report = @(ok, name, detail) fprintf('%s %-38s %s\n', ...
    subsref({'FAIL', 'PASS'}, substruct('{}', {ok + 1})), name, detail);

shadow = nan(numel(grids), 3);      % pEDFM totals, for the cross-grid check

for gi = 1:numel(grids)
    Gm = grids(gi).G;
    Gm.rock = makeRock(Gm, 100*milli*darcy, 0.2);
    frac = struct('points', fracPts, 'aperture', 5*milli*meter, ...
        'poro', 0.4, 'perm', 1*nano*darcy);
    [Gb, frac] = EDFMshalegrid(Gm, frac, 'Tolerance', tol, ...
        'plotgrid', false, 'fracturelist', 1); %#ok<ASGLU>
    Gb = fracturematrixShaleNNC3D(Gb, tol);

    methods = struct('name', {'pEDFM', 'CPEDFM'}, 'G', {[], []});
    methods(1).G = preparePEDFMExport(pMatFracNNCs3D(Gb, tol));
    methods(2).G = pMatFracNNCsCPEDFM_CPG(Gb, tol);

    for mi = 1:numel(methods)
        Gp = methods(mi).G;
        label = sprintf('%s / %s', grids(gi).name, methods(mi).name);

        f = Gp.nnc.pMMfaces(:);
        a = Gp.nnc.pMMarea(:);
        assert(~isempty(f), 'no projection records for %s', label);

        N = Gp.Matrix.faces.normals(f, :);
        [~, dir] = max(abs(N), [], 2);
        got = arrayfun(@(d) sum(a(dir == d)), 1:3);

        % No face may receive more shadow than it has area - this is the
        % failure the multiplier clamp would hide.
        [uf, ~, grp] = unique(f);
        perFace = accumarray(grp, a);
        ratio = max(perFace ./ Gp.Matrix.faces.areas(uf));

        eta = getEmbeddedFractureMultipliers(Gp, tol);
        etaOK = all(eta >= -1e-12 & eta <= 1 + 1e-12);

        if strcmp(methods(mi).name, 'pEDFM')
            shadow(gi, :) = got;
            relErr = abs(got - expect) ./ max(expect, 1);
            ok = all(relErr < 1e-9) && ratio <= 1 + 1e-9 && etaOK;
            report(ok, label, sprintf( ...
                'x %.6f y %.6f (analytic), max area/face %.4f', got(1), got(2), ratio));
            if ~ok
                failures{end+1} = label; %#ok<AGROW>
                fprintf('        rel err [%.3e %.3e %.3e], eta in range %d\n', relErr, etaOK);
            end
        else
            % Continuous sealing: more faces than the pure shadow, by design.
            ok = ratio <= 1 + 1e-9 && etaOK && sum(got) >= sum(shadow(gi, :)) - 1e-9;
            report(ok, label, sprintf( ...
                'sealed %.1f m2 over %d faces (shadow %.1f), max area/face %.4f', ...
                sum(got), numel(uf), sum(shadow(gi, :)), ratio));
            if ~ok
                failures{end+1} = label; %#ok<AGROW>
                fprintf('        over-assigned %d, eta in range %d\n', ratio > 1, etaOK);
            end
        end
    end
end

%% ---- the projection must not depend on the grid ------------------------
okGrid = max(abs(shadow(1,:) - shadow(2,:))) < 1e-9;
report(okGrid, 'pEDFM is grid independent', sprintf( ...
    'uniform [%.6f %.6f] vs non-uniform [%.6f %.6f]', ...
    shadow(1,1), shadow(1,2), shadow(2,1), shadow(2,2)));
if ~okGrid, failures{end+1} = 'grid independence'; end

%% ---- summary -----------------------------------------------------------
fprintf('\n');
if isempty(failures)
    fprintf('=== test37: ALL CHECKS PASSED ===\n');
else
    fprintf('=== test37: %d FAILURES: %s ===\n', numel(failures), ...
        strjoin(failures, ', '));
end
