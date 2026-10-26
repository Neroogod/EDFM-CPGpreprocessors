function G = fracturematrixNNC_U(G, tol, varargin)
% Fracture-matrix NNCs for every fracture plane of an unstructured EDFM grid.
%
% SYNOPSIS:
%   G = fracturematrixNNC_U(G, tol)
%   G = fracturematrixNNC_U(G, tol, 'Verbose', false)
%
% PARAMETERS:
%   G   - global EDFM grid from EDFMgrid_U (has G.Matrix, G.FracGrid, G.rock).
%   tol - geometric tolerance.
%
% RETURNS:
%   G - with G.nnc extended by one row per fracture-matrix connection:
%         cells      Nx2 [matrix cell, fracture cell]
%         T          Nx1 transmissibility
%         type       Nx1 'fracmat boundary' | 'fracmat interior'
%         area       Nx1 intersection area
%         normal     Nx3 fracture plane normal
%         planepoint Nx3 point on the fracture plane
%         CI         Nx1 connectivity index (diagnostic)
%         davg       Nx1 average normal distance (diagnostic)
%
% TRANSMISSIBILITY:
%   T = CI * k_eff, with k_eff the pore-volume-weighted harmonic mean of
%   the fracture permeability and the matrix permeability RESOLVED ALONG
%   THE FRACTURE NORMAL (permNormal_U) - flow into a fracture is normal to
%   its plane, so an anisotropic matrix must be sampled in that direction.
%   For isotropic rock this reduces exactly to the Cartesian
%   fracturematrixShaleNNC3D expression.
%
% SEE ALSO:
%   EDFMgrid_U, fracMatNNC_U, calcDavg_U, fracturefractureNNC_U

opt = struct('Verbose', true);
opt = merge_options(opt, varargin{:});
verbose = opt.Verbose;

assert(isfield(G, 'Matrix') && isfield(G, 'FracGrid'), ...
    'G must be a global EDFM grid (run EDFMgrid_U first).');

t1   = tic;
Gmat = G.Matrix;

if ~isfield(G.nnc, 'cells')
    G.nnc            = struct;
    G.nnc.cells      = zeros(0,2);
    G.nnc.T          = [];
    G.nnc.type       = {};
    G.nnc.area       = [];
    G.nnc.normal     = zeros(0,3);
    G.nnc.planepoint = zeros(0,3);
end
if ~isfield(G.nnc, 'CI'),   G.nnc.CI   = []; end
if ~isfield(G.nnc, 'davg'), G.nnc.davg = []; end

if isfield(G.rock, 'poro')
    pv = poreVolume(G, G.rock);
else
    pv = G.cells.volumes;
end

fracNames = fieldnames(G.FracGrid);
dispif(verbose, 'Processing fracture-matrix NNCs (%d planes)...\n', numel(fracNames));

for i = 1:numel(fracNames)
    fn = fracNames{i};
    F  = G.FracGrid.(fn);

    if F.matrixnnc
        dispif(verbose, '  %s already processed; skipped.\n', fn);
        continue
    end

    [cellsi, CIi, areai, typei, davgi] = fracMatNNC_U(F, Gmat, tol);

    % permeability of each side, resolved along the fracture normal
    nrm = F.matrix_connection.normal;
    k1  = permNormal_U(G.rock.perm(cellsi(:,1), :), nrm);   % matrix
    k2  = permNormal_U(G.rock.perm(cellsi(:,2), :), nrm);   % fracture

    w1 = pv(cellsi(:,1)) ./ k1;
    w2 = pv(cellsi(:,2)) ./ k2;
    wt = pv(cellsi(:,1)) +  pv(cellsi(:,2));
    Ti = CIi .* (wt ./ (w1 + w2));

    G.nnc.cells      = [G.nnc.cells;      cellsi];
    G.nnc.type       = [G.nnc.type;       typei];
    G.nnc.T          = [G.nnc.T;          Ti];
    G.nnc.area       = [G.nnc.area;       areai];
    G.nnc.normal     = [G.nnc.normal;     nrm];
    G.nnc.planepoint = [G.nnc.planepoint; F.matrix_connection.planepoint];
    G.nnc.CI         = [G.nnc.CI;         CIi];
    G.nnc.davg       = [G.nnc.davg;       davgi];

    G.FracGrid.(fn).matrixnnc = true;

    dispif(verbose, '  %s: %d NNCs, <d> in [%.4g, %.4g] m\n', ...
        fn, numel(Ti), min(davgi), max(davgi));
end

dispif(verbose, 'Fracture-matrix NNCs done in %.3g s.\n', toc(t1));
end
