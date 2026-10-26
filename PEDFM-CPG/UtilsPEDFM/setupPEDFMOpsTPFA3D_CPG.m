function s = setupPEDFMOpsTPFA3D_CPG(G, rock, tol)
%SETUPPEDFMOPSTPFA3D_CPG TPFA operators for 3D pEDFM on corner-point grids.
%
% SYNOPSIS:
%   s = setupPEDFMOpsTPFA3D_CPG(G, rock, tol)
%
% Builds the standard EDFM TPFA operators, then applies the pEDFM
% matrix-matrix transmissibility multipliers of Eq. 27
% (Olorode, Wang & Rashid 2020, SPE-201243-PA) before assembling the
% internal-connection transmissibility vector.
%
% Structurally a clone of setupPEDFMOpsTPFA.m with one substitution: the
% multiplier comes from transmultPEDFM3D_CPG (face-indexed, fault-safe)
% instead of transmultpEDFM (cell-pair-indexed, Cartesian-only). Everything
% else - the NNC-aware neighbour list, T assembly, and the R2019b+
% Div/AccDiv operators - matches the original so results stay comparable.
%
% PARAMETERS:
%   G    - grid after pMatFracNNCsPEDFM3D_CPG.
%   rock - rock structure (G.rock is used for the base operators).
%   tol  - geometric tolerance.

    if nargin < 3 || isempty(tol), tol = 1e-5; end

    assert(isfield(G, 'nnc') && isfield(G.nnc, 'pMMfaces'), ...
        'setupPEDFMOpsTPFA3D_CPG:noProjections', ...
        'Run pMatFracNNCsPEDFM3D_CPG before building pEDFM operators.');

    s = setupEDFMOperatorsTPFA(G, G.rock, tol);

    transmultpMM = transmultPEDFM3D_CPG(G, tol);
    T = s.T_all .* transmultpMM;      % Eq. 27
    s.T_all = T;
    T = [T; G.nnc.T];
    T = T(s.internalConn);
    s.T = T;

    % Adsorption support, as in the original.
    if isfield(G.rock, 'shaleMechanisms') && isfield(G.rock.shaleMechanisms, 'sorption')
        s.isSorbed = rock.isMatrix;
        s.gv = s.isSorbed .* (1 - rock.poro) .* G.cells.volumes;
    end
    if isfield(rock, 'isMatrix')
        s.isMatrix = rock.isMatrix;
    end

    % Migrate hfm operators from 2018b to 2019b conventions (as in original).
    n = size(s.N, 1);
    C = (s.C)';
    if n == 0
        s.AccDiv = @(acc, flux) acc;
        s.Div = @(x) zeros(G.cells.num, 1);
    else
        s.AccDiv = @(acc, flux) acc + C*flux;
        s.Div = @(x) C*x;
    end
end
