function s = setupCPEDFMOpsTPFA_CPG(G, rock, tol, varargin)
% TPFA operators for CPEDFM on corner-point grids.
%
% SYNOPSIS:
%   s = setupCPEDFMOpsTPFA_CPG(G, rock, tol)
%
% Clone of setupPEDFMOpsTPFA.m (CPEDFM/UpdateMay27_2023/pEDFM_nnc) with a
% single substitution: the matrix-matrix transmissibility multiplier comes
% from transmultpEDFM_CPG (face-indexed, fault-safe) instead of
% transmultpEDFM (cell-pair-indexed, Cartesian-only). Everything else -
% the NNC-aware neighbor list, T assembly, and the R2019b+ Div/AccDiv
% operators - matches the original so results are comparable line for
% line.

s = setupEDFMOperatorsTPFA(G, G.rock, tol);

transmultpMM = transmultpEDFM_CPG(G, tol);
T = s.T_all .* transmultpMM; % apply pMM transmissibility multiplier
s.T_all = T;
T = [T; G.nnc.T];
T = T(s.internalConn);
s.T = T;

% Include adsorption (as in the original)
if isfield(G.rock, 'shaleMechanisms') && isfield(G.rock.shaleMechanisms, 'sorption')
    s.isSorbed = rock.isMatrix;
    s.gv = s.isSorbed .* (1 - rock.poro) .* G.cells.volumes;
end
if isfield(rock, 'isMatrix')
    s.isMatrix = rock.isMatrix;
end

% Migrate hfm operators from 2018b to 2019b conventions (as in original)
n = size(s.N, 1);
C = (s.C)'; % transpose once for speed
if n == 0
    s.AccDiv = @(acc, flux) acc;
    s.Div = @(x) zeros(G.cells.num, 1);
else
    s.AccDiv = @(acc, flux) acc + C*flux;
    s.Div = @(x) C*x;
end
end
