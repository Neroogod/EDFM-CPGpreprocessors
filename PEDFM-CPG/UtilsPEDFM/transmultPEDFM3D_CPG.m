function transmultp = transmultPEDFM3D_CPG(G, tol)
%TRANSMULTPEDFM3D_CPG Matrix-matrix transmissibility multipliers for pEDFM on CPG.
%
% SYNOPSIS:
%   transmultp = transmultPEDFM3D_CPG(G, tol)
%
% Implements Eq. 27 of Olorode, Wang & Rashid (2020), SPE-201243-PA,
%
%     T_pM-M = T_M-M * (A_ii - A_if_perp_x) / A_ii,
%
% as a per-face multiplier.
%
% RETURNS:
%   transmultp - one multiplier per face of the global grid. Faces carrying
%                fracture projections get max(A_face - A_proj, 0)/A_face;
%                every other face gets 1.
%
% WHY NOT transmultpEDFM.m
%   The Cartesian original recovers the face of each projected pair by
%   matching cell pairs against G.faces.neighbors, which assumes exactly
%   one shared face per pair of cells. That is false across a fault, where
%   a logical interface is split into several faces: the lookup then
%   returns the wrong face, or several, and the multiplier is applied in
%   the wrong place. pMatFracNNCsPEDFM3D_CPG records the global face index
%   of every projection explicitly in G.nnc.pMMfaces, so the multiplier is
%   a direct accumulation over faces.
%
%   Areas accumulate additively across fracture planes; colocated
%   projections were already union-corrected (Eq. 29) at the NNC stage. A
%   projected area exceeding the face area clamps the face shut, which is
%   the sealing behaviour a low-conductivity fracture requires.

    if nargin < 2 || isempty(tol), tol = 1e-5; end

    transmultp = ones(G.faces.num, 1);

    if ~isfield(G, 'nnc') || ~isfield(G.nnc, 'pMMfaces') || isempty(G.nnc.pMMfaces)
        warning('transmultPEDFM3D_CPG:noProjections', ...
            'G.nnc.pMMfaces missing or empty; returning unit multipliers.');
        return
    end

    [uf, ~, ic] = unique(G.nnc.pMMfaces);
    projArea = accumarray(ic, G.nnc.pMMarea);

    oldArea = G.faces.areas(uf);
    newArea = oldArea - projArea;
    newArea(newArea < tol) = 0;

    transmultp(uf) = newArea ./ oldArea;
end
