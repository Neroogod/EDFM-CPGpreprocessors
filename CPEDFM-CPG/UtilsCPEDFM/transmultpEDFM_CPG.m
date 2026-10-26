function transmultp = transmultpEDFM_CPG(G, tol)
% Matrix-matrix transmissibility multipliers from CPEDFM projected areas,
% face-indexed (corner-point safe).
%
% SYNOPSIS:
%   transmultp = transmultpEDFM_CPG(G, tol)
%
% RETURNS:
%   transmultp - one multiplier per face of the global grid. For every
%                face carrying fracture projections the remaining open
%                area fraction is  max(A_face - A_projected, 0) / A_face;
%                all other faces get 1.
%
% Replaces transmultpEDFM.m for corner-point grids. The original recovers
% the face for each projected pair by matching cell pairs against
% G.Matrix.faces.neighbors, which assumes exactly one shared face per
% pair - false across faults, where an interface is split into several
% faces. Here pMatFracNNCsCPEDFM_CPG records the global face index of
% every projection explicitly (G.nnc.pMMfaces), so the multiplier is a
% direct accumulation:
%   - several fractures projecting onto one face compose additively
%     (their areas were already union-corrected at the NNC stage),
%   - a projected area exceeding the face area (a neighborhood collapsing
%     several partial areas onto one face) clamps the face shut, which is
%     exactly the sealing behavior the stair-step projection requires.

if nargin < 2, tol = 1e-5; end

transmultp = ones(G.faces.num, 1);

if ~isfield(G.nnc, 'pMMfaces') || isempty(G.nnc.pMMfaces)
    warning('transmultpEDFM_CPG:noProjections', ...
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
