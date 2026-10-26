function coarse = buildCoarseCPEDFM(G, ops, rock, p)
% Build a coarse connection-list model of a CPEDFM system from a
% (barrier-aware) partition.
%
% SYNOPSIS:
%   coarse = buildCoarseCPEDFM(Gp, ops, rock, p)
%
% PARAMETERS:
%   G    - global EDFM grid after the projection stage.
%   ops  - operators from setupCPEDFMOpsTPFA_CPG (connection list ops.N
%          with transmissibilities ops.T, multipliers applied, NNCs
%          included).
%   rock - rock of the global grid.
%   p    - partition vector over ALL cells (matrix + fracture), e.g.
%          from partitionCPEDFM_MS. With a barrier-aware partition, no
%          block straddles a sealed surface, so the coarse model
%          preserves the compartmentalization of the fine model exactly:
%          sealed fine faces carry zero transmissibility and therefore
%          contribute nothing to any coarse connection.
%
% RETURNS:
%   coarse - struct with
%     .nBlocks   number of coarse blocks
%     .pv        pore volume per block [m^3] (sum of fine pore volumes;
%                fracture-cell volumes are folded into their host block)
%     .depth     pore-volume-weighted centroid depth per block [m]
%     .pairs     nPairs x 2 coarse connections (I < J)
%     .T         aggregated transmissibility per connection [SI], the
%                sum of the fine transmissibilities of all connections
%                crossing the block interface (exact for flow parallel
%                across the interface; the standard aggregation for
%                connection-list coarsening)
%     .partition copy of p
%
% NOTE: aggregation-based transmissibilities are appropriate for
% screening models; flow-based upscaling can replace .T without
% changing the export format.

p = p(:);
assert(numel(p) == G.cells.num, 'partition must cover all cells of G');
nBlocks = max(p);

pvFine = poreVolume(G, rock);
pv = accumarray(p, pvFine, [nBlocks, 1]);
depth = accumarray(p, pvFine .* G.cells.centroids(:,3), [nBlocks, 1]) ./ max(pv, eps);

% aggregate fine connections across block interfaces
I = p(ops.N(:,1));
J = p(ops.N(:,2));
cross = I ~= J;
lo = min(I(cross), J(cross));
hi = max(I(cross), J(cross));
[pairs, ~, ic] = unique([lo, hi], 'rows');
T = accumarray(ic, ops.T(cross));

% drop numerically dead connections (fully sealed interfaces)
keep = T > 0;
coarse = struct('nBlocks', nBlocks, 'pv', pv, 'depth', depth, ...
    'pairs', pairs(keep,:), 'T', T(keep), 'partition', p);
end
