function A = cellNeighbourhood_U(G)
% Node-sharing cell adjacency, for any grid.
%
% SYNOPSIS:
%   A = cellNeighbourhood_U(G)
%
% RETURNS:
%   A - G.cells.num x G.cells.num sparse logical matrix; A(a,b) is true
%       when cells a and b share at least one node (the diagonal is true).
%       find(A(:,c)) is then the cloud of cells touching cell c.
%
% This is the unstructured generalisation of the 26-cell 3x3x3 cloud that
% hfm/findneighbours builds on a Cartesian grid (that routine asserts
% isfield(G,'cartDims') and cannot run here). Sharing a node is the right
% criterion for the fracture-fracture candidate filter: a fracture cell is
% confined to its host matrix cell, so two fracture cells can only touch if
% their host cells touch.

nc = G.cells.num;
[cn, cpos] = gridCellNodes(G, 1:nc);
cellNo = rldecode((1:nc)', diff(cpos));

M = sparse(cn, cellNo, 1, G.nodes.num, nc);   % node-cell incidence
A = (M' * M) > 0;
end
