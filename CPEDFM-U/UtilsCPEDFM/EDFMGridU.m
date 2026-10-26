function [G, fracplanes] = EDFMGridU(G, fracplanes, varargin)
% EDFM fracture gridding on general unstructured grids.
%
% SYNOPSIS:
%   [G, fracplanes] = EDFMGridU(G, fracplanes, 'Tolerance', tol, ...)
%
% DEPRECATED - kept so existing scripts keep running. The implementation
% now lives in EDFM-U/UtilsEDFM/EDFMgrid_U, which this simply calls; prefer
% that name in new code. Options are unchanged ('Tolerance', 'fracturelist',
% 'Verbose').
%
% WHY IT MOVED - the version that used to live here had two defects in the
% face-coplanar ('fracmat boundary') path, both fixed in EDFMgrid_U:
%
%   1. It clipped with the old fracPolyCellIntersect, which abandoned its
%      clip loop on hitting the coplanar face and returned an untrimmed
%      polygon - see that file for the reproduction. This alone made the
%      coplanar case throw in fractureplanegeneralgrid.
%
%   2. It registered a face-coplanar fracture cell against only ONE of the
%      two matrix cells that border it:
%
%          if Gm.faces.neighbors(onFace,1) ~= c, continue; end
%
%      The cell on the far side was left with no connection to the fracture
%      at all. Worse, the single surviving row then fooled the surface-area
%      correction downstream: adjustT infers the factor by counting how
%      many rows share a fracture cell, saw one, concluded "fracture inside
%      a cell", and applied the interior x2 rule to what is actually a
%      boundary fracture. On a 4x4x3 grid with a fracture on an interior
%      face plane that gave 12 connections instead of 24, with 12 matrix
%      cells east of the fracture entirely disconnected from it and every
%      surviving transmissibility exactly 2x too large. The aggregate
%      transfer capacity happened to match, so the error was topological:
%      the fracture only talked to one side.
%
%      EDFMgrid_U creates one fracture cell per coplanar face, registers it
%      against BOTH neighbours, and records the fracture-cell index of each
%      connection explicitly in matrix_connection.fraccells, so the
%      surface-area factor can be taken from the recorded intersection type
%      instead of by counting rows.
%
% Because of (2) the returned matrix_connection now has one row per
% CONNECTION rather than per fracture cell; for an all-interior fracture
% (the common case) the two are identical, so existing code that reads
% matrix_connection.cells/area/type is unaffected.
%
% NOTE - pair this with fracturematrixNNC_U, not fracturematrixShaleNNC3D:
% the latter routes through hfm/calcdavg, which integrates the average
% fracture-matrix distance <d> over the cell's AXIS-ALIGNED BOUNDING BOX.
% On the PEBI and tetrahedral grids this module targets the box is not the
% cell, and since CI = A/<d> that error passes straight into every
% fracture-matrix transmissibility (~16% median on a PEBI grid, ~44% on
% tetrahedra). fracturematrixNNC_U uses calcDavg_U, which integrates over
% the true cell polyhedron exactly.
%
% SEE ALSO:
%   EDFMgrid_U, fracturematrixNNC_U, pMatFracNNCsCPEDFM_U

[G, fracplanes] = EDFMgrid_U(G, fracplanes, varargin{:});
end
