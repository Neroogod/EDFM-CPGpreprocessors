function [xarea, polyOut, onFace] = fracPolyCellIntersect(poly, G, c, tol)
% Intersection of a planar fracture polygon with a convex grid cell.
%
% SYNOPSIS:
%   [xarea, polyOut, onFace] = fracPolyCellIntersect(poly, G, c, tol)
%
% DEPRECATED - kept so existing scripts keep running. The implementation
% now lives in EDFM-U/UtilsEDFM/clipPolyByCell_U, which this simply calls;
% prefer that name in new code.
%
% WHY IT MOVED
%   The version that used to live here returned a WRONG polygon whenever
%   the fracture was coplanar with one of the cell's faces. On detecting
%   the coplanar face it did `break`, abandoning the clip loop, so the
%   polygon was only ever trimmed by the faces that happened to be visited
%   BEFORE that face - the remaining faces never trimmed it at all.
%
%   With a fracture on the x = 20 plane of a 4x4x3 grid over 40x40x30, every
%   one of the 24 cells touching the plane was handed back the ENTIRE
%   936 m^2 fracture polygon instead of its own 64-100 m^2 face footprint:
%   22464 m^2 of intersection area in total against a true 1872 m^2, a
%   factor of 12. Those identical, overlapping polygons then made the
%   fracture grid degenerate and EDFMGridU threw in
%   fractureplanegeneralgrid, so the face-coplanar case did not merely lose
%   accuracy - it did not run.
%
%   clipPolyByCell_U records the coplanar face and CONTINUES the loop.
%   Clipping against a coplanar face is a no-op (every vertex has
%   |d| < tol, so all are kept), so carrying on costs nothing and leaves
%   polyOut equal to (fracture polygon) INTERSECT (that face) - the correct
%   boundary footprint.
%
% SEE ALSO:
%   clipPolyByCell_U, EDFMGridU

[xarea, polyOut, onFace] = clipPolyByCell_U(poly, G, c, tol);
end
