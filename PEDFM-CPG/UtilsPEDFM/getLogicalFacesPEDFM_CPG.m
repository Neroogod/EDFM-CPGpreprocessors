function faceSets = getLogicalFacesPEDFM_CPG(G, c)
%GETLOGICALFACESPEDFM_CPG Group the faces of one cell by logical direction tag.
%
% SYNOPSIS:
%   faceSets = getLogicalFacesPEDFM_CPG(G, c)
%
% PARAMETERS:
%   G - grid whose G.cells.faces carries tags in column 2 (1=left/West,
%       2=right/East, 3=front/South, 4=back/North, 5=bottom, 6=top), as
%       produced by cartGrid, tensorGrid and processGRDECL.
%   c - single cell index.
%
% RETURNS:
%   faceSets - 1x6 cell array. faceSets{t} holds the global face indices of
%              cell c carrying tag t. On a Cartesian grid every set is a
%              singleton. On a faulted corner-point grid a set can hold
%              several faces (a split logical side) and can be empty
%              (pinch-out / degenerate side).
%
% This function exists so that no caller ever indexes cell faces
% positionally. The Cartesian pEDFM reference (pMatFracNNCs3D.m, lines
% 44-49) reads faceCentroids(1..6,:) as left/right/front/back/bottom/top;
% corner-point grids neither keep cell faces sorted by tag nor guarantee
% six of them. See pedfm_implementation.md section 2.1.

    fpos = G.cells.facePos(c):(G.cells.facePos(c+1) - 1);
    f = G.cells.faces(fpos, 1);

    assert(size(G.cells.faces, 2) >= 2, ...
        'getLogicalFacesPEDFM_CPG:noTags', ...
        ['G.cells.faces has no tag column; pEDFM face selection needs ', ...
         'logical directions. Grid type not supported.']);
    t = G.cells.faces(fpos, 2);

    faceSets = cell(1, 6);
    for tag = 1:6
        faceSets{tag} = f(t == tag);
    end
end
