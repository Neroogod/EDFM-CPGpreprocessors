function info = cmgDetectFractureModel(G)
% Determine which embedded-fracture method a grid carries, if any.
%
% SYNOPSIS:
%   info = cmgDetectFractureModel(G)
%
% RETURNS:
%   info - struct with
%     .model      'matrix' | 'EDFM' | 'pEDFM' | 'CPEDFM'
%     .hasFrac    true when G carries an assembled fracture grid
%     .nMatrix    number of matrix cells
%     .nFrac      number of fracture cells
%     .nNNC       number of non-neighbour connections
%     .nProjFace  face-indexed projections (CPEDFM contract)
%     .nProjPair  pair-indexed projections (standard pEDFM contract)
%     .reason     human-readable basis for the classification
%
% NOTE ON THE TEST USED:
%   Membership is decided on NON-EMPTY projection arrays, never on isfield.
%   prepareEDFMExport deliberately writes empty pMMfaces and pMMneighs onto
%   a normalized EDFM grid, so an isfield test would report that grid as
%   CPEDFM. This distinction is the whole reason this function exists.
%
% SEE ALSO:
%   cmgDetectGrid, buildCMGModel, writeCMG

info = struct('model', 'matrix', 'hasFrac', false, ...
    'nMatrix', G.cells.num, 'nFrac', 0, 'nNNC', 0, ...
    'nProjFace', 0, 'nProjPair', 0, 'reason', '');

hasFrac = isfield(G, 'FracGrid') && isfield(G, 'Matrix');
if ~hasFrac
    info.reason = 'no FracGrid/Matrix fields: plain matrix grid';
    return
end

info.hasFrac = true;
info.nMatrix = G.Matrix.cells.num;
info.nFrac = G.cells.num - G.Matrix.cells.num;

if ~isfield(G, 'nnc') || ~isfield(G.nnc, 'cells') || isempty(G.nnc.cells)
    % Fracture cells exist but were never connected to anything.
    info.reason = 'FracGrid present but G.nnc.cells is empty';
    info.model = 'EDFM';
    return
end

info.nNNC = size(G.nnc.cells, 1);
info.nProjFace = localCount(G.nnc, 'pMMfaces');
info.nProjPair = localCount(G.nnc, 'pMMneighs');

if info.nProjFace > 0
    info.model = 'CPEDFM';
    info.reason = sprintf('%d face-indexed projections (pMMfaces)', info.nProjFace);
elseif info.nProjPair > 0
    info.model = 'pEDFM';
    info.reason = sprintf('%d pair-indexed projections (pMMneighs)', info.nProjPair);
else
    info.model = 'EDFM';
    info.reason = 'NNCs present, no projection records';
end

% An explicit stamp from a normalization adapter wins over inference: it
% records what the caller asked for, not just what the arrays look like.
if isfield(G.nnc, 'exportModelType') && ~isempty(G.nnc.exportModelType)
    stamped = char(string(G.nnc.exportModelType));
    if ~strcmpi(stamped, info.model)
        info.reason = sprintf('%s; overridden by G.nnc.exportModelType = %s', ...
            info.reason, stamped);
        info.model = stamped;
    end
end
end

%--------------------------------------------------------------------------
function n = localCount(nnc, field)
n = 0;
if isfield(nnc, field) && ~isempty(nnc.(field))
    n = size(nnc.(field), 1);
end
end
