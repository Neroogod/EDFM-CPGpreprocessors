function v = deckBlockArray(res, cellValues, fillValue)
%DECKBLOCKARRAY Scatter per-cell MRST values onto the exported block grid.
%
% SYNOPSIS:
%   v = deckBlockArray(res, cellValues, fillValue)
%
% PARAMETERS:
%   res        - reservoir description from buildCMGModel.
%   cellValues - one value per MRST cell of the global grid, matrix cells
%                first and fracture cells after, which is how MRST numbers a
%                grid assembled by EDFMshalegrid / EDFMGridCPG / EDFMGridU.
%   fillValue  - value given to blocks that carry no MRST cell (inactive
%                slots and unused padding).
%
% RETURNS:
%   v - nTot x 1 array in the exported grid's block order, ready to be
%       written as PRESSURE / SWAT / *PRES / *SW.
%
% WHY THIS EXISTS
%   On a fracture export buildCMGModel does not simply append the fracture
%   cells: it places each one in the padding layers underneath the column of
%   the matrix cell it intersects, leaving the unused slots inactive. For the
%   8x8x2 example in the tests that turns 148 MRST cells into a 8x8x4 = 256
%   block grid with 108 inactive blocks, and the fracture cells land at
%   scattered addresses rather than at blocks 129..148.
%
%   Writing an initial state by appending values to the end of the matrix
%   array would therefore put fracture pressures and saturations on the wrong
%   blocks - silently, since the deck would still be well formed and still
%   run. buildCMGModel records the addresses it used in grid.matIJK and
%   grid.fracIJK; this function is the only place that mapping is applied.
%
% SEE ALSO: buildFullDeck, writeCMGDeck, writeEclipseDeck.

    nTot = res.props.nTot;
    dims = reshape(double(res.grid.dims), 1, []);
    cellValues = double(cellValues(:));

    v = repmat(double(fillValue), nTot, 1);

    hasMat  = isfield(res.grid, 'matIJK')  && ~isempty(res.grid.matIJK);
    hasFrac = isfield(res.grid, 'fracIJK') && ~isempty(res.grid.fracIJK);

    if ~hasMat
        % Matrix-only export on a grid whose block order is MRST's own: the
        % cell index is the block index. Corner-point grids with inactive
        % cells are covered by the matIJK branch above.
        n = min(numel(cellValues), nTot);
        v(1:n) = cellValues(1:n);
        return
    end

    matIJK = double(res.grid.matIJK);
    nMat = size(matIJK, 1);
    nUse = min(nMat, numel(cellValues));
    linM = sub2ind(dims, matIJK(1:nUse,1), matIJK(1:nUse,2), matIJK(1:nUse,3));
    v(linM) = cellValues(1:nUse);

    if hasFrac && numel(cellValues) > nMat
        fracIJK = double(res.grid.fracIJK);
        nFrac = size(fracIJK, 1);
        avail = numel(cellValues) - nMat;
        nUseF = min(nFrac, avail);
        if nUseF > 0
            linF = sub2ind(dims, fracIJK(1:nUseF,1), fracIJK(1:nUseF,2), ...
                                 fracIJK(1:nUseF,3));
            v(linF) = cellValues(nMat + (1:nUseF));
        end
    end
end
