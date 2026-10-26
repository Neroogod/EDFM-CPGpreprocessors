function info = writeEclipseCPEDFM(G, rock, grdecl, wells, outFile, varargin)
% Export a CPEDFM/pEDFM model to a self-contained Eclipse-format deck.
%
% SYNOPSIS:
%   info = writeEclipseCPEDFM(G, rock, grdecl, wells, 'CASE.DATA')
%
% PARAMETERS:
%   G       - global EDFM grid AFTER the projection stage
%             (pMatFracNNCsCPEDFM_CPG or pMatFracNNCsCPEDFM_U): must carry
%             G.Matrix, G.FracGrid, and G.nnc with cells/T/pMMfaces/pMMarea.
%   rock    - rock of the global grid (matrix + fracture cells).
%   grdecl  - corner-point description of the MATRIX grid (cartDims,
%             COORD, ZCORN, optional ACTNUM). For Cartesian grids, pass
%             the output of simpleGrdecl/cartGrid-derived grdecl.
%   wells   - struct array: .name, .type ('rate'|'bhp'), .val (SI),
%             .comp ('WATER'|'OIL'), .cells (matrix cell indices).
%   outFile - path of the .DATA file to write.
%
% OPTIONAL ('key', value):
%   'tol'    - tolerance for the multiplier assembly (default 1e-5)
%   'p0'     - initial pressure [Pa] (default 200 bar)
%   'sw0'    - initial water saturation (default 0)
%   'dt'     - schedule step lengths [s] (default 10 x 30 days)
%   'title'  - deck title string
%
% RETURNS:
%   info - struct with the fracture-cell address map (.fracIJK), extended
%          dimensions (.dims), and counts of NNC / MULT / EDITNNC records.
%
% METHOD (the non-intrusive EDFM representation):
%   * Fracture cells are appended as extra k-layers of padding cells below
%     the matrix grid. Their geometry is placeholder; their physics is set
%     through EDIT overrides: PORV (true fracture pore volume) and DEPTH
%     (true fracture-cell centroid depth, for correct gravity head).
%   * All geometric transmissibilities of the padding region, and the
%     interface between the matrix and the padding, are zeroed, so the
%     fracture cells communicate exclusively through explicit NNCs.
%   * Every EDFM/CPEDFM connection (matrix-fracture, fracture-fracture,
%     and projection-matrix-fracture) becomes one NNC record with the
%     transmissibility computed by the MRST pipeline, converted to
%     METRIC units (cP.m^3/day/bar).
%   * The CPEDFM matrix-matrix seals become MULTX/MULTY/MULTZ records for
%     faces between logical neighbors, and EDITNNC records for
%     fault-juxtaposed pairs (which Eclipse connects through internally
%     generated NNCs).
%
% The deck is verified by tests/test20_eclipse_roundtrip.m, which reads
% it back with MRST's deckformat parser and reproduces the CPEDFM
% operators from the deck data alone.

opt = struct('tol', 1e-5, 'p0', 200*barsa, 'sw0', 0.0, ...
             'dt', repmat(30*day, 10, 1), ...
             'title', 'CPEDFM_EXPORT');
opt = merge_options(opt, varargin{:});

assert(all(isfield(G, {'Matrix', 'FracGrid'})), 'G must be a global EDFM grid.');
assert(isfield(G.nnc, 'cells') && ~isempty(G.nnc.cells), 'G.nnc is empty.');

nMat   = G.Matrix.cells.num;
nFrac  = G.cells.num - nMat;
NX = grdecl.cartDims(1); NY = grdecl.cartDims(2); NZ = grdecl.cartDims(3);

% METRIC transmissibility unit
transUnit = centi*poise*meter^3/(day*barsa);

%% ---- addressing -------------------------------------------------------
% matrix cells -> (i,j,k) in the original grid
[im, jm, km] = ind2sub([NX NY NZ], G.Matrix.cells.indexMap);
matIJK = [im(:), jm(:), km(:)];

% fracture cells -> slots in the padding layers
nzExtra = ceil(nFrac/(NX*NY));
dims = [NX, NY, NZ + nzExtra];
slot = (1:nFrac)';
kf = NZ + ceil(slot/(NX*NY));
r  = slot - (kf - NZ - 1)*NX*NY;
jf = ceil(r/NX);
if_ = r - (jf-1)*NX;
fracIJK = [if_, jf, kf];

ijkOf = @(c) localIJK(c, nMat, matIJK, fracIJK);

%% ---- extended geometry -------------------------------------------------
zcorn = reshape(grdecl.ZCORN, 2*NX, 2*NY, 2*NZ);
zbot  = zcorn(:,:,end);          % bottom surface of the matrix grid
dzPad = 1.0;                     % placeholder thickness [m]
zcornExt = zeros(2*NX, 2*NY, 2*dims(3));
zcornExt(:,:,1:2*NZ) = zcorn;
for L = 1:nzExtra
    zcornExt(:,:,2*NZ + 2*L - 1) = zbot + (L-1)*dzPad;
    zcornExt(:,:,2*NZ + 2*L    ) = zbot +  L   *dzPad;
end

% ACTNUM: matrix actnum + active fracture slots
if isfield(grdecl, 'ACTNUM')
    actMat = grdecl.ACTNUM(:);
else
    actMat = ones(NX*NY*NZ, 1);
end
actPad = zeros(NX*NY*nzExtra, 1);
actPad(1:nFrac) = 1;
actnum = [actMat; actPad];

%% ---- property arrays over the extended grid ---------------------------
nTotExt = prod(dims);
permx = zeros(nTotExt,1); permy = permx; permz = permx; poro = permx;

matLin = sub2ind(dims, matIJK(:,1), matIJK(:,2), matIJK(:,3));
pm = rock.perm(1:nMat, :);
if size(pm,2) == 1, pm = repmat(pm, 1, 3); end
permx(matLin) = pm(:,1); permy(matLin) = pm(:,2); permz(matLin) = pm(:,min(3,end));
poro(matLin)  = rock.poro(1:nMat);

fracLin = sub2ind(dims, fracIJK(:,1), fracIJK(:,2), fracIJK(:,3));
pf = rock.perm(nMat+1:end, :);
if size(pf,2) == 1, pf = repmat(pf, 1, 3); end
permx(fracLin) = pf(:,1); permy(fracLin) = pf(:,min(2,end)); permz(fracLin) = pf(:,min(3,end));
poro(fracLin)  = max(rock.poro(nMat+1:end), 1e-6);

%% ---- NNC records -------------------------------------------------------
nncCells = G.nnc.cells;
nncT     = convertTo(G.nnc.T, transUnit);
nNNC = size(nncCells, 1);
nncIJK = zeros(nNNC, 6);
for q = 1:nNNC
    nncIJK(q, 1:3) = ijkOf(nncCells(q,1));
    nncIJK(q, 4:6) = ijkOf(nncCells(q,2));
end

%% ---- matrix-matrix seals: MULT records and EDITNNC records -------------
mult = getEmbeddedFractureMultipliers(G, opt.tol);
sealedFaces = find(mult < 1 - 1e-12);
sealedFaces = sealedFaces(sealedFaces <= G.Matrix.faces.num);

multRec = zeros(0, 5);   % [dir iLow jLow kLow value] dir: 1=X,2=Y,3=Z
editnncRec = zeros(0, 7);% [i1 j1 k1 i2 j2 k2 value]
for q = 1:numel(sealedFaces)
    f = sealedFaces(q);
    nb = G.Matrix.faces.neighbors(f, :);
    if any(nb == 0), continue; end
    a = matIJK(nb(1), :); b = matIJK(nb(2), :);
    d = b - a;
    if isequal(abs(d), [1 0 0]) || isequal(abs(d), [0 1 0]) || isequal(abs(d), [0 0 1])
        % logical neighbors: MULTX/Y/Z on the lower-index cell
        dir = find(d ~= 0);
        if d(dir) > 0, low = a; else, low = b; end
        multRec(end+1, :) = [dir, low, mult(f)]; %#ok<AGROW>
    else
        % fault-juxtaposed pair: Eclipse connects these through an
        % internally generated NNC -> scale it with EDITNNC
        editnncRec(end+1, :) = [a, b, mult(f)]; %#ok<AGROW>
    end
end

%% ---- fracture-cell overrides ------------------------------------------
pvG = poreVolume(G, rock);
fracPV    = pvG(nMat+1:end);
fracDepth = G.cells.centroids(nMat+1:end, 3);

%% ---- write the deck ----------------------------------------------------
fid = fopen(outFile, 'w');
assert(fid > 0, 'Could not open %s for writing.', outFile);
c = onCleanup(@() fclose(fid));
w = @(varargin) fprintf(fid, varargin{:});

w('-- Generated by writeEclipseCPEDFM (CPEDFM-export module)\n');
w('-- Fracture cells: %d (padding layers %d-%d); NNCs: %d; MULT seals: %d; EDITNNC seals: %d\n\n', ...
    nFrac, NZ+1, dims(3), nNNC, size(multRec,1), size(editnncRec,1));

% ------------------------------------------------ RUNSPEC
w('RUNSPEC\n\nTITLE\n %s /\n\n', opt.title);
w('DIMENS\n %d %d %d /\n\n', dims);
w('OIL\nWATER\n\nMETRIC\n\n');
w('TABDIMS\n 1 1 20 20 /\n\n');
w('WELLDIMS\n %d %d 2 %d /\n\n', numel(wells), ...
    max(arrayfun(@(x) numel(x.cells), wells)), numel(wells));
w('START\n 1 JAN 2026 /\n\nUNIFOUT\n\n');

% ------------------------------------------------ GRID
w('GRID\n\n');
w('SPECGRID\n %d %d %d 1 F /\n\n', dims);
writeVector(fid, 'COORD', grdecl.COORD, '%.6f');
writeVector(fid, 'ZCORN', zcornExt(:), '%.6f');
writeVector(fid, 'ACTNUM', actnum, '%d');
writeVector(fid, 'PERMX', convertTo(permx, milli*darcy), '%.6e');
writeVector(fid, 'PERMY', convertTo(permy, milli*darcy), '%.6e');
writeVector(fid, 'PERMZ', convertTo(permz, milli*darcy), '%.6e');
writeVector(fid, 'PORO', poro, '%.6f');

% CPEDFM seals between logical neighbors
if ~isempty(multRec)
    w('-- CPEDFM matrix-matrix seals (logical-neighbor faces)\n');
    w('EQUALS\n');
    dirName = {'MULTX', 'MULTY', 'MULTZ'};
    for q = 1:size(multRec, 1)
        w(' ''%s'' %.8e %d %d %d %d %d %d /\n', dirName{multRec(q,1)}, ...
            multRec(q,5), multRec(q,2), multRec(q,2), ...
            multRec(q,3), multRec(q,3), multRec(q,4), multRec(q,4));
    end
    w('/\n\n');
end

% EDFM / CPEDFM non-neighbor connections
w('-- EDFM + CPEDFM projection NNCs (trans in cP.m3/day/bar)\n');
w('NNC\n');
for q = 1:nNNC
    w(' %d %d %d  %d %d %d  %.8e /\n', nncIJK(q,:), nncT(q));
end
w('/\n\n');

% ------------------------------------------------ EDIT
w('EDIT\n\n');
w('-- isolate the fracture padding region: geometry plays no role there\n');
w('BOX\n 1 %d 1 %d %d %d /\n', NX, NY, NZ, NZ);
w('TRANZ\n %d*0.0 /\n', NX*NY);
w('ENDBOX\n\n');
w('BOX\n 1 %d 1 %d %d %d /\n', NX, NY, NZ+1, dims(3));
w('TRANX\n %d*0.0 /\n', NX*NY*nzExtra);
w('TRANY\n %d*0.0 /\n', NX*NY*nzExtra);
w('TRANZ\n %d*0.0 /\n', NX*NY*nzExtra);
w('ENDBOX\n\n');

w('-- fracture-cell pore volumes [m3] and centroid depths [m],\n');
w('-- written as BOX arrays over the padding region (slot order)\n');
padPV = zeros(NX*NY*nzExtra, 1);  padPV(1:nFrac) = fracPV;
padZ  = ones(NX*NY*nzExtra, 1);   padZ(1:nFrac)  = fracDepth;
w('BOX\n 1 %d 1 %d %d %d /\n', NX, NY, NZ+1, dims(3));
writeVector(fid, 'PORV', padPV, '%.8e');
writeVector(fid, 'DEPTH', padZ, '%.6f');
w('ENDBOX\n\n');

if ~isempty(editnncRec)
    w('-- CPEDFM seals on fault-juxtaposed pairs\n');
    w('EDITNNC\n');
    for q = 1:size(editnncRec, 1)
        w(' %d %d %d  %d %d %d  %.8e /\n', editnncRec(q,1:6), editnncRec(q,7));
    end
    w('/\n\n');
end

% ------------------------------------------------ PROPS
w('PROPS\n\n');
w('SWOF\n');
sw = linspace(0, 1, 11)';
for q = 1:numel(sw)
    w(' %.4f %.6f %.6f 0.0\n', sw(q), sw(q)^2, (1-sw(q))^2);
end
w('/\n\n');
w('PVTW\n 200.0 1.0 4.0E-5 0.5 0.0 /\n\n');
w('PVCDO\n 200.0 1.0 1.0E-5 5.0 0.0 /\n\n');
w('ROCK\n 200.0 4.35E-5 /\n\n');
w('DENSITY\n 800.0 1000.0 1.2 /\n\n');

% ------------------------------------------------ SOLUTION
% Equilibrium initialization at the datum pressure; the oil-water
% contact is placed below the model so the reservoir initializes at
% the connate water saturation (parser- and simulator-friendly, unlike
% explicit PRESSURE/SWAT arrays whose expected lengths differ between
% readers).
w('SOLUTION\n\n');
zDatum = mean(G.cells.centroids(1:nMat, 3));
zDeep  = max(zcornExt(:)) + 1000;
w('EQUIL\n %.4f %.4f %.4f 0.0 1* 0.0 1* 1* 0 /\n\n', ...
    zDatum, convertTo(opt.p0, barsa), zDeep);

% ------------------------------------------------ SUMMARY
w('SUMMARY\n\nFPR\nWBHP\n/\nWOPR\n/\nWWIR\n/\n\n');

% ------------------------------------------------ SCHEDULE
w('SCHEDULE\n\n');
w('WELSPECS\n');
for q = 1:numel(wells)
    hijk = ijkOf(wells(q).cells(1));
    w(' ''%s'' ''G'' %d %d 1* ''%s'' /\n', wells(q).name, hijk(1), hijk(2), ...
        upper(wells(q).comp));
end
w('/\n\n');
w('COMPDAT\n');
for q = 1:numel(wells)
    for cc = wells(q).cells(:)'
        hijk = ijkOf(cc);
        w(' ''%s'' %d %d %d %d ''OPEN'' 2* 0.2 /\n', wells(q).name, ...
            hijk(1), hijk(2), hijk(3), hijk(3));
    end
end
w('/\n\n');
for q = 1:numel(wells)
    if strcmpi(wells(q).type, 'rate') && strcmpi(wells(q).comp, 'WATER')
        w('WCONINJE\n ''%s'' ''WATER'' ''OPEN'' ''RATE'' %.6e 1* 1000.0 /\n/\n\n', ...
            wells(q).name, convertTo(wells(q).val, meter^3/day));
    elseif strcmpi(wells(q).type, 'bhp')
        w('WCONPROD\n ''%s'' ''OPEN'' ''BHP'' 5* %.4f /\n/\n\n', ...
            wells(q).name, convertTo(wells(q).val, barsa));
    end
end
w('TSTEP\n');
w(' %.4f', convertTo(opt.dt, day));
w(' /\n\nEND\n');

info = struct('dims', dims, 'fracIJK', fracIJK, 'matIJK', matIJK, ...
    'nNNC', nNNC, 'nMult', size(multRec,1), 'nEditnnc', size(editnncRec,1), ...
    'multRec', multRec, 'editnncRec', editnncRec, ...
    'nncIJK', nncIJK, 'nncT', nncT, 'fracPV', fracPV, 'fracDepth', fracDepth);
end

%% -------------------------------------------------------------------------
function ijk = localIJK(c, nMat, matIJK, fracIJK)
if c <= nMat
    ijk = matIJK(c, :);
else
    ijk = fracIJK(c - nMat, :);
end
end

function writeVector(fid, kw, v, fmt)
fprintf(fid, '%s\n', kw);
fprintf(fid, [' ' fmt ' ' fmt ' ' fmt ' ' fmt ' ' fmt ' ' fmt '\n'], ...
    v(1:6*floor(numel(v)/6)));
rem = v(6*floor(numel(v)/6)+1:end);
if ~isempty(rem)
    fprintf(fid, [' ' fmt], rem);
    fprintf(fid, '\n');
end
fprintf(fid, '/\n\n');
end
