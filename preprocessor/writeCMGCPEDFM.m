function info = writeCMGCPEDFM(G, rock, grdecl, outFile, varargin)
% Export CPEDFM/pEDFM data as CMG-style include files.
%
% SYNOPSIS:
%   info = writeCMGCPEDFM(G, rock, grdecl, 'cpedfm_cmg.inc')
%
% STATUS: the DATA CONTENT written by this function (fracture-cell
% addresses, pore volumes, depths, NNC transmissibilities, and seal
% multipliers) is identical to what the verified Eclipse exporter
% writes (see writeEclipseCPEDFM and test20). The CMG KEYWORD SYNTAX,
% however, differs between CMG simulators (IMEX/GEM/STARS) and
% versions, and could not be machine-verified here. Check the emitted
% cards against the keyword manual of your CMG version before use, or
% import the tables through Builder. The mapping is:
%
%   MRST quantity                    -> CMG concept
%   ------------------------------------------------------------------
%   fracture cells (padding layers)  -> extra k-layers of blocks; make
%                                       unused slots *NULL
%   fracture pore volume [m3]        -> pore-volume override of the
%                                       padding blocks (*VOLMOD or
%                                       porosity override)
%   fracture depth [m]               -> block depth (via grid geometry
%                                       or depth-modification cards)
%   NNC list (pairs + trans)         -> explicit non-neighbor
%                                       connections; GEM/IMEX accept
%                                       user NNCs in recent versions
%                                       (trans written in md.m,
%                                       CMG's transmissibility unit
%                                       for SI dataset)
%   seal multipliers (per face)      -> *TRANSI / *TRANSJ / *TRANSK
%                                       block multipliers via *MOD
%                                       (apply to the +direction face
%                                       of the addressed block)
%   fault-pair seals                 -> multipliers on the fault
%                                       connections (Builder's
%                                       faults/NNC facility)
%
% The function writes one include file with clearly banner-separated
% sections in the order above, using *MOD-card syntax for the
% multipliers and comment-tagged tables for the NNCs and overrides.
%
% PARAMETERS: as in writeEclipseCPEDFM (wells/schedule are left to the
% host dataset).
%
% OPTIONS:
%   simulator   - 'IMEX' (default), 'GEM', or 'STARS'.
%   cmgVersion  - release label recorded in the output manifest.
%   writeTables - emit simulator-tagged CSV tables and JSON manifest
%                 alongside the include file (default true).
%
% RETURNS:
%   info - same bookkeeping struct as writeEclipseCPEDFM (addresses,
%          counts, records).

opt = struct('tol', 1e-5, 'simulator', 'IMEX', ...
    'cmgVersion', 'unspecified', 'writeTables', true);
opt = merge_options(opt, varargin{:});
profile = getCMGSimulatorProfile(opt.simulator, opt.cmgVersion);
modelType = 'CPEDFM';
if isfield(G, 'nnc') && isfield(G.nnc, 'exportModelType')
    modelType = G.nnc.exportModelType;
end

assert(all(isfield(G, {'Matrix', 'FracGrid'})), 'G must be a global EDFM grid.');
nMat  = G.Matrix.cells.num;
nFrac = G.cells.num - nMat;
NX = grdecl.cartDims(1); NY = grdecl.cartDims(2); NZ = grdecl.cartDims(3);

% CMG SI transmissibility unit: md.m
transUnit = milli*darcy*meter;

% ---- addressing (same layout as the Eclipse exporter) ----
[im, jm, km] = ind2sub([NX NY NZ], G.Matrix.cells.indexMap);
matIJK = [im(:), jm(:), km(:)];
nzExtra = ceil(nFrac/(NX*NY));
dims = [NX, NY, NZ + nzExtra];
slot = (1:nFrac)';
kf = NZ + ceil(slot/(NX*NY));
r  = slot - (kf - NZ - 1)*NX*NY;
jf = ceil(r/NX);
if_ = r - (jf-1)*NX;
fracIJK = [if_, jf, kf];
ijkOf = @(c) localIJK(c, nMat, matIJK, fracIJK);

% ---- NNCs ----
nNNC = size(G.nnc.cells, 1);
nncIJK = zeros(nNNC, 6);
for q = 1:nNNC
    nncIJK(q,1:3) = ijkOf(G.nnc.cells(q,1));
    nncIJK(q,4:6) = ijkOf(G.nnc.cells(q,2));
end
nncT = convertTo(G.nnc.T, transUnit);

% ---- seals ----
mult = getEmbeddedFractureMultipliers(G, opt.tol);
sealedFaces = find(mult < 1 - 1e-12);
sealedFaces = sealedFaces(sealedFaces <= G.Matrix.faces.num);
multRec = zeros(0,5); faultRec = zeros(0,7);
for q = 1:numel(sealedFaces)
    f = sealedFaces(q);
    nb = G.Matrix.faces.neighbors(f,:);
    if any(nb == 0), continue; end
    a = matIJK(nb(1),:); b = matIJK(nb(2),:);
    d = b - a;
    if nnz(d) == 1 && max(abs(d)) == 1
        dir = find(d ~= 0);
        if d(dir) > 0, low = a; else, low = b; end
        multRec(end+1,:) = [dir, low, mult(f)]; %#ok<AGROW>
    else
        faultRec(end+1,:) = [a, b, mult(f)];    %#ok<AGROW>
    end
end

% ---- overrides ----
pvG = poreVolume(G, rock);
fracPV = pvG(nMat+1:end);
fracDepth = G.cells.centroids(nMat+1:end, 3);
pf = rock.perm(nMat+1:end, :);
if size(pf,2) == 1, pf = repmat(pf,1,3); end

% ---- write ----
fid = fopen(outFile, 'w');
assert(fid > 0, 'Could not open %s.', outFile);
c = onCleanup(@() fclose(fid));
w = @(varargin) fprintf(fid, varargin{:});

w('** =====================================================================\n');
w('** CPEDFM export for CMG %s (generated by writeCMGCPEDFM)\n', profile.simulator);
w('** TARGET_SIMULATOR: %s\n', profile.simulator);
w('** TARGET_VERSION: %s\n', profile.version);
w('** FRACTURE_MODEL: %s\n', modelType);
w('** PHYSICS_PROFILE: %s\n', profile.physics);
w('** HOST_MODEL: %s\n', profile.hostModel);
w('** %s\n', profile.phaseNote);
w('** %s\n', profile.propertyPlacement);
w('** %s\n', profile.connectionPlacement);
w('** DATA verified against MRST; KEYWORD SYNTAX to be checked against the\n');
w('** manual of the selected %s release. See the mapping table in the\n', profile.simulator);
w('** function header. Grid must be extended to %d x %d x %d, with layers\n', dims);
w('** %d-%d as fracture padding (unused slots NULL).\n', NZ+1, dims(3));
w('** =====================================================================\n\n');

w('** ---- 1. Null out unused padding slots -------------------------------\n');
w('*NULL *IJK\n');
w('** used fracture slots stay live; everything else in the padding is null\n');
used = false(NX, NY, nzExtra);
used(sub2ind([NX NY nzExtra], fracIJK(:,1), fracIJK(:,2), fracIJK(:,3)-NZ)) = true;
for L = 1:nzExtra
    for jj = 1:NY
        row = squeeze(used(:, jj, L));
        idx = find(~row)';
        for ii = idx
            w(' %d %d %d 0\n', ii, jj, NZ+L);
        end
    end
end
w('\n');

w('** ---- 2. Fracture block properties -----------------------------------\n');
w('** columns: i j k  PV[m3]  depth[m]  permI[mD] permJ[mD] permK[mD]\n');
w('** map PV via pore-volume override, depth via depth card of your version\n');
for q = 1:nFrac
    w('** FRACBLK %d %d %d  %.8e  %.4f  %.6e %.6e %.6e\n', fracIJK(q,:), ...
        fracPV(q), fracDepth(q), convertTo(pf(q,1), milli*darcy), ...
        convertTo(pf(q,min(2,end)), milli*darcy), convertTo(pf(q,min(3,end)), milli*darcy));
end
w('\n');

w('** ---- 3. Seal multipliers on logical-neighbor faces ------------------\n');
w('** *TRANSI/J/K multiply the +direction face of the addressed block\n');
dirKw = {'*TRANSI', '*TRANSJ', '*TRANSK'};
for d = 1:3
    rows = multRec(multRec(:,1) == d, :);
    if isempty(rows), continue; end
    w('%s *MOD\n', dirKw{d});
    for q = 1:size(rows,1)
        w(' %d %d %d = %.8e\n', rows(q,2:4), rows(q,5));
    end
    w('\n');
end

if ~isempty(faultRec)
    w('** ---- 3b. Seals on fault-juxtaposed pairs (apply to fault NNCs) ---\n');
    w('** columns: i1 j1 k1  i2 j2 k2  multiplier\n');
    for q = 1:size(faultRec,1)
        w('** FAULTSEAL %d %d %d  %d %d %d  %.8e\n', faultRec(q,1:6), faultRec(q,7));
    end
    w('\n');
end

w('** ---- 4. Non-neighbor connections ------------------------------------\n');
w('** columns: i1 j1 k1  i2 j2 k2  trans[%s]\n', profile.transmissibilityUnit);
w('** %s\n', profile.connectionPlacement);
for q = 1:nNNC
    w('** NNC %d %d %d  %d %d %d  %.8e\n', nncIJK(q,:), nncT(q));
end

info = struct('dims', dims, 'fracIJK', fracIJK, 'nNNC', nNNC, ...
    'nMult', size(multRec,1), 'nFault', size(faultRec,1), ...
    'multRec', multRec, 'faultRec', faultRec, 'nncIJK', nncIJK, 'nncT', nncT, ...
    'simulator', profile.simulator, 'cmgVersion', profile.version, ...
    'modelType', modelType, 'profile', profile, 'tableFiles', struct);

if opt.writeTables
    info.tableFiles = writeProfileTables(outFile, profile, modelType, dims, NZ, ...
        fracIJK, fracPV, fracDepth, pf, nncIJK, nncT, multRec, faultRec);
end
end

function ijk = localIJK(c, nMat, matIJK, fracIJK)
if c <= nMat
    ijk = matIJK(c, :);
else
    ijk = fracIJK(c - nMat, :);
end
end

function files = writeProfileTables(outFile, profile, modelType, dims, matrixLayers, ...
    fracIJK, fracPV, fracDepth, pf, nncIJK, nncT, multRec, faultRec)
[outDir, baseName] = fileparts(outFile);
if isempty(outDir), outDir = pwd; end
prefix = fullfile(outDir, [baseName, '_', lower(profile.simulator)]);

files = struct;
files.fractureBlocks = [prefix, '_fracture_blocks.csv'];
files.nnc = [prefix, '_nnc.csv'];
files.multipliers = [prefix, '_multipliers.csv'];
files.faultSeals = [prefix, '_fault_seals.csv'];
files.manifest = [prefix, '_manifest.json'];

fracTable = [fracIJK, fracPV, fracDepth, ...
    convertTo(pf(:,1), milli*darcy), ...
    convertTo(pf(:,min(2,end)), milli*darcy), ...
    convertTo(pf(:,min(3,end)), milli*darcy)];
writeNumericCSV(files.fractureBlocks, ...
    'i,j,k,pore_volume_m3,depth_m,perm_i_mD,perm_j_mD,perm_k_mD', fracTable);
writeNumericCSV(files.nnc, ...
    'i1,j1,k1,i2,j2,k2,transmissibility_md_m', [nncIJK, nncT]);
writeNumericCSV(files.multipliers, ...
    'direction_1I_2J_3K,i,j,k,multiplier', multRec);
writeNumericCSV(files.faultSeals, ...
    'i1,j1,k1,i2,j2,k2,multiplier', faultRec);

manifest = struct('schema', 'CPEDFM-CMG-export-v1', ...
    'simulator', profile.simulator, 'cmgVersion', profile.version, ...
    'fractureModel', modelType, ...
    'physics', profile.physics, 'hostModel', profile.hostModel, ...
    'unitSystem', profile.unitSystem, 'dims', dims, ...
    'matrixLayers', matrixLayers, ...
    'fracturePaddingLayers', [matrixLayers+1, dims(3)], ...
    'counts', struct('fractureCells', size(fracIJK,1), ...
        'nnc', size(nncIJK,1), 'multipliers', size(multRec,1), ...
        'faultSeals', size(faultRec,1)), ...
    'files', files, 'integration', struct( ...
        'properties', profile.propertyPlacement, ...
        'connections', profile.connectionPlacement, ...
        'phaseModel', profile.phaseNote), 'caveat', profile.caveat);
writeJSON(files.manifest, manifest);
end

function writeNumericCSV(path, header, data)
fid = fopen(path, 'w');
assert(fid > 0, 'Could not open %s.', path);
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', header);
if isempty(data), return; end
format = [repmat('%.15g,', 1, size(data,2)-1), '%.15g\n'];
fprintf(fid, format, data');
end

function writeJSON(path, value)
fid = fopen(path, 'w');
assert(fid > 0, 'Could not open %s.', path);
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', jsonencode(value, 'PrettyPrint', true));
end
