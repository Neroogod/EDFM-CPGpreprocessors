function out = webStage1Grid(outFile, opts)
% Stage 1: build the matrix grid and rock, and checkpoint them.
%
% opts.source is one of
%   'cartesian'  dims, physdim
%   'tensor'     dx, dy, dz  (vectors of block sizes)
%   'grdecl'     file        (path to a .GRDECL)
%   'builtin'    name        ('simple' | 'johansen' | 'norne')
%
% opts.perm [mD], opts.poro, opts.datum [m] supply rock and depth where the
% source does not carry them.

source = lower(localField(opts, 'source', 'cartesian'));
permMD = localField(opts, 'perm', 100);
poro   = localField(opts, 'poro', 0.2);
datum  = localField(opts, 'datum', []);

grdecl = [];
warnings = {};

switch source
    case 'cartesian'
        dims = reshape(localField(opts, 'dims', [10 10 5]), 1, []);
        phys = reshape(localField(opts, 'physdim', [100 100 50]), 1, []);
        G = computeGeometry(cartGrid(dims, phys));

    case 'tensor'
        dx = localField(opts, 'dx', []);
        dy = localField(opts, 'dy', []);
        dz = localField(opts, 'dz', []);
        assert(~isempty(dx) && ~isempty(dy) && ~isempty(dz), ...
            'webStage1Grid:spacing', ...
            'A tensor grid needs dx, dy and dz block-size vectors.');
        G = computeGeometry(tensorGrid([0 cumsum(dx(:)')], ...
            [0 cumsum(dy(:)')], [0 cumsum(dz(:)')]));

    case 'grdecl'
        f = localField(opts, 'file', '');
        assert(exist(f, 'file') == 2, 'webStage1Grid:noFile', ...
            'GRDECL not found: %s', f);
        grdecl = readGRDECL(f);
        Gs = processGRDECL(grdecl);
        if numel(Gs) > 1
            warnings{end+1} = sprintf( ...
                ['processGRDECL produced %d disconnected components; the ' ...
                 'largest is used, as MRST does by default.'], numel(Gs));
        end
        G = computeGeometry(Gs(1));

    case 'builtin'
        [G, grdecl, extra] = localBuiltin(localField(opts, 'name', 'simple'));
        warnings = [warnings, extra];

    otherwise
        error('webStage1Grid:source', 'Unknown grid source "%s".', source);
end

% ---- rock -------------------------------------------------------------
if ~isempty(grdecl) && isfield(grdecl, 'PERMX') && isfield(grdecl, 'PORO')
    rock = grdecl2Rock(grdecl, G.cells.indexMap);
    rock.perm = max(rock.perm, 1e-6*milli*darcy);
    rock.poro = max(rock.poro, 1e-4);
else
    rock = makeRock(G, permMD*milli*darcy, poro);
end
G.rock = rock;

gridInfo = cmgDetectGrid(G);
pv = sum(poreVolume(G, rock));

save(outFile, 'G', 'rock', 'grdecl', 'datum', '-v7.3');

out.warnings = warnings;
out.stats = struct( ...
    'source', source, ...
    'cells', G.cells.num, ...
    'cartDims', localDims(G), ...
    'family', gridInfo.family, ...
    'poreVolume', pv, ...
    'depthTop', gridInfo.ztop, ...
    'depthBottom', gridInfo.zbot, ...
    'permMinMD', convertTo(min(rock.perm(:,1)), milli*darcy), ...
    'permMaxMD', convertTo(max(rock.perm(:,1)), milli*darcy), ...
    'poroMin', min(rock.poro), ...
    'poroMax', max(rock.poro));
out.files = struct('checkpoint', outFile);
end

%--------------------------------------------------------------------------
function [G, grdecl, warnings] = localBuiltin(name)
warnings = {};
switch lower(name)
    case 'simple'
        dims = [20 10 8]; phys = [200 100 40];
        grdecl = simpleGrdecl(dims, 0.12);
        grdecl.COORD(1:3:end) = grdecl.COORD(1:3:end)*phys(1);
        grdecl.COORD(2:3:end) = grdecl.COORD(2:3:end)*phys(2);
        grdecl.COORD(3:3:end) = grdecl.COORD(3:3:end)*phys(3) + 2000;
        grdecl.ZCORN = grdecl.ZCORN*phys(3) + 2000;
        grdecl.ACTNUM = ones(prod(dims), 1);
        G = computeGeometry(processGRDECL(grdecl));

    case 'johansen'
        jdir = fullfile(mrstPath('co2lab-common'), 'data', 'johansen');
        sector = fullfile(jdir, 'NPD5');
        grdecl = readGRDECL([sector '.grdecl']);
        K = reshape(load([sector, '_Permeability.txt'])', [], 1);
        p = reshape(load([sector, '_Porosity.txt'])',     [], 1);
        grdecl.PERMX = K; grdecl.PERMY = K; grdecl.PERMZ = 0.1*K;
        grdecl.PORO  = p;
        Gs = processGRDECL(grdecl);
        G = computeGeometry(Gs(1));
        warnings{end+1} = ...
            'Johansen: the largest connected component is used, as makeJohansenVEgridFull does.';

    case 'norne'
        grdecl = readGRDECL(fullfile(getDatasetPath('norne'), 'NORNE.GRDECL'));
        grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));
        grdecl = cutGrdecl(grdecl, [10 25; 35 55; 1 22]);
        G = computeGeometry(subsref(processGRDECL(grdecl), substruct('()', {1})));
        warnings{end+1} = 'Norne: sector cut [10 25; 35 55; 1 22].';

    otherwise
        error('webStage1Grid:builtin', 'Unknown built-in case "%s".', name);
end
end

%--------------------------------------------------------------------------
function d = localDims(G)
if isfield(G, 'cartDims') && ~isempty(G.cartDims)
    d = reshape(double(G.cartDims), 1, []);
else
    d = [G.cells.num, 1, 1];
end
end

%--------------------------------------------------------------------------
function v = localField(s, name, dflt)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = dflt;
end
end
