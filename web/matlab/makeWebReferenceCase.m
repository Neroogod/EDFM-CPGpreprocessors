%% makeWebReferenceCase
% Produce a matched pair for the browser writer's regression test:
%
%   web/static/testdata/ref.grdecl   input the JS parses
%   web/static/testdata/ref.inc      what writeCMGModel.m produces
%   web/static/testdata/ref.dat      what writeCMG.m's stub deck produces
%
% The browser must reproduce ref.inc and ref.dat byte for byte from
% ref.grdecl alone. For that comparison to mean anything the two sides must
% agree on the active set, so this case is built with every block active and
% asserts that processGRDECL kept them all - MRST also discards disconnected
% components and degenerate cells, which a GRDECL parser cannot know about.

clc; clear;
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
webSetupPaths();

outDir = fullfile(webRoot(), 'web', 'static', 'testdata');
if ~isfolder(outDir), mkdir(outDir); end

%% ---- a small faulted corner-point grid --------------------------------
dims    = [8 6 4];
physdim = [400 300 40];
datum   = 2500;

grdecl = simpleGrdecl(dims, 0.12);
grdecl.COORD(1:3:end) = grdecl.COORD(1:3:end)*physdim(1);
grdecl.COORD(2:3:end) = grdecl.COORD(2:3:end)*physdim(2);
grdecl.COORD(3:3:end) = grdecl.COORD(3:3:end)*physdim(3) + datum;
grdecl.ZCORN          = grdecl.ZCORN*physdim(3) + datum;
grdecl.ACTNUM         = ones(prod(dims), 1);

% Varying properties, so the writer takes the *ALL path rather than *CON.
nCells = prod(dims);
[ii, jj, kk] = ind2sub(dims, (1:nCells)');
grdecl.PERMX = 100 + 5*ii + 3*jj + 20*kk;          % mD
grdecl.PERMY = grdecl.PERMX;
grdecl.PERMZ = 0.1*grdecl.PERMX;
grdecl.PORO  = 0.18 + 0.001*ii + 0.002*kk;

G = processGRDECL(grdecl);
G = computeGeometry(G);

assert(G.cells.num == nCells, ...
    ['processGRDECL kept %d of %d cells. The reference case must have ' ...
     'every block active, otherwise the browser (which only sees ACTNUM) ' ...
     'cannot match it.'], G.cells.num, nCells);

% Permeability in SI, so writeCMG converts it back to exactly the mD values
% that are written into the GRDECL.
rock = makeRock(G, [grdecl.PERMX, grdecl.PERMY, grdecl.PERMZ]*milli*darcy, ...
    grdecl.PORO);

%% ---- write the GRDECL the browser will read ---------------------------
grdeclFile = fullfile(outDir, 'ref.grdecl');
fid = fopen(grdeclFile, 'w');
assert(fid > 0, 'cannot write %s', grdeclFile);
fprintf(fid, '-- Reference case for the browser writer regression test\n');
fprintf(fid, '-- %d x %d x %d corner-point grid, every block active\n\n', dims);
fprintf(fid, 'SPECGRID\n %d %d %d 1 F /\n\n', dims);
localArray(fid, 'COORD',  grdecl.COORD,  '%.6f');
localArray(fid, 'ZCORN',  grdecl.ZCORN,  '%.6f');
localArray(fid, 'ACTNUM', grdecl.ACTNUM, '%d');
localArray(fid, 'PORO',   grdecl.PORO,   '%.6f');
localArray(fid, 'PERMX',  grdecl.PERMX,  '%.6f');
localArray(fid, 'PERMY',  grdecl.PERMY,  '%.6f');
localArray(fid, 'PERMZ',  grdecl.PERMZ,  '%.6f');
fclose(fid);

%% ---- write the reference deck -----------------------------------------
info = writeCMG(G, fullfile(outDir, 'ref.dat'), ...
    'rock', rock, 'grdecl', grdecl, ...
    'simulator', 'IMEX', 'cmgVersion', '2024.30', ...
    'title', 'REF');

fprintf('grdecl  : %s\n', grdeclFile);
fprintf('dataset : %s\n', info.dataset);
fprintf('blocks  : %d (%d inactive)\n', info.counts.blocks, info.counts.inactive);

%% ---- uniform Cartesian, the *CART path --------------------------------
cdims = [10 7 5];
cphys = [500 350 50];
Gc = computeGeometry(cartGrid(cdims, cphys));
nc = Gc.cells.num;
[ci, cj, ck] = ind2sub(cdims, (1:nc)');
rockC = makeRock(Gc, ...
    [80 + 4*ci + 2*cj + 15*ck, 80 + 4*ci + 2*cj + 15*ck, ...
     8 + 0.4*ci + 0.2*cj + 1.5*ck]*milli*darcy, ...
    0.15 + 0.002*ci + 0.003*ck);

infoC = writeCMG(Gc, fullfile(outDir, 'ref_cart.dat'), ...
    'rock', rockC, 'simulator', 'IMEX', 'cmgVersion', '2024.30', ...
    'datum', 1800, 'title', 'REF_CART');
fprintf('cart    : %s  (%d blocks)\n', infoC.dataset, infoC.counts.blocks);

localWriteSpec(fullfile(outDir, 'ref_cart.json'), struct( ...
    'kind', 'cart', 'dims', cdims, 'physdim', cphys, 'datum', 1800, ...
    'title', 'REF_CART', ...
    'permi', convertTo(rockC.perm(:,1), milli*darcy), ...
    'permj', convertTo(rockC.perm(:,2), milli*darcy), ...
    'permk', convertTo(rockC.perm(:,3), milli*darcy), ...
    'poro',  rockC.poro));

%% ---- variable Cartesian, the *VARI path -------------------------------
xs = [0 12 30 40 70 85 130];
ys = [0 20 45 100];
zs = [0 4 12 18];
Gv = computeGeometry(tensorGrid(xs, ys, zs));
nv = Gv.cells.num;
[vi, vj, vk] = ind2sub(Gv.cartDims, (1:nv)');
rockV = makeRock(Gv, ...
    [120 + 6*vi - 3*vj, 120 + 6*vi - 3*vj, 12 + 0.6*vi]*milli*darcy, ...
    0.21 + 0.004*vk);

infoV = writeCMG(Gv, fullfile(outDir, 'ref_vari.dat'), ...
    'rock', rockV, 'simulator', 'IMEX', 'cmgVersion', '2024.30', ...
    'datum', 2200, 'title', 'REF_VARI');
fprintf('vari    : %s  (%d blocks)\n', infoV.dataset, infoV.counts.blocks);

localWriteSpec(fullfile(outDir, 'ref_vari.json'), struct( ...
    'kind', 'vari', 'dims', Gv.cartDims, 'datum', 2200, 'title', 'REF_VARI', ...
    'dx', diff(xs), 'dy', diff(ys), 'dz', diff(zs), ...
    'permi', convertTo(rockV.perm(:,1), milli*darcy), ...
    'permj', convertTo(rockV.perm(:,2), milli*darcy), ...
    'permk', convertTo(rockV.perm(:,3), milli*darcy), ...
    'poro',  rockV.poro));

%% ---- Eclipse counterparts of all three ---------------------------------
ea = writeEclipse(G,  fullfile(outDir, 'ref_ecl.DATA'), ...
    'rock', rock,  'grdecl', grdecl, 'title', 'REF');
eb = writeEclipse(Gc, fullfile(outDir, 'ref_cart_ecl.DATA'), ...
    'rock', rockC, 'datum', 1800, 'title', 'REF_CART');
ec = writeEclipse(Gv, fullfile(outDir, 'ref_vari_ecl.DATA'), ...
    'rock', rockV, 'datum', 2200, 'title', 'REF_VARI');
fprintf('eclipse : %s\n', ea.dataset);
fprintf('          %s\n', eb.dataset);
fprintf('          %s\n', ec.dataset);

fprintf('\n=== reference cases written ===\n');

%%-------------------------------------------------------------------------
function localWriteSpec(path, s)
% The grid description the browser is given, so both sides start from
% identical inputs rather than from separately-typed numbers.
fid = fopen(path, 'w');
assert(fid > 0, 'cannot write %s', path);
c = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', jsonencode(s));
end

%%-------------------------------------------------------------------------
function localArray(fid, kw, v, fmt)
% Eclipse free format: keyword, values, terminating slash.
fprintf(fid, '%s\n', kw);
v = v(:);
perLine = 6;
full = perLine*floor(numel(v)/perLine);
if full > 0
    rowFmt = [repmat([' ' fmt], 1, perLine) '\n'];
    fprintf(fid, rowFmt, reshape(v(1:full), perLine, []));
end
if full < numel(v)
    fprintf(fid, [' ' fmt], v(full+1:end));
    fprintf(fid, '\n');
end
fprintf(fid, '/\n\n');
end
