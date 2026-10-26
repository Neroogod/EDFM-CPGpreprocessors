function mcp_export_driver(modelType, simulator, outDir)
%MCP_EXPORT_DRIVER Non-interactive CMG export for the cmg-mrst MCP.
%
%   mcp_export_driver('pEDFM','IMEX','C:/path/to/out')
%
% Builds the same embedded-fracture test geometry as test22, exports via
% writeCMGEmbeddedFractureModel, and writes mcp_export_result.json for the
% Python MCP to consume.

if nargin < 1 || isempty(modelType), modelType = 'pEDFM'; end
if nargin < 2 || isempty(simulator), simulator = 'IMEX'; end
if nargin < 3 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'mcp_work', lower(modelType));
end

modelType = char(string(modelType));
simulator = upper(char(string(simulator)));
outDir = char(string(outDir));
if ~isfolder(outDir), mkdir(outDir); end

result = struct('ok', false, 'errors', {{}}, 'include', '', ...
    'manifest', '', 'table_files', struct(), 'dims', [], 'counts', struct());

try
    thisDir = fileparts(mfilename('fullpath'));
    exportDir = fullfile(thisDir, '..');
    cpgDir = fullfile(thisDir, '..', '..', 'CPEDFM-CPG');
    mrstRoot = fullfile(thisDir, '..', '..', '..');

    if exist(fullfile(mrstRoot, 'startup.m'), 'file') == 2
        run(fullfile(mrstRoot, 'startup.m'));
    end
    addpath(cpgDir);
    if exist('startupCPEDFMCPG', 'file') == 2
        startupCPEDFMCPG();
    end
    addpath(exportDir);
    mrstModule add deckformat

    tol = 1e-5;
    dims = [12, 10, 2];
    physdim = [120, 100, 20];

    G0 = cartGrid(dims, physdim);
    G0 = computeGeometry(G0);
    G0.rock = makeRock(G0, 100*milli*darcy, 0.2);

    frac = struct;
    frac(1).points = [35, 0, 0; 85, 100, 0; 85, 100, 20; 35, 0, 20];
    frac(1).aperture = 5*milli*meter;
    frac(1).poro = 0.4;
    frac(1).perm = 1*nano*darcy;

    [Gbase, frac] = EDFMshalegrid(G0, frac, 'Tolerance', tol, ...
        'plotgrid', false, 'fracturelist', 1); %#ok<ASGLU>
    Gbase = fracturematrixShaleNNC3D(Gbase, tol);

    switch lower(modelType)
        case 'edfm'
            Gexp = prepareEDFMExport(Gbase);
        case 'pedfm'
            Gexp = pMatFracNNCs3D(Gbase, tol);
        case 'cpedfm'
            Gexp = pMatFracNNCsCPEDFM_CPG(Gbase, tol);
        otherwise
            error('mcp_export_driver:badModel', 'Unknown modelType %s', modelType);
    end

    grdecl = simpleGrdecl(dims, 0);
    grdecl.COORD(1:3:end) = grdecl.COORD(1:3:end)*physdim(1);
    grdecl.COORD(2:3:end) = grdecl.COORD(2:3:end)*physdim(2);
    grdecl.COORD(3:3:end) = grdecl.COORD(3:3:end)*physdim(3);
    grdecl.ZCORN = grdecl.ZCORN*physdim(3);

    outInc = fullfile(outDir, sprintf('%s_%s.inc', lower(modelType), lower(simulator)));
    info = writeCMGEmbeddedFractureModel(modelType, Gexp, Gexp.rock, grdecl, ...
        outInc, 'simulator', simulator, 'cmgVersion', '2024.30', 'tol', tol);

    result.ok = true;
    result.include = outInc;
    result.dims = info.dims;
    result.counts = struct('nnc', info.nNNC, 'multipliers', info.nMult, ...
        'faultSeals', info.nFault, 'fractureCells', size(info.fracIJK, 1));
    if isfield(info, 'tableFiles')
        result.table_files = info.tableFiles;
        if isfield(info.tableFiles, 'manifest')
            result.manifest = info.tableFiles.manifest;
        end
    end
    result.model_type = info.modelType;
    result.simulator = info.simulator;
catch ME
    result.ok = false;
    result.errors = {ME.message};
    if ~isempty(ME.stack)
        result.errors{end+1} = sprintf('%s:%d', ME.stack(1).name, ME.stack(1).line);
    end
end

outJson = fullfile(outDir, 'mcp_export_result.json');
fid = fopen(outJson, 'w');
assert(fid > 0, 'Could not write %s', outJson);
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', jsonencode(result, 'PrettyPrint', true));
clear cleaner;

if ~result.ok
    error('mcp_export_driver:failed', '%s', strjoin(result.errors, ' | '));
end
end
