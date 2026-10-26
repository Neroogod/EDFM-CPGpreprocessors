function startupCPEDFMCPG()
% Add CPEDFM-CPG and its dependencies to the MATLAB path.
%
% Dependencies (all read-only, none of their files are modified):
%   1. CPEDFM/UpdateMay27_2023  - Cartesian CPEDFM release:
%      fracturematrixShaleNNC3D, setupPEDFMOpsTPFA (regression oracle),
%      EDFMshalegrid, assembleShaleGlobalGrid, ADFNE15.
%   2. CPEDFMcornerpointaug24/Updates_8_23 - corner-point gridding:
%      EDFMGridCPG, pebiAABBintersectCPGv5, assembleShaleGlobalGrid.
%   3. hfm - EDFM primitives: setupEDFMOperatorsTPFA, pointplanedistance,
%      fractureplanegeneralgrid, geometry utilities.
%
% Path order matters: this module's UtilsCPEDFM/ is added LAST so its functions
% shadow nothing and are shadowed by nothing (all names are new).

thisDir = fileparts(mfilename('fullpath'));
modRoot = fileparts(thisDir); % .../modules

% --- MRST core modules ---
mrstModule add hfm ad-core ad-props ad-blackoil mrst-gui

% --- Cartesian CPEDFM release (regression oracle + shared NNC stage) ---
cartCPEDFM = fullfile(modRoot, 'CPEDFM', 'UpdateMay27_2023');
assert(isfolder(cartCPEDFM), 'CPEDFM/UpdateMay27_2023 not found at %s', cartCPEDFM);
addpath(fullfile(cartCPEDFM, 'hfmFns'));
addpath(fullfile(cartCPEDFM, 'pEDFM_nnc'));
addpath(fullfile(cartCPEDFM, 'CPEDFM', 'UtilsCPEDFM'));
addpath(fullfile(cartCPEDFM, 'utils'));
addpath(genpath(fullfile(cartCPEDFM, 'ADFNE15')));

% --- Data used by the replicated Eagle Ford paper examples ---
shaleData = fullfile(modRoot, 'shale', 'examples', 'data');
if isfolder(shaleData)
    addpath(shaleData);
end

% --- Corner-point gridding (EDFMGridCPG + intersection engine) ---
cpgDir = fullfile(modRoot, 'CPEDFMcornerpointaug24', 'Updates_8_23');
assert(isfolder(cpgDir), 'CPEDFMcornerpointaug24/Updates_8_23 not found at %s', cpgDir);
addpath(cpgDir);

% --- This module ---
addpath(fullfile(thisDir, 'UtilsCPEDFM'));
addpath(fullfile(thisDir, 'examples'));

% Bundled ADFNE only needs the two-argument exponential CDF/inverse CDF.
% Supply lightweight equivalents when Statistics Toolbox is unavailable,
% but never shadow MATLAB's licensed implementations.
if exist('expcdf', 'file') ~= 2 || exist('expinv', 'file') ~= 2
    addpath(fullfile(thisDir, 'UtilsCPEDFM', 'compat', 'statistics'));
end

% --- Guard against wrong shadow copies ---
% Both dependency folders carry assembleShaleGlobalGrid; EDFMGridCPG must
% resolve to the corner-point version (only location on path by design).
assert(contains(which('EDFMGridCPG'), 'Updates_8_23'), ...
    'EDFMGridCPG resolves to %s, expected the Updates_8_23 copy.', which('EDFMGridCPG'));
assert(contains(which('pMatFracNNCs3DCPEDFM'), 'UtilsCPEDFM'), ...
    'pMatFracNNCs3DCPEDFM resolves to unexpected copy: %s', which('pMatFracNNCs3DCPEDFM'));

fprintf('CPEDFM-CPG paths configured.\n');
end
