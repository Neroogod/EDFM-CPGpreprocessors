function startupPEDFMU()
%STARTUPPEDFMU Add PEDFM-U and its dependencies to the MATLAB path.
%
% Reuses (all read-only; no file outside this module is modified):
%   CPEDFM-U/UtilsCPEDFM      - EDFMGridU (unstructured fracture gridding),
%                               fracPolyCellIntersect
%   PEDFM-CPG/UtilsPEDFM      - transmultPEDFM3D_CPG, setupPEDFMOpsTPFA3D_CPG
%                               and projectedAreaPEDFM_CPG. All three are
%                               grid-agnostic: the multiplier and the Eq. 29
%                               union are face-indexed and never touch
%                               logical directions.
%   CPEDFM/UpdateMay27_2023   - fracturematrixShaleNNC3D,
%                               setupShaleEDFMOpsTPFA, assembleShaleGlobalGrid
%   hfm                       - setupEDFMOperatorsTPFA, geometry utilities
%   upr                       - PEBI grid generation (tests)
%
% Uses `hfm`, not `hfm_v2`: this module grids fractures with EDFMGridU and
% never calls the corner-point intersection engine, so the two-output
% convexpolygonarea that PEDFM-CPG requires is irrelevant here. CPEDFM-U
% loads `hfm` for the same reason.

    thisDir = fileparts(mfilename('fullpath'));
    modRoot = fileparts(thisDir);                    % .../modules

    mrstModule add hfm ad-core ad-props ad-blackoil mrst-gui upr

    % --- Cartesian release: shared NNC stage + regression oracle ---
    cartRoot = fullfile(modRoot, 'CPEDFM', 'UpdateMay27_2023');
    assert(isfolder(cartRoot), 'startupPEDFMU:missingDep', ...
        'CPEDFM/UpdateMay27_2023 not found at %s', cartRoot);
    addpath(fullfile(cartRoot, 'hfmFns'));
    addpath(fullfile(cartRoot, 'pEDFM_nnc'));
    addpath(fullfile(cartRoot, 'utils'));
    addpath(genpath(fullfile(cartRoot, 'ADFNE15')));

    % --- Unstructured fracture gridding (from CPEDFM-U) ---
    uUtils = fullfile(modRoot, 'CPEDFM-U', 'UtilsCPEDFM');
    assert(isfolder(uUtils), 'startupPEDFMU:missingDep', ...
        'CPEDFM-U/UtilsCPEDFM not found at %s', uUtils);
    addpath(uUtils);

    % --- Grid-agnostic operator stage (from PEDFM-CPG) ---
    pcpgUtils = fullfile(modRoot, 'PEDFM-CPG', 'UtilsPEDFM');
    assert(isfolder(pcpgUtils), 'startupPEDFMU:missingDep', ...
        'PEDFM-CPG/UtilsPEDFM not found at %s', pcpgUtils);
    addpath(pcpgUtils);

    % --- This module ---
    addpath(fullfile(thisDir, 'UtilsPEDFM'));

    % --- Guards ---
    assert(contains(which('EDFMGridU'), fullfile('CPEDFM-U', 'UtilsCPEDFM')), ...
        'startupPEDFMU:shadow', ...
        'EDFMGridU resolves to %s, expected the CPEDFM-U copy.', which('EDFMGridU'));

    borrowed = {'transmultPEDFM3D_CPG', 'setupPEDFMOpsTPFA3D_CPG', ...
                'projectedAreaPEDFM_CPG'};
    for k = 1:numel(borrowed)
        w = which(borrowed{k});
        assert(contains(w, fullfile('PEDFM-CPG', 'UtilsPEDFM')), ...
            'startupPEDFMU:shadow', ...
            '%s resolves to %s, expected the PEDFM-CPG copy.', borrowed{k}, w);
    end

    exported = {'pMatFracNNCsPEDFM_U', 'calculateProjectionPEDFM_U', ...
                'selectProjectionFacesPEDFM_U'};
    for k = 1:numel(exported)
        w = which(exported{k});
        assert(contains(w, fullfile('PEDFM-U', 'UtilsPEDFM')), ...
            'startupPEDFMU:shadow', ...
            '%s resolves to %s, expected this module''s copy.', exported{k}, w);
    end

    fprintf('PEDFM-U paths configured.\n');
end
