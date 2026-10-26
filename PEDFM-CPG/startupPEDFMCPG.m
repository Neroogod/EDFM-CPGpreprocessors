function startupPEDFMCPG()
%STARTUPPEDFMCPG Add PEDFM-CPG and its dependencies to the MATLAB path.
%
% Dependencies (all read-only; no file outside this module is modified):
%   1. hfm_v2 - EDFM primitives: setupEDFMOperatorsTPFA, pointplanedistance,
%      geometry utilities. NOT hfm: the corner-point intersection engine
%      calls convexpolygonarea with two outputs ([area, points], see
%      pebiAABBintersectCPGv5 line 648) and only the hfm_v2 copy returns
%      them. Loading hfm instead fails inside EDFMGridCPG's parfor with
%      "Too many output arguments". The guard below asserts the right copy
%      wins.
%   2. CPEDFM/UpdateMay27_2023 - fracturematrixShaleNNC3D,
%      setupShaleEDFMOpsTPFA, ADFNE, and pMatFracNNCs3D (the Cartesian
%      pEDFM used as the regression oracle in
%      tests/test00_regression_cartesian.m).
%   3. CPEDFMcornerpointaug24/Updates_8_23 - corner-point fracture gridding:
%      EDFMGridCPG, pebiAABBintersectCPGv5, assembleShaleGlobalGrid.
%
% This module's exported names all carry a PEDFM3D_CPG / PEDFM_CPG suffix so
% that PEDFM-CPG and CPEDFM-CPG can be on the path at the same time without
% shadowing each other. The guards at the end assert that.

    thisDir = fileparts(mfilename('fullpath'));
    modRoot = fileparts(thisDir);                    % .../modules

    % --- MRST core modules ---
    mrstModule add hfm_v2 ad-core ad-props ad-blackoil mrst-gui

    % --- Cartesian release: shared NNC stage + regression oracle ---
    cartRoot = fullfile(modRoot, 'CPEDFM', 'UpdateMay27_2023');
    assert(isfolder(cartRoot), ...
        'startupPEDFMCPG:missingDep', ...
        'CPEDFM/UpdateMay27_2023 not found at %s', cartRoot);
    addpath(fullfile(cartRoot, 'hfmFns'));
    addpath(fullfile(cartRoot, 'pEDFM_nnc'));
    addpath(fullfile(cartRoot, 'utils'));
    addpath(genpath(fullfile(cartRoot, 'ADFNE15')));

    % --- Corner-point fracture gridding ---
    cpgDir = fullfile(modRoot, 'CPEDFMcornerpointaug24', 'Updates_8_23');
    assert(isfolder(cpgDir), ...
        'startupPEDFMCPG:missingDep', ...
        'CPEDFMcornerpointaug24/Updates_8_23 not found at %s', cpgDir);
    addpath(cpgDir);

    % --- This module ---
    addpath(fullfile(thisDir, 'UtilsPEDFM'));
    addpath(fullfile(thisDir, 'examples'));

    % ADFNE needs the two-argument exponential CDF/inverse CDF. Borrow the
    % lightweight equivalents from CPEDFM-CPG when the Statistics Toolbox is
    % unavailable, but never shadow MATLAB's licensed implementations.
    if exist('expcdf', 'file') ~= 2 || exist('expinv', 'file') ~= 2
        compatDir = fullfile(modRoot, 'CPEDFM-CPG', 'UtilsCPEDFM', 'compat', 'statistics');
        if isfolder(compatDir)
            addpath(compatDir);
        end
    end

    % --- Guards ---
    assert(contains(which('EDFMGridCPG'), 'Updates_8_23'), ...
        'startupPEDFMCPG:shadow', ...
        'EDFMGridCPG resolves to %s, expected the Updates_8_23 copy.', ...
        which('EDFMGridCPG'));

    % convexpolygonarea must be the two-output copy (see the note above).
    % Check the signature rather than the path, so the guard keeps working
    % if the file moves.
    cpa = which('convexpolygonarea');
    assert(~isempty(cpa), 'startupPEDFMCPG:missingDep', ...
        'convexpolygonarea not found on the path.');
    sigLine = '';
    fid = fopen(cpa, 'r');
    if fid > 0
        c = onCleanup(@() fclose(fid));
        while ischar(sigLine) && ~startsWith(strtrim(sigLine), 'function')
            sigLine = fgetl(fid);
        end
    end
    assert(ischar(sigLine) && contains(sigLine, '['), ...
        'startupPEDFMCPG:shadow', ...
        ['convexpolygonarea resolves to the single-output copy at %s. ', ...
         'The corner-point intersection engine needs the two-output ', ...
         'version from hfm_v2; loading the hfm module instead shadows it ', ...
         'and EDFMGridCPG then fails with "Too many output arguments".'], cpa);

    % Every exported name of this module must resolve inside UtilsPEDFM.
    exported = {'pMatFracNNCsPEDFM3D_CPG', 'selectProjectionFacesPEDFM_CPG', ...
                'getLogicalFacesPEDFM_CPG', 'projectedAreaPEDFM_CPG', ...
                'transmultPEDFM3D_CPG', 'setupPEDFMOpsTPFA3D_CPG'};
    for k = 1:numel(exported)
        w = which(exported{k});
        assert(contains(w, fullfile('PEDFM-CPG', 'UtilsPEDFM')), ...
            'startupPEDFMCPG:shadow', ...
            '%s resolves to %s, expected this module''s copy.', exported{k}, w);
    end

    fprintf('PEDFM-CPG paths configured.\n');
end
