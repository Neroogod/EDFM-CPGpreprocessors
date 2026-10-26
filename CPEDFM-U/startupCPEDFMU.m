function startupCPEDFMU()
% Add CPEDFM-U and its dependencies to the MATLAB path.
%
% Reuses (read-only):
%   EDFM-U/UtilsEDFM           - EDFMgrid_U (gridding) and
%                                fracturematrixNNC_U (fracture-matrix NNCs
%                                with the exact <d> of calcDavg_U). CPEDFM-U
%                                supplies only the projection stage on top.
%   CPEDFM-CPG/UtilsCPEDFM     - transmultpEDFM_CPG, setupCPEDFMOpsTPFA_CPG
%                                (grid-agnostic operator stage)
%   CPEDFM/UpdateMay27_2023    - assembleShaleGlobalGrid,
%                                setupShaleEDFMOpsTPFA
%   hfm                        - fractureplanegeneralgrid,
%                                setupEDFMOperatorsTPFA, geometry utils
%   upr                        - PEBI grid generation (tests)

thisDir = fileparts(mfilename('fullpath'));
modRoot = fileparts(thisDir);

mrstModule add hfm ad-core ad-props ad-blackoil mrst-gui upr

cartCPEDFM = fullfile(modRoot, 'CPEDFM', 'UpdateMay27_2023');
assert(isfolder(cartCPEDFM), 'CPEDFM/UpdateMay27_2023 not found at %s', cartCPEDFM);
addpath(fullfile(cartCPEDFM, 'hfmFns'));
addpath(fullfile(cartCPEDFM, 'pEDFM_nnc'));
addpath(fullfile(cartCPEDFM, 'utils'));

cpgUtils = fullfile(modRoot, 'CPEDFM-CPG', 'UtilsCPEDFM');
assert(isfolder(cpgUtils), 'CPEDFM-CPG/UtilsCPEDFM not found at %s', cpgUtils);
addpath(cpgUtils);

edfmUtils = fullfile(modRoot, 'EDFM-U', 'UtilsEDFM');
assert(isfolder(edfmUtils), ...
    ['EDFM-U/UtilsEDFM not found at %s. CPEDFM-U''s gridding and ', ...
     'fracture-matrix NNC stages live there.'], edfmUtils);
addpath(edfmUtils);

addpath(fullfile(thisDir, 'UtilsCPEDFM'));

assert(contains(which('transmultpEDFM_CPG'), 'CPEDFM-CPG'), ...
    'transmultpEDFM_CPG resolves to unexpected copy: %s', which('transmultpEDFM_CPG'));
for f = {'EDFMgrid_U', 'fracturematrixNNC_U', 'calcDavg_U', 'clipPolyByCell_U'}
    w = which(f{1});
    assert(contains(w, fullfile('EDFM-U', 'UtilsEDFM')), ...
        '%s resolves to an unexpected copy: %s', f{1}, w);
end

fprintf('CPEDFM-U paths configured.\n');
end
