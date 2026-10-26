function startupEDFMU()
% Add EDFM-U and its dependencies to the MATLAB path.
%
% Reuses (read-only, nothing is modified):
%   CPEDFM/UpdateMay27_2023/hfmFns - assembleShaleGlobalGrid,
%                                    setupShaleEDFMOpsTPFA,
%                                    fracturefractureShaleNNCs3D
%   hfm                            - fractureplanegeneralgrid,
%                                    setupEDFMOperatorsTPFA,
%                                    PEBIPEBIintersection3D, geometry utils
%   ad-core / ad-props / ad-blackoil - simulation stack
%   upr                            - PEBI grid generation (tests, examples)

thisDir = fileparts(mfilename('fullpath'));
modRoot = fileparts(thisDir);

mrstModule add hfm ad-core ad-props ad-blackoil mrst-gui upr

shaleFns = fullfile(modRoot, 'CPEDFM', 'UpdateMay27_2023', 'hfmFns');
assert(isfolder(shaleFns), ...
    'CPEDFM/UpdateMay27_2023/hfmFns not found at %s', shaleFns);
addpath(shaleFns);

addpath(fullfile(thisDir, 'UtilsEDFM'));

% EDFM-U must own these names; nothing else in the tree defines them.
for f = {'EDFMgrid_U', 'calcDavg_U', 'clipCellByPlane_U', 'fracMatNNC_U'}
    w = which(f{1});
    assert(contains(w, fullfile('EDFM-U', 'UtilsEDFM')), ...
        '%s resolves to an unexpected copy: %s', f{1}, w);
end

fprintf('EDFM-U paths configured.\n');
end
