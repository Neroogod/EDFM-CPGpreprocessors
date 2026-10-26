function info = webSetupPaths()
% Put MRST, the preprocessor and the EDFM modules on the path, in the one
% order that works, and report what actually resolved.
%
% ORDER IS LOAD-BEARING. startupCPEDFMCPG internally runs "mrstModule add
% hfm", and hfm's convexpolygonarea returns a single output while
% EDFMGridCPG calls it with two. Adding hfm_v2 afterwards puts the
% two-output copy back in front. Reverse the order and the corner-point
% gridder dies with "Too many output arguments" at
% pebiAABBintersectCPGv5:648.
%
% The workspace holds six copies of EDFMGridCPG and three of
% pebiAABBintersectCPGv5, so which() is recorded for the functions whose
% identity changes results. A run log that names them turns a wrong pairing
% into something visible instead of a strange number.

root = webRoot();

if exist('mrstPath', 'file') ~= 2
    run(fullfile(root, 'startup.m'));
end

cpgDir = fullfile(root, 'modules', 'CPEDFM-CPG');
if isfolder(cpgDir)
    addpath(cpgDir);
    if exist('startupCPEDFMCPG', 'file') == 2
        evalc('startupCPEDFMCPG()');      % it is chatty; we log our own
    end
end

% After the CPG startup, never before.
mrstModule add hfm_v2 deckformat ad-core ad-props

preDir = fullfile(root, 'modules', 'preprocessor');
addpath(preDir);

shaleUtils = '';
if ~isempty(mrstPath('shale'))
    shaleUtils = fullfile(mrstPath('shale'), 'utils');
    if isfolder(shaleUtils), addpath(shaleUtils); end
end

info = struct( ...
    'root', root, ...
    'matlab', version(), ...
    'preprocessor', preDir, ...
    'shaleUtils', shaleUtils, ...
    'resolved', struct( ...
        'EDFMGridCPG',              localWhich('EDFMGridCPG'), ...
        'EDFMshalegrid',            localWhich('EDFMshalegrid'), ...
        'convexpolygonarea',        localWhich('convexpolygonarea'), ...
        'calcdavg',                 localWhich('calcdavg'), ...
        'fracturematrixShaleNNC3D', localWhich('fracturematrixShaleNNC3D'), ...
        'writeCMG',                 localWhich('writeCMG')));
end

%--------------------------------------------------------------------------
function p = localWhich(name)
p = which(name);
if isempty(p), p = '(not found)'; end
end
