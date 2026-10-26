function info = writeEclipsePEDFM(G, rock, grdecl, wells, outFile, varargin)
% Export a standard pEDFM model to a self-contained Eclipse deck.
%
% Exporter-specific option 'ambiguousFacePolicy' is consumed here; all
% remaining options are passed to writeEclipseCPEDFM.

[policy, writerArgs] = localTakeOption(varargin, 'ambiguousFacePolicy', 'error');
Gp = preparePEDFMExport(G, 'ambiguousFacePolicy', policy);
info = writeEclipseCPEDFM(Gp, rock, grdecl, wells, outFile, ...
    'title', 'PEDFM_EXPORT', writerArgs{:});
info.modelType = 'pEDFM';
end

function [value, args] = localTakeOption(args, key, defaultValue)
value = defaultValue;
names = args(1:2:end);
idx = find(cellfun(@(x) ischar(x) || isstring(x), names) & ...
    cellfun(@(x) strcmpi(string(x), key), names), 1);
if ~isempty(idx)
    pos = 2*idx - 1;
    value = args{pos+1};
    args(pos:pos+1) = [];
end
end
