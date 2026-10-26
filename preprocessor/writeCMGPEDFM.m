function info = writeCMGPEDFM(G, rock, grdecl, outFile, varargin)
% Export standard pEDFM data as a CMG-style include file.

[policy, writerArgs] = localTakeOption(varargin, 'ambiguousFacePolicy', 'error');
Gp = preparePEDFMExport(G, 'ambiguousFacePolicy', policy);
info = writeCMGCPEDFM(Gp, rock, grdecl, outFile, writerArgs{:});
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
