function info = writeCMGEDFM(G, rock, grdecl, outFile, varargin)
% Export non-intrusive EDFM data as a CMG-style include file.

Ge = prepareEDFMExport(G);
info = writeCMGCPEDFM(Ge, rock, grdecl, outFile, varargin{:});
info.modelType = 'EDFM';
assert(info.nMult == 0 && info.nFault == 0, ...
    'EDFM export unexpectedly produced projection multipliers.');
end
