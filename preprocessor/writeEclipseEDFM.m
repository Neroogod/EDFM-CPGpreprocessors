function info = writeEclipseEDFM(G, rock, grdecl, wells, outFile, varargin)
% Export a non-intrusive EDFM model to a self-contained Eclipse deck.

Ge = prepareEDFMExport(G);
info = writeEclipseCPEDFM(Ge, rock, grdecl, wells, outFile, ...
    'title', 'EDFM_EXPORT', varargin{:});
info.modelType = 'EDFM';
assert(info.nMult == 0 && info.nEditnnc == 0, ...
    'EDFM export unexpectedly produced projection multipliers.');
end
