function info = writeEclipseEmbeddedFractureModel(modelType, G, rock, grdecl, wells, outFile, varargin)
% Unified Eclipse entry point for EDFM, pEDFM, and CPEDFM models.

modelType = validatestring(modelType, {'EDFM', 'pEDFM', 'CPEDFM'});
switch lower(modelType)
    case 'edfm'
        info = writeEclipseEDFM(G, rock, grdecl, wells, outFile, varargin{:});
    case 'pedfm'
        info = writeEclipsePEDFM(G, rock, grdecl, wells, outFile, varargin{:});
    case 'cpedfm'
        info = writeEclipseCPEDFM(G, rock, grdecl, wells, outFile, varargin{:});
        info.modelType = 'CPEDFM';
end
end
