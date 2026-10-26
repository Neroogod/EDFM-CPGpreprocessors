function info = writeCMGEmbeddedFractureModel(modelType, G, rock, grdecl, outFile, varargin)
% Unified CMG entry point for EDFM, pEDFM, and CPEDFM models.

modelType = validatestring(modelType, {'EDFM', 'pEDFM', 'CPEDFM'});
switch lower(modelType)
    case 'edfm'
        info = writeCMGEDFM(G, rock, grdecl, outFile, varargin{:});
    case 'pedfm'
        info = writeCMGPEDFM(G, rock, grdecl, outFile, varargin{:});
    case 'cpedfm'
        info = writeCMGCPEDFM(G, rock, grdecl, outFile, varargin{:});
        info.modelType = 'CPEDFM';
end
end
