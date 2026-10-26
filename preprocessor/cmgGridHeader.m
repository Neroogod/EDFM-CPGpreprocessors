function lines = cmgGridHeader(model)
% The *GRID / *KDIR cards that declare the fundamental grid.
%
% SYNOPSIS:
%   lines = cmgGridHeader(model)   % cellstr, one card per entry
%
% WHY THIS IS SEPARATE:
%   Builder will not accept the fundamental grid definition from inside an
%   *INCLUDE file. Given
%
%       *INCLUDE 'case.inc'        with *GRID *CORNER ... inside it
%
%   Builder reads the dataset, reports an error at the end of the *ZCORN
%   block, and shows an empty grid ("Plane 1 of 1"), while IMEX accepts the
%   same file without complaint. Moving these two cards into the .dat and
%   leaving the arrays in the include makes Builder render the grid.
%
%   So the include never carries *GRID/*KDIR; the dataset does. If you are
%   including the file into a deck of your own, that deck needs these cards
%   ahead of the *INCLUDE line - they are repeated in the include banner for
%   exactly that reason.
%
% SEE ALSO:
%   writeCMGModel, writeCMG

dims = model.grid.dims;

switch model.grid.kind
    case 'cart'
        kw = '*CART';
    case {'vari', 'connlist'}
        kw = '*VARI';
    case 'corner'
        kw = '*CORNER';
    otherwise
        error('cmgGridHeader:kind', 'Unhandled grid kind %s', model.grid.kind);
end

lines = { ...
    sprintf('*GRID %s %d %d %d', kw, dims(1), dims(2), dims(3)), ...
    '*KDIR *DOWN'};
end
