function k = permNormal_U(perm, n)
% Normal-projected permeability n'*K*n, row-wise.
%
% SYNOPSIS:
%   k = permNormal_U(perm, n)
%
% PARAMETERS:
%   perm - Nx1, Nx3, Nx6 or Nx9 MRST permeability rows.
%   n    - 1x3 unit direction, or Nx3 (one direction per row).
%
% RETURNS:
%   k - Nx1 permeability seen along n.
%
% Flow between a fracture and its host matrix cell is normal to the
% fracture plane, so an anisotropic matrix must be sampled along the plane
% normal rather than by an arbitrary principal component.

N = size(perm, 1);
if size(n,1) == 1
    n = repmat(n, N, 1);
end
n = n ./ sqrt(sum(n.^2, 2));

switch size(perm, 2)
    case 1
        k = perm;
    case 3
        k = perm(:,1).*n(:,1).^2 + perm(:,2).*n(:,2).^2 + perm(:,3).*n(:,3).^2;
    case 6
        % [k11 k12 k13 k22 k23 k33]
        k = perm(:,1).*n(:,1).^2 + perm(:,4).*n(:,2).^2 + perm(:,6).*n(:,3).^2 ...
          + 2*(perm(:,2).*n(:,1).*n(:,2) + perm(:,3).*n(:,1).*n(:,3) ...
             + perm(:,5).*n(:,2).*n(:,3));
    case 9
        k = perm(:,1).*n(:,1).^2 + perm(:,5).*n(:,2).^2 + perm(:,9).*n(:,3).^2 ...
          + (perm(:,2)+perm(:,4)).*n(:,1).*n(:,2) ...
          + (perm(:,3)+perm(:,7)).*n(:,1).*n(:,3) ...
          + (perm(:,6)+perm(:,8)).*n(:,2).*n(:,3);
    otherwise
        error('permNormal_U:format', ...
            'Unsupported permeability format (%d components).', size(perm,2));
end
end
