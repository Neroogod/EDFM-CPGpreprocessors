function P = dedupPoly_U(P, tol)
% Drop vertices that coincide with their cyclic successor to within tol.
%
% SYNOPSIS:
%   P = dedupPoly_U(P, tol)
%
% Clipping a polygon against several planes routinely produces repeated
% vertices where a plane grazes an existing corner. Those duplicates are
% harmless for area but break the loop ordering used when a polygon is
% turned into grid faces, so they are removed here.

n = size(P,1);
if n < 2
    return
end

keep = true(n,1);
for i = 1:n
    j = mod(i, n) + 1;
    if norm(P(i,:) - P(j,:)) < tol
        keep(j) = false;
    end
end
P = P(keep,:);
end
