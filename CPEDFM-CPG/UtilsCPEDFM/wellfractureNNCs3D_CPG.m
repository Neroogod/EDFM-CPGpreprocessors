function [G, wells] = wellfractureNNCs3D_CPG(G, fracplanes, wells, tol)
%WELLFRACTURENNCS3D_CPG Locate well/fracture intersections on CPG grids.
%
% The legacy routine maps an intersection point through Pts2CellIdx and
% assumes exactly one matrix cell and one matching fracture NNC. A point on
% a distorted or fault-split CPG interface can legitimately return several
% matrix cells. This implementation first computes the geometric
% well/polygon intersection and then selects the nearest lower-dimensional
% cell belonging to that specific fracture plane.

if nargin < 4, tol = 1e-5; end
if isempty(wells) || isempty(fieldnames(wells)), return; end

for wi = 1:numel(wells)
    points = wells(wi).points;
    lines = [points(1:end-1,:), points(2:end,:)];
    hitPts = zeros(0,3);
    hitIDs = zeros(0,1);

    for fi = 1:numel(fracplanes)
        pts = linePolygonIntersections(lines, fracplanes(fi).points, tol);
        if ~isempty(pts)
            hitPts = [hitPts; pts]; %#ok<AGROW>
            hitIDs = [hitIDs; repmat(fi, size(pts,1), 1)]; %#ok<AGROW>
        end
    end

    wells(wi).WellXFracID = hitIDs;
    wells(wi).WellXFracPts = hitPts;
    xFrac = zeros(numel(hitIDs), 1);
    for h = 1:numel(hitIDs)
        fi = hitIDs(h);
        name = ['Frac', num2str(fi)];
        assert(isfield(G.FracGrid, name), 'Missing G.FracGrid.%s.', name);
        F = G.FracGrid.(name);
        ids = (F.cells.start:F.cells.start + F.cells.num - 1)';
        delta = G.cells.centroids(ids,:) - hitPts(h,:);
        [~, nearest] = min(sum(delta.^2, 2));
        xFrac(h) = ids(nearest);
    end
    wells(wi).XFracCellIDs = [xFrac, hitIDs];
    fprintf('Total number of fractures that intersect well #%d is %d\n', ...
        wi, numel(hitIDs));
end

for fi = 1:numel(fracplanes)
    name = ['Frac', num2str(fi)];
    if isfield(G.FracGrid, name) && isfield(fracplanes(fi), 'aperture')
        G.FracGrid.(name).aperture = fracplanes(fi).aperture;
    end
end
end

function hits = linePolygonIntersections(lines, polygon, tol)
p0 = polygon(1,:);
u = polygon(2,:) - p0;
v = polygon(3,:) - p0;
nrm = cross(u, v);
nrm = nrm/norm(nrm);

hits = zeros(0,3);
for q = 1:size(lines,1)
    a = lines(q,1:3);
    d = lines(q,4:6) - a;
    denom = dot(nrm, d);
    if abs(denom) < tol, continue; end
    t = dot(nrm, p0 - a)/denom;
    if t < -tol || t > 1 + tol, continue; end
    p = a + t*d;
    if pointInPolygon3D(p, polygon, nrm, tol)
        hits(end+1,:) = p; %#ok<AGROW>
    end
end
if ~isempty(hits)
    hits = uniquetol(hits, max(tol, 1e-12), 'ByRows', true);
end
end

function inside = pointInPolygon3D(point, polygon, normal, tol)
% Drop the coordinate most parallel to the plane normal and use inpolygon.
[~, drop] = max(abs(normal));
keep = setdiff(1:3, drop);
[inside, on] = inpolygon(point(keep(1)), point(keep(2)), ...
    polygon(:,keep(1)), polygon(:,keep(2)));
inside = inside || on || min(vecnorm(polygon - point, 2, 2)) <= tol;
end
