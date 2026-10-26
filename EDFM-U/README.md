# EDFM-U

Embedded Discrete Fracture Model (EDFM) for **general unstructured grids**:
PEBI/Voronoi, tetrahedral, extruded prisms, corner-point — any grid whose
cells are convex polyhedra. Companion to the Cartesian release in
`CPEDFM/UpdateMay27_2023` and to the projection-based modules `PEDFM-*` /
`CPEDFM-*`. No existing file is modified.

## Usage

```matlab
run modules/EDFM-U/startupEDFMU.m

G = computeGeometry(makeLayeredGrid(pebiGrid2D(h, dims), nz));
G.rock = makeRock(G, perm, poro);

[G, fracplanes] = EDFMgrid_U(G, fracplanes, 'Tolerance', tol);  % gridding
G = fracturematrixNNC_U(G, tol);                                % frac-matrix NNCs
[G, fracplanes] = fracturefractureNNC_U(G, fracplanes, tol);    % frac-frac NNCs
ops = setupShaleEDFMOpsTPFA(G, G.rock, tol);                    % reused operators
```

`fracplanes(i)` needs `points` (Nx3 ordered planar polygon), `aperture`,
`poro` and `perm`.

## What is actually different from the Cartesian code

Running the existing Cartesian EDFM on an unstructured grid is not just a
matter of calling it with a different `G`. Three things break, and this
module fixes each:

### 1. `<d>` was integrated over the cell's bounding box

EDFM couples a fracture to its host cell through the connectivity index
`CI = A / <d>`, where `<d>` is the average normal distance from the cell to
the fracture plane. `hfm/calcdavg` computes it by averaging `|n·(x-x0)|` at
the centroids of a 10×10×10 sample grid laid over the cell's **axis-aligned
bounding box**. On a Cartesian cell the box *is* the cell. On a PEBI or
tetrahedral cell it is not.

`calcDavg_U` integrates over the true cell polyhedron, exactly and without
sampling. The signed distance is affine, so over any region on one side of
the plane `∫|d| dV = V·|d(C)|`; splitting the cell at the plane gives

```
<d> = ( V⁻·|d(C⁻)| + V⁺·|d(C⁺)| ) / V
```

with `(V⁻, C⁻)` from one half-space clip and the other side by complement,
so the two parts reproduce `G.cells.volumes` and `G.cells.centroids` by
construction.

Measured against Monte-Carlo ground truth (`tests/test03`), relative error
in `<d>` — and therefore in every fracture-matrix transmissibility:

| grid | `calcDavg_U` | bounding-box `calcdavg` |
|---|---|---|
| Cartesian (control) | 3e-14 | 4.6e-3 |
| extruded PEBI | 1.0e-4 (MC noise) | **1.8e-2** median, 8.4e-2 worst |
| tetrahedral | 8.1e-4 (MC noise) | **4.4e-1** median, 7.8e-1 worst |

On tetrahedra the legacy estimator is wrong by ~44%. Even on Cartesian
cells its 10³ midpoint sampling carries ~0.5% bias, because `|d|` kinks at
the fracture plane; it converges onto `calcDavg_U` as the sampling refines
(7.7e-3 at 10³ samples, 1.2e-6 at 640³).

### 2. Face-coplanar fractures lost a neighbour

A fracture lying in the plane of a matrix face borders **two** cells. The
Cartesian gridder registers such intersections for both; the unstructured
gridder in `CPEDFM-U` keeps only the first neighbour, after which `adjustT`
— which infers the surface-area factor by counting how many rows share a
fracture cell — sees one row and doubles the transmissibility, which is the
rule for an *interior* fracture, not a boundary one.

`EDFMgrid_U` creates one fracture cell per coplanar face and registers it
against both neighbours, recording the mapping explicitly in
`matrix_connection.fraccells`. `fracMatNNC_U` then takes the surface-area
factor from the recorded intersection **type** rather than by counting
rows, so a fracture that is face-coplanar in one part of the domain and
interior elsewhere is handled correctly — a case the Cartesian code cannot
express at all (it asserts `Nm == Nf || Nm == 2*Nf`).

### 3. Fracture-fracture NNCs were Cartesian-only

`fracturefractureShaleNNCs3D` builds its candidate filter with
`hfm/findneighbours`, which opens with `assert(isfield(G,'cartDims'))` and
constructs a 26-cell 3×3×3 cloud by logical indexing. It throws on a PEBI
grid. Its inner routine also assumes row `i` of `matrix_connection.cells`
belongs to fracture cell `mod(i,Nf)`, which conflicts with the explicit
`fraccells` map needed for (2).

`fracturefractureNNC_U` replaces those two grid-dependent pieces — the
cloud by the node-sharing adjacency of `cellNeighbourhood_U`, the implicit
ordering by `fraccells` — and keeps the geometry identical to the Cartesian
release (`PEBIPEBIintersection3D` plus the Moinfar harmonic
transmissibility). On a Cartesian grid it reproduces the original to 1.2e-4
(`tests/test02` check F).

### Also

* Fracture-matrix transmissibilities resolve an anisotropic matrix
  permeability **along the fracture normal** (`permNormal_U`, supporting
  1/3/6/9-component perm), since flow into a fracture is normal to its
  plane. For isotropic rock this reduces exactly to the Cartesian formula.
* `G.nnc` carries diagnostic `CI` and `davg` columns, kept aligned with
  `G.nnc.cells` through the frac-frac stage.

## Files

```
startupEDFMU.m                    paths (reuses hfm + CPEDFM operator stage)
UtilsEDFM/
  clipPolyByCell_U.m              fracture polygon x convex cell (Sutherland-Hodgman)
  clipCellByPlane_U.m             convex cell x half-space -> volume + centroid
  calcDavg_U.m                    exact <d> over the true cell polyhedron
  EDFMgrid_U.m                    fracture gridding, incl. face-coplanar handling
  fracMatNNC_U.m                  per-plane CI = A_surface / <d>
  fracturematrixNNC_U.m           loops planes, assembles G.nnc
  fracturefractureNNC_U.m         unstructured fracture-fracture NNCs
  cellNeighbourhood_U.m           node-sharing cell adjacency
  permNormal_U.m                  n'*K*n for 1/3/6/9-component perm
  polyArea3D_U.m, dedupPoly_U.m, orderPlanarLoop_U.m    geometry helpers
tests/
  test01_geometry_U.m             geometry kernel vs analytic
  test02_cartesian_crosscheck_U.m chain-back vs the validated Cartesian code
  test03_davg_accuracy_U.m        <d> vs Monte-Carlo on PEBI and tetrahedra
  test04_pebi_dfm_validation_U.m  flow vs a fracture-conforming DFM reference
  boxMeanAbsDist_U.m              exact box integral (test reference only)
examples/
  examplePebiFracturedWaterflood_U.m
```

## Test status (MATLAB R2021b, all passing)

* **test01 — geometry kernel.** `<d>` matches the closed-form box integral
  to 7e-13 over seven plane positions and orientations. Clipped volume and
  centroid close against the grid's own to 1e-15 (Cartesian, PEBI) and
  3e-14 (tetrahedral). Clipped fracture areas sum to the polygon area to
  machine precision on all three grids. A face-coplanar fracture produces
  12 fracture cells and 24 connections, 12 host cells on each side.
  Fracture grid volume = area × aperture to 4e-14.

* **test02 — Cartesian chain-back.** Against `EDFMshalegrid` +
  `fracturematrixShaleNNC3D` on the same grid: identical host cells,
  intersection areas to 8.6e-16, NNC transmissibilities to 4.5e-3, matrix
  pressure fields to 6e-7. The residual in T is the *reference's* `<d>`
  quadrature error, not EDFM-U's — check D grades both against the exact
  box integral (EDFM-U 3.5e-11, reference 5.5e-3). Repeated for two
  intersecting fractures; fracture-fracture NNCs match in count and cell
  pairs, with T to 1.2e-4.

* **test03 — `<d>` accuracy.** The table above. `calcDavg_U` sits within
  2.0 Monte-Carlo standard errors on every graded cell, i.e. it is
  indistinguishable from exact at the resolution of the reference.

* **test04 — flow vs a conforming DFM reference.** Non-conforming extruded
  PEBI + EDFM-U against a fracture-conforming PEBI grid (3.4k cells) with
  matched fracture conductivity `k·w`.
  - *Conductive throughflow*: EDFM-U within **0.54%** of the reference;
    probe-line pressure rel. L2 error 5.3%.
  - *Transfer-limited drainage* (closed domain, fracture pinned at low
    pressure, comparison masked to the 7–88% decline window): EDFM-U rel.
    L2 **0.104** vs **0.118** for the bounding-box `<d>`.

## When the exact `<d>` actually changes an answer

Worth being explicit, because it is easy to oversell. Both `test04` case 1
and `examples/` are *transport-dominated*: the fractures are conductive and
matrix-to-fracture transfer is not rate-limiting, so a 7–29% error in `<d>`
moves cumulative oil by ~0% and throughflow by ~0.6%. `<d>` governs the
answer when transfer *is* the bottleneck — tight matrix, low-conductivity
or partially sealing features, early-time transients — which is what
`test04` case 2 isolates.

That residual 10% in case 2 is the intrinsic EDFM transfer error, not a
resolution artefact: refining the EDFM grid from 608 to ~1200 cells and
narrowing the reference band from 2.5 m to 1.5 m moved it only from 10.8%
to 10.4%. EDFM assumes a single average distance per host cell, so with
metre-scale cells it resolves the near-fracture gradient only to first
order.

## Limitations

* Cells must be **convex** (PEBI/Voronoi, tetrahedra, extruded prisms and
  Cartesian always are; corner-point hexahedra usually are). For strongly
  non-convex or warped-face cells the clip is an approximation of the same
  order as MRST's own geometry.
* Fracture polygons are planar per plane. Curved faults are represented as
  several `fracplanes`.
* Plain EDFM cannot represent **flow barriers** — a sealing fracture leaks
  regardless of `<d>`, because the matrix-matrix connection across it is
  never removed. Use `PEDFM-*` / `CPEDFM-*` (projection-based) for those.
  EDFM-U is the conductive-fracture method.
* `fracturefractureNNC_U` assumes each fracture grid is one layer thick
  (as produced by `fractureplanegeneralgrid`).

## Reused, unmodified

`fractureplanegeneralgrid`, `PEBIPEBIintersection3D`, `gridCellNodes` and
the geometry utilities from `hfm`; `assembleShaleGlobalGrid` and
`setupShaleEDFMOpsTPFA` from `CPEDFM/UpdateMay27_2023/hfmFns`. `startupEDFMU`
asserts that EDFM-U's own function names resolve to this module.
