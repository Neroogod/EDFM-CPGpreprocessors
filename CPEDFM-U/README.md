# CPEDFM-U

Continuous projection-based EDFM (CPEDFM) for **general unstructured
grids**: PEBI/Voronoi, tetrahedra, extruded prisms — any grid whose
cells are convex polyhedra. Companion module to `CPEDFM-CPG`
(corner-point) and the Cartesian release in `CPEDFM/UpdateMay27_2023`
(Rashid & Olorode, SPE-217469-PA). No existing file is modified.

## Usage

```matlab
run modules/CPEDFM-U/startupCPEDFMU.m

G = computeGeometry(makeLayeredGrid(pebiGrid2D(h, dims), nz)); G.rock = ...;
[G, fracplanes] = EDFMgrid_U(G, fracplanes, 'Tolerance', tol);  % from EDFM-U
G  = fracturematrixNNC_U(G, tol);                               % from EDFM-U
G  = pMatFracNNCsCPEDFM_U(G, tol);                              % projection (here)
ops = setupCPEDFMOpsTPFA_CPG(G, G.rock, tol);                   % reused from CPEDFM-CPG
```

The gridding and fracture-matrix NNC stages now come from the **EDFM-U**
module; CPEDFM-U supplies the projection stage on top. The former local
names `EDFMGridU` and `fracPolyCellIntersect` still work as deprecated
wrappers, but `fracturematrixShaleNNC3D` should no longer be used here —
see *Corrections* below.

## Method: the projection reduces to a centroid-side cut

The paper's structured-grid machinery — the Fig. 8 orientation branch
table, the Fig. 9 IJK local neighborhoods, the first/last-member
collapse — has no meaning on unstructured grids (no face tags, no
logical directions, arbitrary face counts). All of it approximates one
property: a *continuous* set of sealed faces. On unstructured grids
that set is computed directly:

> Seal every face of a fracture-hosting cell whose two neighboring
> cell centroids lie on opposite sides of the fracture surface,
> restricted to the fracture's lateral extent (tip guard).

Properties:

* The sign-flip set is by construction a **cut of the matrix
  connectivity graph**, so continuity holds provably for any cell
  shape — the property CPEDFM constructs combinatorially on Cartesian
  grids.
* For the classic single-cell case it selects the face nearest the
  fracture, recovering standard pEDFM selection but from a globally
  consistent criterion (no Jiang & Younis inconsistency, hence no
  leakage).
* The **tip guard** (face centroid projected into the fracture polygon,
  with a buffer of half the face's characteristic length) prevents
  sealing beyond the fracture's extent.
* Sealed faces get the complete face area (default
  `'areaMode','fullface'`, per Figs. 10/13 of the paper) or the clipped
  projected fracture-cell polygon (`'partialsum'`). Eq. 10
  transmissibilities use normal-projected permeabilities n'*K*n
  (supports 1-, 3- and 6-component perm).

## Corrections (2026-08-01)

Three defects were found and fixed. The first two only fired for fractures
lying **exactly** in the plane of a grid face, which is why none of
test10/11/12 — all of which use slanted fractures — caught them;
`test13_coplanar_fracture_U` now covers that path.

1. **The coplanar clip returned an untrimmed polygon.**
   `fracPolyCellIntersect` did `break` on detecting the coplanar face, so
   the polygon was only trimmed by the faces visited *before* it. With a
   fracture on the x=20 plane of a 4×4×3 grid every one of the 24 cells
   along the plane was handed the entire 936 m² polygon instead of its own
   64–100 m² footprint — 22464 m² of total intersection area against a
   true 1872 m². The duplicate polygons then made the fracture grid
   degenerate and `EDFMGridU` threw in `fractureplanegeneralgrid`, so the
   case did not run at all.

2. **Only one of the two bordering cells was registered.** `EDFMGridU`
   kept just the first neighbour, leaving the far side with no connection
   to the fracture; the lone surviving row then fooled the row-counting
   surface-area rule (`adjustT`), which applied the interior ×2 factor to
   a boundary fracture. 12 connections instead of 24, 12 matrix cells
   disconnected, every surviving T exactly 2× too large. The aggregate
   transfer capacity happened to match, so the error was *topological*:
   the fracture only talked to one side.

3. **`<d>` was integrated over the cell's bounding box.** Routing through
   `fracturematrixShaleNNC3D` → `fracmatnnc` → `hfm/calcdavg` estimates the
   average fracture-matrix distance over the cell's axis-aligned bounding
   box, which is the cell only for Cartesian grids. Since `CI = A/<d>`,
   that biased every fracture-matrix transmissibility by ~16% (median) on
   a PEBI grid and ~44% on tetrahedra. The pipeline now uses
   `fracturematrixNNC_U`, whose `calcDavg_U` integrates over the true cell
   polyhedron exactly.

Fixing (1) and (2) made the coplanar path functional and exposed a fourth
gap, also fixed: `pMatFracNNCsCPEDFM_U` projected only *interior* NNCs and
returned early when there were none, so a coplanar fracture had nothing
sealed — and an unsealed sealing fracture on a face plane is simply
bypassed by matrix-matrix flow across that same face. Such faces are now
closed directly (`sealCoplanarFaces`); no projection is needed, because the
face to close is the one the fracture occupies and the boundary NNCs to
both neighbours already carry the flow.

Nothing outside CPEDFM-U was modified: `fracturematrixShaleNNC3D`,
`fracmatnnc`, `calcdavg` and `adjustT` are shared with other modules and
are left untouched.

## Files

```
startupCPEDFMU.m                paths (reuses EDFM-U + CPEDFM-CPG stages)
UtilsCPEDFM/fracPolyCellIntersect.m    DEPRECATED -> EDFM-U/clipPolyByCell_U
UtilsCPEDFM/EDFMGridU.m                DEPRECATED -> EDFM-U/EDFMgrid_U
UtilsCPEDFM/calculateProjectionCPEDFM_U.m
                                       per-plane flip-cut projection + tip guard
                                       + Eq. 10 NNCs; holds the local helpers
                                       permNormal / dist2poly /
                                       clippedProjArea / orderCCW
UtilsCPEDFM/pMatFracNNCsCPEDFM_U.m     loops planes, assembles G.nnc, and seals
                                       face-coplanar fractures
tests/runAllTestsCPEDFMU.m             runs the whole suite
tests/test10_regression_cartesian_U.m
tests/test11_pebi_seal_negcontrol.m
tests/test12_pebi_dfm_validation.m
tests/test13_coplanar_fracture_U.m     face-coplanar path (the regression the
                                       other three did not cover)
examples/examplePebiWaterfloodCPEDFM_U.m
```

## Example (executed end to end)

`examples/examplePebiWaterfloodCPEDFM_U.m`: 10-year waterflood in a
500x500x50 m extruded-PEBI reservoir (1,992 cells) with a sealing fault
slanted in plan and dip, two conductive (10 D) natural fractures near
the injector, and three wells - injector + producer inside the western
fault block, observation producer east of the fault. Results (CPEDFM-U
vs plain EDFM on the identical grid):

* PROD-E (across the fault): 82 m^3 vs 1.29e6 m^3 cumulative oil - EDFM
  overestimates across-fault recovery by ~4 orders of magnitude.
* Eastern block end pressure: 100 bar (isolated, at producer BHP) vs
  152 bar (spurious support).
* PROD-W (inside the block): 8.2e5 vs 3.5e5 m^3 - with the fault
  sealed, the injected water sweeps the western block instead of
  bleeding across.

Figures written to `examples/output/`: fracture system + sealed-face
projection, top-layer water saturation comparison, producer-rate
comparison.

Reused verbatim from `CPEDFM-CPG/UtilsCPEDFM`: `transmultpEDFM_CPG` (face-
indexed matrix-matrix multipliers) and `setupCPEDFMOpsTPFA_CPG`
(operators) — both grid-agnostic by design.

## Test status (MATLAB R2021b, 2026-08-01, all passing in 25.5 s)

`run tests/runAllTestsCPEDFMU.m`. test10/11/12 below are unchanged by the
corrections above (their numbers moved only in the last digit, e.g. test12's
west influx 6.385e-4 -> 6.386e-4 from the exact `<d>`).

* **test13** (face-coplanar path, new): clipped areas sum to exactly
  2 × the polygon area with no cell exceeding its own face footprint;
  every fracture cell shared by exactly 2 hosts, 12 on each side of the
  plane; boundary `CI·<d>/A = 1.000000` (no spurious interior doubling);
  a full-span sealing coplanar fracture splits the matrix into 2
  components where plain EDFM leaves 1; and a mixed model with one
  coplanar and one slanted fracture assembles — a configuration the
  Cartesian implementation cannot represent at all, since `fracmatnnc`
  asserts `Nm == Nf || Nm == 2*Nf`. PASS.

## Earlier test status (MATLAB R2021b, 2026-07-10)

* **test10** (Cartesian chain-back): `EDFMGridU` reproduces the
  validated `EDFMshalegrid` intersections exactly (200/200 host cells,
  area rel. diff 2e-15); the flip-cut compartmentalizes the Fig. 11a
  fracture into 2 components vs 1 for EDFM; every sealed face straddles
  the fracture plane. PASS.
* **test11** (extruded PEBI, 1062 cells, 6-9 faces/cell): full slanted
  sealing fracture -> exactly 2 macro compartments vs 1 for EDFM;
  partial central fracture -> 1 component (tip guard); multipliers and
  face-pair consistency clean; projection figure in `tests/output/`.
  PASS.
* **test12** (DFM validation): west-east flow across a slanted fracture,
  CPEDFM-U on a non-conforming PEBI grid vs an explicit
  fracture-conforming PEBI reference (3426 cells, fracture cells
  gridded).
  - Sealing fracture (1 nD): probe-line rel. L2 error 0.000 vs 0.494
    for EDFM; cross-fracture dP 100.00 bar vs reference 100.00 (EDFM:
    58.07); residual through-flow matches the reference to < 0.1% of
    the EDFM leak (which is ~6 orders of magnitude larger).
  - Finite-conductivity barrier (0.02 mD, aperture matched to the
    reference fracture-cell width so w/k resistances are identical):
    probe-line rel. L2 error 0.001; cross-fracture dP 99.74 bar vs
    reference 99.86 - the Eq. 10 projected transmissibilities are
    quantitatively validated, not just the binary seal.
  Figure: `tests/output/test12_validation.png` (pressure fields for
  reference / CPEDFM-U / EDFM + probe-line overlay). PASS.

## Limitations / notes

* Cells must be convex (true for PEBI/Voronoi, tets, extruded prisms,
  Cartesian; corner-point hexahedra are usually convex — for strongly
  non-convex CPG cells prefer the CPEDFM-CPG module).
* Fracture polygons are assumed planar per plane (as in the paper);
  the flip test is evaluated against each plane's own geometry, so
  piecewise-planar (curved) fault surfaces can be represented as
  multiple fracplanes.
* Fractures coplanar with a cell face are registered as
  'fracmat boundary' against BOTH bordering cells and sealed directly on
  the face they occupy rather than projected (see *Corrections*). This
  departs from the Cartesian release's convention, which leaves them
  unsealed.
* Fracture-fracture intersection NNCs: use the existing
  `fracturefractureShaleNNCs3D` if intersecting fracture networks are
  needed; not exercised by the tests here.
