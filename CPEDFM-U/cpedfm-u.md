# CPEDFM-U: Implementation Log

A step-by-step record of how the unstructured-grid CPEDFM was designed,
implemented, debugged, and verified (2026-07-10, MATLAB R2021b).
Companion to `CPEDFM-CPG/cpedfm-cpg.md` (which documents the
corner-point module this one builds on) and to the manuscript in
`SPEJ_CPEDFM_for_CPG_and_unstructed_grids/`.

---

## Step 0 — The design insight

The corner-point work (CPEDFM-CPG) ended with a discovery that defined
this module: the **flip-face rule** — seal every host-adjacent face
whose two neighbor centroids straddle the fracture plane — added there
as a robustness fix, is by construction a cut of the matrix
connectivity graph, and is expressible with only `faces.neighbors`,
cell centroids, and the plane. The paper's Cartesian machinery
(orientation branch table, local neighborhoods, first/last-member
collapse) exists to *construct* continuity combinatorially; the flip
criterion *computes* it directly. On fully unstructured grids
(PEBI/Voronoi, tets, extruded prisms) — no face tags, no (i,j,k), 6–9+
faces per cell — the Cartesian machinery is undefined, but the flip
criterion survives unchanged. The design decision: discard the
orientation table and neighborhoods entirely and make the flip
criterion (plus a lateral tip guard) the whole selection algorithm.
Bonus property established up front: for the classical single-host-cell
fracture, the sign flip occurs across the face nearest the fracture —
recovering standard pEDFM selection, but from a globally consistent
rule.

## Step 1 — Reuse assessment and contract reading

Before writing new code, two questions were settled by reading the
existing pipeline:

- **What is grid-agnostic already?** `transmultpEDFM_CPG` (face-indexed
  multipliers) and `setupCPEDFMOpsTPFA_CPG` from CPEDFM-CPG consume
  only `G.nnc.pMMfaces/pMMarea` and the operator structures — reused
  verbatim, zero new operator code.
- **What must the gridding produce?** `fracturematrixShaleNNC3D` and
  `fracmatnnc` were read line by line: the contract is
  `F.matrix_connection` with `cells/area/type/normal/planepoint`,
  plus FracGrid bookkeeping (`cells.start`, `matrixnnc`,
  `fracgridnnc`, rock fields). Reproducing this contract exactly lets the
  standard EDFM NNC stage run on unstructured grids unmodified.

  > **This conclusion was wrong, and was corrected on 2026-08-01.** The
  > original note here claimed the connectivity-index calculation
  > (`calcdavg`/`calcfracmatCI`) was "already node-based and hence valid
  > for arbitrary polyhedra". It is node-based, but only to take
  > `min`/`max` of the nodes: `calcdavg` lays a 10×10×10 sample grid over
  > the cell's axis-aligned **bounding box**. That is the cell only on a
  > Cartesian grid. On PEBI and tetrahedral cells the box is a different
  > shape with a different centroid, so `<d>` — and through `CI = A/<d>`
  > every fracture-matrix transmissibility — was biased by ~16% (median)
  > and ~44% respectively. The pipeline now uses `fracturematrixNNC_U`
  > from EDFM-U, whose `calcDavg_U` integrates `|d|` over the true cell
  > polyhedron in closed form. Reusing an existing routine is only safe
  > once its *assumptions*, not just its interface, have been checked
  > against the new setting.

> **Steps 2 and 3 below describe the original implementation, which is
> superseded.** Both files are now deprecated wrappers over
> `EDFM-U/UtilsEDFM/clipPolyByCell_U` and `EDFM-U/UtilsEDFM/EDFMgrid_U`,
> which fix two defects in the face-coplanar path: the clip abandoned its
> loop on the coplanar face and returned an untrimmed polygon (a 12×
> area overestimate that also made the fracture grid degenerate), and the
> gridder registered such a fracture against only one of the two cells
> bordering it. See the *Corrections* section of `README.md` and
> `tests/test13_coplanar_fracture_U.m`.

## Step 2 — Geometry engine: `fracPolyCellIntersect.m`

Fracture-polygon × cell intersection by successive half-space clipping
(3D Sutherland--Hodgman): for each cell face, clip the polygon against
the plane (face centroid, outward unit normal — sign resolved from
`faces.neighbors` orientation). Exact for convex cells, which covers
every target grid family; ~100 lines versus the 995-line
hexahedra-only `pebiAABBintersectCPGv5`. Includes: detection of the
on-face (coplanar) case, returned as a boundary flag; near-duplicate
vertex cleanup; fan-triangulation area. One cleanup pass during review
removed a redundant clamped-crossing variable.

## Step 3 — Gridding driver: `EDFMGridU.m`

Modeled on `EDFMGridCPG` but general: restructure G
(Matrix/FracGrid/nnc), per-fracture loop with a bounding-box coarse
filter, per-cell clipping, fracture-plane grid via the existing
`fractureplanegeneralgrid`, `matrix_connection` fields per the
contract (normal/planepoint replicated per row), global assembly via
the existing `assembleShaleGlobalGrid`. On-face intersections are
registered once (for the first neighbor) with type 1
('fracmat boundary'), matching the pipeline convention; all others are
type 2 ('interior').

## Step 4 — Projection: `pMatFracNNCsCPEDFM_U.m`

Per fracture plane: signed distances of all matrix centroids to the
plane; for each host cell, each face with an interior neighbor is
sealed iff (a) the neighbor centroids straddle the plane and (b) the
**tip guard** holds — the face centroid, projected onto the plane,
lies inside the fracture polygon within a buffer of half the face's
characteristic length (sqrt(A_f)/2), computed in an in-plane basis
with a point-to-polygon distance helper. The guard is what keeps
partial fractures from sealing beyond their extent, and is sharper
than the host-adjacency proxy used on the CPG side. Per sealed face:
full face area (paper's Figs. 10/13 rule; `'partialsum'` mode clips
the fracture polygon against the face polygon via
`getUnionProjectAreaCPG`-style projection), Eq.-10 transmissibility
with **normal-projected permeability** n'·K·n (helper supporting 1-,
3-, and 6-component perm — required once faces are not axis-aligned),
and the deterministic pM/host pairing (non-host neighbor if one
exists; otherwise the positive-side cell). Output appends the same
`G.nnc` fields as the CPG module, with `pDir` redefined as the
dominant axis of the face normal (plot coloring only).

## Step 5 — test10: chain back to the validated stack

On a cartGrid (a valid unstructured grid), `EDFMGridU` was compared
against the existing `EDFMshalegrid`: identical host-cell set (200/200)
and intersection areas to 2e-15 — the new clipping engine reproduces
the validated geometry. The flip-cut then compartmentalizes the
Fig.-11a fracture (2 components, 250 sealed faces; EDFM stays at 1),
and every sealed face was verified to straddle the fracture plane.
Passed on the first run.

## Step 6 — test11: genuine unstructured grids

Extruded `pebiGrid2D` (upr module) + `makeLayeredGrid`: 1,062 Voronoi
prisms, 6–9 faces per cell, no tags, no logical indices. A sealing
fracture inclined to all axes and spanning the domain: exactly **2
macro compartments** vs 1 for EDFM (111 host cells, 181 sealed faces).
Negative control — a fracture covering only the central patch — stays
at **1 component** (53 sealed faces, correctly not a cut), validating
the tip guard. Multiplier/face-pair sanity checks clean; Fig.-11-style
projection figure written. Passed on the first run.

## Step 7 — test12: quantitative DFM validation

The paper-style validation (west→east flow, fixed boundary pressures)
against an explicit fracture-conforming PEBI reference
(`pebiGrid2D` with the fracture line as a cell constraint, 3,426
cells, 840 fracture cells; the 2D→3D constraint-tag replication was
asserted geometrically). Three debugging episodes are worth recording:

1. **Flux criterion**: the sealing case's influx check originally
   compared two near-zero numbers relatively (reference 2.9e-4 vs
   CPEDFM 6.4e-4 m³/d — both negligible leaks) and failed; the
   meaningful scale is the EDFM leak (1,701 m³/d). Criterion changed
   to |ΔQ| < 0.1% of the EDFM leak.
2. **Reference fracture width**: for the finite-conductivity case, the
   conforming grid's fracture cells are ~6 m wide, not
   aperture-width; for a like-for-like w/k cross-resistance the EDFM
   aperture was set to the *measured* mean width (from 3D fracture-cell
   volumes — the 2D grid object carries no geometry, the first
   attempt at `G2D.cells.volumes` errored).
3. **The processed-flag gotcha**: reusing the fracplane struct returned
   by a previous `EDFMGridU` call (`fp2 = fp`) silently skips the
   fracture (`processed=true`), and case 2 initially ran
   *fracture-free* — detected because its dP exactly matched the
   no-barrier profile. Always build fresh fracplane structs.

Final results: sealing case — probe-line rel. L2 error 0.000 vs 0.494
for EDFM, cross-fracture dP 100.00 bar vs reference 100.00 (EDFM
58.07); finite-conductivity barrier (0.02 mD) — probe error 0.001, dP
99.74 vs 99.86 bar (0.12%), quantitatively validating the Eq.-10
projected transmissibilities, not just the binary seal. The four-panel
validation figure (`tests/output/test12_validation.png`) is
manuscript-ready.

## Step 8 — Hardening inherited from the Johansen debugging

The zero-permeability/zero-distance guards added to the CPG projection
(0/0 → NaN on degenerate cells) were applied identically here:
kF+kM ≤ 0 → zero exchange, distances floored at tolerance.

## Final verified state

```
startupCPEDFMU.m               (reuses CPEDFM-CPG operator stage)
UtilsCPEDFM/fracPolyCellIntersect.m   UtilsCPEDFM/EDFMGridU.m
UtilsCPEDFM/pMatFracNNCsCPEDFM_U.m
tests/test10_regression_cartesian_U.m   PASS
tests/test11_pebi_seal_negcontrol.m     PASS
tests/test12_pebi_dfm_validation.m      PASS
examples/examplePebiWaterfloodCPEDFM_U.m (executed: 3-well waterflood,
  sealing fault + 2 conductive fractures; across-fault producer 82 m³
  vs 1.29e6 m³ under EDFM; BHP injector chosen after a rate-controlled
  variant demonstrated unbounded pressurization of the sealed block)
```

## Lessons

1. **Simplicity was the result, not the goal**: removing the
   structured-grid machinery made the unstructured algorithm shorter
   than its ancestors *and* gave it the provable continuity property
   (graph cut) the Cartesian construction only approximates.
2. **Match resistances, not just permeabilities**, when validating
   embedded models against conforming references (the measured-width
   aperture matching in test12).
3. **Fresh fracplane structs every pipeline run** — the `processed`
   flag is a silent skip.
4. The tip guard's projected-centroid-in-polygon test is the right
   general form of lateral restriction; the CPG module's
   host-adjacency proxy is the special case worth backporting if tip
   over-sealing is ever observed there.
