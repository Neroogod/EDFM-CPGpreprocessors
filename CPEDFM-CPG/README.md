# CPEDFM-CPG

Continuous projection-based EDFM (CPEDFM) for **corner-point grids**.

Extends the Cartesian CPEDFM implementation of Rashid & Olorode (2024),
*"A Continuous Projection-Based EDFM Model for Flow in Fractured
Reservoirs"*, SPE Journal 29(01), SPE-217469-PA, to general corner-point
geometry: non-axis-aligned faces, fault-split interfaces, pinch-outs.

Nothing in the existing modules is modified. This module only *reads*
from:

| Dependency | Used for |
|---|---|
| `modules/CPEDFM/UpdateMay27_2023` | `fracturematrixShaleNNC3D`, regression oracle (`pMatFracNNCs3DCPEDFM`, `transmultpEDFM`), ADFNE |
| `modules/CPEDFMcornerpointaug24/Updates_8_23` | `EDFMGridCPG`, `pebiAABBintersectCPGv5` (CPG fracture gridding) |
| `modules/hfm` | `setupEDFMOperatorsTPFA`, geometry utilities |

## Usage

```matlab
run modules/CPEDFM-CPG/startupCPEDFMCPG.m

G = processGRDECL(grdecl);  G = computeGeometry(G);  G.rock = ...;
[G, fracplanes] = EDFMGridCPG(G, fracplanes, 'Tolerance', tol, ...);  % existing
G  = fracturematrixShaleNNC3D(G, tol);                                % existing
G  = pMatFracNNCsCPEDFM_CPG(G, tol);                                  % NEW
ops = setupCPEDFMOpsTPFA_CPG(G, G.rock, tol);                         % NEW
model.operators = ops;
```

## What breaks on corner-point grids, and how it is fixed

The Cartesian release (`CPEDFM/UpdateMay27_2023/CPEDFM/UtilsCPEDFM/`)
assumes:

1. **Six faces per cell, positionally ordered by tag.** CPG cells keep
   the tags (`G.cells.faces(:,2)` = 1..6) but not the ordering, faulted
   cells have 7+ faces (split logical sides), pinched cells fewer.
   -> `getLogicalFaces` returns *face sets per tag*; nothing indexes
   faces positionally.
2. **Neighborhoods from exact centroid-coordinate equality**
   (`calculateProjectionCPEDFM.m` lines 127-139). On any distorted grid
   every neighborhood silently degenerates to one cell.
   -> `buildNeighborhoods` groups host cells by logical IJK indices
   (`ind2sub(G.cartDims, G.cells.indexMap(...))`), equivalent on
   Cartesian grids, valid under distortion and fault throw.
3. **One shared face per matrix-cell pair** (`transmultpEDFM.m`
   pair-matching). False across faults.
   -> `pMatFracNNCsCPEDFM_CPG` records the global face index of every
   projection (`G.nnc.pMMfaces`); `transmultpEDFM_CPG` accumulates
   per-face directly.
4. **Axis-aligned projected-area formula.** -> areas and Eq. 10
   distances use the actual `G.faces.normals`/centroids.

## Deliberate deviations from the Cartesian release

* **Full-face sealing (default `'areaMode','fullface'`).** The paper
  (Fig. 10 "complete or continuous area", Fig. 13 "the area ... equals
  the area of the blue-colored interface between Cells 9 and 10")
  assigns the complete interface to the kept face. The Cartesian code
  instead *sums the partial projected areas*; the two coincide for
  fractures parallel to a grid axis (all of its validated flow cases)
  but the partial sum covers only a fraction of each face for fractures
  inclined to all three axes, so no face seals and the barrier leaks
  (verified: 15/134 crossing faces partially reduced, 0 sealed, domain
  stays connected). `'areaMode','partialsum'` reproduces the
  original behavior if wanted.
* **Flip-face sealing.** The paper's collapse ("first or last face in
  the neighborhood") seals the host-slab *envelope*. When the host slab
  is 1-2 cells thick (the paper's cases) the envelope is a watertight
  cut. On coarse/steep/faulted CPG columns the slab can be several cells
  thick and the envelope provably leaks through side-steps (verified on
  a 5-cell-thick column). This module additionally seals every
  host-adjacent face whose two neighbor cells lie on opposite sides of
  the fracture plane. That set is by construction a graph cut, so the
  barrier is watertight for any slab thickness; on thin slabs it adds
  faces the fracture geometrically crosses (a strict improvement).
* **Extreme-member collapse per orientation.** The original's x/y
  collapse loops keep the same member/face regardless of the orientation
  branch (and return zero merged transmissibility for positive-side
  selections). This module follows the paper: first member's face for
  negative sides, last member's for positive sides.
* Directional permeability column (`min(d, size(perm,2))`) instead of
  always column 1; identical for isotropic rock.

## Files

```
startupCPEDFMCPG.m              path setup + shadow guards
UtilsCPEDFM/calculateProjectionCPEDFM_CPG.m
                                       per-plane projection; holds the four
                                       helpers as local functions:
                                         fracOrientation     Fig. 8 branch table
                                         getLogicalFaces     cell -> face sets per tag
                                         buildNeighborhoods  IJK neighborhoods (Fig. 9)
                                         getUnionProjectArea polygon-union area
UtilsCPEDFM/pMatFracNNCsCPEDFM_CPG.m   main: loops planes, assembles G.nnc
UtilsCPEDFM/pMatFracNNCsPEDFM_CPG.m    conventional per-cell pEDFM on CPG
UtilsCPEDFM/transmultpEDFM_CPG.m       face-indexed M-M multipliers
UtilsCPEDFM/setupCPEDFMOpsTPFA_CPG.m   TPFA operators (clone of setupPEDFMOpsTPFA)
UtilsCPEDFM/makeExampleCornerPointGrid.m GRDECL generator for example grids
UtilsCPEDFM/wellfractureNNCs3D_CPG.m CPG-safe well/fracture intersection mapping
tests/test00_regression_cartesian.m
tests/test01_skewedCPG_singleFrac.m
tests/test02_multifrac_negcontrol.m
examples/exampleJohansenFaultCPEDFM_CPG.m
```

The `examples/Fig*.m` scripts adapt the complete Cartesian CPEDFM paper
suite to genuine GRDECL corner-point reservoir grids. Their embedded-model
grids are created by `makeExampleCornerPointGrid`, which writes COORD and
ZCORN, applies skewed pillars and warped layers, and calls `processGRDECL`
with preserved corner-point nodes. The original fracture systems, controls,
and figure workflows remain in place, and the scripts use `EDFMGridCPG`,
`pMatFracNNCsCPEDFM_CPG`, and `setupCPEDFMOpsTPFA_CPG` for every CPEDFM
branch. Plain EDFM/pEDFM branches and fine conforming-grid reference models
are retained where they form part of a published comparison; those reference
meshes are not the embedded reservoir grid being tested. The Eagle Ford
scripts retain `meshHFsystem` only to determine their original logical
refinement dimensions, then rebuild the simulation grid as GRDECL CPG.
Each script configures the module automatically through
`startupCPEDFMCPG.m`.

For controlled sensitivity studies, the shared constructor accepts
`pillarSkew`, `layerWarp`, `faultDrop`, and `actnum` options. Setting the
first two to zero gives an orthogonal grid represented in corner-point
format; nonzero defaults ensure that the examples exercise actual distorted
CPG geometry rather than Cartesian geometry with a different function name.

## Test status (MATLAB R2021b, 2026-07-10)

* **test00** (Cartesian regression, Fig. 11a case): all original pMM
  pairs, per-face areas (rel. diff 2e-15), multipliers, and pM-F
  transmissibilities (rel. diff 3e-20) reproduced; 200 flip-face
  additions, all host-adjacent and fully sealed. PASS.
* **test01** (skewed faulted CPG, slanted sealing fracture): matrix
  splits into >= 2 components under CPEDFM vs 1 under EDFM; the sealed
  side sustains the injected pressure while EDFM leaks (dP ratio >> 1).
  Projection figure written to `tests/output/`. PASS.
* **test02** (negative control + multi-fracture): a partial fracture
  does NOT compartmentalize (1 component); two opposing full fractures
  give exactly 3 macro compartments; fault-split faces consistent. PASS.

Known limitations / notes:

* Fractures lying exactly on cell faces are handled by the standard
  EDFM path ('fracmat boundary'), not projected - same as the Cartesian
  release.
* Host cells pinched between the flip layer and the envelope layer can
  become matrix-graph slivers; they stay coupled through their
  fracture NNCs. See test02 check B.
* With rate-controlled injection into a fully sealed compartment,
  pressure grows unbounded (test01 shows this deliberately); use BHP
  controls for realistic studies.
* `exampleJohansenFaultCPEDFM_CPG.m` requires the co2lab Johansen
  dataset and was not part of the automated test run.
