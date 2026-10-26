# PEDFM-CPG

3D projection-based EDFM (pEDFM) for **corner-point grids**.

Extends the Cartesian 3D pEDFM of Olorode, Wang & Rashid (2020),
*"Three-Dimensional Projection-Based Embedded Discrete-Fracture Model for
Compositional Simulation of Fractured Reservoirs"*, SPE Journal 25(04)
2143-2161, SPE-201243-PA, to general corner-point geometry: non-axis-
aligned faces, fault-split interfaces, pinch-outs, and cells with other
than six faces.

The paper restricts itself to Cartesian grids and defers this explicitly
(p. 2147): *"Although our algorithm could be adapted for simple
corner-point grids with six faces, we leave this (as well as more
complicated corner-point and unstructured grids) to a future study."*

The design, the five Cartesian assumptions that had to be replaced, and
the reformulation of the paper's Table 1 are documented in
[pedfm_implementation.md](pedfm_implementation.md).

## Usage

```matlab
run modules/PEDFM-CPG/startupPEDFMCPG.m

G = processGRDECL(grdecl);  G = computeGeometry(G);  G.rock = ...;
[G, fracplanes] = EDFMGridCPG(G, fracplanes, 'Tolerance', tol);  % existing
G   = fracturematrixShaleNNC3D(G, tol);                          % existing
G   = pMatFracNNCsPEDFM3D_CPG(G, tol);                           % NEW
ops = setupPEDFMOpsTPFA3D_CPG(G, G.rock, tol);                   % NEW
model.operators = ops;
```

Options on `pMatFracNNCsPEDFM3D_CPG`:

| Option | Default | Meaning |
|---|---|---|
| `defaultFace` | `'bottom'` | Arbitrary pick for Case 3 (all three pairs equidistant), as in the Cartesian release |
| `relTol` | `1e-6` | Equidistance tolerance, as a fraction of the host cell's `V^(1/3)` |
| `distanceMode` | `'centroid'` | `'centroid'` is the paper's rule; `'plane'` measures along the fracture normal |
| `unionArea` | `true` | Apply the Eq. 29 union when several fracture cells share a face |

## pEDFM vs CPEDFM

This module is **not** `modules/CPEDFM-CPG`. pEDFM selects projection
faces independently for each host cell. CPEDFM (Rashid & Olorode 2024)
collapses host cells into logical neighbourhoods and seals a connected
envelope. The practical consequence is visible in test01: pEDFM seals 91%
of the faces the fracture crosses but leaves a few gaps mid-fracture, so a
sealing fracture throttles flow without fully compartmentalizing the
domain. Producing a watertight cut is what the continuous variant was
introduced for. Use CPEDFM-CPG when you need guaranteed
compartmentalization; use this module when you want the published pEDFM
algorithm.

`CPEDFM-CPG/UtilsCPEDFM/pMatFracNNCsPEDFM_CPG.m` contains a *simplified*
pEDFM-on-CPG used inside that module as a comparison baseline: it selects
faces from the fracture-plane orientation table alone and never performs
the paper's distance test or tie-break. This module supersedes it. All
names here carry a `PEDFM3D_CPG` / `PEDFM_CPG` suffix so both modules can
share a path; `startupPEDFMCPG.m` asserts no shadowing.

## Dependencies

Read-only; no file outside this module is modified.

| Dependency | Used for |
|---|---|
| `hfm_v2` | `setupEDFMOperatorsTPFA`, geometry utilities |
| `CPEDFM/UpdateMay27_2023` | `fracturematrixShaleNNC3D`, `setupShaleEDFMOpsTPFA`, ADFNE, and `pMatFracNNCs3D` as the regression oracle |
| `CPEDFMcornerpointaug24/Updates_8_23` | `EDFMGridCPG`, `pebiAABBintersectCPGv5` |

**`hfm_v2`, not `hfm`.** The corner-point intersection engine calls
`convexpolygonarea` with two outputs (`pebiAABBintersectCPGv5` line 648)
and only the `hfm_v2` copy returns them; loading `hfm` instead makes
`EDFMGridCPG` fail inside its `parfor` with *"Too many output
arguments"*. `startupPEDFMCPG.m` checks the resolved signature and fails
loudly rather than letting this surface later.

## Files

```
startupPEDFMCPG.m                    path setup + anti-shadow guards
pedfm_implementation.md              design document
UtilsPEDFM/
  pMatFracNNCsPEDFM3D_CPG.m          main driver; assembles G.nnc
  selectProjectionFacesPEDFM_CPG.m   Fig. 2 flow chart: distance test,
                                     Cases 0-3, node-containment tie-break
  getLogicalFacesPEDFM_CPG.m         cell -> face sets per tag
  projectedAreaPEDFM_CPG.m           Eq. 29 union on the real face plane
  transmultPEDFM3D_CPG.m             Eq. 27 face-indexed M-M multipliers
  setupPEDFMOpsTPFA3D_CPG.m          TPFA operators
tests/
  test00_regression_cartesian.m
  test01_skewedCPG_sealingFrac.m
  test02_faults_and_ties.m
```

## Test status (MATLAB R2021b, 2026-08-01)

All three pass.

* **test00 — Cartesian regression.** 9x9x3 grid, one fracture inclined to
  all three axes (normal `[0.744, -0.496, -0.447]`), 60 interior NNCs.
  Against `pMatFracNNCs3D`: pM-M pair sets **identical** (101 pairs);
  projected areas match to **0.000e+00** on the 49 singleton pairs and
  **2.1e-15** on the 52 colocated ones; pM-F transmissibilities match as a
  multiset to **0.000e+00** (153 values); face multipliers match to
  **2.8e-16**. Case counts C0=57, C3=3, zero tie-break fallbacks.
* **test01 — skewed CPG, sealing fracture.** `simpleGrdecl` with skewed
  pillars: 1600/1600 distinct centroid x-coordinates and 3069/5337
  non-axis-aligned faces, so every Cartesian assumption is violated. 358
  projections, all Case 0, zero fallbacks. 91% of the 123 faces the
  fracture crosses are sealed (median multiplier 0.0000); cross-fracture
  pressure drop 5.05 bar vs 3.43 bar for plain EDFM.
* **test02 — faults and tie-breaking.** `makeModel3`, 516 cells of which
  216 have more than six faces; 36 of 48 host cells sit on fault-split
  sides. Cases 0, 1 and 3 all observed across the three configurations,
  with zero fallbacks. The corner invariant that validates the Table 1
  reformulation — the three selected faces meet at exactly one node —
  holds on **69/69** six-faced hosts.

## Known limitations

* **The barrier is not watertight.** test01 leaves 11 of 123 crossed faces
  open, mid-fracture. This is inherent to pEDFM's per-cell selection, not
  a defect of the port; see the pEDFM vs CPEDFM section.
* **Case 2 is untested against real geometry.** Cases 0, 1 and 3 occur in
  the test suite; a configuration producing exactly two equidistant pairs
  was not found. The branch is implemented and shares the node-containment
  helper with Case 3, but it has not been exercised.
* **Exact ties are measure-zero on warped grids.** On a genuinely
  distorted CPG the equidistance branches essentially never fire on their
  own (test01 and test02 Part A are 100% Case 0). Cases 1-3 are artifacts
  of Cartesian symmetry. test02 widens `relTol` deliberately to reach them.
* Fractures lying exactly on a cell face take the standard EDFM path
  (`'fracmat boundary'`) and are not projected — same as the Cartesian
  release. The paper notes pEDFM degenerates to DFM in that limit.
* The Case 3 default direction is a genuine arbitrary choice, as in the
  paper. It is exposed as `defaultFace` rather than hidden.

## A bug in the Cartesian reference

`pMatFracNNCs3D.m` appends the **pre-collapse** transmissibility vector to
`G.nnc.T` (line 519) but the **post-collapse** `cells`/`area`/`type`
arrays (lines 570-572). When any two host cells project onto the same
interface, `G.nnc.T` ends up longer than `G.nnc.cells` — in test00's case
213 versus 161, a 52-entry mismatch. Downstream, `setupPEDFMOpsTPFA` does
`T = [T; G.nnc.T]` and then indexes with `s.internalConn`, so the fracture
NNC transmissibilities are misaligned against the connection list. This
module appends matching lengths; test00 check F asserts it. The reference
file was not modified.
