# CPEDFM-CPG: Implementation Log

A step-by-step record of how the corner-point-grid extension of CPEDFM
was designed, implemented, debugged, and verified (2026-07-10/11,
MATLAB R2021b). The README documents *what* the module is; this file
documents *how it came to be*, including the dead ends and the
discoveries that changed the design. Companion to the manuscript in
`SPEJ_CPEDFM_for_CPG_and_unstructed_grids/`.

---

## Step 0 — Understand the source material

Before writing any code, three bodies of material were studied:

1. **The paper** (Rashid & Olorode 2024, SPE-217469-PA,
   `knowledge/journals/CPEDFM.pdf`). Extracted the full text of the
   "CPEDFM Algorithm" section (pdfminer; the PDF's CID-encoded body
   text defeats naive extraction). Key algorithmic content: the Fig. 8
   flow chart (orientation angles dip / dip direction / normal-z select
   one face tag per direction, uniformly for all host cells of a
   plane), the Fig. 9 local neighborhoods (host cells grouped by shared
   transverse coordinates), the collapse onto the first/last interface
   of each neighborhood, and the Fig. 13 transmissibility construction
   (the projected NNC uses the **complete** area of the selected
   interface).

2. **The Cartesian implementation**
   (`modules/CPEDFM/UpdateMay27_2023/CPEDFM/UtilsCPEDFM/`):
   `pMatFracNNCs3DCPEDFM.m` (1,096 lines) and
   `calculateProjectionCPEDFM.m` (819 lines) were read in full. This
   established the exact data contract consumed downstream
   (`G.nnc.cells/T/type/area/normal/planepoint`, `pMMneighs`, and the
   multiplier stage `transmultpEDFM.m` / `setupPEDFMOpsTPFA.m`), the
   ADFNE-based orientation computation, the 8+2-branch face table, and
   the per-branch zeroing patterns of the collapse loops.

3. **The existing corner-point modules**
   (`CPEDFMcornerpointaug24`, `CPEDFMcornerpoint-main`). Finding: they
   contain the *geometry* half only — `EDFMGridCPG.m` +
   `pebiAABBintersectCPGv5.m` grid fractures on corner-point cells and
   store the same `matrix_connection` fields as the Cartesian gridding
   — but their case scripts run plain EDFM operators and handle
   sealing faults by hand-zeroing face transmissibilities
   (`JohansenTEDFMNoFaultCase.m`). The automated continuous projection
   did not exist on corner-point grids. That gap defined this module.

## Step 1 — Identify the Cartesian assumptions that break

Line-level analysis of `calculateProjectionCPEDFM.m` identified four
hard assumptions:

1. **Positional face indexing** — the branch table sets
   `faceSelected(t) = 1` assuming face *position* t of a cell is its
   logical face t (left/right/front/back/bottom/top).
2. **Neighborhoods by exact centroid-coordinate equality** (lines
   127–139: `localcoord(:,1)==localXCord(i) & ...`). On any distorted
   grid, no two centroids share coordinates, so every neighborhood
   silently degenerates to one cell.
3. **Axis-aligned projected areas** (direction cosines against x/y/z).
4. **One face per cell pair** — `transmultpEDFM.m` recovers faces by
   `ismember` on cell-pair rows, ill-posed across fault-split
   interfaces.

Two further, subtler limitations surfaced during testing (Steps 6–7).

## Step 2 — Probe the actual corner-point data structures

A probe script on `simpleGrdecl` grids established the facts the design
had to respect:

- `processGRDECL` **preserves logical tags** in `G.cells.faces(:,2)`
  (1..6), but cell faces are **not stored in tag order** even on
  unfaulted grids (a probed cell listed tags 2,1,3,4,5,6).
- Faulted cells carry **7–8 faces** (e.g., three tag-1 faces on one
  cell); pinch-outs remove faces.
- `G.cells.indexMap` + `G.cartDims` give logical (i,j,k) for every
  cell, unaffected by geometric distortion.

Conclusion: make all continuity-critical decisions in logical space
(tags + ijk), evaluate all areas/distances/permeabilities from true
geometry.

## Step 3 — Scaffold the module

Created `modules/CPEDFM-CPG/` with `UtilsCPEDFM/`, `tests/`, and
`examples/`, matching the Cartesian CPEDFM module layout while retaining
the CPG-specific regression tests.
Hard constraint honored throughout: **no existing file in any other
module is modified**. `startupCPEDFMCPG.m` wires the dependencies
(`CPEDFM/UpdateMay27_2023` for the NNC stage and ADFNE;
`CPEDFMcornerpointaug24/Updates_8_23` for `EDFMGridCPG`; `hfm`) and
asserts against wrong shadow copies of duplicated function names
(several functions exist in 2–3 places across the repo).

A useful discovery here: `CPEDFM/UpdateMay27_2023` is self-contained
(own `hfmFns/`, `pEDFM_nnc/`, `ADFNE15/`), so the missing `shale`
module was not needed.

## Step 4 — Implement the core functions

Written in dependency order:

All four helpers below are local functions of
**`calculateProjectionCPEDFM_CPG.m`** (mirroring how the Cartesian release
keeps its helpers inside `calculateProjectionCPEDFM.m`). Callers that need
them directly fetch handles with
`h = calculateProjectionCPEDFM_CPG('helpers')`.

- **`getLogicalFaces`** — cell → six *face sets* keyed by tag
  (multi-face sets across faults, empty sets at pinch-outs). No caller
  ever indexes faces positionally.
- **`fracOrientation`** — dip / dip direction / signed
  inclination computed from the plane normal *stored by the gridding*
  (`G.nnc.normal`, keeping one sign convention through the pipeline),
  plus the paper's branch table ported verbatim → tag triple
  (t1,t2,t3).
- **`buildNeighborhoods`** — neighborhoods from logical ijk:
  direction-d groups share the two transverse indices, members sorted
  by the direction-d index. Equivalent to the coordinate rule on
  Cartesian grids; valid under distortion; keeps cells across a fault
  throw in the same column.
- **`pMatFracNNCsCPEDFM_CPG.m`** — the driver: per fracture plane,
  select tags, compute per-host partial areas
  (A_fm·|cos(faceNormal, fracNormal)| with area-weighted set normals),
  form neighborhoods, seal interfaces, emit projected NNCs
  (Eq.-10-style T = k_harm·A/dist with true centroids), and append the
  face-indexed records `pMMfaces`/`pMMarea` that everything downstream
  consumes.
- **`transmultpEDFM_CPG.m`** — multipliers by direct per-face
  accumulation of `pMMarea` (accumarray + clamp), replacing pair
  matching; composes across multiple fractures and fault-split faces.
- **`setupCPEDFMOpsTPFA_CPG.m`** — 1:1 clone of `setupPEDFMOpsTPFA`
  with the multiplier call swapped.

## Step 5 — Regression against the Cartesian oracle (test00)

`test00_regression_cartesian.m` runs the paper's Fig. 11a case through
the **original** `pMatFracNNCs3DCPEDFM` and the new code on the same
gridding, and compares: neighborhoods (ijk vs coordinate grouping),
selected face tags (vs the ADFNE-based classification), sealed
cell-pair sets, per-face projected areas, multipliers, and projected
transmissibilities. First run: **exact agreement** — 245/245 pairs,
areas to 2e-15, multipliers to 0, T to 3e-20. (A cartGrid is a valid
corner-point grid, so this chains the CPG code directly to the
published implementation.)

## Step 6 — Discovery #1: partial-sum areas do not seal (test01)

`test01_skewedCPG_singleFrac.m` (skewed *faulted* `simpleGrdecl`,
fracture inclined to **all three** axes, graph-cut + flow checks)
initially FAILED: only a 1.5× pressure-jump ratio vs EDFM and one
connected component. A diagnostic script showed **0 of 134
fracture-crossing faces sealed**; projected areas covered as little as
13% of their faces.

Root cause: the Cartesian implementation assigns each collapsed
interface the *sum of partial* projected areas — which equals the full
face area only for fractures parallel to a grid axis (all of the
paper's validated flow cases). For triple-inclined fractures no
multiplier reaches zero. The paper's own Figs. 10/13, however, assign
the **complete interface area** ("the area … equals the area of the
blue-colored interface between Cells 9 and 10").

Fix: default `'areaMode','fullface'` (complete area per kept face →
multiplier exactly zero), with `'partialsum'` retained as an option
reproducing the Cartesian numerics. test00 still passes under the new
default because sum == full face for its vertical fracture. test01 then
passed: 2 components, 5.8× dP ratio.

## Step 7 — Discovery #2: the envelope collapse leaks on thick host
slabs (test02)

`test02_multifrac_negcontrol.m` (negative control + two opposing
fractures) FAILED its two-fracture case: 1 component instead of 3 —
even a single fracture that sealed in test01 failed on this
configuration. A shortest-path trace through the transmissibility graph
pinpointed the leak exactly: a (j,k) column with **five** host cells
(I = 5..9), the plane crossing between I = 6 and 7 (open, interior to
the slab), and the collapse sealing only the *envelope* face at
I = 9|10 — flow crossed mid-slab and escaped sideways through a
neighboring column with a thinner slab.

Conclusion: the paper's "first/last interface of the neighborhood"
rule seals the host-slab envelope, which coincides with the fracture
surface only for 1–2-cell-thick slabs (fine Cartesian grids). On
coarse/skewed/faulted corner-point grids the slab is locally thick and
the envelope provably leaks.

Fix: the **flip-face rule** — additionally seal every direction-d face
of a host cell whose two neighbor centroids lie on opposite sides of
the fracture plane. The flip set is, by construction, a cut of the
matrix connectivity graph within the fracture's lateral extent, so the
barrier is watertight for any slab thickness. Envelope faces are kept
(paper compliance + regression continuity); duplicates de-duplicated
per face; the pM/host pairing for both-host flip faces follows the
grid-direction side convention.

Follow-on adjustments:
- test00 check C became a **superset check**: all original pairs
  reproduced, plus 200 flip additions (every one host-adjacent, every
  one fully sealed) — a documented robustness improvement, not a
  regression.
- test02's compartment counts moved to **macro components** (>5% of
  cells): host cells pinched between the flip and envelope layers form
  tiny matrix-graph slivers that remain coupled through their fracture
  NNCs, which the matrix-only diagnostic ignores.

Result: full suite green — test00 (regression), test01 (2 components,
sealed), test02 (negative control 1 macro; two fractures → 3 macro;
fault-split faces consistent).

## Step 8 — Application example and robustness guards

`examples/exampleJohansenFaultCPEDFM_CPG.m` recreates the group's
Johansen CO2 fault case with the manual face-zeroing replaced by the
automatic projection. Running it surfaced two hardening needs:

- **Johansen zero-pore-volume pinched cells** make the standard EDFM
  NNC transmissibilities evaluate 0/0 → NaN, tripping
  `assert(all(isfinite(T)))` in ad-core. Same quirk the original TEDFM
  script patches; the example now sanitizes non-finite NNC values
  (reporting the count — 6) and the core projection functions gained
  guards: zero-permeability pairs → zero exchange (not 0/0), distances
  floored at tolerance.
- `plotToolbar` guarded behind `usejava('desktop')` for headless runs.

End-to-end run: 8/8 CO2 timesteps, ~75 s.

## Step 9 — The latent geometry bug (found later, fixed everywhere)

While building the multiscale layer (test32), `addWell` failed with a
negative Peaceman well index on east-edge cells. Diagnosis: the
`simpleGrdecl`-based **test grids** scaled `COORD` x/y and `ZCORN` to
field dimensions but never `COORD`'s z-components — pillars spanned
[0,1] m while corners sat at up to 40 m, and the 40× extrapolation
along tilted pillars **inverted cells** (negative volumes) at the east
taper. The bug never affected earlier conclusions (no check touched
those volumes) but was wrong.

Fix: `grdecl.COORD(3:3:end) = grdecl.COORD(3:3:end) * physdim(3)` in
all five affected tests (01, 02, 20, 30, 32; a double-application in
test01 from the batch edit was caught and removed). All suites
re-verified on the corrected geometry. Side effects that *strengthened*
the tests: test02's corrected geometry produces the benign all-host
slivers in its negative control too (criterion updated to macro
components, matching case B), and the Eclipse-export test began
exercising the previously-dormant EDITNNC path (47 fault-juxtaposed
seal records).

## Step 10 — Final verified state

Files:

```
startupCPEDFMCPG.m
UtilsCPEDFM/calculateProjectionCPEDFM_CPG.m   (holds the four helpers)
UtilsCPEDFM/pMatFracNNCsCPEDFM_CPG.m    UtilsCPEDFM/transmultpEDFM_CPG.m
UtilsCPEDFM/setupCPEDFMOpsTPFA_CPG.m
tests/test00_regression_cartesian.m
tests/test01_skewedCPG_singleFrac.m
tests/test02_multifrac_negcontrol.m
examples/exampleJohansenFaultCPEDFM_CPG.m
```

The Cartesian release's `Fig*.m` paper examples and their shared data/helper
files are also replicated under `examples/`. CPEDFM branches use the CPG
gridding, projection, and operator functions; deliberately comparative EDFM
and pEDFM branches continue to use their original formulations.

Results (corrected geometry, MATLAB R2021b):

| Test | Outcome |
|---|---|
| test00 | Original pairs/areas/multipliers/T reproduced exactly; +flip additions, all host-adjacent, all sealed |
| test01 | Matrix splits into compartments under CPEDFM only; sealed side holds injected pressure; projection figure written |
| test02 | Partial fracture: 1 macro component (no over-seal); two opposing fractures: 3 macro compartments; fault-split faces consistent |
| Johansen example | Runs end to end (6 NaN NNCs sanitized — known data quirk) |

Deliberate deviations from the Cartesian release (all documented in
the README): full-face areas by default; flip-face augmentation;
orientation-consistent extreme-member collapse (the original's x/y
loops keep a fixed member/face regardless of branch and return zero
merged T for positive sides); directional permeability column.

## Lessons that shaped the later modules

1. **Face-indexed projection records** (`pMMfaces`/`pMMarea`) —
   adopted for multiplier robustness — turned out to be the enabling
   design for everything downstream: the unstructured variant, the
   Eclipse/CMG preprocessor, and the barrier-aware multiscale
   partition all consume the sealed-face list directly.
2. **The flip-face criterion generalizes**: on unstructured grids it
   *replaces* the orientation table and neighborhoods entirely
   (CPEDFM-U), with continuity provable as a graph cut.
3. **Test the negative control**: the over-sealing guard (partial
   fracture must not compartmentalize) caught design errors that
   sealing tests alone never would.
4. **Regression against the published implementation first**: chaining
   every generalization back to the Cartesian oracle on a cartGrid
   made it possible to distinguish deliberate deviations from bugs.
