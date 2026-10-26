# Implementing 3D pEDFM on Corner-Point Grids

Design document for `modules/PEDFM-CPG`.

Reference: Olorode, Wang & Rashid (2020), *"Three-Dimensional
Projection-Based Embedded Discrete-Fracture Model for Compositional
Simulation of Fractured Reservoirs"*, SPE Journal 25(04) 2143-2161,
SPE-201243-PA. Cited below as **the paper**; equation and figure numbers
refer to it.

---

## 1. What this module is, and what it is not

The paper presents the first fully 3D pEDFM algorithm, but restricts it
to **structured Cartesian grids**. It says so explicitly (p. 2147):

> "Although our algorithm could be adapted for simple corner-point grids
> with six faces, we leave this (as well as more complicated corner-point
> and unstructured grids) to a future study."

This module is that future study for corner-point geometry. It implements
the paper's algorithm — Fig. 2's flow chart, Cases 0-3, Table 1's face
selection, Eqs. 24-29 — on general CPG: non-axis-aligned faces,
fault-split interfaces, pinch-outs, and cells with more or fewer than six
faces.

**pEDFM is not CPEDFM.** `modules/CPEDFM-CPG` implements the *continuous*
projection variant of Rashid & Olorode (2024), which collapses host cells
into logical neighborhoods and seals a connected envelope. pEDFM makes an
**independent per-host-cell face selection** with no neighborhood
collapse and no continuity completion. The two produce different
projection sets by design; this module deliberately does not borrow
CPEDFM's neighborhood or flip-face machinery.

`CPEDFM-CPG/UtilsCPEDFM/pMatFracNNCsPEDFM_CPG.m` already provides a
*simplified* pEDFM-on-CPG used inside that module as a comparison
baseline. It selects faces purely from the fracture-plane orientation
table and never performs the paper's distance test or equidistance
tie-breaking. It is not a faithful pEDFM and is not reused here; this
module supersedes it. All exported names carry a `PEDFM3D_CPG` suffix so
that both modules can sit on the path simultaneously without shadowing
(the startup script asserts this).

---

## 2. The algorithm, and the five places Cartesian assumptions are baked in

The reference Cartesian implementation is
`modules/CPEDFM/UpdateMay27_2023/pEDFM_nnc/pMatFracNNCs3D.m`. Reading it
against `processGRDECL` output, five assumptions fail on CPG. Each drives
one design decision below.

| # | Cartesian assumption | Where | Why it fails on CPG |
|---|---|---|---|
| 1 | Every cell has exactly 6 faces, ordered L,R,F,B,Bot,Top | `faceCentroids(1,:)`…`(6,:)`, lines 44-49 | Faulted cells have 7+ faces (split logical sides); pinched cells fewer. CPG does not keep `cells.faces` sorted by tag. |
| 2 | `G.faces.neighbors(f,1)` is the "negative-side" cell | `projectionCells`, "row1 is Left,Me" comment, line 38 | MRST orders `faces.neighbors` by face-normal direction, not by tag. Holds incidentally on `cartGrid`, not on CPG. |
| 3 | A face normal has exactly one nonzero component | `projDir = find(G.faces.normals(...))`, line 599 | CPG normals are non-axis-aligned; `find` returns 2-3 indices and `p.Dir` becomes garbage or errors. |
| 4 | One shared face per matrix-cell pair | `ismember(G.faces.neighbors, pMMneighs, 'rows')`, line 585 | Across a fault an interface is split into several faces; the lookup returns the wrong face or several. |
| 5 | Projection area = polygon area on a global YZ/XZ/XY plane | `getUnionProjectArea`, lines 652-660 | Correct plane is the actual (tilted) face plane. |

### 2.1 Face sets instead of positional indices

Everything is driven by *tag sets*, never by position:

```matlab
faceSets = getLogicalFaces(G, c);   % 1x6 cell, faceSets{t} = faces with tag t
```

Tags come from `G.cells.faces(:,2)` (1=left, 2=right, 3=front, 4=back,
5=bottom, 6=top), which `processGRDECL` populates. On a Cartesian grid
every set is a singleton and the code reduces exactly to the original. On
CPG a set may hold several faces (fault split) or be empty (pinch-out).
Direction `d` owns the tag pair `(2d-1, 2d)`.

This kills assumption 1. Assumption 2 dies with it: the projection cell
is resolved as *"the neighbor of face f that is not the host"*, never by
column index:

```matlab
nb   = G.faces.neighbors(f,:);
proj = nb(nb ~= host);
```

Assumption 3 dies too — the direction `d` is already known from the tag,
so `p.Dir` is never recovered from the normal.

### 2.2 The distance test (Cases 0-3)

The paper compares the distance from the fracture-cell centroid to the
two opposing face centroids in each direction. With tag *sets* the
representative point is the **area-weighted centroid** of the set:

```
x_tag = Σ_f (A_f · x_f) / Σ_f A_f
```

which is the true centroid of a fault-split interface and reduces to the
single face centroid on Cartesian grids. Then

```
dNeg = ||x_{2d-1} − x_frac||,   dPos = ||x_{2d} − x_frac||
```

- one set empty → the other side is selected, no tie;
- both empty → direction degenerate, skipped;
- otherwise the paper's comparison, with the tolerance made **relative**.

The original uses an absolute test, `(distLeft + tol) < distRight`. With
`tol = 1e-5` m and field-scale cells that is fine, but on a CPG with cells
spanning orders of magnitude in size an absolute tolerance either never
or always fires. Here the comparison is

```
equidistant  ⟺  |dNeg − dPos| ≤ max(relTol · L_c, tol)
```

with `L_c = V_c^{1/3}` the cell's characteristic length and
`relTol = 1e-6` by default. On uniform grids this reproduces the original
classification.

`Case 0` (no ties) needs nothing further. Cases 1-3 need the tie-break.

### 2.3 Table 1 is a node-containment rule — this is the key step

Cases 1-3 resolve ties through Table 1, which maps *(direction of the
selected face, index of the node of that face farthest from the fracture
plane)* to the two remaining faces. The Cartesian code hard-codes it as a
36-branch `if/elseif` ladder (lines 199-337), and it depends entirely on
the node-ordering scheme of Fig. 7 — an ordering `processGRDECL` does not
guarantee.

Reading Table 1 as geometry rather than as a lookup:

| Selected dir | Node idx | Remaining faces |
|---|---|---|
| X (left/right) | 1 | Front, Bottom |
| X | 2 | Back, Bottom |
| X | 3 | Back, Top |
| X | 4 | Front, Top |
| Y (front/back) | 1 | Left, Bottom |
| Y | 2 | Left, Top |
| … | … | … |

With the Fig. 7 node ordering, node 1 of the left face is the
front-bottom-left corner; the faces listed for it are Front and Bottom.
Node 2 is back-bottom-left; the faces listed are Back and Bottom. In
every row, **the two selected faces are exactly the two other faces of
the cell that contain the chosen node.** Table 1 is a hexahedron
incidence table written out longhand.

So the CPG-safe formulation of Cases 1-3 is:

> Pick the node `n*` (of the already-selected face, or of the shared edge
> in Case 1) that is farthest from the fracture plane. In each still-
> unresolved direction, select the tag whose face set contains `n*`.

This is *identical* to Table 1 on a hexahedron, requires no node ordering,
no fixed face count, and extends unchanged to fault-split sides (search
both faces of a split tag for `n*`). It also preserves the property the
table exists for — the chosen faces do not intersect the fracture path,
because they meet at the corner pushed farthest from the fracture plane.

Concretely:

- **Case 1 (one tie).** Two faces already selected. Their shared nodes
  (an edge, 2 nodes on a hexahedron) are candidates; take the one with
  max `|pointplanedistance|` to the fracture plane; the third face is the
  face of the tied direction containing it.
- **Case 2 (two ties).** One face selected. Take its own nodes, pick the
  farthest, resolve both tied directions by containment.
- **Case 3 (three ties).** No face selected. Pick a default direction
  (`'defaultFace'`, default `'bottom'`, matching the Cartesian release),
  select the nearer/either face there, then run Case 2.

**Fallback.** On a degenerate CPG cell the containment search can come up
empty (a fault-split side may not touch `n*` at all; a pinched side has
no face). The code then falls back to the tag whose area-weighted
centroid lies farther from the fracture plane — the same "push away from
the fracture" intent, decided on centroids instead of nodes. Every
fallback is counted and reported in `G.nnc.pStats.nFallback` so a user can
tell how much of a run rested on it.

### 2.4 Areas and transmissibilities (Eqs. 24-29)

Projected area onto a selected face (assumption 5 fixed by using the real
face normal):

```
A_proj = A_fm · |n_f · n_frac| / (|n_f| · |n_frac|)
```

`A_fm` is the fracture-matrix contact area already stored in
`G.nnc.area`. For a fault-split tag set, `A_proj` is distributed over the
member faces in proportion to their areas. `A_proj` is clamped to the face
area so the multiplier of Eq. 27 cannot go negative.

Projection matrix-fracture transmissibility, Eqs. 24-25:

```
k_pMF = 2 k_pM k_f / (k_pM + k_f)          (harmonic; 0 if both are 0)
d_pMF = ||x_projCell − x_fracCell||         (floored at tol)
T_pMF = k_pMF · A_proj / d_pMF
```

Permeability is read from column `min(d, size(perm,2))` so anisotropic
rock uses the directional component; identical to column 1 for isotropic
rock. The Cartesian original always uses column 1 and computes the
harmonic mean with `length(pMFcells)/sum(1./perm)`, which divides by zero
for a sealing fracture with `k_f = 0`; the guarded form above returns 0
instead, keeping the seal without an `Inf`.

Projection matrix-matrix multiplier, Eq. 27:

```
T_pMM = T_MM · (A_ii − A_proj) / A_ii
```

implemented as a per-face multiplier. Assumption 4 is fixed by recording
the **global face index** of every projection in `G.nnc.pMMfaces` at
generation time, so the multiplier is a direct `accumarray` over faces
rather than a cell-pair lookup. This is the same fix CPEDFM-CPG applies,
and it is what makes fault-split interfaces work.

Eq. 29 (multiple fracture segments in one host cell) takes the **union**
of projected areas, not the sum. The union is computed by projecting each
fracture cell's mid-plane polygon onto the actual face plane using an
orthonormal basis built from the face normal, then `union` of
`polyshape`s. Falls back to the summed area if `polyshape` fails
(non-simple polygon) — flagged in `pStats.nUnionFallback`.

---

## 3. Module layout

```
PEDFM-CPG/
  startupPEDFMCPG.m                     path setup + anti-shadow guards
  pedfm_implementation.md               this document
  README.md                             usage
  UtilsPEDFM/
    pMatFracNNCsPEDFM3D_CPG.m           main driver: loops interior fracture
                                        NNCs, assembles G.nnc
    selectProjectionFacesPEDFM_CPG.m    Fig. 2 flow chart: distance test,
                                        Cases 0-3, node-containment tie-break
    getLogicalFacesPEDFM_CPG.m          cell -> face sets per tag
    projectedAreaPEDFM_CPG.m            Eq. 24 area + Eq. 29 union
    transmultPEDFM3D_CPG.m              Eq. 27 face-indexed M-M multipliers
    setupPEDFMOpsTPFA3D_CPG.m           TPFA operators
  tests/
    test00_regression_cartesian.m       must reproduce pMatFracNNCs3D exactly
    test01_skewedCPG_sealingFrac.m      distorted CPG, sealing fracture
    test02_faulted_equidistant.m        fault splits + forced Cases 1-3
  examples/
    exampleTeneValidation_CPG.m         paper Figs. 9-13 on a CPG
```

Pipeline:

```matlab
run modules/PEDFM-CPG/startupPEDFMCPG.m

G = processGRDECL(grdecl);  G = computeGeometry(G);  G.rock = ...;
[G, fracplanes] = EDFMGridCPG(G, fracplanes, 'Tolerance', tol);  % existing
G   = fracturematrixShaleNNC3D(G, tol);                          % existing
G   = pMatFracNNCsPEDFM3D_CPG(G, tol);                           % this module
ops = setupPEDFMOpsTPFA3D_CPG(G, G.rock, tol);                   % this module
model.operators = ops;
```

Dependencies are read-only: `hfm` (`setupEDFMOperatorsTPFA`,
`pointplanedistance`), `CPEDFM/UpdateMay27_2023` (`fracturematrixShaleNNC3D`,
and `pMatFracNNCs3D` as the regression oracle),
`CPEDFMcornerpointaug24/Updates_8_23` (`EDFMGridCPG`). No file outside
this module is modified.

## 4. Data written to `G.nnc`

| Field | Meaning |
|---|---|
| `cells`, `T`, `type`, `area` | appended pM-F NNCs, `type = 'fracpmat'` |
| `pMMneighs` | `[hostCell projCell]` in side order (neg side first) |
| `pMMfaces` | global face index carrying each projection |
| `pMMarea` | projected area on that face (Eq. 24, union-corrected) |
| `pDir` | logical direction 1..3 |
| `pFracPlane` | fracture-grid number |
| `pStats` | `nCase0..nCase3`, `nFallback`, `nUnionFallback`, `nDropped` |

`pStats` exists because the tie-break fallbacks are the parts of this
port most likely to be wrong on a grid nobody has tried yet; a run that
silently took 4 000 fallbacks should be visible without instrumenting the
code.

## 5. Validation — results

All three tests pass on MATLAB R2021b (2026-08-01). Numbers below are
measured, not planned.

1. **Cartesian regression (test00).** The strongest available check: on
   `cartGrid` this module must reproduce `pMatFracNNCs3D`. A single
   fracture inclined to all three axes (normal `[0.744, -0.496, -0.447]`,
   60 interior NNCs) gives **identical pM-M pair sets** (101 pairs),
   areas matching to **0.000e+00** on singleton pairs and **2.1e-15** on
   colocated ones, transmissibilities matching as a multiset to
   **0.000e+00**, and face multipliers to **2.8e-16**. C0=57, C3=3, zero
   fallbacks.

   Two differences are expected and asserted rather than ignored. The
   reference emits zero-area NNCs where a projection direction is
   orthogonal to the fracture normal; this module skips them, and the test
   asserts every reference-only pair has reference area exactly 0. And the
   reference collapses colocated projections to one NNC per cell pair
   while this module keeps one per fracture cell (section 2.4).

2. **Distorted CPG (test01).** `simpleGrdecl` with skewed pillars:
   1600/1600 distinct centroid x-coordinates, 3069/5337 non-axis-aligned
   faces. 358 projections, **all Case 0, zero fallbacks**. 91% of the 123
   faces the fracture crosses are sealed (median multiplier 0.0000);
   cross-fracture pressure drop 5.05 bar vs 3.43 bar for EDFM.

   The original expectation here — full compartmentalization — was wrong
   for pEDFM and was corrected. 11 crossed faces stay open mid-fracture
   because per-host-cell selection does not produce a connected surface.
   That is the limitation CPEDFM was later introduced to remove, so the
   test asserts a materially larger pressure drop than EDFM and records
   the residual leak explicitly instead.

3. **Faults and ties (test02).** `makeModel3`: 516 cells, 216 with more
   than six faces; 36 of 48 hosts on fault-split sides. Cases 0, 1 and 3
   observed, zero fallbacks. The **corner invariant holds 69/69**: on
   every six-faced host the three selected faces meet at exactly one node,
   which is the property Table 1 encodes (section 2.3). This validates the
   node-containment reformulation without reference to any node ordering.

   Two findings worth recording. On a genuinely warped CPG exact
   equidistance is measure-zero, so Cases 1-3 essentially never fire on
   their own — test01 and test02 Part A are 100% Case 0. The tie-break is
   reached by widening `relTol` deliberately. And no configuration
   producing exactly **two** equidistant pairs (Case 2) was found; that
   branch remains unexercised.

Not done: the Hajibeygi et al. (2011) crossing-fracture validation of the
paper's Figs. 9-13 on a corner-point representation, which would compare
`k_f = 10 nd` (barrier) and `k_f = 1000 d` (conduit) against a fine
explicit-fracture reference. The `examples/` directory is empty.

## 6. Known limitations

- **The barrier is not watertight** (test01: 11 of 123 crossed faces stay
  open). Inherent to pEDFM's per-cell selection; use CPEDFM-CPG when
  guaranteed compartmentalization is needed.
- **Case 2 is implemented but unexercised** — no test configuration
  produced exactly two equidistant pairs. It shares the node-containment
  helper with Case 3, which is covered.
- Fractures lying exactly on a cell face are handled by the standard EDFM
  path (`'fracmat boundary'`), not projected — same as the Cartesian
  release. The paper notes pEDFM degenerates to DFM in that limit.
- The tie-break needs cell nodes; fully degenerate cells (zero-volume
  pinch-outs) are skipped and counted in `pStats.nDropped`.
- Case 3's default direction is a genuine arbitrary choice, as in the
  paper ("we randomly select one out of the cell faces"). It is exposed
  as an option rather than hidden.
- Projections onto domain-boundary faces are dropped (no projection cell
  exists), matching `'ProjMatrix outside'` in the original.

## 7. A bug found in the Cartesian reference

`pMatFracNNCs3D.m` appends the **pre-collapse** transmissibility vector to
`G.nnc.T` (line 519) but the **post-collapse** `cells`/`area`/`type`
arrays (lines 570-572). Whenever two host cells project onto the same
interface, `G.nnc.T` ends up longer than `G.nnc.cells` — 213 versus 161 in
test00's case, a 52-entry mismatch. `setupPEDFMOpsTPFA` then does
`T = [T; G.nnc.T]` and indexes with `s.internalConn`, misaligning the
fracture NNC transmissibilities against the connection list.

This module appends matching lengths and test00 check F asserts it. The
reference file was left unmodified — it is a dependency of other modules
and fixing it is a separate decision.
