# PEDFM-U

3D projection-based EDFM (pEDFM) for **general unstructured grids**:
PEBI/Voronoi, tetrahedra, extruded prisms — any grid whose cells are
convex polyhedra.

Extends the Cartesian 3D pEDFM of Olorode, Wang & Rashid (2020),
*"Three-Dimensional Projection-Based Embedded Discrete-Fracture Model for
Compositional Simulation of Fractured Reservoirs"*, SPE Journal 25(04)
2143-2161, SPE-201243-PA. The paper defers everything beyond Cartesian
grids explicitly (p. 2147): *"Although our algorithm could be adapted for
simple corner-point grids with six faces, we leave this (as well as more
complicated corner-point and unstructured grids) to a future study."*

Companion to `PEDFM-CPG` (corner-point) and `CPEDFM-U` (the continuous
variant on unstructured grids). No file outside this module is modified.

## Usage

```matlab
run modules/PEDFM-U/startupPEDFMU.m

G = computeGeometry(makeLayeredGrid(pebiGrid2D(h, dims), nz)); G.rock = ...;
[G, fracplanes] = EDFMGridU(G, fracplanes, 'Tolerance', tol);  % from CPEDFM-U
G   = fracturematrixShaleNNC3D(G, tol);                        % existing
G   = pMatFracNNCsPEDFM_U(G, tol);                             % NEW
ops = setupPEDFMOpsTPFA3D_CPG(G, G.rock, tol);                 % from PEDFM-CPG
model.operators = ops;
```

Options on `pMatFracNNCsPEDFM_U`:

| Option | Default | Meaning |
|---|---|---|
| `relTol` | `1e-6` | Equidistance tolerance, as a fraction of the host cell's `V^(1/3)` |
| `defaultFace` | `'bottom'` | Seeds the paper's Case 3, where every antipodal pair is equidistant and the first face is arbitrary |
| `unionArea` | `true` | Apply the Eq. 29 union when several fracture cells share a face |

## Method: antipodal face pairs replace the three axes

The paper's algorithm rests on the host cell having exactly six faces
forming three opposing pairs, one per grid direction, and picking the
nearer face of each pair (Fig. 2, Case 0). Unstructured cells have no
logical directions, no face tags, and a varying number of faces, so *"the
pair of faces in the X-direction"* is undefined.

The observation this module is built on:

> The three axis pairs of a hexahedron are just its three **antipodal face
> pairs**, and that notion survives on any convex polyhedron.

So:

> For each face `f` of the host cell, let its partner `g` be the face whose
> outward unit normal is most anti-parallel to `f`'s. Select `f` when the
> fracture is closer to `f`'s plane than to `g`'s, measuring perpendicular
> distance from the fracture cell centroid.

On a hexahedron the partner of *left* is *right* (dot product −1), and the
perpendicular-distance test is **identical** to the paper's
centroid-distance test: both faces of an axis pair share the cell's
mid-plane in the other two coordinates, so those terms cancel from the
comparison. The rule therefore reproduces the published selection exactly
— which `tests/test20_regression_cartesian_U.m` verifies against
`pMatFracNNCs3D` rather than leaving it asserted.

On a PEBI cell with *n* faces the rule yields about *n*/2 selections, one
per antipodal direction — the natural reading of "project the fracture
once per direction" for that cell shape.

**Equidistant pairs** (the paper's Cases 1–3) are resolved by the same
node-containment rule used in `PEDFM-CPG`: Table 1 of the paper is a
hexahedron incidence table, in which the two listed faces are exactly the
other faces of the cell containing the node farthest from the fracture
plane. That reading needs no node ordering and transfers unchanged here.
Case 3 (every pair tied) is seeded by `defaultFace`, matching the
Cartesian release's arbitrary pick so the two agree.

**Tip guard.** Selection sees only the infinite fracture *plane*, so a
selected face is dropped unless its centroid, projected onto that plane,
lies inside the fracture polygon with a buffer of half the face's
characteristic length. Without it a fracture spanning part of the domain
would project beyond its own extent. Same construction as CPEDFM-U.

Areas use Eq. 24 with the real face normal; Eq. 25 uses normal-projected
permeabilities `n'*K*n`, so anisotropic rock is handled on faces of
arbitrary orientation (1-, 3- and 6-component perm supported).

## pEDFM vs CPEDFM on unstructured grids

`CPEDFM-U` seals every face whose two neighbouring cell centroids straddle
the fracture. That set is by construction a **cut** of the matrix
connectivity graph, so the barrier is watertight for any cell shape.

pEDFM selects faces **independently per host cell**, so the sealed set need
not be connected: a sealing fracture throttles flow without necessarily
compartmentalizing the domain. That difference is the reason the
continuous variant exists. Use CPEDFM-U when you need guaranteed
compartmentalization; use this module when you want the published pEDFM
algorithm.

## Files

```
startupPEDFMU.m                       path setup + anti-shadow guards
UtilsPEDFM/
  selectProjectionFacesPEDFM_U.m      antipodal-pair rule + Cases 1-3
  calculateProjectionPEDFM_U.m        per-plane projection, tip guard, Eqs. 24-25
  pMatFracNNCsPEDFM_U.m               driver; assembles G.nnc
tests/
  test20_regression_cartesian_U.m
  test21_pebi_sealing.m
```

Reused rather than duplicated: `EDFMGridU` and the polygon utilities from
`CPEDFM-U/UtilsCPEDFM`, and `transmultPEDFM3D_CPG`,
`setupPEDFMOpsTPFA3D_CPG`, `projectedAreaPEDFM_CPG` from
`PEDFM-CPG/UtilsPEDFM` — all three are grid-agnostic, being face-indexed
and never touching logical directions. `startupPEDFMU.m` asserts each
resolves to the intended copy.

Uses `hfm`, not `hfm_v2`: this module grids fractures with `EDFMGridU` and
never calls the corner-point intersection engine, so the two-output
`convexpolygonarea` that PEDFM-CPG requires is irrelevant here.

## Test status (MATLAB R2021b, 2026-08-01)

Both pass.

* **test20 — Cartesian regression.** 9x9x3 grid, one fracture inclined to
  all three axes, 60 interior NNCs. Against `pMatFracNNCs3D`: pM-M pair
  sets **identical** (101 pairs); projected areas match to **2.2e-16**
  (49 singleton pairs) and **2.0e-15** (52 colocated); pM-F
  transmissibilities match as a multiset to **2.7e-20** (153 values); face
  multipliers to **4.2e-16**. This is the check that turns the
  "antipodal pairs = axis pairs" argument into a verified claim.
* **test21 — extruded PEBI.** 1062 cells, **6–9 faces per cell**, 111 host
  cells, 407 projections over 284 faces, zero ties and zero fallbacks.
  100% of the 181 faces the fracture crosses within its extent are
  throttled (median multiplier 0.0000); no host ever selects more than one
  face of an antipodal pair (max 0.57 of a cell's faces). Tip guard
  confirmed: a 10x30 m central fracture patch projects only onto faces
  within ~9 m of its footprint. Cross-fracture pressure drop 1.29 bar vs
  0.90 bar for plain EDFM.

## Known limitations

* **Cells must be convex** (true for PEBI/Voronoi, tets, extruded prisms,
  Cartesian). The antipodal-partner search and the perpendicular-distance
  test both assume it.
* **The barrier is not watertight** — inherent to pEDFM's per-cell
  selection, not a defect of this port. See the pEDFM vs CPEDFM section.
* **Pairing need not be mutual** on a general polyhedron: `f`'s partner may
  be `g` while `g`'s partner is `h`, so both `f` and `g` can be selected.
  This is harmless — each face carries at most its own projected area and
  the Eq. 27 multiplier is clamped at zero — and cannot occur on a
  hexahedron. It did not occur in test21 (check C).
* **Exact ties are measure-zero** on unstructured grids: test21 recorded
  zero. Cases 1–3 are Cartesian-symmetry artifacts; they are covered by
  test20, where 3 of 60 hosts land exactly on Case 3.
* Fractures coplanar with a cell face are registered as
  `'fracmat boundary'` and handled by the standard EDFM path, not
  projected — same convention as the Cartesian release.
* Fracture-fracture intersection NNCs are not exercised here; use the
  existing `fracturefractureShaleNNCs3D` for intersecting networks.
