# MRST → CMG Preprocessor: Implementation Plan

Plan of record for turning the existing export code into one function that
takes an MRST grid — any grid family, with or without fractures — and writes a
CMG dataset that runs. Companion to `README.md` (reference) and
`cpedfm_export.md` (implementation log of the Eclipse/CMG work to date).

Drafted 2026-07-27 against CMG 2024.30 and MATLAB R2021b on this machine.

## Status

| Phase | State | Evidence |
|---|---|---|
| 0 — contract and dispatch | **done** | `test30` checks A–E |
| 1 — structured grid writers, matrix only | **done** | `test30` F.1–F.3, `test31`, 5/5 decks run in IMEX |
| 2 — unstructured connection list | **done** | `test30` D and F.4 |
| 3 — EDFM fracture cells | **done** (all grid families) | `test32`–`test34`, `test38` |
| 4 — projection on non-uniform spacing | **done** | `test37`; R1 closed |
| 5 — pEDFM and CPEDFM | **done** (structured grids) | `test32` |
| 6 — wells and fluid templates | **wells done** | `examples/exportNorneWaterfloodCMG.m` |
| — Eclipse matrix export | **done** | `test35` |
| — Eclipse fracture export | **done** | `test36` |

Implemented so far: `cmgDetectGrid.m`, `cmgDetectFractureModel.m`,
`buildCMGModel.m`, `writeCMGModel.m`, `writeCMG.m`, and
`tests/test30_cmg_grid_export.m`, `tests/test31_eagleford_matrix_cmg.m`.

Measured on the Eagle Ford matrix grid (`test31`): 1976 cells, `cartDims`
[26 19 4], dx ranging 2.5–7.222 m, correctly classified `vari`. IMEX reports
1976 total and 1976 active blocks. Bulk and pore volume agree with MRST to
~4e-15 relative on the MRST side; IMEX's reported pore volume of 205.75 M rm³
against a geometric 205.74 M rm³ is the rock-compressibility response to the
hydrostatic gradient across the section, not an export error.

Three CMG constraints were discovered by running rather than by reading, and
are now encoded in `writeCMGModel.m`:

1. `*DEPTH` addresses a **single block**; the array form is `*DEPTH-TOP`.
2. `*DEPTH-TOP` under `*GRID *CART` requires every block of a layer to share
   one depth. The connection-list pseudo-grid puts all blocks in layer K=1 at
   differing depths, so it must be `*GRID *VARI` — only `*VARI` permits blocks
   of one layer to lie at different depths.
3. `*TRANSI *MOD` is rejected unless the fundamental grid already carries a
   transmissibility assignment, so `*TRANSI/J/K *CON 1.0` is emitted first.
4. **Builder will not accept `*GRID` from inside an `*INCLUDE`.** IMEX reads
   such a deck without complaint, but Builder reports an error at the end of
   the `*ZCORN` block and renders an empty grid ("Plane 1 of 1"). The
   fundamental grid must be declared in the `.dat`; only the arrays belong in
   the include. `cmgGridHeader.m` emits those two cards and both deck writers
   use it. Verified by opening the dataset in Builder 2024.30: with the header
   moved out, the Norne sector renders as 7392 grid blocks / 5204 exterior
   faces.

### Fracture export, and proving the connections are live

`test32` exports EDFM, pEDFM and CPEDFM from one shared grid and reproduces
the counts `test22` records on the Eclipse side — EDFM 20 NNC / 0 seals,
pEDFM 58 / 28, CPEDFM 68 / 48 — with pore volume conserved to 1.5e-16.

Deck acceptance proves nothing about the connection list. **IMEX runs to
Normal Termination even when a `*SCONNECT` record addresses a null block** —
verified by deliberately corrupting one. A silently ignored connection list
would therefore look exactly like a working one, which is the same trap the
earlier CMG writer fell into by emitting its NNCs as comments.

`test33` settles it with a flow test. A producer is completed in a **fracture
cell only**; every geometric connection of the padding region is zeroed, so
its sole path to the reservoir is the exported `*SCONNECT` list. Two decks are
run, identical but for the connection block being commented out:

| | live | severed |
|---|---|---|
| Cumulative oil | 3 887.7 m³ | 0.0196 m³ |
| Average reservoir pressure | 18 058 → 5 229 kPa | 18 058 kPa, unchanged |

The well drained 870× the fracture's own pore volume (4.47 m³) out of the
matrix, and depleted the reservoir by 12 830 kPa. Severed, it produced
essentially nothing and the reservoir was untouched — which also confirms the
padding isolation is complete. `tools/cmg-mrst-mcp/verify_sconnect.py`
automates the comparison.

Two implementation notes. Fracture grids are always emitted as `*GRID *VARI`
even when the matrix grid is uniform, because padding blocks sit at depths
unrelated to the layer they occupy and only `*VARI` permits that. And
`cmgDetectGrid` classifies using `G.Matrix` when the grid carries fractures:
the assembled grid's nodes include the fracture cells, so it is no longer an
axis-aligned tensor product even when the reservoir grid is.

### Corner-point fractures

`test34` extends the same machinery to corner-point grids. Padding layers are
appended to ZCORN, and each used padding block is given corners at
`z_fracture ± dz/2`, so its centre **is** the fracture cell's true depth.
Depth then comes from the geometry and the matrix blocks keep their own
corner-derived depths untouched — no `*DEPTH-TOP` override is needed, which
matters because that keyword would otherwise overwrite every matrix depth
too. Pillar bottoms in COORD are extended to span the deepest padding corner.
Padding blocks need not be conforming: every connection among them is
removed.

Measured (20×10×8 faulted CPG, one fracture): ZCORN grows to exactly
8·NX·NY·(NZ+pad); pore volume transfers with **zero** relative error; and the
fracture block centres match the MRST centroids to **0.000e+00 m**.

The full Johansen case now runs end to end — 88 775 matrix cells + 1 141
fracture cells, 1 155 NNCs, exported as 120 000 blocks and rendered in
Builder with the fracture cells sitting at their true depths.

### Unstructured grids with fractures

The connection-list representation turns out to suit fractures *better* than
the structured one. Padding layers exist only because a structured grid has no
spare block addresses; a connection-list model gives every cell its own
address already, so matrix and fracture cells are numbered straight through
and the whole connectivity — internal matrix faces plus every NNC — is one
`*SCONNECT` list. No padding, no `*NULL`, no addressing scheme.

Projection seals fold in rather than being emitted separately. In the
structured path a sealed matrix face becomes a `*TRANS?*MOD` multiplier
because the face is only implied by the grid; here the face **is** an explicit
connection, so the multiplier multiplies its transmissibility directly. Same
physics, one mechanism fewer.

`test38` (PEBI footprint extruded to three layers, 180 matrix + 36 fracture
cells) checks that every cell becomes a block with zero padding layers, that
pore volume transfers exactly, and that the connection count equals internal
matrix faces minus sealed ones plus NNCs — 624 = 588 − 0 + 36. The deck runs
in IMEX, and the same model exports to Eclipse as an `Nc x 1 x 1` strip with
the connectivity as `NNC` records.

Two things worth recording:

* **Gridder.** `EDFMshalegrid` and `EDFMGridCPG` both assume 8-corner cells
  and fail inside `checkplaneAABBintersect` on a PEBI grid. The unstructured
  gridder is `EDFMGridU` from CPEDFM-U, which uses half-space clipping and is
  exact for any convex cell. Its startup also runs `mrstModule add hfm`, so
  `hfm_v2` must be loaded after it — the same trap as `startupCPEDFMCPG`.
* **Tolerances.** The in-memory pore-volume check is exact (0.00e+00); the
  Eclipse check of the same quantity is asserted only to 1e-8, because `PORV`
  is written as `%.8e` — nine significant figures — so a value that has been
  through the file cannot agree more closely. Tightening it would be testing
  the format, not the code.

### Where fracture cells are placed

A fracture cell **always** takes the host matrix column it intersects, read
from the matrix-fracture NNC that created it. Cells sharing a column stack
downwards through successive padding layers, so the padding-layer count is
the deepest such stack.

The rule that makes this work: a block's depth must lie inside its own
column's vertical extent. Two alternatives were tried and both break it.

**Sequential slots** — `(i,j)` from a running counter, depth from the
fracture. Numerically identical (padding connections are all zeroed; the model
talks only through `*SCONNECT`) but the blocks land at map positions unrelated
to their fracture and float above and below the reservoir. A Builder render of
the Johansen export is how the defect was found.

**Nearest free slot** — host column when free, otherwise the nearest free
column in the same padding layer. This keeps the padding to one layer and
looks like a clean win on block count, but it reintroduces the same floaters:
a displaced block still carries its fracture's true depth while sitting in a
column that lies at a different depth. Rejected for the same reason as
sequential slots.

| placement | blocks | grid | include | placement correct? |
|---|---|---|---|---|
| sequential | 120 000 | 100×100×12 | 19.5 MB | no — scattered |
| nearest free slot | 120 000 | 100×100×12 | 19.5 MB | no — displaced cells float |
| **host column, stacked** | 330 000 | 100×100×**33** | 52 MB | **yes** |

The cost is real: property arrays span the full `NX·NY·NZ` index space, so
240 084 of the 330 000 Johansen blocks are inactive padding. Johansen is the
bad case — four near-vertical fracture planes spanning all 11 layers, two of
them crossing the same columns, hence 22 padding layers. It is accepted
deliberately. An export that cannot be inspected cannot be trusted, and the
inactive blocks cost file size rather than simulation time.

Worth revisiting only if file size becomes a practical limit; the lead would
be CMG's sparse `*IJK` array notation for the padding region, which addresses
values by block rather than requiring a full-grid array.

### Pore volume: *PORVOL, not *VOLMOD

Fracture pore volume was originally carried by `*VOLMOD`, a bulk-volume
multiplier, which meant computing the bulk volume of a placeholder padding
block — awkward for corner-point, where that depends on pillar geometry.
`*PORVOL` sets pore volume directly and overrides the porosity-times-bulk
calculation, so MRST's pore volume transfers exactly for every block on every
grid family. It also removes the `*VOLMOD` accumulation hazard. Relative
error went from 1.5e-16 to exactly zero.

### Which EDFM module: hfm_v2

**`hfm_v2` is the canonical module for this workspace.** Load it *after*
`startupCPEDFMCPG`, which internally runs `mrstModule add hfm` and would
otherwise shadow it.

`hfm` and `hfm_v2` hold the same 180 files, and 173 of them differ only in
line endings. Seven differ in substance, two of which matter:

* `convexpolygonarea` — `hfm` returns one output, `hfm_v2` two. `EDFMGridCPG`
  calls it with two, so loading `hfm` kills the corner-point gridder with
  "Too many output arguments" at `pebiAABBintersectCPGv5:648`.
* `calcdavg` — sets the fracture-matrix connectivity index (Moinfar et al.,
  2013). `hfm` subdivides the matrix cell 10×10×10 uniformly; `hfm_v2` uses
  10×10×3 with a seeded z-perturbation, evidently an adaptation for
  corner-point cells, which are not boxes.

Measured on the Cartesian test case, the two give **identical** NNC
transmissibilities — `sum(T) = 1.328887254593e-15`, matching to every printed
digit — so standardising on `hfm_v2` changes none of the results recorded
here. Divergence is expected only on non-box cells, and has not been measured.

`tests/test32`, `test33` and `test34` all load `hfm_v2`.

### Eclipse export

`writeEclipse.m` and `writeEclipseModel.m` are siblings of the CMG pair,
consuming the same `cmgModel` struct. Corner-point grids emit
`SPECGRID`/`COORD`/`ZCORN`; Cartesian and variable grids emit
`DXV`/`DYV`/`DZV`/`TOPS`. `test35` reads each deck back with
`readEclipseDeck` and rebuilds the grid from the file alone — pore volume
recovered to 3.8e-10.

Two choices worth recording:

* **`TOPS` carries every block**, not just layer 1. Eclipse accepts the
  `NX*NY` form, but a reader then has to infer which convention was used, and
  MRST's own `readGRID` expects one value per block. One extra array buys
  unambiguity.
* **No run-length compression.** Eclipse's `n*value` form would shrink the
  files considerably, but the browser port has to reproduce the output byte
  for byte, and two independent implementations of a compression heuristic is
  a bug waiting to happen.

Fracture models export too. The padding cells become ordinary blocks and
their physics enters through an `EDIT` section following the grid keywords in
the same file, mirroring the CMG path keyword for keyword:

| model field | CMG | Eclipse |
|---|---|---|
| connections | `*SCONNECT` (md·m) | `NNC` (cP·m³/day/bar) |
| fracture pore volume | `*PORVOL` | `PORV` |
| block depth | `*DEPTH-TOP` | `DEPTH` |
| padding isolation | `*TRANS?*MOD` | `EQUALS` on `TRANX/Y/Z` |
| projection seals | `*TRANS?*MOD` | `EQUALS` on `MULTX/Y/Z` |

Transmissibility is converted through SI with MRST's own unit constants
rather than a hand-derived factor: the model stores md·m for the CMG path, and
Eclipse METRIC wants cP·m³/(day·bar).

**Wells** are written as a commented `SCHEDULE` block — `WELSPECS`,
`COMPDAT` and the controls — because they belong after PROPS/REGIONS/SOLUTION
and this deck stops at the end of GRID. Connection factors go in `COMPDAT`
item 8 explicitly rather than being left to Eclipse's radial model: on a
padding-cell or connection-list grid the block geometry is a placeholder, so a
computed well index would be meaningless. `test39` reads the deck as text
(a commented block is invisible to a parser, which is what makes it safe and
also what makes it easy to get silently wrong), checks one `COMPDAT` row per
perforation, and confirms the factors are MRST's `WI` converted, to 1.1e-09.
It also asserts the deck still parses with `SCHEDULE.control` empty, proving
the block really is inert until uncommented.

`test36` reads each deck back and checks the NNC count, the converted
transmissibilities and the pore volume — EDFM 20, pEDFM 58, CPEDFM 68
connections, `PORV` to 1e-13. It also asserts the two writers agree on the
same model: 720 blocks, 68 connections and 48 seals in both. That agreement
is the point of keeping `buildCMGModel` separate from the writers.

### On writing properties

`writeCMG` writes **no invented data by default**. `'hostDeck'` selects what
accompanies the include: `'stub'` (default) is grid only, for opening in
Builder and populating there; `'none'` is the include alone; `'runnable'`
adds a generic placeholder fluid model so the deck can be executed by the
IMEX harness, and must be asked for explicitly. The earlier default emitted a
placeholder PVT table, relative-permeability curves and contacts that looked
like reservoir data but came from nowhere — an easy thing to mistake for
derived output.

Note that a `'stub'` deck is deliberately incomplete and IMEX will not read
it; it reports "End-of-File before *RUN" for any deck lacking a fluid model,
however the section markers are arranged (checked across four variants).
Builder is the intended consumer. The include itself is independently valid
and runs unchanged inside a complete deck.

### Field-case validation (Norne sector)

`examples/exportNorneWaterfloodCMG.m` builds a real faulted corner-point case
in MRST — the same `[10 25; 35 55; 1 22]` cut the hand-written
`NorneOut2CMGOMO.m` produced the committed CMG includes from — adds a water
injector and a BHP producer, solves it, and exports grid and wells in one
`writeCMG` call. The result runs a **360-day waterflood in IMEX to Normal
Termination**, 27 timesteps, material balances 1.000/1.000/1.000.

| Quantity | MRST | IMEX | Agreement |
|---|---|---|---|
| Active blocks | 7056 | 7056 | exact |
| Total blocks | 7392 | 7392 | exact |
| Well completion layers | 42 | 42 | exact |
| Water injection rate | 2000 m³/day | 2000.000 m³/day | exact |
| Pore volume (`cpor = 0`) | 135 265.921 M rm³ | 135 265 M rm³ | 6.8e-6, IMEX display precision |

The pore-volume check is run with rock compressibility zeroed on purpose.
With the default `cpor`, IMEX reports 135 316 M rm³ — a 0.037% difference that
is entirely the compressibility response to initialization pressure sitting
above `*PRPOR`, not an export error. `writeCMG` exposes `'cpor'` so this
reconciliation can be made exact on demand.

Well indices transfer through `*PERF *WI` in md·m. CMG defines
`q = wi · λ · Δp` with `wi` the geometric part only, which is exactly MRST's
`W.WI`, so the values move across without a model change. Letting CMG
recompute them from `*GEOMETRY` would apply a radial-flow model to the
exported block geometry rather than the geometry MRST derived them from.

A fourth CMG constraint surfaced here: `cmgDetectGrid` originally returned
early for corner-point grids and never computed the depth range, so the
initialization datum defaulted to 0 m and produced a negative `*DGOC`. Depth
range is now computed for every family before any family-specific return, and
the contacts are derived from it (`*DGOC` at the top, `*DWOC` below the base)
rather than from a fixed offset.

---

## 1. Objective

```matlab
mrstModule add preprocessor
% ... build G in MRST exactly as you do today ...
info = writeCMG(G, 'case.dat', 'simulator', 'IMEX', 'wells', wells);
```

One entry point. The grid family, the presence of fractures, and which
embedded-fracture method was used are all **detected from `G`**, not passed in
as flags. Supported inputs:

| Grid family | MRST source | Fractures |
|---|---|---|
| Uniform Cartesian | `cartGrid` | none / EDFM / pEDFM / CPEDFM |
| Variable Cartesian | `tensorGrid`, `meshHFsystem` | none / EDFM / pEDFM / CPEDFM |
| Corner-point | `processGRDECL` | none / EDFM / pEDFM / CPEDFM |
| Unstructured | `pebi`, `pebiGrid2D`, `triangleGrid`, tets | none / EDFM / pEDFM / CPEDFM |

---

## 2. Verified foundations

Everything below was checked on this machine rather than assumed. This matters
because the module's current documentation states that CMG keyword syntax
"could not be machine-verified here" — that is no longer true, and the plan
depends on it being false.

### 2.1 CMG has a native keyword for every ingredient

Verified against the installed HTML keyword manual at
`C:\Program Files (x86)\CMG\Manuals\2024.30\{IMEX,GEM,STARS}`. Use
`tools/cmg_manual.py find|show <keyword> [simulator]` to re-read any of these.

| Keyword | Form | Why it matters |
|---|---|---|
| `*SCONNECT` | `uba1 uba2 trans` | **The enabler.** Arbitrary connection with transmissibility given directly in **md·m**. Negative value = multiplier. Can also override an existing regular connection. Present in IMEX, GEM and STARS. |
| `*IRCONNECT` | `uba1 uba2 L1 L2 A12 *I|*J|*K` | Geometric alternative; CMG derives T from permeability. Required instead of `*SCONNECT` only under `*DYNAGRID`. |
| `*VOLMOD` | array | Bulk-volume **multiplier** (default 1.0), not an absolute volume. Zero makes a block null. **Accumulates** if a block is specified twice. |
| `*NULL` | array | Deactivates unused padding slots. |
| `*MOD` | `i1(:i2) j1(:j2) k1(:k2) op value`, op ∈ `+ - * / =` | I-J-K ranges supported, so `*TRANSI *MOD 1:20 1:20 6:15 = 0` isolates a padding region compactly. |
| `*GRID` | `*CART | *VARI | *CORNER | *RADIAL  ni nj nk` | **No unstructured grid type exists.** See §5.4. |
| `*CORNERS` | array | Complete per-block corner locations — lets a corner-point grid be written from `G.nodes` without a `grdecl`. Must not be combined with other corner-point array keywords. |

The critical consequence: `*SCONNECT` takes exactly what an EDFM NNC list
already is — a cell pair and a transmissibility — in exactly the unit
`writeCMGCPEDFM` already converts to (`milli*darcy*meter`). The NNC records
currently emitted as `**` comments become live cards with no numerical change.

### 2.2 The module already loads as an MRST module

`startup.m` calls `build_mrst_path_tree`, which auto-scans `modules/`. Confirmed
in a live session:

```
REGISTERED: yes -> ...\modules\preprocessor
  writeCMGEmbeddedFractureModel  exist=2
  writeCMGCPEDFM                 exist=2
  preparePEDFMExport             exist=2
LOADED_MODULES: preprocessor
```

No packaging work is needed. `mrstModule add preprocessor` works today.

### 2.3 Grid family is detectable from `G`

Measured stamps:

| Built with | `G.type` | `cartDims`? |
|---|---|---|
| `cartGrid` | `tensorGrid, cartGrid, computeGeometry` | yes |
| `tensorGrid` | `tensorGrid, computeGeometry` | yes |
| `processGRDECL` | `processGRDECL, computeGeometry` | yes |
| `triangleGrid` | `triangleGrid, computeGeometry` | **no** |
| `pebi` | `pebi, computeGeometry` | **no** |
| `pebiGrid2D` (upr) | `clippedPebi2D, computeGeometry` | **no** |

Unstructured constructors are numerous and their type strings vary, so the
robust discriminator is `isfield(G,'cartDims')` for structured-vs-unstructured,
with `G.type` used only to split the structured cases. Do not enumerate
unstructured type names.

### 2.4 Corner-point geometry is *not* on `G`

After `processGRDECL`, the grid's top-level fields are only
`type, griddim, nodes, cartDims, cells, faces`. `COORD` and `ZCORN` are
consumed and discarded. So corner-point export either takes an optional
`grdecl` argument, or rebuilds corners from `G.nodes` and writes `*CORNERS`.
A probe confirmed each cell of a `processGRDECL` grid resolves to 6 faces and
8 unique nodes, so the rebuild is feasible for well-formed hexahedra.

---

## 3. Dispatch

```
writeCMG(G, ...)
│
├─ grid family
│   ├─ ~isfield(G,'cartDims')            → UNSTRUCTURED  → connection-list (§5.4)
│   ├─ any(strcmp(G.type,'processGRDECL')) → CORNER-POINT → *GRID *CORNER
│   ├─ any(strcmp(G.type,'cartGrid'))    → UNIFORM       → *GRID *CART + *DI/*DJ/*DK *CON
│   └─ any(strcmp(G.type,'tensorGrid'))  → VARIABLE      → *GRID *VARI + *DI *IVAR ...
│
├─ fractures present?
│   └─ isfield(G,'FracGrid') && isfield(G,'Matrix')
│         no  → matrix-only deck, done
│
└─ method (only if fractures)
    ├─ ~isempty(G.nnc.pMMfaces)   → CPEDFM
    ├─ ~isempty(G.nnc.pMMneighs)  → pEDFM
    └─ otherwise                  → EDFM
```

**Trap:** test *non-empty*, never `isfield`. `prepareEDFMExport` deliberately
writes empty `pMMfaces` and `pMMneighs` onto a normalized EDFM grid, so an
`isfield` check misclassifies EDFM as CPEDFM.

The detected combination must be **reported** — in the returned `info`, in the
deck banner, and in the JSON manifest — never applied silently. A `'model'`
option forces a specific method and errors if the grid cannot supply it.

---

## 4. The mapping

| MRST concept | CMG expression |
|---|---|
| matrix cells | grid blocks in the fundamental grid |
| cell dimensions | `*DI/*DJ/*DK` (`*CON` uniform, `*IVAR/*JVAR/*KVAR` variable) |
| corner-point geometry | `*COORD`/`*ZCORN` (from `grdecl`) or `*CORNERS` (from `G.nodes`) |
| `rock.perm` | `*PERMI` / `*PERMJ` / `*PERMK`, mD |
| `rock.poro` | `*POR` |
| inactive cells | `*NULL` |
| fracture cells | extra padding k-layers; unused slots `*NULL` |
| fracture pore volume | `*VOLMOD` = `PV_target / (bulk × poro)` |
| fracture cell depth | per-block depth override |
| **any NNC** (matrix–frac, frac–frac, well–frac, projection) | `*SCONNECT uba1 uba2 trans` [md·m] |
| pEDFM/CPEDFM seal, logical neighbours | `*TRANSI/J/K *MOD i j k = η` |
| pEDFM/CPEDFM seal, irregular pair | `*SCONNECT uba1 uba2 -η` (negative = multiplier) |
| padding isolation | `*TRANSI/J/K *MOD <range> = 0` |
| wells | `*WELL` / `*PERF` with **explicit** well index |
| fluid model, schedule | host-deck template per simulator (not derived from `G`) |

All three fracture methods reduce to the same two mechanisms. Only the record
counts differ:

| Method | Connections | Matrix-face modification |
|---|---|---|
| EDFM | matrix–frac, frac–frac | none |
| pEDFM | + projection NNCs | partial multipliers |
| CPEDFM | + continuous projection NNCs | continuous / full-face seals |

---

## 5. Design decisions

### 5.1 Target API

```matlab
info = writeCMG(G, outFile, ...
    'simulator', 'IMEX', ...   % IMEX | GEM | STARS
    'wells',     wells, ...    % not derivable from G
    'rock',      G.rock, ...   % defaults to G.rock
    'grdecl',    grdecl, ...   % optional; corner-point fast path
    'model',     'auto', ...   % auto | matrix | EDFM | pEDFM | CPEDFM
    'template',  '', ...       % host-deck template; default per simulator
    'tol',       1e-5);
```

Only `G` and `outFile` are required. `writeCMGEmbeddedFractureModel` stays as a
thin explicit-method wrapper for backward compatibility.

### 5.2 One intermediate contract

Every grid family and every fracture method funnels through a single struct, and
one writer consumes it. This is what keeps "fractures optional" from becoming a
second code path.

```matlab
cmgModel.meta      % simulator, version, modelType, gridFamily, units, counts
cmgModel.grid      % kind: 'cart'|'vari'|'corner'|'connlist'
                   % dims, di/dj/dk | coord/zcorn | corners, null mask
cmgModel.props     % permi, permj, permk, por, volmod, depth
cmgModel.conn      % nConn x 7: [i1 j1 k1 i2 j2 k2 trans]   (md·m)
cmgModel.mult      % logical-neighbour seals: [dir i j k value]
cmgModel.wells     % name, type, control, completions [uba, WI]
```

Builders fill `.grid`/`.props`; adapters fill `.conn`/`.mult`; the writer is
pure formatting. No fractures simply means `.conn` and `.mult` are empty.

### 5.3 Fracture cells

Keep the padding-layer scheme already implemented in `writeCMGCPEDFM`:
`nzExtra = ceil(nFrac/(NX·NY))` extra k-layers, unused slots `*NULL`.

Two changes from the current code. Pore volume moves from a comment to
`*VOLMOD`, computed as `PV_target / (bulk_volume × porosity)` — and because
`*VOLMOD` accumulates across data lines, each block must appear exactly once.
Padding blocks must also have their automatic connections killed with
`*TRANSI/J/K *MOD` over the padding range, so fracture cells communicate only
through `*SCONNECT`.

### 5.4 Unstructured grids

CMG has no unstructured grid type, so there is exactly one viable route:
a **connection-list pseudo-grid**. This is the same pattern
`writeEclipseCoarseCPEDFM.m` already uses successfully on the Eclipse side.

- `*GRID *CART Nc 1 1` — placeholder geometry, physically meaningless
- every geometric connection killed with `*TRANSI *MOD 1:Nc 1 1 = 0`
- `*VOLMOD` carries each cell's true pore volume
- per-block depth carries the true centroid depth, for gravity
- **the entire connectivity** — every internal matrix face plus every fracture
  NNC — becomes `*SCONNECT` records

This is exact rather than approximate: a TPFA simulator consumes only
`(cell pair, transmissibility)`, so a PEBI model is reproduced faithfully by a
pseudo-grid whose whole connectivity is the connection list.

Two consequences to plan for. **Size**: internal faces run about 2.6× cell count
for a structured 3D grid (measured) and higher for PEBI, so expect roughly 3–4
`*SCONNECT` records per cell; a 100k-cell model means ~350k records, which needs
its own `*INCLUDE` file and a written record-count budget. **Visualization**:
I-J-K meaning is gone, so CMG Results will render a 1-D string of blocks — keep
MRST for spatial QC and treat the CMG run as a numerical check.

The fracture side needs no new work here: `CPEDFM-U` already produces the same
`G.nnc.pMMfaces`/`pMMarea` contract as the structured modules, so the adapters
and the connection writer are reused unchanged.

---

## 6. Phases

Each phase has one acceptance test. Do not start the next until it passes.

### Phase 0 — Contract and dispatch
Define `cmgModel` (§5.2) and `detectGridFamily`/`detectFractureModel` (§3).
No deck writing yet.
**Accept:** a unit test classifies `cartGrid`, `tensorGrid`, `processGRDECL`,
`pebi`, and Eagle Ford's `G` correctly, and detects EDFM / pEDFM / CPEDFM /
none on the test22 grids — including the empty-`pMMfaces` trap.

### Phase 1 — Structured grid writers, matrix only
`*CART`, `*VARI`, `*CORNER`. Lift the corner-point path from
`NorneOut2CMGOMO.m`, fixing the defects in §8.
**Accept:** Eagle Ford's `G_matrix` (no fractures) exports and reaches
`Normal Termination` in IMEX; total pore volume and active block count match
MRST to machine precision.

### Phase 2 — Unstructured connection-list writer
Per §5.4.
**Accept:** a small PEBI case exports and runs; pore volume matches MRST;
`*SCONNECT` record count equals the internal-face count.

### Phase 3 — EDFM fracture cells
Padding layers, `*NULL`, `*VOLMOD`, depth, padding isolation, and `*SCONNECT`
for matrix–fracture and fracture–fracture NNCs. **One** hydraulic fracture, not
153.
**Accept:** single-phase incompressible pressure from CMG matches the MRST
solution cell by cell to a stated tolerance. This is the first test that proves
transmissibilities landed correctly; deck acceptance alone does not.

### Phase 4 — Projection validation on non-uniform spacing
Before trusting pEDFM/CPEDFM on Eagle Ford, verify projected areas and
transmissibilities on a deliberately uneven tensor grid. See Risk R1.
**Accept:** projected areas match hand-computed values on a case with a 3:1
spacing contrast.

### Phase 5 — pEDFM and CPEDFM
Add seals: `*TRANS**MOD` for logical neighbours, negative `*SCONNECT` for
irregular pairs. `getEmbeddedFractureMultipliers` already produces the numbers.
**Accept:** all three methods export from one shared grid, with the record
counts test22 reports; the CPEDFM deck reproduces the MRST pressure solution;
a sealing fracture demonstrably blocks flow (negative control).

### Phase 6 — Wells and host-deck templates
Map wells to `*PERF` with explicit well indices — never let CMG compute a
Peaceman factor on placeholder geometry. Handle the well–fracture NNCs from
`wellfractureNNCs3D`. Add IMEX black-oil / GEM compositional / STARS thermal
templates.
**Accept:** full Eagle Ford case runs with its horizontal well; well rates
tracked against MRST.

---

## 7. Risks

**R1 — Projection functions on non-uniform grids. CLOSED, no defect found.**
`test37` checks the projection against a value computable by hand: for a
planar fracture of area A and unit normal n, the shadow gathered on faces
normal to direction d is exactly `A·|n_d|`, whatever the grid.

On a 4:1 dx-contrast, 2:1 dz-contrast tensor grid and on a uniform grid, the
pEDFM projection reproduces the analytic totals to **1.6e-16** and the two
grids agree exactly (720.000000 / 360.000000 on both), despite spreading the
shadow over different numbers of faces. No face is ever over-assigned, so the
`max(faceArea - projArea, 0)` clamp in `getEmbeddedFractureMultipliers` is
never the thing hiding an error.

One trap the first version of this test fell into: with the fracture touching
the domain boundary, the y-total came out 800 against an expected 1000 on
*both* grids. That is correct behaviour — shadow reaching a boundary face has
no cell pair, so no projection record exists — but it looks like a defect. The
test now insets the fracture on every side.

CPEDFM cannot be held to the same invariant. Its flip-face rule seals every
host-adjacent face whose neighbour centroids straddle the plane, which is a
cut of the connectivity graph rather than a shadow, so the sealed area
legitimately depends on the grid (2400 m² over 96 faces uniform, 2700 m² over
105 non-uniform). What is asserted instead: no over-assignment, multipliers in
[0, 1], and at least as much sealing as the pure shadow.

**R2 — `preparePEDFMExport` transmissibility reconstruction.** It rebuilds
projected `T` from harmonic permeability and centroid distance to work around a
legacy defect. That reconstruction should be re-derived under variable spacing.
*Mitigation: folded into Phase 4.*

**R3 — Scale.** Eagle Ford has ~153 fracture planes. Fracture-cell count,
`*SCONNECT` count, and seal count all need measuring before the full case is
attempted. *Mitigation: Phase 3 uses one fracture; report counts before scaling.*

**R4 — `*VOLMOD` accumulation.** Specifying a block twice multiplies rather than
replaces. *Mitigation: assert each block appears exactly once before writing.*

**R5 — Depth and datum.** Eagle Ford's `physdim` starts at z=0, which is not a
physical reservoir depth and will affect initialization and gravity.
*Mitigation: require an explicit datum depth for fracture cases.*

**R6 — Non-hexahedral corner-point cells.** The `*CORNERS` rebuild assumes 8
unique nodes per cell; pinched or degenerate cells break that.
*Mitigation: detect and fall back to requiring `grdecl`.*

---

## 8. Defects to fix on the way in

Both existing grid scripts write `*PERMI` as the keyword in all three
permeability files — [`NorneOut2CMGOMO.m:51,56`](NorneOut2CMGOMO.m) and
[`writeCMGjohansenFullValidation.m:50,55`](writeCMGjohansenFullValidation.m).
The Johansen script additionally writes `PERMX` *values* into the PERMJ and
PERMK files, so its `PERMZ = 0.1*K` anisotropy is silently lost.

The committed include files in `modules/CMG input files/` have correct
`*PERMJ` / `*PERMK` headers and a real 0.2 K-ratio, which means **they cannot be
reproduced by the scripts in this repo** — they were corrected by hand at some
point. Reconcile before either script becomes the basis for the new writer.

Also: `writeCMGCPEDFM` returns no `modelType` field when called directly (only
the wrappers set it), and `getEmbeddedFractureMultipliers` compares an absolute
area against `tol`, which is not scale-invariant.

---

## 9. Verification harness

The MCP bridge in `tools/cmg-mrst-mcp` closes the loop automatically: one call
exports from MATLAB, assembles a host deck, runs IMEX, and returns parsed
errors. `tools/matlab-mrst-mcp` runs MATLAB and MRST tests headlessly. Wire both
in from Phase 1, not at the end, so every keyword change is checked against the
real simulator immediately.

Three levels, in increasing strength:

1. **Deck acceptance** — IMEX reaches `Normal Termination`. Catches syntax and
   addressing errors. Already automated. *Necessary, far from sufficient.*
2. **Static reconciliation** — pore volume, block count, connection count and
   transmissibility values read back and compared against MRST.
3. **Solution match** — single-phase incompressible pressure, CMG vs MRST, cell
   by cell. **This is the only test that proves the model was transferred.**

A note carried over from the Eclipse work: the current CMG include writes NNCs
and fracture properties as comments, so the IMEX run recorded in
`tests/mcp_work/` exercised only `*NULL` and `*TRANS**MOD`. Level-1 passes there
did *not* validate connectivity. Do not let that pattern recur — a passing run
means nothing about cards that are commented out.

---

## 10. Open questions

1. Do IMEX, GEM and STARS impose a practical ceiling on `*SCONNECT` record
   count? Relevant to R3 and to unstructured models at scale.
2. Preferred seal mechanism when a pair is *both* a logical neighbour and needs
   a multiplier — `*TRANS**MOD` or negative `*SCONNECT`? Both are valid;
   pick one and apply it consistently.
3. Should the connection-list route be offered for *structured* grids too, as a
   fallback when padding-layer addressing is awkward? It would be uniform, at
   the cost of losing I-J-K meaning.
4. Fracture cells in a corner-point model: padding layers need ZCORN extension
   below the true bottom surface, as the Eclipse writer does. Confirm CMG
   accepts degenerate padding cells whose volume is overridden by `*VOLMOD`.
