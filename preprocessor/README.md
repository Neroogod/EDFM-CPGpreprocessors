# preprocessor

Preprocessor that exports EDFM, pEDFM, and CPEDFM models to commercial
reservoir simulators. All three techniques share one non-intrusive
embedded-fracture representation and the same Eclipse/CMG writing path.

## writeCMG: one entry point for any MRST grid

`writeCMG` takes an MRST grid and writes a CMG dataset. The grid family and
the fracture method are both detected from `G` and reported in the output,
so the same call works whether or not fractures are present.

```matlab
mrstModule add preprocessor
info = writeCMG(G, 'case.dat', 'simulator', 'IMEX', 'datum', 3000);
info = writeCMG(G, 'case.dat', 'grdecl', grdecl);   % corner-point grids
```

| MRST grid | Detected as | CMG grid |
|---|---|---|
| `cartGrid` | `cart` | `*GRID *CART` + `*DI/*DJ/*DK *CON` |
| `tensorGrid` (e.g. `meshHFsystem`) | `vari` | `*GRID *VARI` + `*DI *IVAR` … |
| `processGRDECL` | `corner` | `*GRID *CORNER` + `*COORD`/`*ZCORN` |
| `pebi`, `triangleGrid`, tets | `unstructured` | connection-list pseudo-grid |

CMG has no unstructured grid type, so unstructured grids export as an
`Nc x 1 x 1` `*VARI` pseudo-grid whose whole connectivity is an `*SCONNECT`
list, with `*PORVOL` restoring true pore volume and `*DEPTH-TOP` the true
depth. This is exact for a TPFA simulator, which consumes only
`(cell pair, transmissibility)`.

Two files are written: `case.inc` holding the reservoir description, and
`case.dat`, a dataset that declares the grid and `*INCLUDE`s it. **No fluid
properties are invented.** `'hostDeck'` selects what accompanies the include:
`'stub'` (default) is grid only, ready to open in Builder and populate there;
`'none'` is the include alone, for dropping into a deck of your own;
`'runnable'` adds a generic placeholder fluid model so the deck can be
executed for testing, and must be requested explicitly.

`*GRID` is written in the `.dat`, never in the include — Builder rejects the
fundamental grid definition when it arrives through an `*INCLUDE` and renders
an empty grid, though IMEX accepts it either way.

Fractures are detected from `G` and exported too. On a **structured or
corner-point** grid, fracture cells become padding k-layers — each placed in
the column of the matrix cell it intersects, with `*NULL` on unused slots,
`*PORVOL` carrying true pore volume and `*DEPTH-TOP` true depth. Every NNC
becomes a `*SCONNECT` record in md·m; pEDFM/CPEDFM seals become `*TRANS?*MOD`
on logical neighbours and negative `*SCONNECT` (a multiplier) on irregular
pairs. The padding region's automatic connections are zeroed, so fracture
cells communicate only through the explicit list.

On an **unstructured** grid there is no padding at all: the connection-list
model already gives every cell its own block, so matrix and fracture cells are
numbered straight through and seals fold into the connection
transmissibilities directly.

`writeEclipse` is the Eclipse sibling, consuming the same model — matrix and
fracture models alike, with connections as `NNC`, pore volume as `PORV`,
depths as `DEPTH` and seals as `MULTX/Y/Z`. Wells go out as a commented
`SCHEDULE` block, since they belong after the sections this deck stops short
of; their `COMPDAT` connection factors are MRST's well indices written
explicitly, never left to a radial model on placeholder geometry.

**Status:** verified against IMEX 2024.30 and MRST's deck reader.

| grid family | matrix | EDFM / pEDFM / CPEDFM |
|---|---|---|
| uniform Cartesian | CMG + Eclipse | CMG + Eclipse |
| variable Cartesian | CMG + Eclipse | CMG + Eclipse |
| corner-point | CMG + Eclipse | CMG + Eclipse |
| unstructured | CMG + Eclipse | CMG + Eclipse (connection list) |

The connection list is proven to carry flow, not merely to parse: see
`tests/test33_sconnect_is_live.m`. IMEX accepts a `*SCONNECT` pointing at a
null block without complaint, so deck acceptance is not evidence on its own.


The key point: after the projection stage, a CPEDFM model is completely
described by three keyword-expressible ingredients --- fracture cells
(pore volume, depth, properties), an NNC list (matrix-fracture,
fracture-fracture, and projection NNCs with precomputed
transmissibilities), and per-face matrix-matrix multipliers (the
seals). All three are emitted from the face-indexed data
(`G.nnc.pMMfaces` / `pMMarea`) that the projection functions already
produce, so no algorithmic work happens at export time.

## Files

```
writeEclipseCPEDFM.m        Eclipse-format .DATA deck writer (VERIFIED)
writeEclipseEDFM.m          EDFM adapter + Eclipse deck writer
writeEclipsePEDFM.m         pEDFM adapter + Eclipse deck writer
writeEclipseEmbeddedFractureModel.m  unified EDFM/pEDFM/CPEDFM entry point
writeCMGCPEDFM.m            CMG include-file writer (data verified,
                            keyword syntax to check against CMG manuals)
writeCMGEDFM.m              EDFM adapter + CMG include writer
writeCMGPEDFM.m             pEDFM adapter + CMG include writer
writeCMGEmbeddedFractureModel.m      unified CMG entry point
prepareEDFMExport.m         removes projection NNCs/multipliers
preparePEDFMExport.m        converts pair-indexed projections to faces
getEmbeddedFractureMultipliers.m     shared face multiplier assembly
buildCoarseCPEDFM.m         coarse connection-list model from a
                            barrier-aware partition (screening models)
writeEclipseCoarseCPEDFM.m  coarse-model Eclipse deck writer (VERIFIED)

writeCMG.m                  unified MRST -> CMG entry point (any grid family)
writeEclipse.m              unified MRST -> Eclipse entry point (any grid)
buildCMGModel.m             G -> simulator-neutral cmgModel contract
writeCMGModel.m             cmgModel -> CMG reservoir-description cards
writeEclipseModel.m         cmgModel -> Eclipse GRID-section keywords
cmgGridHeader.m             the *GRID/*KDIR cards (must live in the .dat)
cmgDetectGrid.m             grid-family detection (cart/vari/corner/unstr)
cmgDetectFractureModel.m    fracture-method detection (none/EDFM/pEDFM/CPEDFM)

tests/test20_eclipse_roundtrip.m
tests/test21_coarse_export.m
tests/test22_edfm_pedfm_export.m
tests/test30_cmg_grid_export.m       detection + matrix export, all families
tests/test32_cmg_fracture_export.m   EDFM/pEDFM/CPEDFM export + conservation
tests/test33_sconnect_is_live.m      proves the connection list carries flow
tests/test34_cornerpoint_fracture_cmg.m  corner-point + fractures
tests/test35_eclipse_matrix_export.m Eclipse matrix export, round-trip
tests/test36_eclipse_fracture_export.m Eclipse EDFM/pEDFM/CPEDFM export
tests/test37_projection_nonuniform.m projection vs analytic, uneven spacing
tests/test38_unstructured_fracture_export.m  PEBI + fractures, no padding
tests/test39_eclipse_wells.m         WELSPECS/COMPDAT with explicit WI
tests/test31_eagleford_matrix_cmg.m  Eagle Ford *VARI export
```

## Coarse-model export ("screening models")

`buildCoarseCPEDFM` + `writeEclipseCoarseCPEDFM` turn a CPEDFM model
into a genuinely coarse Eclipse deck: the coarse blocks of a
barrier-aware partition (CPEDFM-MS/partitionCPEDFM_MS) become the
cells of a connection-list model - block pore volumes and depths via
EDIT overrides, aggregated interface transmissibilities as NNCs, and
wells with explicitly summed connection factors (no Peaceman on the
placeholder geometry). Because no barrier-aware block straddles a
sealed surface, sealed fine faces contribute zero to every coarse
connection and the compartmentalization survives coarsening exactly;
the tiny fracture-mediated leakage paths of the fine model (matrix ->
sealing fracture -> matrix) are retained at their true magnitude
(~1e-9 of matrix transmissibility in the test case).

test21 verifies (skewed faulted CPG, sealing fracture, 1,742 fine
cells -> 46 blocks): pore volume conserved to machine precision; the
coarse matrix-connectivity graph reproduces the fine model's 2 macro
compartments with zero significant cross-compartment connections;
under the fine Dirichlet drive the coarse model reproduces the fine
compartment plateau pressures to 1e-4; and the exported deck
round-trips through the MRST parser with all PORV/DEPTH/NNC values
matching. Aggregated (sum-of-fine-face) transmissibilities are the
standard screening-level choice; flow-based upscaled values can be
substituted into `coarse.T` without changing the export format.

## Usage

```matlab
% after: EDFMGridCPG -> fracturematrixShaleNNC3D -> pMatFracNNCsCPEDFM_CPG
wells = struct('name', {'INJ','PROD'}, 'type', {'rate','bhp'}, ...
    'val', {50*meter^3/day, 100*barsa}, 'comp', {'WATER','OIL'}, ...
    'cells', {injCells, prodCells});
info = writeEclipseCPEDFM(Gp, Gp.rock, grdecl, wells, 'CASE.DATA');
info = writeCMGCPEDFM(Gp, Gp.rock, grdecl, 'cpedfm_cmg.inc', ...
    'simulator', 'IMEX', 'cmgVersion', '2024.10');

% EDFM: original matrix-fracture/fracture-fracture NNCs, no face seals
info = writeEclipseEDFM(Ge, Ge.rock, grdecl, wells, 'EDFM.DATA');
info = writeCMGEDFM(Ge, Ge.rock, grdecl, 'edfm.inc');

% pEDFM: projected NNCs plus partial projection multipliers
info = writeEclipsePEDFM(Gp, Gp.rock, grdecl, wells, 'PEDFM.DATA');
info = writeCMGPEDFM(Gp, Gp.rock, grdecl, 'pedfm.inc');

% Unified entry points
info = writeEclipseEmbeddedFractureModel('EDFM', Ge, Ge.rock, ...
    grdecl, wells, 'EDFM.DATA');
info = writeCMGEmbeddedFractureModel('pEDFM', Gp, Gp.rock, ...
    grdecl, 'pedfm.inc', 'simulator', 'GEM');
```

For every CMG writer, `simulator` accepts `IMEX`, `GEM`, or `STARS`
(case-insensitive). The selected profile changes the target metadata and
integration guidance: IMEX is black-oil/conventional multiphase, GEM is
compositional/EOS, and STARS is thermal/reactive. Alongside the include,
the writer emits simulator-tagged fracture-block, NNC, multiplier, and
fault-seal CSV files plus a JSON manifest. These tables contain the
authoritative numerical preprocessor output and can be imported through
Builder or mapped to the exact connection cards supported by the selected
CMG release. Set `'writeTables', false` to suppress the sidecar package.

## Model-specific behavior

| Model | Exported fracture connections | Matrix-face modification |
|---|---|---|
| EDFM | matrix-fracture and fracture-fracture NNCs | none |
| pEDFM | EDFM NNCs plus projected matrix-fracture NNCs | partial projection multipliers |
| CPEDFM | EDFM NNCs plus continuous projected NNCs | continuous/full-face multipliers |

Standard Cartesian pEDFM stores projected matrix pairs rather than physical
face IDs. `preparePEDFMExport` resolves unique pairs and records
`G.nnc.pMMfaces/pMMarea`. On a faulted grid a pair can share multiple split
faces; the default is a hard ambiguity error. The optional
`'ambiguousFacePolicy','areaWeighted'` distributes the projected area in
proportion to the candidate face areas. Recording face IDs during pEDFM
projection remains the preferred production workflow.

## Representation (the non-intrusive EDFM pattern)

* Fracture cells are appended as padding k-layers below the matrix
  grid (placeholder geometry). Their physics enters through EDIT
  overrides: `PORV` (true fracture pore volume) and `DEPTH` (true
  centroid depth, for gravity). All geometric transmissibilities of the
  padding, and the matrix-padding interface, are zeroed (`TRANX/Y/Z`),
  so fracture cells communicate only through NNCs.
* Every MRST NNC becomes one Eclipse `NNC` record, with the
  transmissibility converted to METRIC units (cP.m3/day/bar).
* CPEDFM seals between logical neighbors become `MULTX/MULTY/MULTZ`
  records (`EQUALS` cards); seals on fault-juxtaposed pairs, which
  Eclipse connects through internally generated NNCs, become `EDITNNC`
  records.
* The deck is self-contained (RUNSPEC through SCHEDULE, oil-water,
  METRIC, EQUIL initialization) and can also serve as an include-file
  source for a larger study.

Unstructured (CPEDFM-U) models export through the same writer in
connection-list form: since TPFA sees only (cell pair, transmissibility),
a PEBI model is reproduced exactly by a padding-style pseudo-grid whose
entire connectivity is the NNC list. (Native unstructured input is
available in INTERSECT/tNavigator; OPM Flow runs the Eclipse deck.)

## Verification status (2026-07-10, MATLAB R2021b)

`test20_eclipse_roundtrip.m` builds the CPEDFM-CPG test01 case (skewed
faulted corner-point grid, slanted sealing fracture), exports the deck,
reads it back with MRST's `deckformat` parser, and reproduces every
CPEDFM data structure from the deck alone:

* grid round-trip: 1742 active cells = matrix + fracture cells. PASS
* NNC records: 378/378, ijk exact, transmissibility rel. diff 2.3e-10
  (unit conversion round-off). PASS
* MULT seals: 236/236 records, exact. PASS
* fracture PORV/DEPTH overrides: rel. diff 4e-9 / 5e-7 m. PASS
* EDITNNC: not exercised by this case (no fault-juxtaposed seals);
  the writer path exists and parses, but add a case with sealed
  fault-pairs before relying on it.

`test22_edfm_pedfm_export.m` builds one shared embedded-fracture grid and
exports EDFM, pEDFM, and CPEDFM variants. It verifies that EDFM has no
projection multipliers, both projected methods add NNCs and face
multipliers, Eclipse decks round-trip their NNC counts, and the Eclipse and
CMG writers consume identical normalized EDFM/pEDFM data. MATLAB R2021b
result: EDFM 20 NNC/0 MULT, pEDFM 58 NNC/28 MULT, and CPEDFM 68 NNC/48
MULT. PASS.

Not yet verified (no simulator available on this machine):
* an actual Eclipse / OPM Flow / CMG run of the exported deck. The
  deck is syntactically parseable by MRST's Eclipse-format reader and
  numerically faithful; run it through OPM Flow (free) or Eclipse at
  LSU and compare against the MRST solution as the final acceptance
  test. The group's earlier CMG cross-validation script
  (CPEDFMcornerpointaug24/Updates_8_23/JohansenCMGval.m) is a template
  for the CMG side.

Known limitations:
* Wells must complete in matrix cells (fracture-cell completions would
  need WELSPECS on padding addresses --- straightforward, not written).
* The schedule writer covers water-rate injectors and BHP producers;
  extend the SCHEDULE section for other controls.
* CMG keyword syntax (\*NULL \*IJK, \*TRANSI/J/K \*MOD, explicit-NNC
  facility) varies by simulator and version --- the include file
  carries the verified numbers with banner comments; check cards
  against your manuals or import through Builder.
