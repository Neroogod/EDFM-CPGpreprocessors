# CPEDFM-export: Implementation Log

A step-by-step record of how the commercial-simulator preprocessor
(fine-scale Eclipse/CMG export + coarse screening-model export) was
designed, implemented, debugged, and verified (2026-07-10/11, MATLAB
R2021b). Companion to `CPEDFM-CPG/cpedfm-cpg.md`,
`CPEDFM-U/cpedfm-u.md`, `CPEDFM-MS/cpedfm-ms.md`, and to the
manuscript in `SPEJ_CPEDFMCPGandUnstructuredgridspreprocessorforeclipseandcmg/`.

---

## Step 0 — Feasibility analysis: three keyword-expressible ingredients

The starting question was whether CPEDFM/pEDFM could reach Eclipse and
CMG at all, given that commercial simulators cannot be modified. The
answer came from decomposing a *fully processed* CPEDFM model: after
the geometric work is done, it reduces to exactly three ingredients,
each already expressible in standard keyword languages:

1. **Fracture cells** (pore volumes, depths, properties) — the
   established non-intrusive EDFM representation: auxiliary grid cells
   with `PORV`/`DEPTH` overrides (Xu & Sepehrnoori's line of work,
   which carried conventional EDFM into commercial practice).
2. **The NNC list** (matrix-fracture, fracture-fracture, and the
   CPEDFM projection connections with precomputed transmissibilities)
   — Eclipse `NNC` records.
3. **The seals** — the new ingredient no published preprocessor
   delivers: per-face multipliers η_f, mapping to `MULTX/MULTY/MULTZ`
   for logical-neighbor faces and `EDITNNC` for fault-juxtaposed pairs
   (which Eclipse connects through internally generated NNCs).

The decisive enabling fact: the CPEDFM implementation records its
projections **face-indexed** (`G.nnc.pMMfaces`/`pMMarea`) — a design
made in CPEDFM-CPG for multiplier robustness — so the export is a
formatting exercise with no algorithmic work remaining.

## Step 1 — Choose the validation strategy (no simulator available)

`flow` (OPM), Docker, and any Eclipse/CMG installation were all absent
from the machine. Decision: validate by **round-trip through an
independent Eclipse-format parser** — MRST's `deckformat` reader —
reproducing every exported data structure from the deck alone. This is
arguably stronger than a simulation comparison for verifying deck
*semantics* (values land in the right keywords with the right units),
while the final acceptance test (an actual Eclipse/OPM/CMG run vs the
MRST solution) is documented as the remaining step for LSU. Parser
keyword support (`NNC`, `EQUALS`, `EDIT`-section arrays) was grepped
in `deckformat/deckinput` before committing to the format.

## Step 2 — `writeEclipseCPEDFM.m` (fine-scale deck writer)

Design decisions, in order:

- **Addressing**: matrix cells keep their logical (i,j,k); the N_f
  fracture cells are assigned slots in ceil(N_f/(NX·NY)) padding
  k-layers appended below the grid. ZCORN is extended from the actual
  bottom surface (per-pillar, avoiding inverted padding cells); COORD
  needs no change; ACTNUM activates exactly the used slots.
- **Isolation**: EDIT `TRANX/TRANY/TRANZ` zeroed over the padding
  region and the matrix→padding interface (BOX-limited arrays), so
  fracture cells communicate only through explicit NNCs.
- **Units**: transmissibilities converted with MRST unit constants to
  METRIC cP·m³/(day·bar) — `convertTo(T, centi*poise*meter^3/(day*barsa))`
  rather than a hand-derived factor.
- **Seals**: per sealed face, classify the neighbor pair by logical
  offset — unit step in one index → `MULTX/Y/Z` (EQUALS record at the
  lower-index cell); anything else (fault juxtaposition) → `EDITNNC`.
- **Completeness**: the deck is self-contained (RUNSPEC…SCHEDULE,
  oil-water, wells, timesteps) so it runs as generated and doubles as
  an include-file source.

## Step 3 — test20 first run: three parser findings

The round-trip test (test01's skewed faulted CPG + sealing fracture)
passed grid/NNC/MULT checks immediately — 378 NNCs with ijk exact and
transmissibility differences 2.3e-10 (unit round-off), 236 MULT seals
exact — but produced three findings that reshaped details:

1. **EDIT-section EQUALS records are not materialized** by MRST's
   reader (the fracture `PORV`/`DEPTH` overrides never appeared).
   Fix: write them as BOX-limited plain arrays over the padding
   region — equally standard Eclipse, parser-friendly.
2. **Explicit `PRESSURE`/`SWAT` solution arrays warn** in the parser
   (expected lengths differ between readers). Fix: `EQUIL`
   initialization with the contact placed below the model — the
   simulator-standard choice anyway.
3. **MRST merges EDIT keywords into `deck.GRID`** (there is no
   `deck.EDIT` struct; `readEDIT`'s get_state/set_state operate on
   `deck.GRID`). The test's field lookups were corrected accordingly.

After these fixes: all five checks green (grid 1,742 active cells
exact; NNC; MULT; PORV rel. 4e-9; DEPTH 5e-7 m).

## Step 4 — `writeCMGCPEDFM.m` (honesty-first CMG writer)

CMG keyword syntax varies by simulator (IMEX/GEM/STARS) and version,
and could not be machine-verified here. Decision: the CMG writer emits
the **identical verified data** in an include file with best-effort
cards (`*NULL *IJK` for unused padding slots, `*TRANSI/J/K *MOD`
multipliers) and banner-commented tables for the fracture-block
properties and NNC list (transmissibilities in CMG's md·m unit), with
the caveat stated prominently in the header and README: check card
syntax against the manuals or import through Builder. The group's
existing `JohansenCMGval.m` workflow is the natural validation
template. Sanity-run on the test case: 378 NNCs, 236 multipliers,
32 kB include.

The writer now accepts `'simulator', 'IMEX'`, `'GEM'`, or `'STARS'`
and an optional `'cmgVersion'` label. `getCMGSimulatorProfile.m` selects
black-oil, compositional/EOS, or thermal/reactive integration guidance,
respectively. Each run also creates simulator-tagged CSV tables for
fracture blocks, NNCs, multipliers, and fault seals, plus a JSON manifest
that records the target, dimensions, units, model type, and record counts.
This keeps the numeric preprocessing common and verified while making the
handoff simulator-specific without claiming unverified proprietary cards.

## Step 5 — The geometry fix wakes up the EDITNNC path

When the project-wide COORD z-scaling bug was fixed (see
`CPEDFM-CPG/cpedfm-cpg.md`, Step 9), the corrected grid produced
fault-juxtaposed seals for the first time — test20 went from 0 to
**47 EDITNNC records**, and its check D promptly failed. Diagnosis: a
fresh grep showed MRST's parser does **not** implement EDITNNC at all
(it lands in `deck.UnhandledKeywords.EDIT`); an earlier interleaved
grep had suggested otherwise. Eclipse and OPM Flow support the
keyword natively. Fix: check D verifies the EDITNNC block by parsing
the deck text directly (regexp + sscanf on our own known format),
and reports that the MRST parser skips the keyword. All checks green
again (344 NNCs, 175 MULT, 47 EDITNNC on the corrected geometry).

## Step 6 — Coarse screening-model export

Motivated by a design question ("should multiscale be applied to the
preprocessor?") whose answer was no — commercial simulators cannot
consume MsRSB bases, and their AMG solvers are already barrier-aware
via strength-of-connection coarsening — but which surfaced the one
genuine crossover: the **barrier-aware partition as a coarsening for
screening models**.

- **`buildCoarseCPEDFM.m`**: from a barrier-aware partition
  (`CPEDFM-MS/partitionCPEDFM_MS`), aggregate block pore volumes,
  pv-weighted depths, and coarse connection transmissibilities as the
  sum of fine connections crossing each block interface. Two
  properties fall out of the construction: sealed fine faces carry
  zero and contribute nothing to any coarse connection, and no
  barrier-aware block straddles a barrier — so **compartmentalization
  survives coarsening exactly** — while the fracture-mediated leakage
  paths of the fine model are retained at their true (negligible)
  magnitude. Aggregated transmissibilities are the screening-level
  choice; flow-based values can replace `coarse.T` without format
  changes.
- **`writeEclipseCoarseCPEDFM.m`**: connection-list deck (Nc×1×1
  placeholder grid; PORV/DEPTH overrides; all geometric trans zeroed;
  NNC list is the entire connectivity). Wells are mapped to blocks
  with **explicitly summed connection factors** written in COMPDAT
  item 8 — a hard requirement, since a Peaceman calculation on the
  placeholder geometry would be meaningless.

## Step 7 — test21: one criterion refinement

First run: pore volume conserved to 1.5e-16, coarse plateaus matched
the fine Dirichlet solution to 1e-4, deck round-trip exact — but the
compartment check FAILED with "7 cross-compartment connections".
Diagnosis: those are the fracture-mediated matrix→sealing-fracture→
matrix paths, correctly exported at ~2e-9 of the matrix
transmissibility scale; an unweighted (and even a 1e-12-relative
thresholded) graph counts them as connectivity. The principled fix:
the compartment-topology check aggregates **matrix-matrix connections
only**, exactly parallel to the fine-side check (which also excludes
NNCs), while the export keeps the weak paths. Final: 46 blocks, 108
connections, 2 macro compartments preserved, all checks green.

## Final verified state

```
writeEclipseCPEDFM.m         fine-scale deck writer      (test20 PASS)
writeCMGCPEDFM.m             CMG include writer          (data verified;
                                                          syntax to check)
buildCoarseCPEDFM.m          coarse model builder        (test21 PASS)
writeEclipseCoarseCPEDFM.m   coarse deck writer          (test21 PASS)
prepareEDFMExport.m          EDFM normalization adapter
preparePEDFMExport.m         pair-to-face pEDFM adapter
getEmbeddedFractureMultipliers.m shared multiplier assembly
writeEclipseEDFM.m           EDFM Eclipse front end
writeEclipsePEDFM.m          pEDFM Eclipse front end
writeCMGEDFM.m               EDFM CMG front end
writeCMGPEDFM.m              pEDFM CMG front end
tests/test20_eclipse_roundtrip.m
tests/test21_coarse_export.m
tests/test22_edfm_pedfm_export.m
tests/output/CPEDFM_TEST20.DATA, CPEDFM_COARSE.DATA, cpedfm_cmg.inc
```

## EDFM and pEDFM exporter extension

The original CPEDFM writer has been generalized through model adapters.
`prepareEDFMExport` removes projected `fracpmat` connections and emits unit
matrix-face multipliers. `preparePEDFMExport` converts the standard
pair-indexed `pMMneighs` and appended projection areas into the face-indexed
contract used by the CPEDFM writers. Unique matrix pairs map directly;
fault-split ambiguity raises an error unless explicit area-weighted
distribution is requested. The Eclipse and CMG model-specific wrappers and
unified dispatch functions all reuse the same padding-cell, NNC, property,
unit-conversion, and seal-record implementations.

`test22_edfm_pedfm_export.m` passes in MATLAB R2021b: the shared case
exports 20 NNC/0 MULT for EDFM, 58 NNC/28 MULT for pEDFM, and 68 NNC/48
MULT for CPEDFM. All three Eclipse decks round-trip their NNC records, and
the EDFM/pEDFM CMG includes contain the same normalized counts and
multipliers. The run also exposed and fixed a legacy Cartesian pEDFM defect:
duplicate projected pairs were collapsed in `cells/area` but not in `T`;
the adapter now rebuilds the retained projected transmissibilities from the
retained union areas, harmonic permeability, and centroid distance.

Verified: grid round-trip exact; 344 NNC transmissibilities to
2.5e-10; 175 MULT seals exact; 47 EDITNNC records (text-verified);
fracture PORV/DEPTH to 4e-9 / 5e-7 m; coarse model PV-conservative to
machine precision, compartment-preserving, plateau-accurate to 1e-4,
and round-trip exact.

Not yet verified (requires software unavailable on this machine): an
actual Eclipse / OPM Flow / CMG run of the exported decks against the
MRST reference — the final acceptance test, and the content of the
manuscript's empty "Simulation Comparison" section.

## Lessons

1. **Separate the model from the solver.** The preprocessor ships
   discretization data; solution methods (multiscale) do not export.
   The one legitimate crossover is the barrier-aware partition reused
   as a *coarsening*, which became the screening-model feature.
2. **Validate against a parser you can run**, but know its gaps: the
   MRST reader's silent non-support of EDIT-EQUALS and EDITNNC would
   have masked real content had the checks trusted `isfield` alone;
   text-level verification of self-generated formats is a legitimate
   fallback.
3. **Never let a simulator compute well indices on placeholder
   geometry** — always write explicit connection factors.
4. **Weak connections are data, not noise**: the coarse export's
   fracture-mediated paths are physically faithful; diagnostics must
   distinguish topology checks (matrix-only) from the exported model
   (everything).
5. A dormant code path (EDITNNC) stayed untested until an unrelated
   bug fix changed the geometry — coverage follows the data, so
   test-case geometry deserves as much scrutiny as test logic.
