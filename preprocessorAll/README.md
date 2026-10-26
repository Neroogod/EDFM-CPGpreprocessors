# preprocessorAll

Export a **complete, runnable** simulator deck from an MRST case.

`modules/preprocessor` exports the reservoir description — geometry, rock,
fractures, connections — and stops there, deliberately inventing no fluid
data. Its own README lists the consequence: *"Still to be supplied, in
Builder or by hand: \*MODEL / component and PVT data, \*ROCKFLUID,
\*INITIAL, \*NUMERICAL, RUN section"*, and its verification section records
*"Not yet verified: an actual Eclipse / OPM Flow / CMG run of the exported
deck."*

This module closes both gaps. Give it the same triple an MRST run takes —
`(model, state0, schedule)` — and it writes a dataset that executes as it
stands. **The CMG path is verified by actually running IMEX 2024.30**, not
just by parsing.

## Usage

```matlab
mrstModule add preprocessor
addpath modules/preprocessorAll

% ... your ordinary MRST script, up to and including:
%     model = GenericBlackOilModel(G, rock, fluid, ...);
%     state0 = initResSol(...);
%     schedule = simpleSchedule(dt, 'W', W);

deck = buildFullDeck(model, state0, schedule);      % simulator-neutral
writeCMGDeck(deck,     'case.dat');                 % runnable IMEX dataset
writeEclipseDeck(deck, 'CASE.DATA');                % runnable Eclipse deck
```

For corner-point grids pass `'grdecl', grdecl`, exactly as in `writeCMG`.

| Option (on `buildFullDeck`) | Default | Meaning |
|---|---|---|
| `grdecl` | `[]` | corner-point description, required for `processGRDECL` grids |
| `simulator` | `'IMEX'` | `IMEX`, `GEM` or `STARS` |
| `pRange`, `nP` | from `state0`, 20 | PVT sampling window and resolution |
| `nSat` | 21 | relative-permeability samples |
| `fitTol` | `1e-3` | when a property must be collapsed to a two-parameter deck form, the largest residual accepted before a warning is recorded |
| `startDate` | `[2026 1 1]` | recurrent-section start |

## What makes this possible: sampling the fluid

An MRST fluid is a set of anonymous functions — `bW(p)`, `muO(p)`,
`krW(s)`. Eclipse and CMG both want tables. `sampleFluidTables` evaluates
the handles on a pressure/saturation grid and, crucially, **measures how
faithfully a deck format can represent the result**.

That matters because some deck entries are two-parameter models. Eclipse's
`PVTW` and CMG's `*CW`/`*CVW` assume `b(p)` is log-linear; `*CO` assumes the
same for oil. `sampleFluidTables` fits that form, reports the worst relative
error over the sampled range, and `buildFullDeck` turns anything above
`fitTol` into an entry in `deck.warnings` — which is then printed into the
deck banner itself. A fluid the format cannot carry is visible in the file,
not silently approximated.

For the `initSimpleADIFluid` exponential fluids used in the tests the
residual is **2.1e-16**, i.e. exact.

## What is and is not carried over

**Carried**: grid geometry, porosity, permeability, fracture cells, NNCs and
transmissibility multipliers (all via the existing `buildCMGModel`), PVT,
relative permeability, capillary pressure, rock compressibility, initial
pressure and saturations cell by cell, wells with MRST's own well indices,
controls, and the report-step schedule.

**Not carried**: anything with no deck equivalent — solver tuning, MRST-side
discretization choices, non-tabulable property functions.

Initialization is written as explicit per-cell `PRESSURE`/`SWAT` (Eclipse)
and `*PRES`/`*SW` (CMG) rather than as an equilibration region, because
`state0` is what the MRST run actually started from; re-deriving it from
contacts would introduce a difference nobody asked for.

## Two things that required care

**Dead oil cannot use a saturated PVT table.** IMEX rejects a multi-row
`*PVT` table whose `Rs` is constant (`SETTBP`: *"entry 2 for Rs is constant
or decreasing but should be increasing"*), and a dead oil's `Bo` *decreases*
with pressure, which the table also forbids. Two-phase oil-water therefore
uses the single-row form at a reference pressure with `*CO`/`*CVO` carrying
the pressure dependence — the form IMEX's own oil-water templates use
(`TPL/frr/mxfrr003.dat`). Live oil uses the multi-row saturated table
(`TPL/cmb/mxcmb001.dat`), with monotonicity checked in MATLAB so a
non-monotone fluid is reported where it can be acted on.

**Fracture cells do not live at the end of the block array.** On a fracture
export `buildCMGModel` places each fracture cell in the padding layers under
the column of the matrix cell it intersects, leaving unused slots inactive —
148 MRST cells become a 8×8×4 = 256 block grid with 108 inactive blocks.
Appending initial-state values to the end of the matrix array would put
fracture pressures on the wrong blocks, silently, since the deck would still
be well formed and still run. `deckBlockArray` is the single place that
applies the `grid.matIJK` / `grid.fracIJK` mapping. (In test52's single-plane
case the two orderings happen to coincide; the mapping is used because it is
the correct general form, not because that test distinguishes them.)

## Single-phase water

Neither simulator has a single-phase option. A water-only MRST model is
exported as oil-water with `So = 0` everywhere and the oil given the
*water's own* properties. With no oil in place and none injected the oil
equations cannot influence the result, so this is a faithful mapping rather
than an invented second phase. The deck says so in a comment at the point it
happens.

## Files

```
buildFullDeck.m        (model, state0, schedule) -> simulator-neutral deck
sampleFluidTables.m    fluid function handles -> PVT / relperm tables + fit residuals
deckBlockArray.m       scatter per-cell values onto the exported block grid
writeCMGDeck.m         -> complete runnable IMEX dataset (.dat + .inc)
writeEclipseDeck.m     -> complete runnable Eclipse deck (.DATA + .grdecl)
tests/test50_cmg_runnable_ow.m       export -> run IMEX -> compare with MRST
tests/test51_eclipse_roundtrip.m     export -> parse -> compare every section
tests/test52_fractured_full_deck.m   EDFM case -> run IMEX
```

The reservoir description is not rewritten: `writeCMGModel`,
`writeEclipseModel`, `cmgGridHeader` and `buildCMGModel` are reused from
`modules/preprocessor`, so every grid family and fracture method it supports
is supported here too.

## Verification (MATLAB R2021b, IMEX 2024.30, 2026-08-01)

All three tests pass.

**test50 — oil-water waterflood, exported and run.** 10×10×3 grid at 2000 m,
two wells, 360 days.

| | result |
|---|---|
| Sampled tables vs MRST fluid functions | **0.00e+00** (Bw, Bo, muO, krw, krow) |
| IMEX run | **Normal Termination, 0 errors** |
| Material balance (o/w/g) | **1.000 / 1.000 / 1.000** |
| Cumulative water injected vs schedule | **exact** (36 000 m³) |
| Cumulative oil: IMEX vs MRST | 30 846 vs 31 463 m³ — **2.0 %** |
| Mean reservoir pressure: IMEX vs MRST | 20 614 vs 20 659 kPa — **0.2 %** |

The residual difference is expected: IMEX and MRST use different timestep
control and implicitness, so this is agreement between two independent
simulators on the same physics, not a round-trip identity.

**test51 — Eclipse deck round-trip.** 8×8×3 case read back with MRST's
`readEclipseDeck`: RUNSPEC phases and dimensions correct; `PVTW`/`PVDO` to
**4.3e-8**; `SWOF` to **1.1e-16**; `PRESSURE` and `SWAT` **exact**; 48
`COMPDAT` rows with connection factors matching MRST's well indices to
**4.2e-10**; `TSTEP` summing to the schedule's 300 days exactly.

**test52 — embedded fractures, exported and run.** 148-cell EDFM grid
(128 matrix + 20 fracture, 20 NNCs) exported as 256 blocks and run:
**Normal Termination, 0 errors**; total pore volume IMEX 25 603.0 m³ vs MRST
25 603.2 m³ (**8.8e-6**), so the fracture pore volume survived; cumulative
injection exact.

## Known limitations

* **Eclipse is verified by parser, not by a simulator run.** No Eclipse or
  OPM Flow binary is on this machine. The deck round-trips through MRST's
  reader with every section intact, and the CMG path — which *is* run —
  shares `buildFullDeck` and `sampleFluidTables`, so the physics content is
  the same. Run it through OPM Flow (free) before relying on it.
* **Three-phase requires dissolved gas.** A three-phase model with
  `disgas = false` has no single-table IMEX form and raises a clear error;
  enable disgas or drop the gas phase.
* **Vaporized oil (`vapoil`) is rejected**, not silently dropped.
* **Wells must complete in matrix cells.** Fracture-cell completions are
  exported at their padded addresses with a warning; verify them in the
  target simulator.
* **Controls covered**: rate and BHP, injectors and producers, with the
  schedule's control periods re-declared when they change. Other MRST
  control types (group control, VFP tables, `lims` beyond the primary
  control) are not translated.
* **Capillary pressure** is exported when the fluid defines `pcOW`/`pcOG`
  and written as zero otherwise; MRST fluids frequently omit it.
* CMG keyword syntax varies by simulator and release. The IMEX path is
  exercised against 2024.30; GEM and STARS select different metadata but
  their fluid sections are not written by this module.
