# EDFM Grid Builder — Web App Implementation Plan

A browser front end for the MRST → CMG/Eclipse preprocessor. The user supplies
a grid and a set of fractures, walks the pipeline one stage at a time,
inspects what each stage produced, and downloads a simulator deck.

Companion to `modules/preprocessor/preprocessor_implementation.md`, which
covers the export engine this app drives.

Drafted 2026-07-27.

**Decisions taken:** local-first now, hosted on merslab.com later, alongside
the Flash Calculation Predictor. The app produces files only — it does **not**
run IMEX.

---

## 1. What the user gets

The workflow mirrors `JohansenGridpre.m` and `ValidationCase.m`, because that
is the pipeline that works:

```
  1  matrix grid            processGRDECL / cartGrid / tensorGrid
  2  fracture planes        fracplanes: points, aperture, poro, perm
  3  EDFM gridding          EDFMGridCPG  |  EDFMshalegrid
  4  fracture-matrix NNCs   fracturematrixShaleNNC3D
  5  fracture-fracture NNCs fracturefractureShaleNNCs3D
  6  projection (optional)  pMatFracNNCs3D | pMatFracNNCsCPEDFM_CPG
  7  export                 writeCMG | writeEclipse* | MRST .mat
```

Every stage is a checkpoint: look at what it produced, download it, continue
or stop. **Stage 1 alone is a valid product** — "spit out the matrix grid
first" is just running the export after stage 1, which `writeCMG` already
supports today (`fractureModel: matrix`).

---

## 2. The hosting target decides the architecture

merslab.com is a **static site**. The Flash Calculation Predictor states
plainly that "all computation runs privately in your browser — your data never
leaves your computer", and a fetch of the page shows no API endpoints, no
POSTs, no compute backend. Model weights ship as static assets and inference
runs in JS/WASM.

**A static host cannot run MATLAB.** So "host it on merslab later" is not a
deployment detail to defer; it decides what the app can be. Two things follow.

**MATLAB licensing is not the blocker — the absence of a server is.** Even
with MATLAB Compiler, a compiled binary needs a machine to run on. merslab has
none. (Worth knowing anyway: MATLAB Compiler cannot `run` or `eval`
user-supplied `.m` files at runtime, so an "upload your MRST script" feature
could never work in a compiled deployment regardless.)

**Therefore the product splits in two**, along a line that already exists in
the code.

---

## 3. The split: what needs MATLAB and what does not

| Needs MATLAB + MRST | Portable to the browser |
|---|---|
| `EDFMGridCPG` / `EDFMshalegrid` — fracture-plane × cell intersection | GRDECL parsing (COORD/ZCORN/PERM/PORO) |
| `fracturematrixShaleNNC3D` — areas, `calcdavg`, connectivity index | Cartesian / tensor grid construction |
| `fracturefractureShaleNNCs3D` — plane-plane intersection | `buildCMGModel` — addressing, unit conversion, array assembly |
| projection stages (pEDFM / CPEDFM) | `writeCMGModel` — card formatting |
| `processGRDECL` geometry processing | `writeCMG` — deck assembly, manifest |

Only the **geometry** genuinely needs MATLAB. Everything downstream of it is
array arithmetic and string formatting.

This is not a lucky accident — `buildCMGModel` was deliberately separated from
`writeCMGModel` so that one produces a simulator-neutral struct and the other
only formats it. **That seam is exactly where the browser/MATLAB boundary
falls.** Porting `writeCMGModel` to JavaScript is a transcription job, not a
redesign.

### Two deliverables

**A — Local app.** Full pipeline, drives the installed MATLAB. For building
fractured grids from scratch. Runs at `localhost`, nothing uploaded.

**B — merslab browser app.** Fully static, no server, no MATLAB. Two
capabilities:

- **B1 — grid → deck.** Upload a GRDECL (or describe a Cartesian grid), get a
  CMG or Eclipse matrix deck. Needs no MATLAB at any point, and is the broader
  win: anyone with a corner-point grid gets a Builder-ready deck with no
  install.
- **B2 — processed model → deck.** Upload the portable JSON that A exports
  (§4) and get the fractured deck, including the choice of simulator, host-deck
  mode and output format — without re-running MATLAB.

**Long term — C.** Port the EDFM geometry itself to JS/WASM, making the full
pipeline browser-native. Large, but there is a strong lead: `CPEDFM-U`'s
`fracPolyCellIntersect.m` already implements fracture-polygon × cell
intersection by half-space clipping in ~100 portable lines, replacing the
995-line hexahedra-only version. It is exact for convex cells, which covers
every target grid family. That is the kernel a WASM port would be built on.

---

## 4. The bridge: a portable model format

The interchange artifact is `cmgModel` — the struct `buildCMGModel` already
produces — serialised as JSON:

```jsonc
{
  "meta":  { "gridFamily": "corner", "fractureModel": "EDFM",
             "unitSystem": "SI", "transUnit": "md.m", "counts": { ... } },
  "grid":  { "kind": "corner", "dims": [100,100,33],
             "coord": [...], "zcorn": [...] },        // or di/dj/dk
  "props": { "permi": [...], "permj": [...], "permk": [...],
             "poro": [...], "porvol": [...], "active": [...] },
  "conn":  [[i1,j1,k1, i2,j2,k2, trans], ...],        // md.m
  "mult":  [[dir,i,j,k, value], ...],
  "isolate": [ { "keyword": "*TRANSK", "range": [...], "value": 0 } ],
  "wells": [ ... ]
}
```

One new MATLAB function, `exportCMGModelJSON(model, file)`, and one JS module
that consumes it. Both A and B then write decks from the same description,
which means **one format specification, two implementations** — and a natural
regression test: the JS writer's output must be byte-identical to the MATLAB
writer's for the same model.

Size is the obvious concern: Johansen's `props` arrays are 330 000 values
each. Mitigations: gzip (these compress hard — mostly zeros), omit inactive
entries via run-length encoding, or ship binary arrays alongside a JSON
header. Decide when B2 is built, not before.

---

## 5. Architecture

**A — local app**

```
browser (static HTML/JS)  ──REST+polling──▶  FastAPI  ──▶  MATLAB worker
                                               │             (runner.py)
                                               ▼
                                       web/runs/<run-id>/
```

The MATLAB execution layer already exists and is tested:
`tools/matlab-mrst-mcp/matlab_mrst_mcp/runner.py` runs MATLAB headlessly,
bootstraps MRST, parses errors, strips JVM shutdown noise, returns structured
JSON. The backend calls it rather than re-inventing subprocess handling.

**B — merslab app**

```
browser (static HTML/JS + deck-writer module)  ──▶  download
```

No server at all. Same UI shell as A, with the MATLAB-only stages disabled and
explained.

**Shared front end.** A and B run the same UI code, differing only in which
stages are enabled and whether a stage executes locally in JS or via the API.
Build the stage runner behind one interface from day one:

```js
runStage(n, input, opts) -> { output, stats, preview }
```

with a `MatlabStageRunner` (A) and a `BrowserStageRunner` (B). Getting this
seam right at the start is what makes B cheap later.

### State between stages

Each stage writes a checkpoint — a `.mat` in A, an in-memory/IndexedDB object
in B. Stateless stages (load previous, run, save new) rather than a long-lived
MATLAB session: slower per stage, but stages become resumable, individually
downloadable and separately cacheable, which is exactly the "look at it, then
continue" behaviour wanted. The `.mat` files are also the app's real currency —
a user can download `stage3_edfm.mat`, load it in their own MATLAB and carry
on by hand.

---

## 6. The stage contract

Every stage is one MATLAB entry point with the same shape, so the backend
treats them uniformly and new stages are cheap:

```matlab
function out = webStage<N>(inFile, outFile, optsJson)
```

These live in `web/matlab/` and are thin — parse options, call the existing
MRST and preprocessor functions, report. **No new numerical code.**

| Stage | Calls | Reports | Browser-only (B)? |
|---|---|---|---|
| 1 grid | `processGRDECL`/`cartGrid`/`tensorGrid`, `makeRock` | cells, cartDims, family, pore volume, ranges | yes (B1) |
| 2 fractures | build/validate `fracplanes` | plane count, areas, aperture/perm | yes |
| 3 EDFM grid | `EDFMGridCPG` \| `EDFMshalegrid` | fracture cells, time | **no** |
| 4 frac–matrix | `fracturematrixShaleNNC3D` | NNC count | **no** |
| 5 frac–frac | `fracturefractureShaleNNCs3D` | NNC count, pairs checked | **no** |
| 6 projection | `pMatFracNNCs3D` \| `pMatFracNNCsCPEDFM_CPG` | projected NNCs, seals | **no** |
| 7 export | `writeCMG` \| `writeEclipse*` | blocks, inactive, padding layers, connections | yes |

Stage 7 runs after **any** of 1, 3–6. In B, stages 3–6 are visibly disabled
with a link to download the local app.

---

## 7. Front end

Static HTML/CSS/JS, no build step, matching the clean scientific look of the
merslab tools. A stepper down the left, the active stage's form in the middle,
results and preview on the right. Completed stages collapse to a summary row
with their key numbers and a download button, so the whole run history stays
visible — the flash-predictor's "metrics panel + downloadable result" pattern,
repeated per stage.

**Preview.** In A, phase 1 renders server-side with MATLAB `plotGrid` to PNG.
In B there is no MATLAB, so preview must be `three.js` over the grid's
exterior faces — which means the browser mesh builder is needed for B anyway
and should be written once, used by both.

**Long stages.** EDFM gridding on Johansen took 317 s. Live status and elapsed
time, never a frozen button; the backend runs the job in a worker thread and
the page polls `GET /api/runs/<id>/status`.

**Validation first.** Cheap client-side checks — planes inside the grid
bounding box, aperture > 0, positive permeability, ≥3 non-collinear points per
plane — so obvious mistakes fail in milliseconds rather than after a
five-minute grid build.

---

## 8. Inputs

**Grid source** — upload `.GRDECL`; upload MRST `.mat` with `G` and `rock`;
built-in cases (Johansen NPD5, Norne sector, `simpleGrdecl` faulted box);
parametric Cartesian / variable Cartesian.

**Rock** — constant, per-direction anisotropic, or arrays from the GRDECL/file.

**Fractures** — upload a table of plane corner points plus aperture, porosity,
permeability; parametric multi-stage hydraulic fractures
(`createMultiStageHFs`); stochastic DFN where ADFNE is installed, with a clear
message where not; or none.

**Method** — EDFM / pEDFM / CPEDFM. The app should *detect* what the grid
carries via `cmgDetectFractureModel` and display it rather than trusting the
user's selection, exactly as `writeCMG` does today.

**Output** — MRST `.mat`; CMG for IMEX/GEM/STARS with `hostDeck` = `stub`
(Builder-ready, grid only), `none` (include only) or `runnable`; Eclipse
`.DATA`.

**Script upload** is an **A-only** feature. It is arbitrary code execution —
fine locally, since the user is running their own code on their own machine,
but it can never ship in B (no server to run it on, and MATLAB Compiler cannot
execute user `.m` at runtime even if there were). The UI must say so rather
than appearing to be a missing feature.

---

## 9. Phases

### Status

Phases 0–5 are built and tested. What was delivered differs from the plan in
three ways, all recorded here rather than quietly:

* **No FastAPI.** The backend is Python standard library only
  (`http.server`). The only Python on the target machine is a bare install
  with `mcp` in it, and the hosted build has no backend at all — asking
  someone to `pip install` a web framework to run the *local* convenience
  tool was the wrong trade. The API shape is unchanged.
* **Phases 4 and 5 came first**, in the order section 11 anticipated: the
  browser writers and the byte-identity harness were built before the
  MATLAB-backed stages, because the static app is the deliverable that gets
  hosted.
* **Ordering within the app is by engine, not by phase.** The page detects a
  backend and shows a mode picker; with none, it is exactly the Phase 5
  static app.

Acceptance evidence, all currently green:

| check | how |
|---|---|
| browser writers match MATLAB byte for byte | `selftest.html`, 12/12 files |
| the seven stages run in order | `webSelfTest`, 3 chains |
| EDFM, pEDFM and CPEDFM all export | `webSelfTest` chains A and C |
| matrix-only export off stage 1 | `webSelfTest` chain B |
| HTTP API, routing and guards | `selftest.py`, 21/21 |
| the page is wired to the API | `browsertest.py --matlab`, 15/15 |

Not built: **Phase 6 script upload.** Running arbitrary user MATLAB is a
different problem from running a fixed pipeline — it is an execution sandbox
question, not a preprocessing one — and the form covers what the scripts in
this repository actually build.

### The original plan

**Phase 0 — skeleton (A).** FastAPI app, run store, `webStage1`, a page that
builds a Cartesian grid and shows statistics. Proves browser → API → MATLAB →
JSON end to end.
*Accept:* a `cartGrid` from the form whose cell count and pore volume match a
hand-run MATLAB session.

**Phase 1 — matrix grid, all sources, export (A).** GRDECL upload, built-in
cases, parametric grids, PNG preview, stage 7 restricted to matrix-only.
*Accept:* Johansen NPD5 uploaded through the browser produces `.dat` + `.inc`
byte-identical to `JohansenGridpre.m`'s matrix export, and opens in Builder.

**Phase 2 — fractures and EDFM stages (A).** Stages 2–5, checkpointing,
per-stage download, progress reporting.
*Accept:* the Johansen 4-plane case reproduces 1141 fracture cells and 1155
NNCs through the browser, and its deck matches the script's.

**Phase 3 — projection and full export (A).** Stage 6, all three methods,
three simulators, Eclipse output.
*Accept:* EDFM/pEDFM/CPEDFM on one shared grid reproduce `test32`'s counts
(20/0, 58/28, 68/48).

**Phase 4 — the portable format.** `exportCMGModelJSON`, the JSON schema, and
the JS deck-writer module.
*Accept:* for the same model, the JS writer's `.inc` is byte-identical to
`writeCMGModel`'s.

**Phase 5 — merslab app (B).** Static build: B1 (GRDECL → deck) and B2 (JSON →
deck), three.js preview, MATLAB-only stages disabled with explanation.
*Accept:* deployable as static files; produces a Builder-openable Johansen
matrix deck with no server and no MATLAB.

**Phase 6 — script upload (A only), polish, packaging.**

Phases 0–3 deliver the tool you can use immediately; 4–5 deliver the thing
that can live on merslab. Phase 4 is the hinge — everything hosted depends on
that format being right.

---

## 10. Risks

**R1 — B's scope disappoints.** The hosted app cannot do EDFM gridding, which
is the interesting part. Mitigation: be explicit in the UI about what runs
where, and make B1 excellent — GRDECL → Builder-ready deck with no install is
genuinely useful on its own.

**R2 — Two writers drift.** A MATLAB writer and a JS writer for the same
format will diverge unless pinned. Mitigation: the byte-identical test in
Phase 4 runs in CI over a fixture set.

**R3 — JSON payload size.** Johansen `props` arrays are 330 000 values each.
Mitigation: gzip, run-length encoding of inactive slots, or binary sidecars.

**R4 — Long stages and browser timeouts.** 317 s gridding exceeds default
proxy and fetch timeouts. Mitigation: job + polling from the start, never a
synchronous request.

**R5 — Duplicate MRST functions on the path.** *Decided: pin `hfm_v2`, loaded
**after** `startupCPEDFMCPG` (which internally runs `mrstModule add hfm` and
would otherwise shadow it).* The two modules differ in seven files, including
`convexpolygonarea` (one output vs two — loading `hfm` makes `EDFMGridCPG`
fail outright) and `calcdavg`, which sets the fracture-matrix connectivity
index. Measured on the Cartesian case both give identical transmissibilities
(`sum(T) = 1.328887254593e-15`); divergence is expected only on non-box
corner-point cells and has not been measured. Six copies of `EDFMGridCPG`
exist. The worker must set an explicit ordered path per run and log `which()`
for `EDFMGridCPG`, `convexpolygonarea`, `calcdavg` and
`fracturematrixShaleNNC3D`, so a wrong pairing shows in the log rather than as
a strange number.

**R6 — Upstream fragility.** `EDFMGridCPG` fails on some geometries (a
concave-intersection branch). Surface the MATLAB error verbatim, never
"something went wrong".

**R7 — Checkpoint size and concurrency (A).** Post-EDFM `G` is large; two
simultaneous runs mean two MATLAB licence seats. Mitigate with `-v7.3`, a disk
budget, eviction, and a queue depth of one.

---

## 11. Open questions

1. ~~Local-first or hosted?~~ **Answered:** local first, merslab later.
2. ~~Should the app run the exported deck in IMEX?~~ **Answered: no.** The app
   produces files. Validation stays in the MATLAB test suite, and the app
   carries no CMG dependency.
3. Is B1 (GRDECL → deck, no MATLAB) worth building before the full local app?
   It is the smallest useful thing, ships to merslab immediately, and needs no
   backend — a case for reordering Phase 5 ahead of Phase 2.
4. How much fracture authoring belongs in the browser? Uploading a table is
   easy; drawing planes on a 3D grid is a large piece of work.
5. Does the merslab page need to match the site's existing template
   (navigation, branding), or is a standalone page under `/apps/` sufficient —
   as the Flash Calculation Predictor appears to be?
