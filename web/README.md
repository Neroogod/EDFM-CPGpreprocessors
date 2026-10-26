# EDFM Grid Builder — web app

A front end for the MRST → CMG/Eclipse preprocessor. See
`web_implementation.md` for the design and the phase plan.

It runs in two modes, and the mode picker only appears when the second one is
actually available:

| | **Browser** | **Local MATLAB** |
|---|---|---|
| needs | nothing | MATLAB + MRST + `python web/server/serve.py` |
| grids | Cartesian, tensor, corner-point | those plus unstructured |
| fractures | — | EDFM, pEDFM, CPEDFM |
| deployable to a static host | yes | no |

The browser mode is the one that gets hosted, so it must keep working with no
backend at all: the health probe simply fails and the page stays as it was.

## Browser mode

**Matrix grid → CMG or Eclipse deck, entirely in the tab.** No MATLAB, no
server, no install.

Two ways to supply the grid:

* **Upload** a corner-point `.GRDECL` → exports `*GRID *CORNER`
* **Describe** one in the form — block counts plus either a total extent
  (uniform, → `*GRID *CART`) or per-direction block sizes (→ `*GRID *VARI`)

Either grid can be written in either format:

* **CMG** — `case.dat` (opens in Builder) plus the `case.inc` it includes
* **Eclipse** — `case.DATA` (RUNSPEC + GRID) plus `case.grdecl`

Both come with a manifest recording what was exported. Neither carries a fluid
model: the CMG deck stops before `*MODEL`, the Eclipse deck at the end of
`GRID`.

```
cd web/static
python -m http.server 8765
```

A plain `file://` open will not work: the page uses ES modules, which browsers
refuse to load over `file://`. Any static server will do.

## Local MATLAB mode

```
python web/server/serve.py --open
```

No dependencies beyond the Python standard library — deliberately, since the
hosted version has no backend at all and making people install a web framework
to get the local one would be a poor trade. The server finds MATLAB from
`MATLAB_EXE`, then `PATH`, then the usual install directories, and says so on
startup.

This unlocks the full pipeline. Each stage writes a checkpoint, so the deck can
be exported at any point — the matrix grid alone, or after any of the fracture
stages:

| stage | does | calls |
|---|---|---|
| 1 | matrix grid and rock | `cartGrid` / `tensorGrid` / `processGRDECL` |
| 2 | fracture planes | `createMultiStageHFs`, or corner points |
| 3 | cut fracture cells in | `EDFMshalegrid` / `EDFMGridCPG` / `EDFMGridU` |
| 4 | fracture–matrix NNCs | `fracturematrixShaleNNC3D` |
| 5 | fracture–fracture NNCs | `fracturefractureShaleNNCs3D` |
| 6 | projection | `pMatFracNNCs3D` / `pMatFracNNCsCPEDFM_CPG` / `_U` |
| 7 | write the deck | `writeCMG` / `writeEclipse` |

Stage 3 picks the gridder from the detected grid family. This is not a
convenience: `EDFMshalegrid` and `EDFMGridCPG` both assume eight corners per
cell and fail inside `checkplaneAABBintersect` on a PEBI grid.

### Environment

| variable | default |
|---|---|
| `MATLAB_EXE` | auto-detected |
| `MRST_ROOT` | the repository root |
| `WEB_WORK_DIR` | `web/work` |
| `WEB_TIMEOUT_S` | `1800` per stage |
| `WEB_MAX_UPLOAD` | 200 MB |

Session directories older than 24 h are deleted at startup.

### API

| | |
|---|---|
| `GET /api/health` | MATLAB present? paths? stage list? |
| `POST /api/session` | → `{sessionId}` |
| `GET /api/session/<sid>` | checkpoint and download state |
| `DELETE /api/session/<sid>` | delete the working directory |
| `POST /api/session/<sid>/stage/<n>` | body = options → `{jobId}` |
| `GET /api/job/<jid>?wait=25` | long-poll until the stage ends |
| `GET /api/session/<sid>/file/<name>` | download |
| `POST /api/session/<sid>/upload/<name>` | store a GRDECL for stage 1 |

`?wait=` is what keeps a five-minute corner-point run from being a five-minute
poll loop: the request hangs on the server until the stage finishes.

## Tests

```
python web/server/selftest.py            # HTTP API, 21 checks
python web/server/selftest.py --quick    # routing and validation only, no MATLAB
python web/server/browsertest.py         # real headless Chrome, browser mode
python web/server/browsertest.py --matlab   # ... and the whole MATLAB pipeline
```

```matlab
webSelfTest        % the seven stages through webStage, as the backend calls them
```

`webSelfTest` runs three chains off one grid: EDFM → pEDFM → CMG, the same
grid CPEDFM → Eclipse, and a matrix-only Eclipse export straight from stage 1.
Nothing in it calls a stage function directly, so a dispatcher that forgets a
case fails there.

`browsertest.py` uses Chrome's `--dump-dom`, and app.js sets
`data-selftest="finished"` on `<body>` when a scripted run completes, so the
dumped DOM is the finished state rather than a snapshot of a page still
working. A worked example can also be linked directly:
`index.html?grid=param&run=1`, plus `&spacing=variable`, `&engine=matlab`,
`&frac=1`, `&export=1`.

### The writer regression test

The browser writer and `writeCMGModel.m` must produce **byte-identical** output
for the same grid. Two implementations of one format drift silently otherwise,
and a deck that is subtly wrong is worse than one that fails.

Open <http://127.0.0.1:8765/selftest.html>. Current status: **PASS**, all
twelve files byte-identical.

| case | grid | CMG `.inc` / `.dat` | Eclipse `.grdecl` / `.DATA` |
|---|---|---|---|
| `ref` | corner-point, 8×6×4 faulted | 33,225 / 1,274 B | 32,779 / 903 B |
| `ref_cart` | uniform Cartesian, 10×7×5 | 18,899 / 1,288 B | 22,016 / 919 B |
| `ref_vari` | variable Cartesian, 6×3×3 | 3,879 / 1,276 B | 3,890 / 907 B |

The parametric cases take their grid description from a `.json` written by the
same MATLAB script, so both sides start from identical inputs rather than from
separately typed numbers.

Regenerate the reference pair after any change to the MATLAB writers:

```matlab
run('web/matlab/makeWebReferenceCase.m')
```

That case is built with every block active on purpose. MRST's `processGRDECL`
can deactivate more cells than `ACTNUM` does — it discards disconnected
components and degenerate cells — and a GRDECL parser cannot know about that.
Byte-identity is only meaningful where the two sides agree on the active set.

## Layout

```
web/
  static/                 the app - this is the deployable artifact
    index.html
    selftest.html         writer regression test
    css/style.css
    js/
      app.js              UI wiring, both engines
      backend.js          the /api client; absent backend is not an error
      grdecl.js           GRDECL reader
      cmgmodel.js         port of buildCMGModel.m (matrix grids)
      cmgwriter.js        port of writeCMGModel.m and writeCMG.m's stub deck
      eclwriter.js        port of writeEclipseModel.m and writeEclipse.m
      fmt.js              MATLAB-compatible printf formatting
      preview.js          top-layer depth map
    testdata/             reference pair, generated by MATLAB
  server/                 local backend, standard library only
    serve.py              entry point
    app.py                routing, static files, path-traversal guards
    store.py              sessions and background jobs
    runner.py             one MATLAB -batch call per stage
    config.py             MATLAB discovery and settings
    selftest.py           HTTP API test
    browsertest.py        headless Chrome test
  matlab/
    webSetupPaths.m       path pinning - hfm_v2 after startupCPEDFMCPG
    webStage.m            stage dispatcher; always writes a result JSON
    webStage1Grid.m       ... stage 1
    webStage2Fractures.m  ... stage 2
    webStage3EDFM.m       ... stage 3
    webStage456.m         ... stages 4, 5 and 6
    webStage7Export.m     ... stage 7
    webSelfTest.m         the seven stages end to end
    makeWebReferenceCase.m  regenerates testdata/
  web_implementation.md   design and phases
```

## Things worth knowing

* **Tolerance is a length, in metres.** The gridders default to 1e-5, which
  suits a shale model a few hundred metres across and fails on a field-scale
  grid with the unhelpful message `Tolerance too tight`. Stage 3 scales the
  starting value with the grid extent and loosens it on retry, reporting what
  it did.
* **A fracture lying exactly on a cell face** produces only boundary
  connections, so pEDFM and CPEDFM have nothing to project and the export is
  EDFM in all but name. Stage 6 says so rather than letting it pass silently.
* **The depth datum is applied at export**, not to the grid, so fracture
  coordinates are in the grid's own system — z from 0, not from 2000 m.

## Known limits

* `INCLUDE` in a GRDECL is not followed; the file must be self-contained.
* In browser mode active blocks come from `ACTNUM` alone (see above).
* Parametric grids have every block active and uniform rock; per-block
  properties come only from an uploaded file.
* Uploading an MRST *script* is not supported; running arbitrary user MATLAB
  is a different problem from running a fixed pipeline, and the form covers
  the cases the scripts in this repository actually build.
* One stage at a time per session. Starting a second returns 409 rather than
  queueing, because the stages share checkpoint files.
