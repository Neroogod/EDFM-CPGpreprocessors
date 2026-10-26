// UI wiring. Everything runs in this tab; nothing is uploaded.

import { parseGRDECL } from "./grdecl.js";
import { buildMatrixModel, buildParametricModel } from "./cmgmodel.js";
import { writeInclude, writeStubDeck, writeManifest } from "./cmgwriter.js";
import { writeEclipseInclude, writeEclipseStub } from "./eclwriter.js";
import { drawTopLayer } from "./preview.js";
import { detect, Session } from "./backend.js";

const $ = (id) => document.getElementById(id);
const state = {
  source: "file",
  engine: "browser",      // "browser" | "matlab"
  grdecl: null,
  model: null,            // browser engine only
  session: null,          // matlab engine only
  fractured: false,
  urls: [],
};

// ---------------------------------------------------------------- source

for (const tab of document.querySelectorAll('.tab[data-source]')) {
  tab.addEventListener("click", () => {
    state.source = tab.dataset.source;
    for (const t of document.querySelectorAll('.tab[data-source]')) {
      t.classList.toggle("active", t === tab);
    }
    for (const p of document.querySelectorAll('.pane[data-pane]')) {
      p.hidden = p.dataset.pane !== state.source;
    }
    $("btn-load").disabled = state.source === "file" && !$("file").files.length;
  });
}

// --------------------------------------------------------------- engine

for (const tab of document.querySelectorAll(".tab[data-engine]")) {
  tab.addEventListener("click", () => {
    setEngine(tab.dataset.engine);
  });
}

function setEngine(engine) {
  state.engine = engine;
  for (const t of document.querySelectorAll(".tab[data-engine]")) {
    t.classList.toggle("active", t.dataset.engine === engine);
  }
  $("engine-note-browser").hidden = engine !== "browser";
  $("engine-note-matlab").hidden = engine === "browser";
  $("fld-from").hidden = engine === "browser";
  $("rail-note-browser").hidden = engine !== "browser";
  $("rail-note-matlab").hidden = engine === "browser";
  // The privacy claim is only true in browser mode; in MATLAB mode the grid
  // does leave the tab, even if only as far as localhost.
  $("foot-privacy").textContent = engine === "browser"
    ? "Generated decks contain only data derived from your file. Nothing is " +
      "uploaded — all processing happens in this browser."
    : "Your grid is sent to the local server on this machine and written to " +
      "its work directory. It does not leave the computer.";
  for (const n of [2, 3, 4, 5, 6]) markStage(n, engine === "matlab" ? "pending" : "locked");
  // Switching engines invalidates whatever the other one built.
  state.model = null;
  state.grdecl = null;
  state.fractured = false;
  $("card-summary").hidden = true;
  $("card-frac").hidden = true;
  $("card-export").hidden = true;
  $("grid-status").hidden = true;
  markStage(1, "active");
  markStage(7, "pending");
}

detect().then((info) => {
  if (!info.available || !info.matlab) return;
  $("card-engine").hidden = false;
  $("engine-workdir").textContent = info.workDir || "the server's work directory";
  $("rail-matlab-path").textContent = info.matlabExe ? `(${info.matlabExe})` : "";
  // Default to the engine that can do everything. Landing in browser mode
  // when MATLAB is right there means the first thing the page says is that
  // fractures are unavailable, which is only true of the mode it picked.
  if (new URLSearchParams(location.search).get("engine") !== "browser") {
    setEngine("matlab");
  }
});

$("spacingMode").addEventListener("change", () => {
  const variable = $("spacingMode").value === "variable";
  $("fld-uniform").hidden = variable;
  $("fld-variable").hidden = !variable;
  // Seed the lists from the uniform grid so the counts always start valid.
  // Leaving stale lists behind is a guaranteed "got 6 values, need 20" error
  // the first time anyone switches modes.
  if (variable) seedSpacingLists();
});

for (const id of ["nx", "ny", "nz", "Lx", "Ly", "Lz"]) {
  $(id).addEventListener("change", () => {
    if ($("spacingMode").value === "variable") seedSpacingLists();
  });
}

function seedSpacingLists() {
  const n = [Number($("nx").value), Number($("ny").value), Number($("nz").value)];
  const L = [Number($("Lx").value), Number($("Ly").value), Number($("Lz").value)];
  const ids = ["dxs", "dys", "dzs"];
  for (let d = 0; d < 3; d++) {
    const count = Math.max(1, Math.round(n[d]));
    const have = $(ids[d]).value.split(/[\s,]+/).filter((s) => s.length);
    if (have.length === count) continue;          // user's own list is fine
    const size = (L[d] > 0 ? L[d] : count) / count;
    $(ids[d]).value = new Array(count).fill(trim(size)).join(" ");
  }
}

const trim = (v) => String(Number(v.toFixed(6)));

$("file").addEventListener("change", () => {
  $("btn-load").disabled = !$("file").files.length;
  const f = $("file").files[0];
  if (f) {
    const base = f.name.replace(/\.[^.]+$/, "").replace(/[^\w.-]/g, "_");
    $("caseName").value = base || "grid_export";
    $("title").value = (base || "grid_export").toUpperCase().slice(0, 40);
  }
});

// The parametric form needs no file, so the button is live from the start.
$("btn-load").disabled = false;
document.querySelector('.tab[data-source="file"]').click();

// ?grid=param[&spacing=variable][&run=1][&engine=matlab][&frac=1][&export=1]
// preselects the form and optionally drives it, so a worked example can be
// linked to directly. It is also how the DOM wiring gets exercised in an
// automated browser test - see web/server/browsertest.py.
{
  const q = new URLSearchParams(location.search);
  const wantMatlab = q.get("engine") === "matlab";

  const drive = async () => {
    if (wantMatlab) {
      // The engine tab only exists once detect() has answered.
      await waitFor(() => !$("card-engine").hidden, 6000);
      document.querySelector('.tab[data-engine="matlab"]').click();
      document.querySelector('.tab[data-source="param"]').click();
    }
    if (!q.get("run")) return;
    $("btn-load").click();
    if (!q.get("frac") && !q.get("export")) return;

    await waitFor(() => !$("card-export").hidden, 600000);
    if (q.get("frac")) {
      $("btn-frac").click();
      await waitFor(() => $("frac-status").classList.contains("ok") ||
                          $("frac-status").classList.contains("err"), 600000);
    }
    if (q.get("export")) {
      $("btn-export").click();
      await waitFor(() => $("export-status").classList.contains("ok") ||
                          $("export-status").classList.contains("err"), 600000);
    }
    document.body.dataset.selftest = "finished";
  };

  if (q.get("grid") === "param") {
    if (!wantMatlab) document.querySelector('.tab[data-source="param"]').click();
    if (q.get("spacing") === "variable") {
      $("spacingMode").value = "variable";
      $("spacingMode").dispatchEvent(new Event("change"));
    }
    // As a microtask, not inline: the click handlers are registered further
    // down this module, so clicking now would dispatch into nothing.
    queueMicrotask(() => drive().catch((err) => {
      document.body.dataset.selftest = `failed: ${err.message}`;
    }));
  }
}

function waitFor(pred, limitMs) {
  return new Promise((resolve, reject) => {
    const t0 = Date.now();
    const tick = () => {
      if (pred()) return resolve();
      if (Date.now() - t0 > limitMs) return reject(new Error("waitFor timed out"));
      setTimeout(tick, 120);
    };
    tick();
  });
}

// ------------------------------------------------------------ read grid

$("btn-load").addEventListener("click", async () => {
  $("btn-load").disabled = true;
  try {
    if (state.engine === "matlab") {
      await loadGridInMatlab();
      return;
    }
    const t0 = performance.now();
    let model, grdecl = null;

    if (state.source === "file") {
      const file = $("file").files[0];
      if (!file) throw new Error("Choose a .GRDECL file first.");
      setStatus("grid-status", "busy",
        `Reading ${file.name} (${fmtBytes(file.size)})…`);
      const text = await file.text();
      // Yield once so the busy state paints before a multi-megabyte parse.
      await nextFrame();
      grdecl = parseGRDECL(text);
      model = buildMatrixModel(grdecl, {
        permUnit: $("permUnit").value,
        perm: Number($("permDefault").value),
        poro: Number($("poroDefault").value),
        simulator: $("simulator").value,
        cmgVersion: $("cmgVersion").value,
      });
    } else {
      setStatus("grid-status", "busy", "Building grid…");
      await nextFrame();
      model = buildParametricModel(readParametricSpec(), {
        simulator: $("simulator").value,
        cmgVersion: $("cmgVersion").value,
      });
    }

    const ms = performance.now() - t0;
    state.grdecl = grdecl;
    state.model = model;

    renderSummary(grdecl, model, ms);
    if (grdecl) {
      $("preview").parentElement.hidden = false;
      drawTopLayer($("preview"), grdecl);
    } else {
      // A parametric grid is a box at a single depth: a depth map would be
      // one flat colour and tell the user nothing.
      $("preview").parentElement.hidden = true;
    }

    markStage(1, "done");
    markStage(7, "active");
    $("card-summary").hidden = false;
    $("card-export").hidden = false;
    setStatus("grid-status", "ok",
      `${model.meta.counts.cells.toLocaleString()} active blocks · ` +
      `exports as *GRID ${model.grid.kind.toUpperCase()} · ${ms.toFixed(0)} ms.`);
    $("card-export").scrollIntoView({ behavior: "smooth", block: "nearest" });
  } catch (err) {
    setStatus("grid-status", "err", `Could not build the grid.\n${err.message}`);
    $("card-summary").hidden = true;
    $("card-export").hidden = true;
  } finally {
    $("btn-load").disabled = false;
  }
});

// ------------------------------------------------------- MATLAB engine

async function loadGridInMatlab() {
  setStatus("grid-status", "busy", "Starting a MATLAB session…");
  state.session = await Session.create();
  state.fractured = false;

  let opts;
  if (state.source === "file") {
    const file = $("file").files[0];
    if (!file) throw new Error("Choose a .GRDECL file first.");
    setStatus("grid-status", "busy",
      `Sending ${file.name} (${fmtBytes(file.size)}) to the local server…`);
    const stored = await state.session.upload("input.grdecl", await file.text());
    opts = {
      source: "grdecl",
      file: stored.path,
      perm: Number($("permDefault").value),
      poro: Number($("poroDefault").value),
    };
  } else {
    const spec = readParametricSpec();
    opts = {
      source: spec.dx ? "tensor" : "cartesian",
      dims: spec.dims,
      physdim: spec.physdim,
      dx: spec.dx, dy: spec.dy, dz: spec.dz,
      perm: spec.permi,
      poro: spec.poro,
      datum: spec.datum,
    };
  }

  const result = await runStage(1, opts, "grid-status", "Building the grid in MRST");
  renderServerSummary(result.stats, result.warnings);
  $("preview").parentElement.hidden = true;   // no client-side geometry to draw

  markStage(1, "done");
  markStage(2, "active");
  markStage(7, "active");
  $("card-summary").hidden = false;
  $("card-frac").hidden = false;
  $("card-export").hidden = false;
  setStatus("grid-status", "ok",
    `${Number(result.stats.cells).toLocaleString()} cells · ` +
    `${result.stats.family} grid · ${result.elapsedSeconds}s. ` +
    `Export now, or add fractures below.`);
}

/** Run one stage with a live clock in the given status box. */
async function runStage(stage, opts, statusId, label) {
  const started = Date.now();
  return state.session.runStage(stage, opts, () => {
    const s = ((Date.now() - started) / 1000).toFixed(0);
    setStatus(statusId, "busy", `${label}… ${s}s`);
  });
}

$("btn-frac-fill").addEventListener("click", () => {
  // Two orthogonal planes through the middle of the grid, inset from the
  // faces. Without a bounding box from the server this is only a starting
  // point, so it is a suggestion the user edits, not a computed answer.
  const s = state.gridBounds;
  if (!s) {
    setStatus("frac-status", "err",
      "Read a grid first — the suggestion is based on its extent.");
    return;
  }
  const [x0, x1, y0, y1, z0, z1] = s;
  const f = (lo, hi, t) => (lo + (hi - lo) * t).toFixed(2);
  const za = f(z0, z1, 0.1), zb = f(z0, z1, 0.9);
  const xm = f(x0, x1, 0.47), ym = f(y0, y1, 0.47);
  $("fracPlanes").value =
    `${xm} ${f(y0, y1, 0.1)} ${za}; ${xm} ${f(y1, y0, 0.1)} ${za}; ` +
    `${xm} ${f(y1, y0, 0.1)} ${zb}; ${xm} ${f(y0, y1, 0.1)} ${zb}\n` +
    `${f(x0, x1, 0.1)} ${ym} ${za}; ${f(x1, x0, 0.1)} ${ym} ${za}; ` +
    `${f(x1, x0, 0.1)} ${ym} ${zb}; ${f(x0, x1, 0.1)} ${ym} ${zb}`;
  setStatus("frac-status", "ok", "Two crossing planes suggested — edit as needed.");
});

for (const tab of document.querySelectorAll(".tab[data-frac]")) {
  tab.addEventListener("click", () => {
    for (const t of document.querySelectorAll(".tab[data-frac]")) {
      t.classList.toggle("active", t === tab);
    }
    for (const p of document.querySelectorAll("[data-fracpane]")) {
      p.hidden = p.dataset.fracpane !== tab.dataset.frac;
    }
  });
}

$("btn-frac").addEventListener("click", async () => {
  if (!state.session) return;
  $("btn-frac").disabled = true;
  $("btn-export").disabled = true;
  const method = $("fracMethod").value;

  // No tolerance is sent: the gridder's tolerance is an absolute length, so
  // the right value depends on the model's size. MATLAB derives it from the
  // grid rather than taking a number the user has no basis to choose.
  const steps = [
    [2, "Defining fracture planes", () => readFractureOpts()],
    [3, "Cutting fracture cells into the grid", () => ({})],
    [4, "Fracture–matrix connections", () => ({})],
    [5, "Fracture–fracture connections", () => ({})],
    [6, "Projection", () => ({ method })],
  ];

  const list = $("frac-progress");
  list.innerHTML = "";
  list.hidden = false;
  const rows = new Map();
  for (const [n, label] of steps) {
    const li = document.createElement("li");
    li.className = "waiting";
    li.innerHTML = `<span class="pn">${n}</span><span class="pl">${label}</span>` +
                   `<span class="pv"></span>`;
    list.appendChild(li);
    rows.set(n, li);
  }

  try {
    for (const [n, label, makeOpts] of steps) {
      const li = rows.get(n);
      li.className = "running";
      markStage(n, "active");
      const result = await runStage(n, makeOpts(), "frac-status", label);
      li.className = "ok";
      li.querySelector(".pv").textContent = describeStats(result.stats);
      markStage(n, "done");
      if (result.warnings && result.warnings.length) {
        const note = document.createElement("small");
        note.className = "pw";
        note.textContent = result.warnings.join(" ");
        li.appendChild(note);
      }
    }
    state.fractured = true;
    setStatus("frac-status", "ok",
      `${method} model complete. Export the deck below — connections are ` +
      `written as *SCONNECT (CMG) or NNC (Eclipse).`);
    markStage(7, "active");
  } catch (err) {
    for (const li of rows.values()) {
      if (li.className === "running") li.className = "err";
    }
    setStatus("frac-status", "err",
      [`Stage failed.`, err.message, ...(err.detail || [])].join("\n"));
  } finally {
    $("btn-frac").disabled = false;
    $("btn-export").disabled = false;
  }
});

function readFractureOpts() {
  const common = {
    aperture: Number($("fracAperture").value) / 1000,   // mm -> m
    poro: Number($("fracPoro").value),
    perm: Number($("fracPerm").value),                  // mD, converted in MATLAB
  };
  const mode = document.querySelector(".tab[data-frac].active").dataset.frac;
  if (mode === "multistage") {
    return {
      ...common,
      mode: "multistage",
      numStages: Number($("msStages").value),
      fracsPerStage: Number($("msPerStage").value),
      spacing: Number($("msSpacing").value),
      clusterSpacing: Number($("msCluster").value),
      halfLength: Number($("msHalf").value),
      height: Number($("msHeight").value),
      centroid1: [Number($("msCx").value), Number($("msCy").value),
                  Number($("msCz").value)],
    };
  }

  const planes = [];
  const lines = $("fracPlanes").value.split(/\n+/).map((s) => s.trim())
    .filter((s) => s.length && !s.startsWith("#"));
  if (!lines.length) throw new Error("No fracture planes were given.");
  lines.forEach((line, i) => {
    const corners = line.split(";").map((s) => s.trim()).filter((s) => s.length);
    if (corners.length < 3) {
      throw new Error(`Plane ${i + 1}: needs at least three corners, got ${corners.length}.`);
    }
    const points = corners.map((c, k) => {
      const v = c.split(/[\s,]+/).map(Number);
      if (v.length !== 3 || v.some((x) => !Number.isFinite(x))) {
        throw new Error(`Plane ${i + 1}, corner ${k + 1}: expected three numbers, got "${c}".`);
      }
      return v;
    });
    planes.push({ points, ...common });
  });
  return { mode: "explicit", planes, ...common };
}

function describeStats(stats) {
  if (!stats) return "";
  const keys = ["fractureCells", "added", "total", "planes", "projectionRecords"];
  return keys.filter((k) => stats[k] !== undefined)
             .map((k) => `${k} ${Number(stats[k]).toLocaleString()}`)
             .join(" · ");
}

function renderServerSummary(stats, warnings) {
  const dims = Array.isArray(stats.cartDims) ? stats.cartDims : [];
  state.gridBounds = null;
  const metrics = [
    ["Dimensions", dims.length === 3 ? dims.join(" × ") : "unstructured",
     `${Number(stats.cells).toLocaleString()} cells`],
    ["Depth range",
     `${Number(stats.depthTop).toFixed(1)} – ${Number(stats.depthBottom).toFixed(1)}`, "m"],
    ["Grid family", stats.family, stats.source],
    ["Pore volume", Number(stats.poreVolume).toPrecision(5), "m³"],
    ["Permeability",
     `${Number(stats.permMinMD).toPrecision(4)} – ${Number(stats.permMaxMD).toPrecision(4)}`,
     "mD"],
    ["Porosity",
     `${Number(stats.poroMin).toFixed(3)} – ${Number(stats.poroMax).toFixed(3)}`, ""],
  ];
  $("summary").innerHTML = metrics.map(([k, v, sub]) =>
    `<div><dt>${k}</dt><dd>${escapeHtml(String(v))} <small>${escapeHtml(String(sub ?? ""))}</small></dd></div>`).join("");

  const box = $("warnings");
  if (warnings && warnings.length) {
    box.innerHTML = `<strong>Notes</strong><ul>${warnings.map((w) =>
      `<li>${escapeHtml(w)}</li>`).join("")}</ul>`;
    box.hidden = false;
  } else {
    box.hidden = true;
  }

  // The suggestion button needs an extent. Parametric grids give one directly;
  // an uploaded GRDECL does not, and guessing one would be worse than saying so.
  if (Array.isArray(stats.cartDims) && state.source === "param") {
    const spec = readParametricSpec();
    const L = spec.physdim || [sum(spec.dx), sum(spec.dy), sum(spec.dz)];
    state.gridBounds = [0, L[0], 0, L[1], 0, L[2]];
  }
}

const sum = (a) => (a || []).reduce((x, y) => x + y, 0);

function readParametricSpec() {
  const dims = [Number($("nx").value), Number($("ny").value), Number($("nz").value)];
  if (dims.some((n) => !Number.isFinite(n) || n < 1)) {
    throw new Error("Block counts must be positive whole numbers.");
  }
  const spec = {
    dims,
    datum: Number($("datum").value) || 0,
    permi: Number($("permParam").value),
    poro: Number($("poroParam").value),
  };
  if ($("spacingMode").value === "variable") {
    spec.dx = parseList($("dxs").value, dims[0], "dx");
    spec.dy = parseList($("dys").value, dims[1], "dy");
    spec.dz = parseList($("dzs").value, dims[2], "dz");
  } else {
    const L = [Number($("Lx").value), Number($("Ly").value), Number($("Lz").value)];
    if (L.some((v) => !Number.isFinite(v) || v <= 0)) {
      throw new Error("Reservoir extent must be positive in every direction.");
    }
    spec.physdim = L;
  }
  return spec;
}

function parseList(text, want, name) {
  const v = text.split(/[\s,]+/).filter((s) => s.length).map(Number);
  if (v.some((x) => !Number.isFinite(x) || x <= 0)) {
    throw new Error(`${name}: every block size must be a positive number.`);
  }
  if (v.length !== want) {
    throw new Error(`${name}: got ${v.length} values but the grid has ${want} blocks in that direction.`);
  }
  return v;
}

$("format").addEventListener("change", () => {
  const ecl = $("format").value === "eclipse";
  $("fld-simulator").hidden = ecl;
  $("fld-version").hidden = ecl;
  $("note-cmg").hidden = ecl;
  $("note-ecl").hidden = !ecl;
  $("lbl-title").textContent = ecl ? "TITLE" : "TITLE1";
});

$("btn-export").addEventListener("click", async () => {
  if (state.engine === "matlab") {
    await exportFromMatlab();
    return;
  }
  if (!state.model) return;
  setStatus("export-status", "busy", "Writing deck…");

  try {
    // Re-read the options that affect the banner, so changing the simulator
    // or format does not require re-reading the grid.
    const m = state.model;
    m.meta.simulator = $("simulator").value;
    m.meta.cmgVersion = $("cmgVersion").value || "2024.30";

    const eclipse = $("format").value === "eclipse";
    const base = ($("caseName").value || "grid_export").replace(/[^\w.-]/g, "_");
    const title = $("title").value || base.toUpperCase();
    const incName = eclipse ? `${base}.grdecl` : `${base}.inc`;
    const datName = eclipse ? `${base}.DATA` : `${base}.dat`;

    const t0 = performance.now();
    const inc = eclipse ? writeEclipseInclude(m) : writeInclude(m);
    const dat = eclipse ? writeEclipseStub(m, incName, title)
                        : writeStubDeck(m, incName, title);
    const man = writeManifest(m, incName, datName);
    const ms = performance.now() - t0;

    showDownloads([
      [datName, dat, eclipse ? "Eclipse deck — RUNSPEC and GRID"
                             : "CMG dataset — open this in Builder"],
      [incName, inc, "grid description (arrays)"],
      [`${base}_manifest.json`, man, "what was exported"],
    ]);

    setStatus("export-status", "ok",
      `Wrote ${fmtBytes(dat.length + inc.length)} in ${ms.toFixed(0)} ms. ` +
      (eclipse
        ? `Add PROPS onwards to ${datName} — keep both files together.`
        : `Open ${datName} in Builder — keep both files in the same folder.`));
    markStage(7, "done");
  } catch (err) {
    setStatus("export-status", "err", `Export failed.\n${err.message}`);
  }
});

async function exportFromMatlab() {
  if (!state.session) return;
  $("btn-export").disabled = true;
  try {
    const eclipse = $("format").value === "eclipse";
    const base = ($("caseName").value || "grid_export").replace(/[^\w.-]/g, "_");
    const opts = {
      format: eclipse ? "eclipse" : "cmg",
      caseName: base,
      title: $("title").value || base.toUpperCase(),
      simulator: $("simulator").value,
      cmgVersion: $("cmgVersion").value || "2024.30",
      hostDeck: "stub",
    };
    const from = Number($("fromStage").value);
    if (from > 0) opts.fromStage = from;
    const result = await runStage(7, opts, "export-status", "Writing the deck");

    // The stage reports absolute paths; turn them into download links.
    const names = Object.values(result.files || {})
      .filter((p) => typeof p === "string")
      .map((p) => p.split(/[\\/]/).pop());
    showServerDownloads(names);

    const s = result.stats || {};
    setStatus("export-status", "ok",
      `${s.fractureModel} model, ${Number(s.blocks).toLocaleString()} blocks, ` +
      `${Number(s.connections || 0).toLocaleString()} connections. ` +
      (eclipse ? "Add PROPS onwards yourself — nothing is invented."
               : "Open the .dat in Builder — keep both files together."));
    markStage(7, "done");
  } catch (err) {
    setStatus("export-status", "err",
      ["Export failed.", err.message, ...(err.detail || [])].join("\n"));
  } finally {
    $("btn-export").disabled = false;
  }
}

function showServerDownloads(names) {
  const ul = $("downloads");
  ul.innerHTML = "";
  for (const name of names) {
    const li = document.createElement("li");
    li.innerHTML = `<span class="name">${escapeHtml(name)}</span>` +
                   `<span class="size">on the local server</span>`;
    const a = document.createElement("a");
    a.href = state.session.fileUrl(name);
    a.download = name;
    a.textContent = "Download";
    li.appendChild(a);
    ul.appendChild(li);
  }
  ul.hidden = false;
}

//------------------------------------------------------------------ render

function renderSummary(g, model, ms) {
  const [nx, ny, nz] = model.grid.dims;
  const c = model.meta.counts;
  const [zt, zb] = [model.grid.ztop, model.grid.zbot];

  const permOn = !!(g && g.PERMX), poroOn = !!(g && g.PORO);
  const metrics = [
    ["Dimensions", `${nx} × ${ny} × ${nz}`, `${c.blocks.toLocaleString()} blocks`],
    ["Active blocks", c.cells.toLocaleString(), `${c.inactive.toLocaleString()} inactive`],
    ["Depth range", `${zt.toFixed(1)} – ${zb.toFixed(1)}`, "m"],
    ["Grid type", `*${model.grid.kind.toUpperCase()}`, model.meta.gridFamily],
    ["Permeability",
      g ? (permOn ? "from file" : "fallback") : "uniform",
      g ? (permOn ? "PERMX/Y/Z" : `${$("permDefault").value} mD`) : `${$("permParam").value} mD`],
    ["Porosity",
      g ? (poroOn ? "from file" : "fallback") : "uniform",
      g ? (poroOn ? "PORO" : $("poroDefault").value) : $("poroParam").value],
  ];

  $("summary").innerHTML = metrics.map(([k, v, sub]) =>
    `<div><dt>${k}</dt><dd>${v} <small>${sub ?? ""}</small></dd></div>`).join("");

  const warns = g ? [...(g.warnings || [])] : [];
  if (g && g.skipped && g.skipped.length) {
    warns.push(`Ignored keywords not needed for the grid: ${g.skipped.slice(0, 12).join(", ")}` +
      (g.skipped.length > 12 ? `, +${g.skipped.length - 12} more` : ""));
  }
  const box = $("warnings");
  if (warns.length) {
    box.innerHTML = `<strong>Notes</strong><ul>${warns.map((w) => `<li>${escapeHtml(w)}</li>`).join("")}</ul>`;
    box.hidden = false;
  } else {
    box.hidden = true;
  }
}

function showDownloads(files) {
  for (const u of state.urls) URL.revokeObjectURL(u);
  state.urls = [];

  const ul = $("downloads");
  ul.innerHTML = "";
  for (const [name, text, note] of files) {
    const blob = new Blob([text], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    state.urls.push(url);

    const li = document.createElement("li");
    li.innerHTML =
      `<span class="name">${escapeHtml(name)}</span>` +
      `<span class="size">${escapeHtml(note)} · ${fmtBytes(blob.size)}</span>`;
    const a = document.createElement("a");
    a.href = url; a.download = name; a.textContent = "Download";
    li.appendChild(a);
    ul.appendChild(li);
  }
  ul.hidden = false;
}

//------------------------------------------------------------------- utils

function markStage(n, cls) {
  const li = document.querySelector(`#rail li[data-stage="${n}"]`);
  if (!li) return;
  li.classList.remove("active", "done", "pending", "locked");
  li.classList.add(cls);
}

function setStatus(id, cls, text) {
  const el = $(id);
  el.className = `status ${cls}`;
  el.textContent = text;
  el.hidden = false;
}

function fmtBytes(n) {
  if (n < 1024) return `${n} B`;
  if (n < 1024 * 1024) return `${(n / 1024).toFixed(1)} kB`;
  return `${(n / 1048576).toFixed(1)} MB`;
}

const nextFrame = () => new Promise((r) => requestAnimationFrame(() => r()));

function escapeHtml(s) {
  return String(s).replace(/[&<>"]/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
}
