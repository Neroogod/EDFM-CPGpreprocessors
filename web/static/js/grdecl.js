// Minimal Eclipse GRDECL reader — enough for the keywords a grid export
// needs, and no more.
//
// Handles: SPECGRID/DIMENS, COORD, ZCORN, ACTNUM, PORO, PERMX/PERMY/PERMZ,
// comments (`--`), the `n*value` repeat form, and `/` terminators.
//
// Does NOT handle INCLUDE (the file must be self-contained), units other
// than METRIC, or the many keywords irrelevant to geometry. Anything
// unrecognised is skipped rather than treated as an error, because real
// decks carry a great deal we do not need.

const NUMERIC_KEYWORDS = new Set([
  "COORD", "ZCORN", "ACTNUM", "PORO",
  "PERMX", "PERMY", "PERMZ", "NTG", "MULTPV",
]);

export function parseGRDECL(text) {
  const out = {
    cartDims: null,
    COORD: null, ZCORN: null, ACTNUM: null, PORO: null,
    PERMX: null, PERMY: null, PERMZ: null,
    keywordsSeen: [], skipped: [], warnings: [],
  };

  // Strip comments first; `--` runs to end of line.
  const lines = text.split(/\r?\n/);
  const cleaned = [];
  for (const line of lines) {
    const i = line.indexOf("--");
    cleaned.push(i >= 0 ? line.slice(0, i) : line);
  }
  const tokens = cleaned.join("\n").split(/\s+/).filter((t) => t.length > 0);

  let i = 0;
  while (i < tokens.length) {
    const tok = tokens[i++];
    if (tok === "/") continue;
    const kw = tok.toUpperCase();

    if (kw === "SPECGRID" || kw === "DIMENS") {
      const vals = readUntilSlash(tokens, i, true);
      i = vals.next;
      const n = vals.data.slice(0, 3).map((v) => Math.round(Number(v)));
      if (n.length === 3 && n.every(Number.isFinite)) out.cartDims = n;
      out.keywordsSeen.push(kw);
      continue;
    }

    if (NUMERIC_KEYWORDS.has(kw)) {
      const vals = readUntilSlash(tokens, i, false);
      i = vals.next;
      out[kw === "PERMX" ? "PERMX" : kw] = vals.data;
      out.keywordsSeen.push(kw);
      continue;
    }

    if (kw === "INCLUDE") {
      out.warnings.push(
        "INCLUDE is not supported: the .GRDECL must be self-contained. " +
        "Concatenate the included files first."
      );
      const vals = readUntilSlash(tokens, i, true);
      i = vals.next;
      continue;
    }

    // Unknown keyword. If a data block follows, skip it; otherwise move on.
    if (/^[A-Za-z]/.test(kw)) {
      if (!out.skipped.includes(kw)) out.skipped.push(kw);
      const vals = readUntilSlash(tokens, i, true);
      i = vals.next;
    }
  }

  validate(out);
  return out;
}

/**
 * Read tokens up to the next `/`. Expands `n*v` repeats. `raw` keeps strings
 * (for SPECGRID's trailing F/T flags); otherwise values are numbers.
 */
function readUntilSlash(tokens, start, raw) {
  const data = [];
  let i = start;
  for (; i < tokens.length; i++) {
    const t = tokens[i];
    if (t === "/") { i++; break; }
    // A bare keyword means the previous block had no terminator.
    if (!raw && /^[A-Za-z]{2,}$/.test(t)) break;
    const star = t.indexOf("*");
    if (star > 0) {
      const count = Number(t.slice(0, star));
      const value = t.slice(star + 1);
      if (Number.isFinite(count) && count > 0 && value.length > 0) {
        const v = raw ? value : Number(value);
        for (let k = 0; k < count; k++) data.push(v);
        continue;
      }
    }
    data.push(raw ? t : Number(t));
  }
  return { data, next: i };
}

function validate(g) {
  if (!g.cartDims) {
    throw new Error("No SPECGRID or DIMENS keyword found; grid size unknown.");
  }
  const [nx, ny, nz] = g.cartDims;
  const nCells = nx * ny * nz;

  if (!g.COORD) throw new Error("No COORD keyword found.");
  if (!g.ZCORN) throw new Error("No ZCORN keyword found.");

  const wantCoord = 6 * (nx + 1) * (ny + 1);
  const wantZcorn = 8 * nCells;
  if (g.COORD.length !== wantCoord) {
    throw new Error(
      `COORD has ${g.COORD.length} values, expected ${wantCoord} for a ` +
      `${nx}x${ny}x${nz} grid.`);
  }
  if (g.ZCORN.length !== wantZcorn) {
    throw new Error(
      `ZCORN has ${g.ZCORN.length} values, expected ${wantZcorn} for a ` +
      `${nx}x${ny}x${nz} grid.`);
  }

  for (const key of ["ACTNUM", "PORO", "PERMX", "PERMY", "PERMZ"]) {
    if (g[key] && g[key].length !== nCells) {
      g.warnings.push(
        `${key} has ${g[key].length} values, expected ${nCells}; ignored.`);
      g[key] = null;
    }
  }
  if (!g.ACTNUM) {
    g.ACTNUM = new Array(nCells).fill(1);
    g.warnings.push("No ACTNUM found; every block treated as active.");
  }
}

/** Depth range from ZCORN, for reporting. */
export function depthRange(zcorn) {
  let lo = Infinity, hi = -Infinity;
  for (let i = 0; i < zcorn.length; i++) {
    const v = zcorn[i];
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  return [lo, hi];
}
