// Browser port of buildCMGModel.m, restricted to matrix-only grids.
//
// The MATLAB version splits cleanly into "assemble a simulator-neutral
// description" (buildCMGModel) and "format it as CMG cards"
// (writeCMGModel). That seam is what makes this port possible: only the
// second half is needed to produce a deck, and neither half touches MRST.
//
// Fracture grids are deliberately out of scope here. Creating fracture cells
// requires plane-cell intersection geometry that lives in MATLAB; this module
// only handles what a GRDECL alone can describe.

const MILLI_DARCY = 9.869232667160128e-16;   // m^2, MRST's darcy constant

export function buildMatrixModel(grdecl, opts = {}) {
  const [nx, ny, nz] = grdecl.cartDims;
  const nTot = nx * ny * nz;

  const simulator = (opts.simulator || "IMEX").toUpperCase();
  const permUnit = opts.permUnit || "mD";     // input unit of PERMX etc.

  // Active mask straight from ACTNUM.
  //
  // NOTE: MRST's processGRDECL can deactivate more than ACTNUM does — it
  // discards disconnected components and degenerate cells. A deck exported
  // here therefore matches MATLAB exactly only when ACTNUM already describes
  // the active set. The UI says so rather than hiding it.
  const active = new Uint8Array(nTot);
  let nInactive = 0;
  for (let i = 0; i < nTot; i++) {
    const on = grdecl.ACTNUM[i] ? 1 : 0;
    active[i] = on;
    if (!on) nInactive++;
  }

  const props = {
    permi: fillArray(grdecl.PERMX, nTot, opts.perm ?? 100, permUnit),
    permj: fillArray(grdecl.PERMY ?? grdecl.PERMX, nTot, opts.perm ?? 100, permUnit),
    permk: fillArray(grdecl.PERMZ ?? grdecl.PERMX, nTot, opts.permz ?? opts.perm ?? 100, permUnit),
    poro:  fillArray(grdecl.PORO, nTot, opts.poro ?? 0.2, null),
    active,
    nTot,
  };

  // Zero out inactive blocks so they cannot contribute a stray value.
  for (let i = 0; i < nTot; i++) {
    if (!active[i]) {
      props.permi[i] = 0; props.permj[i] = 0; props.permk[i] = 0;
      props.poro[i] = 0;
    }
  }

  const [ztop, zbot] = zRange(grdecl.ZCORN);

  return {
    meta: {
      simulator,
      cmgVersion: opts.cmgVersion || "2024.30",
      unitSystem: "SI",
      transUnit: "md.m",
      gridFamily: "corner",
      // These two strings describe the grid, not the tool, and are matched to
      // what MATLAB emits for the same GRDECL so that a deck written here is
      // byte-identical to one written by writeCMGModel.m. Which tool produced
      // the file is recorded in the manifest instead, where it belongs.
      gridType: "processGRDECL,computeGeometry",
      fractureModel: "matrix",
      detection: "no FracGrid/Matrix fields: plain matrix grid",
      counts: {
        cells: nTot - nInactive,
        blocks: nTot,
        inactive: nInactive,
        connections: 0,
        multipliers: 0,
      },
    },
    grid: {
      kind: "corner",
      dims: [nx, ny, nz],
      coord: grdecl.COORD,
      zcorn: grdecl.ZCORN,
      ztop, zbot,
    },
    props,
    conn: [],
    mult: [],
    isolate: [],
    wells: [],
  };
}

/**
 * A grid described by numbers rather than a file: uniform Cartesian, or
 * variable spacing per direction.
 *
 * spec = { dims:[nx,ny,nz], physdim:[Lx,Ly,Lz] }        uniform  -> *CART
 *        { dims:[nx,ny,nz], dx:[], dy:[], dz:[] }       variable -> *VARI
 *        datum, permi/permj/permk/poro (arrays or scalars)
 *
 * The family follows cmgDetectGrid.m: constant spacing in every direction is
 * *CART, anything else is *VARI. MATLAB decides this from the grid's node
 * coordinates; here it is decided from the spacings, which is the same test.
 */
export function buildParametricModel(spec, opts = {}) {
  const [nx, ny, nz] = spec.dims.map((v) => Math.round(v));
  const nTot = nx * ny * nz;

  let dx, dy, dz;
  if (spec.dx && spec.dy && spec.dz) {
    dx = Array.from(spec.dx); dy = Array.from(spec.dy); dz = Array.from(spec.dz);
  } else {
    const [Lx, Ly, Lz] = spec.physdim;
    dx = new Array(nx).fill(Lx / nx);
    dy = new Array(ny).fill(Ly / ny);
    dz = new Array(nz).fill(Lz / nz);
  }
  if (dx.length !== nx || dy.length !== ny || dz.length !== nz) {
    throw new Error(
      `Block sizes do not match the dimensions: got ${dx.length}/${dy.length}/` +
      `${dz.length} for a ${nx}x${ny}x${nz} grid.`);
  }

  const uniform = isUniform(dx) && isUniform(dy) && isUniform(dz);
  const kind = uniform ? "cart" : "vari";
  const gridType = uniform
    ? "tensorGrid,cartGrid,computeGeometry"
    : "tensorGrid,computeGeometry";

  const ztop = spec.datum ?? 0;
  const zbot = ztop + dz.reduce((a, b) => a + b, 0);

  const props = {
    permi: expand(spec.permi, nTot, opts.perm ?? 100),
    permj: expand(spec.permj ?? spec.permi, nTot, opts.perm ?? 100),
    permk: expand(spec.permk ?? spec.permi, nTot, opts.permz ?? opts.perm ?? 100),
    poro:  expand(spec.poro,  nTot, opts.poro ?? 0.2),
    active: new Uint8Array(nTot).fill(1),
    nTot,
  };

  return {
    meta: {
      simulator: (opts.simulator || "IMEX").toUpperCase(),
      cmgVersion: opts.cmgVersion || "2024.30",
      unitSystem: "SI",
      transUnit: "md.m",
      gridFamily: kind,
      gridType,
      fractureModel: "matrix",
      detection: "no FracGrid/Matrix fields: plain matrix grid",
      counts: { cells: nTot, blocks: nTot, inactive: 0,
                connections: 0, multipliers: 0 },
    },
    grid: { kind, dims: [nx, ny, nz], di: dx, dj: dy, dk: dz, ztop, zbot },
    props,
    conn: [], mult: [], isolate: [], wells: [],
  };
}

function isUniform(d) {
  if (d.length <= 1) return true;
  let lo = d[0], hi = d[0];
  for (const v of d) { if (v < lo) lo = v; if (v > hi) hi = v; }
  return (hi - lo) <= 1e-9 * Math.max(1, Math.abs(hi));
}

function expand(src, n, fallback) {
  const out = new Float64Array(n);
  if (Array.isArray(src) && src.length === n) {
    for (let i = 0; i < n; i++) out[i] = src[i];
  } else if (typeof src === "number") {
    out.fill(src);
  } else {
    out.fill(fallback);
  }
  return out;
}

function fillArray(src, n, fallback, unit) {
  const out = new Float64Array(n);
  if (src && src.length === n) {
    if (unit === "m2") {
      for (let i = 0; i < n; i++) out[i] = src[i] / MILLI_DARCY;
    } else {
      for (let i = 0; i < n; i++) out[i] = src[i];
    }
  } else {
    out.fill(fallback);
  }
  return out;
}

function zRange(zcorn) {
  let lo = Infinity, hi = -Infinity;
  for (let i = 0; i < zcorn.length; i++) {
    const v = zcorn[i];
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  return [lo, hi];
}

/** cmgGridHeader.m — the two cards that must sit in the .dat, not the .inc. */
export function gridHeader(model) {
  const d = model.grid.dims;
  const kw = { cart: "*CART", vari: "*VARI", corner: "*CORNER",
               connlist: "*VARI" }[model.grid.kind];
  return [`*GRID ${kw} ${d[0]} ${d[1]} ${d[2]}`, "*KDIR *DOWN"];
}
