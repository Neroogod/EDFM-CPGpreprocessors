// MATLAB-compatible number formatting.
//
// The deck this module writes must be byte-identical to writeCMGModel.m's
// output, so the printf conversions have to match exactly. JavaScript gets
// two of them wrong by default:
//
//   %.6e   MATLAB gives 6.292206e+02; toExponential(6) gives 6.292206e+2
//          (exponent is not zero-padded to two digits)
//   %d     MATLAB prints a non-integer in %d as though it were %e, which is
//          how the legacy NorneOut2CMGOMO output came to be in e-notation
//
// Everything here is deliberately allocation-light: a Johansen ZCORN array
// is 880,000 values, and naive string concatenation is quadratic.

export function fmtF(v, prec) {
  // C "%.<prec>f". toFixed matches, including the -0 case.
  if (!Number.isFinite(v)) return v > 0 ? "Inf" : (v < 0 ? "-Inf" : "NaN");
  let s = v.toFixed(prec);
  if (s === "-" + (0).toFixed(prec)) s = (0).toFixed(prec); // avoid "-0.000000"
  return s;
}

export function fmtE(v, prec) {
  // C "%.<prec>e": mantissa, 'e', sign, at least two exponent digits.
  if (!Number.isFinite(v)) return v > 0 ? "Inf" : (v < 0 ? "-Inf" : "NaN");
  const s = v.toExponential(prec);
  const m = /^(-?[\d.]+)e([+-])(\d+)$/.exec(s);
  if (!m) return s;
  const exp = m[3].length < 2 ? "0" + m[3] : m[3];
  return `${m[1]}e${m[2]}${exp}`;
}

export function fmtD(v) {
  return String(Math.round(v));
}

export function fmtG(v, prec) {
  // C "%.<prec>g": choose %f or %e by exponent, then strip trailing zeros.
  // Used for the *DI/*DJ/*DK block sizes, where MATLAB writes 50 rather than
  // 50.000000 — so a naive toFixed would not match.
  if (!Number.isFinite(v)) return v > 0 ? "Inf" : (v < 0 ? "-Inf" : "NaN");
  if (v === 0) return "0";
  const P = prec > 0 ? prec : 1;
  const X = Number(v.toExponential(P - 1).split("e")[1]);
  let s;
  if (X >= -4 && X < P) {
    s = v.toFixed(Math.max(0, P - 1 - X));
    if (s.indexOf(".") >= 0) s = s.replace(/\.?0+$/, "");
  } else {
    s = v.toExponential(P - 1);
    let [m, e] = s.split("e");
    if (m.indexOf(".") >= 0) m = m.replace(/\.?0+$/, "");
    const sign = e[0];
    let digits = e.slice(1);
    if (digits.length < 2) digits = "0" + digits;
    s = `${m}e${sign}${digits}`;
  }
  return s;
}

/** Pick the formatter for a printf-style spec used by writeCMGModel. */
export function formatter(spec) {
  switch (spec) {
    case "%.6f": return (v) => fmtF(v, 6);
    case "%.4f": return (v) => fmtF(v, 4);
    case "%.6e": return (v) => fmtE(v, 6);
    case "%.8e": return (v) => fmtE(v, 8);
    case "%.8g": return (v) => fmtG(v, 8);
    case "%d":   return fmtD;
    default: throw new Error(`unsupported format ${spec}`);
  }
}

/**
 * writeCMGModel's localWrap: `perLine` values per line, each preceded by one
 * space, then one blank line after the block.
 */
export function wrapValues(out, values, spec, perLine = 6) {
  const f = formatter(spec);
  const n = values.length;
  const full = perLine * Math.floor(n / perLine);
  let row = new Array(perLine);
  for (let i = 0; i < full; i += perLine) {
    for (let k = 0; k < perLine; k++) row[k] = f(values[i + k]);
    out.push(" " + row.join(" "));
  }
  if (full < n) {
    const tail = [];
    for (let i = full; i < n; i++) tail.push(f(values[i]));
    out.push(" " + tail.join(" "));
  }
  out.push("");
}

/** True when every entry equals the first, to within writeCMGModel's test. */
export function isConstant(values) {
  if (values.length === 0) return false;
  let lo = values[0], hi = values[0];
  for (let i = 1; i < values.length; i++) {
    const v = values[i];
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  return (hi - lo) <= 1e-12 * Math.max(1, Math.abs(hi));
}
