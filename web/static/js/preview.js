// Top-layer depth map, drawn straight from ZCORN.
//
// Not a 3D view — deliberately. A depth map of layer 1 is enough to answer
// the question the user actually has at this point ("did my grid load, and
// does it look like my field?"), costs no dependency, and cannot be fooled by
// a camera that failed to fit.
//
// ZCORN is indexed (2*nx, 2*ny, 2*nz); the four top corners of block (i,j,k)
// live at (2i-1:2i, 2j-1:2j, 2k-1) in 1-based MATLAB terms.

export function drawTopLayer(canvas, g) {
  const [nx, ny] = g.cartDims;
  const ctx = canvas.getContext("2d");
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);

  const depth = new Float64Array(nx * ny);
  const shown = new Uint8Array(nx * ny);
  let lo = Infinity, hi = -Infinity;

  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const cell = i + nx * j;                 // layer k = 0
      if (!g.ACTNUM[cell]) continue;
      const z = topFaceDepth(g.ZCORN, nx, ny, i, j, 0);
      depth[cell] = z;
      shown[cell] = 1;
      if (z < lo) lo = z;
      if (z > hi) hi = z;
    }
  }

  if (!shown.some(Boolean)) {
    ctx.fillStyle = "#8a94a0";
    ctx.font = "13px sans-serif";
    ctx.fillText("Top layer has no active blocks.", 14, 24);
    return;
  }

  // Fit the map into the canvas, preserving aspect ratio.
  const pad = 8;
  const s = Math.min((W - 2 * pad) / nx, (H - 2 * pad) / ny);
  const ox = (W - s * nx) / 2, oy = (H - s * ny) / 2;
  const span = hi - lo || 1;

  const img = ctx.createImageData(W, H);
  const data = img.data;
  for (let py = 0; py < H; py++) {
    const jy = Math.floor((py - oy) / s);
    for (let px = 0; px < W; px++) {
      const ix = Math.floor((px - ox) / s);
      const o = 4 * (py * W + px);
      if (ix < 0 || iy_out(jy, ny) || ix >= nx) { data[o + 3] = 0; continue; }
      const cell = ix + nx * jy;
      if (!shown[cell]) { data[o + 3] = 0; continue; }
      const [r, gg, b] = ramp((depth[cell] - lo) / span);
      data[o] = r; data[o + 1] = gg; data[o + 2] = b; data[o + 3] = 255;
    }
  }
  ctx.putImageData(img, 0, 0);

  drawLegend(ctx, W, H, lo, hi);
}

function iy_out(j, ny) { return j < 0 || j >= ny; }

function topFaceDepth(zcorn, nx, ny, i, j, k) {
  // Linear index into ZCORN laid out (2nx, 2ny, 2nz), column-major.
  const NX2 = 2 * nx, NY2 = 2 * ny;
  const kz = 2 * k;                       // top plane of the block
  let sum = 0;
  for (let dj = 0; dj < 2; dj++) {
    for (let di = 0; di < 2; di++) {
      const ii = 2 * i + di, jj = 2 * j + dj;
      sum += zcorn[ii + NX2 * (jj + NY2 * kz)];
    }
  }
  return sum / 4;
}

/** Blue (shallow) to red (deep) — the convention Builder uses for Grid Top. */
function ramp(t) {
  t = Math.min(1, Math.max(0, t));
  const stops = [
    [0.00, [ 33,  70, 160]],
    [0.25, [ 40, 170, 200]],
    [0.50, [ 60, 180,  85]],
    [0.75, [235, 190,  50]],
    [1.00, [190,  40,  40]],
  ];
  for (let s = 1; s < stops.length; s++) {
    if (t <= stops[s][0]) {
      const [t0, c0] = stops[s - 1], [t1, c1] = stops[s];
      const f = (t - t0) / (t1 - t0);
      return [0, 1, 2].map((q) => Math.round(c0[q] + f * (c1[q] - c0[q])));
    }
  }
  return stops[stops.length - 1][1];
}

function drawLegend(ctx, W, H, lo, hi) {
  const w = 14, h = 110, x = W - w - 12, y = H - h - 26;
  const grad = ctx.createLinearGradient(0, y + h, 0, y);
  for (let s = 0; s <= 10; s++) {
    const [r, g, b] = ramp(s / 10);
    grad.addColorStop(s / 10, `rgb(${r},${g},${b})`);
  }
  ctx.fillStyle = grad;
  ctx.fillRect(x, y, w, h);
  ctx.strokeStyle = "#c9d1d9";
  ctx.strokeRect(x + 0.5, y + 0.5, w, h);

  ctx.fillStyle = "#3b4652";
  ctx.font = "11px ui-monospace, Consolas, monospace";
  ctx.textAlign = "right";
  ctx.fillText(`${hi.toFixed(0)} m`, x - 4, y + 9);
  ctx.fillText(`${lo.toFixed(0)} m`, x - 4, y + h);
  ctx.textAlign = "left";
}
