// Talk to the local MATLAB backend, when there is one.
//
// Every call here is optional. Served as a static site there is no /api, the
// probe fails, and the page stays in browser-only mode - which is the mode
// that has to keep working, since that is how it will be hosted.

const TIMEOUT_MS = 4000;

async function json(method, path, body) {
  const opts = { method, headers: {} };
  if (body !== undefined) {
    opts.headers["Content-Type"] = "application/json";
    opts.body = JSON.stringify(body);
  }
  const resp = await fetch(path, opts);
  let payload;
  try {
    payload = await resp.json();
  } catch {
    throw new Error(`${method} ${path} returned ${resp.status} (not JSON).`);
  }
  if (!resp.ok) {
    throw new Error(payload.error || `${method} ${path} failed (${resp.status}).`);
  }
  return payload;
}

/** Is a backend there, and does it have MATLAB? Never throws. */
export async function detect() {
  try {
    const ctl = new AbortController();
    const timer = setTimeout(() => ctl.abort(), TIMEOUT_MS);
    const resp = await fetch("api/health", { signal: ctl.signal });
    clearTimeout(timer);
    if (!resp.ok) return { available: false };
    const health = await resp.json();
    return {
      available: true,
      matlab: !!health.matlabAvailable,
      matlabExe: health.matlabExe,
      workDir: health.workDir,
      stages: health.stages || [],
    };
  } catch {
    return { available: false };
  }
}

export class Session {
  constructor(id) {
    this.id = id;
  }

  static async create() {
    const s = await json("POST", "api/session");
    return new Session(s.sessionId);
  }

  state() {
    return json("GET", `api/session/${this.id}`);
  }

  destroy() {
    return json("DELETE", `api/session/${this.id}`);
  }

  async upload(name, text) {
    const resp = await fetch(`api/session/${this.id}/upload/${name}`, {
      method: "POST",
      headers: { "Content-Type": "application/octet-stream" },
      body: text,
    });
    if (!resp.ok) {
      const err = await resp.json().catch(() => ({}));
      throw new Error(err.error || `Upload failed (${resp.status}).`);
    }
    return resp.json();
  }

  fileUrl(name) {
    return `api/session/${this.id}/file/${encodeURIComponent(name)}`;
  }

  /**
   * Start a stage and wait for it to finish.
   *
   * State comes from a long poll - the request hangs on the server until the
   * stage ends or 25 s pass - so the result arrives the moment it is ready
   * rather than up to a poll interval later. The clock is driven separately
   * by a local timer, because a corner-point EDFM run takes minutes and a
   * spinner with no clock looks hung.
   */
  async runStage(stage, opts, onTick) {
    const started = await json("POST", `api/session/${this.id}/stage/${stage}`,
                               opts || {});
    let job = started;
    const ticker = onTick ? setInterval(() => onTick(job), 1000) : null;
    try {
      while (job.state === "queued" || job.state === "running") {
        job = await json("GET", `api/job/${started.jobId}?wait=25`);
      }
    } finally {
      if (ticker) clearInterval(ticker);
    }
    if (job.state !== "done") {
      const errs = (job.result && job.result.errors) || ["The stage failed."];
      const err = new Error(errs[0]);
      err.detail = errs.slice(1);
      err.job = job;
      throw err;
    }
    return job.result;
  }
}
