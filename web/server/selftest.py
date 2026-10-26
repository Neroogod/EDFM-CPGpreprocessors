"""End-to-end check of the HTTP API.

Starts the real server on an ephemeral port and drives it with urllib, so
routing, JSON handling, polling and downloads are all exercised the way the
browser will exercise them.

    python web/server/selftest.py            full run, needs MATLAB
    python web/server/selftest.py --quick    skip the MATLAB stages
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from server.app import make_server
    from server.config import Settings
else:
    from .app import make_server
    from .config import Settings

PASS, FAIL = "PASS", "FAIL"
_results: list[tuple[str, bool, str]] = []


def check(label: str, ok: bool, note: str = "") -> bool:
    _results.append((label, ok, note))
    print(f"  [{PASS if ok else FAIL}] {label:<28} {note}")
    return ok


def call(base: str, method: str, path: str, payload=None, raw: bool = False):
    data = None
    headers = {}
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(base + path, data=data, headers=headers,
                                 method=method)
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            body = resp.read()
            return resp.status, (body if raw else json.loads(body))
    except urllib.error.HTTPError as exc:
        body = exc.read()
        try:
            return exc.code, json.loads(body)
        except json.JSONDecodeError:
            return exc.code, {"error": body.decode("utf-8", "replace")}


def wait(base: str, job_id: str, limit_s: float = 1800) -> dict:
    t0 = time.time()
    while time.time() - t0 < limit_s:
        _, job = call(base, "GET", f"/api/job/{job_id}")
        if job.get("state") in ("done", "error"):
            return job
        time.sleep(1.0)
    return {"state": "timeout", "result": {"errors": ["poll timed out"]}}


def run_stage(base: str, sid: str, stage: int, opts: dict) -> dict:
    status, started = call(base, "POST", f"/api/session/{sid}/stage/{stage}", opts)
    if status != 200:
        return {"state": "error", "result": {"errors": [str(started)]}}
    return wait(base, started["jobId"])


def stage_note(job: dict) -> str:
    result = job.get("result") or {}
    if job.get("state") == "done":
        stats = result.get("stats") or {}
        bits = []
        for key in list(stats)[:4]:
            value = stats[key]
            if isinstance(value, (int, float, str)):
                bits.append(f"{key}={value}")
        return " ".join(bits)
    return "; ".join((result.get("errors") or ["no detail"])[:2])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true",
                        help="skip stages that need MATLAB")
    args = parser.parse_args(argv)

    settings = Settings.from_env()
    server = make_server(settings, "127.0.0.1", 0)
    port = server.server_address[1]
    base = f"http://127.0.0.1:{port}"
    threading.Thread(target=server.serve_forever, daemon=True).start()
    print(f"web API self-test\n  server: {base}\n")

    try:
        # ---- health and static ----------------------------------------
        status, health = call(base, "GET", "/api/health")
        check("health", status == 200 and health["ok"],
              f"matlab={health.get('matlabAvailable')}")
        check("stage list", len(health.get("stages", [])) == 7)

        status, page = call(base, "GET", "/", raw=True)
        check("static index", status == 200 and b"<" in page,
              f"{len(page)} bytes")

        # ---- traversal and validation must be refused ------------------
        status, _ = call(base, "GET", "/../server/config.py", raw=True)
        check("path traversal blocked", status == 404)

        status, sess = call(base, "POST", "/api/session")
        sid = sess["sessionId"]
        check("create session", status == 200 and bool(sid), sid)

        status, body = call(base, "POST", f"/api/session/{sid}/stage/9", {})
        check("bad stage refused", status == 400, str(body.get("error"))[:48])

        status, body = call(base, "POST", "/api/session/nosuch/stage/1", {})
        check("bad session refused", status == 404)

        status, body = call(base, "POST", f"/api/session/{sid}/stage/4", {})
        check("out-of-order refused", status == 409,
              str(body.get("error"))[:60])

        status, body = call(base, "GET", f"/api/session/{sid}/file/..%2Fconfig.py")
        check("download traversal blocked", status in (400, 404))

        if args.quick or not health.get("matlabAvailable"):
            print("\n  (skipping MATLAB stages)")
            return report()

        # ---- the real pipeline ----------------------------------------
        job = run_stage(base, sid, 1, {
            "source": "cartesian", "dims": [12, 10, 6],
            "physdim": [120, 100, 30], "perm": 100, "poro": 0.2,
            "datum": 2000,
        })
        check("stage 1 grid", job["state"] == "done", stage_note(job))

        # Export straight off stage 1: the matrix-only branch.
        job = run_stage(base, sid, 7, {
            "format": "cmg", "caseName": "matrixonly", "fromStage": 1,
            "hostDeck": "stub",
        })
        ok = job["state"] == "done"
        check("stage 7 matrix only", ok, stage_note(job))

        status, state = call(base, "GET", f"/api/session/{sid}")
        names = state.get("downloads", [])
        check("downloads listed", any(n.endswith(".dat") for n in names),
              ", ".join(names[:4]))

        target = next((n for n in names if n.endswith(".inc")), None)
        if target:
            status, blob = call(base, "GET",
                                f"/api/session/{sid}/file/{target}", raw=True)
            check("file download", status == 200 and len(blob) > 100,
                  f"{target} {len(blob)} bytes")

        job = run_stage(base, sid, 2, {
            "mode": "explicit",
            "planes": [
                {"points": [[35, 10, 5], [35, 90, 5], [35, 90, 25], [35, 10, 25]],
                 "aperture": 1e-3, "poro": 0.5, "perm": 1e4},
                {"points": [[10, 55, 5], [110, 55, 5], [110, 55, 25], [10, 55, 25]],
                 "aperture": 1e-3, "poro": 0.5, "perm": 1e4},
            ],
        })
        check("stage 2 fractures", job["state"] == "done", stage_note(job))

        job = run_stage(base, sid, 3, {"tol": 1e-5})
        check("stage 3 EDFM grid", job["state"] == "done", stage_note(job))

        job = run_stage(base, sid, 4, {"tol": 1e-5})
        check("stage 4 frac-matrix", job["state"] == "done", stage_note(job))

        job = run_stage(base, sid, 5, {"tol": 1e-5})
        check("stage 5 frac-frac", job["state"] == "done", stage_note(job))

        job = run_stage(base, sid, 6, {"method": "pEDFM", "tol": 1e-5})
        check("stage 6 projection", job["state"] == "done", stage_note(job))

        job = run_stage(base, sid, 7, {
            "format": "cmg", "caseName": "fractured", "hostDeck": "stub",
        })
        stats = (job.get("result") or {}).get("stats") or {}
        check("stage 7 fractured", job["state"] == "done", stage_note(job))
        check("export is pEDFM",
              str(stats.get("fractureModel", "")).lower() == "pedfm",
              f"connections={stats.get('connections')}")

        status, _ = call(base, "DELETE", f"/api/session/{sid}")
        check("delete session", status == 200)

        return report()
    finally:
        server.shutdown()
        server.server_close()


def report() -> int:
    failed = [label for label, ok, _ in _results if not ok]
    print(f"\n  {len(_results) - len(failed)}/{len(_results)} checks passed")
    if failed:
        print("  FAILED: " + ", ".join(failed))
        return 1
    print("  ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
