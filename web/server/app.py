"""HTTP layer: static files plus a small JSON API over the MATLAB stages.

Standard library only, on purpose. The hosted version of this app is pure
browser JavaScript with no backend at all; this server exists so that someone
with MATLAB can run the geometry stages locally. Making them pip-install a
web framework first would be a poor trade for that.

API
    GET  /api/health                     MATLAB present? paths? stage list?
    POST /api/session                    -> {sessionId}
    GET  /api/session/<sid>              stage/checkpoint state
    DELETE /api/session/<sid>            delete the working directory
    POST /api/session/<sid>/stage/<n>    body = stage options -> {jobId}
    GET  /api/job/<jid>                  poll: queued|running|done|error
    GET  /api/session/<sid>/file/<name>  download a produced file
    POST /api/session/<sid>/upload/<name>  store a GRDECL for stage 1
"""

from __future__ import annotations

import json
import mimetypes
import re
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

from .config import Settings
from .store import Store, stage_summary

SAFE_NAME = re.compile(r"^[A-Za-z0-9._-]+$")


class Handler(BaseHTTPRequestHandler):
    server_version = "MRSTPreprocessor/1.0"
    settings: Settings
    store: Store

    # ---- plumbing ------------------------------------------------------
    def log_message(self, fmt: str, *args) -> None:      # noqa: A002
        print(f"  {self.address_string()} {fmt % args}")

    def _json(self, payload: dict, status: int = 200) -> None:
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _error(self, message: str, status: int = 400) -> None:
        self._json({"ok": False, "error": message}, status)

    def _body(self) -> dict:
        length = int(self.headers.get("Content-Length") or 0)
        if length <= 0:
            return {}
        if length > self.settings.max_upload_bytes:
            raise ValueError("Request body is too large.")
        raw = self.rfile.read(length)
        try:
            return json.loads(raw.decode("utf-8"))
        except json.JSONDecodeError as exc:
            raise ValueError(f"Body was not valid JSON: {exc}") from exc

    # ---- routing -------------------------------------------------------
    def do_GET(self) -> None:                            # noqa: N802
        parsed = urlparse(self.path)
        path = unquote(parsed.path)
        if path.startswith("/api/"):
            self._api_get(path, parse_qs(parsed.query))
        else:
            self._static(path)

    def do_POST(self) -> None:                           # noqa: N802
        path = unquote(urlparse(self.path).path)
        parts = [p for p in path.split("/") if p]
        try:
            if path == "/api/session":
                self._json(self.store.create_session())
                return

            # /api/session/<sid>/stage/<n>
            if len(parts) == 5 and parts[:2] == ["api", "session"] \
                    and parts[3] == "stage":
                self._start_stage(parts[2], parts[4])
                return

            # /api/session/<sid>/upload/<name>
            if len(parts) == 5 and parts[:2] == ["api", "session"] \
                    and parts[3] == "upload":
                self._upload(parts[2], parts[4])
                return

            self._error("Unknown endpoint.", 404)
        except KeyError:
            self._error("No such session.", 404)
        except FileNotFoundError as exc:
            self._error(str(exc), 409)
        except RuntimeError as exc:
            self._error(str(exc), 409)
        except ValueError as exc:
            self._error(str(exc), 400)

    def do_DELETE(self) -> None:                         # noqa: N802
        parts = [p for p in unquote(urlparse(self.path).path).split("/") if p]
        if len(parts) == 3 and parts[:2] == ["api", "session"]:
            try:
                self.store.delete_session(parts[2])
            except KeyError:
                self._error("No such session.", 404)
                return
            self._json({"ok": True})
        else:
            self._error("Unknown endpoint.", 404)

    # ---- API -----------------------------------------------------------
    def _api_get(self, path: str, query: dict[str, list[str]] | None = None) -> None:
        parts = [p for p in path.split("/") if p]
        query = query or {}
        try:
            if path == "/api/health":
                self._json(
                    {
                        "ok": True,
                        **self.settings.describe(),
                        "stages": stage_summary(),
                    }
                )
                return

            if len(parts) == 3 and parts[:2] == ["api", "job"]:
                job = self.store.job(parts[2])
                # ?wait=N blocks until the job ends or N seconds pass. The cap
                # keeps the socket well inside any proxy or browser idle limit.
                wait = _float(query.get("wait", ["0"])[0])
                if wait > 0 and job.state in ("queued", "running"):
                    job.wait(min(wait, 30.0))
                self._json(job.public())
                return

            if len(parts) == 3 and parts[:2] == ["api", "session"]:
                self._json(self.store.session_state(parts[2]))
                return

            # /api/session/<sid>/file/<name>
            if len(parts) == 5 and parts[:2] == ["api", "session"] \
                    and parts[3] == "file":
                self._download(parts[2], parts[4])
                return

            self._error("Unknown endpoint.", 404)
        except KeyError:
            self._error("No such session or job.", 404)

    def _start_stage(self, sid: str, stage_text: str) -> None:
        try:
            stage = int(stage_text)
        except ValueError:
            raise ValueError(f"Stage must be a number, got {stage_text!r}.")
        if not 1 <= stage <= 7:
            raise ValueError(f"Stage {stage} is out of range; expected 1 to 7.")
        opts = self._body()
        if not isinstance(opts, dict):
            raise ValueError("Stage options must be a JSON object.")
        self._json(self.store.start_stage(sid, stage, opts))

    def _upload(self, sid: str, name: str) -> None:
        if not SAFE_NAME.match(name):
            raise ValueError(f"Unsafe file name: {name!r}")
        target = self.store.session_dir(sid) / name
        length = int(self.headers.get("Content-Length") or 0)
        if length <= 0:
            raise ValueError("The upload was empty.")
        if length > self.settings.max_upload_bytes:
            raise ValueError(
                f"Upload is {length} bytes; the limit is "
                f"{self.settings.max_upload_bytes}. Raise WEB_MAX_UPLOAD."
            )
        remaining = length
        with open(target, "wb") as fh:
            while remaining > 0:
                chunk = self.rfile.read(min(1 << 20, remaining))
                if not chunk:
                    break
                fh.write(chunk)
                remaining -= len(chunk)
        self._json({"ok": True, "path": str(target), "bytes": length})

    def _download(self, sid: str, name: str) -> None:
        if not SAFE_NAME.match(name):
            self._error("Unsafe file name.", 400)
            return
        target = self.store.session_dir(sid) / name
        if not target.is_file():
            self._error(f"No such file: {name}", 404)
            return
        data = target.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(len(data)))
        self.send_header(
            "Content-Disposition", f'attachment; filename="{name}"'
        )
        self.end_headers()
        self.wfile.write(data)

    # ---- static --------------------------------------------------------
    def _static(self, path: str) -> None:
        if path in ("", "/"):
            path = "/index.html"
        root = self.settings.static_dir.resolve()
        target = (root / path.lstrip("/")).resolve()
        # resolve() then check containment: this is what stops ../ traversal.
        if not str(target).startswith(str(root)) or not target.is_file():
            self.send_error(404, "Not found")
            return
        ctype, _ = mimetypes.guess_type(str(target))
        data = target.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", ctype or "application/octet-stream")
        self.send_header("Content-Length", str(len(data)))
        # no-store, not no-cache. This server is edited while it is running,
        # and a stale app.js that merely *looks* like a bug costs far more
        # than re-sending a few hundred kilobytes over loopback.
        self.send_header("Cache-Control", "no-store, must-revalidate")
        self.end_headers()
        self.wfile.write(data)


def _float(text: str) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return 0.0


def make_server(settings: Settings, host: str, port: int) -> ThreadingHTTPServer:
    store = Store(settings)

    class Bound(Handler):
        pass

    Bound.settings = settings
    Bound.store = store
    server = ThreadingHTTPServer((host, port), Bound)
    server.daemon_threads = True
    return server
