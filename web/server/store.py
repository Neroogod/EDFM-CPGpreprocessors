"""Sessions and background jobs.

A session is a directory holding the stage checkpoints (s1.mat .. s6.mat) and
the exported deck. A job is one stage running in a worker thread. State lives
in memory and on disk; nothing is shared between sessions, so two users of a
local server cannot collide.
"""

from __future__ import annotations

import shutil
import threading
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from .config import Settings
from .runner import run_stage

STAGE_NAMES = {
    1: "Matrix grid",
    2: "Fracture planes",
    3: "EDFM gridding",
    4: "Fracture-matrix NNCs",
    5: "Fracture-fracture NNCs",
    6: "Projection",
    7: "Export",
}

# Stage 7 reads whatever the latest checkpoint is, so it has no fixed input.
STAGE_INPUT = {2: 1, 3: 2, 4: 3, 5: 4, 6: 5}


@dataclass
class Job:
    id: str
    session_id: str
    stage: int
    state: str = "queued"           # queued | running | done | error
    created: float = field(default_factory=time.time)
    started: float | None = None
    finished: float | None = None
    result: dict[str, Any] | None = None
    # Set once the job reaches a terminal state, so a client can long-poll
    # instead of asking again every second for the several minutes a
    # corner-point EDFM run takes.
    done_event: threading.Event = field(default_factory=threading.Event)

    def wait(self, timeout_s: float) -> None:
        self.done_event.wait(timeout_s)

    def public(self) -> dict[str, Any]:
        return {
            "jobId": self.id,
            "sessionId": self.session_id,
            "stage": self.stage,
            "stageName": STAGE_NAMES.get(self.stage, f"Stage {self.stage}"),
            "state": self.state,
            "queuedSeconds": round((self.started or time.time()) - self.created, 1),
            "runningSeconds": round(
                (self.finished or time.time()) - (self.started or time.time()), 1
            ),
            "result": self.result,
        }


class Store:
    def __init__(self, settings: Settings) -> None:
        self.settings = settings
        self._lock = threading.Lock()
        self._jobs: dict[str, Job] = {}
        self._sessions: dict[str, dict[str, Any]] = {}

    # ---- sessions ------------------------------------------------------
    def create_session(self) -> dict[str, Any]:
        sid = uuid.uuid4().hex[:12]
        path = self.settings.work_dir / sid
        path.mkdir(parents=True, exist_ok=True)
        with self._lock:
            self._sessions[sid] = {
                "id": sid,
                "created": time.time(),
                "completed": [],
                "activeJob": None,
            }
        return self.session_state(sid)

    def session_dir(self, sid: str) -> Path:
        # uuid hex only; reject anything that could climb out of work_dir.
        if not sid.isalnum():
            raise KeyError(sid)
        path = self.settings.work_dir / sid
        if not path.is_dir():
            raise KeyError(sid)
        return path

    def session_state(self, sid: str) -> dict[str, Any]:
        path = self.session_dir(sid)
        with self._lock:
            meta = self._sessions.get(sid, {"id": sid, "completed": []})
            active = meta.get("activeJob")
        checkpoints = {
            n: (path / f"s{n}.mat").is_file() for n in range(1, 7)
        }
        return {
            "sessionId": sid,
            "checkpoints": checkpoints,
            "completed": sorted(n for n, ok in checkpoints.items() if ok),
            "activeJob": active,
            "canExport": any(checkpoints.values()),
            "downloads": [
                f.name
                for f in sorted(path.glob("*"))
                if f.suffix.lower()
                in {".dat", ".inc", ".data", ".grdecl", ".json", ".log"}
            ],
        }

    def delete_session(self, sid: str) -> None:
        path = self.session_dir(sid)
        shutil.rmtree(path, ignore_errors=True)
        with self._lock:
            self._sessions.pop(sid, None)

    def sweep(self, max_age_hours: float = 24.0) -> int:
        """Drop session directories older than max_age_hours."""
        cutoff = time.time() - max_age_hours * 3600
        removed = 0
        for path in self.settings.work_dir.glob("*"):
            if path.is_dir() and path.stat().st_mtime < cutoff:
                shutil.rmtree(path, ignore_errors=True)
                removed += 1
        return removed

    # ---- jobs ----------------------------------------------------------
    def start_stage(self, sid: str, stage: int, opts: dict) -> dict[str, Any]:
        path = self.session_dir(sid)

        with self._lock:
            meta = self._sessions.setdefault(
                sid, {"id": sid, "created": time.time(), "completed": []}
            )
            active_id = meta.get("activeJob")
            if active_id and self._jobs[active_id].state in ("queued", "running"):
                raise RuntimeError(
                    f"Stage {self._jobs[active_id].stage} is still running in "
                    "this session. Wait for it or start a new session."
                )

        in_file, out_file = self._stage_files(path, stage, opts)
        if stage in STAGE_INPUT and not in_file.is_file():
            raise FileNotFoundError(
                f"Stage {stage} needs the output of stage {STAGE_INPUT[stage]}, "
                "which has not run in this session."
            )

        job = Job(id=uuid.uuid4().hex[:12], session_id=sid, stage=stage)
        with self._lock:
            self._jobs[job.id] = job
            self._sessions[sid]["activeJob"] = job.id

        thread = threading.Thread(
            target=self._run,
            args=(job, path, stage, opts, in_file, out_file),
            daemon=True,
        )
        thread.start()
        return job.public()

    def _stage_files(
        self, path: Path, stage: int, opts: dict
    ) -> tuple[Path | None, Path | None]:
        if stage == 1:
            return None, path / "s1.mat"
        if stage == 7:
            # Export from an explicit stage if asked, else the latest one.
            wanted = opts.get("fromStage")
            if wanted:
                return path / f"s{int(wanted)}.mat", None
            for n in range(6, 0, -1):
                if (path / f"s{n}.mat").is_file():
                    return path / f"s{n}.mat", None
            return path / "s1.mat", None
        return path / f"s{STAGE_INPUT[stage]}.mat", path / f"s{stage}.mat"

    def _run(
        self,
        job: Job,
        path: Path,
        stage: int,
        opts: dict,
        in_file: Path | None,
        out_file: Path | None,
    ) -> None:
        job.started = time.time()
        job.state = "running"
        try:
            result = run_stage(
                self.settings,
                stage,
                path,
                opts,
                in_file,
                out_file,
                log_path=path / f"stage{stage}.log",
            )
        except Exception as exc:                      # noqa: BLE001
            result = {
                "ok": False,
                "stage": stage,
                "errors": [f"The server failed to run the stage: {exc}"],
                "warnings": [],
                "stats": {},
                "files": {},
            }

        job.result = result
        job.finished = time.time()
        job.state = "done" if result.get("ok") else "error"
        job.done_event.set()
        with self._lock:
            meta = self._sessions.get(job.session_id)
            if meta is not None:
                meta["activeJob"] = None
                if result.get("ok") and stage not in meta["completed"]:
                    meta["completed"].append(stage)

    def job(self, jid: str) -> Job:
        with self._lock:
            return self._jobs[jid]

    def cancel_marker(self, jid: str) -> None:
        """Forget a job. MATLAB is not killed - it is a batch run that ends."""
        with self._lock:
            self._jobs.pop(jid, None)


def stage_summary() -> list[dict[str, Any]]:
    return [
        {"stage": n, "name": STAGE_NAMES[n], "needs": STAGE_INPUT.get(n)}
        for n in sorted(STAGE_NAMES)
    ]
