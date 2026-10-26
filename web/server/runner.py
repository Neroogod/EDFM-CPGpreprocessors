"""Invoke one pipeline stage in MATLAB and read back its result JSON.

The stage functions already write a result JSON on both success and failure,
so this module does not parse console text to decide what happened - it reads
the file. Console output is kept only as a log for the user to look at.
"""

from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path
from typing import Any

from .config import Settings


def ml_str(value: str | Path) -> str:
    """Quote a value as a MATLAB single-quoted literal."""
    return "'" + str(value).replace("\\", "/").replace("'", "''") + "'"


def build_command(
    settings: Settings,
    stage: int,
    in_file: Path | None,
    out_file: Path | None,
    opts_file: Path,
    result_file: Path,
) -> str:
    return "; ".join(
        [
            f"addpath({ml_str(settings.matlab_dir)})",
            "webStage({}, {}, {}, {}, {})".format(
                stage,
                ml_str(in_file or ""),
                ml_str(out_file or ""),
                ml_str(opts_file),
                ml_str(result_file),
            ),
        ]
    )


def run_stage(
    settings: Settings,
    stage: int,
    session_dir: Path,
    opts: dict,
    in_file: Path | None,
    out_file: Path | None,
    log_path: Path | None = None,
) -> dict[str, Any]:
    if settings.matlab_exe is None:
        return {
            "ok": False,
            "stage": stage,
            "errors": [
                "MATLAB was not found. Set MATLAB_EXE to the full path of "
                "matlab.exe, or use the browser-only pages, which need no MATLAB."
            ],
            "warnings": [],
            "stats": {},
            "files": {},
        }

    session_dir.mkdir(parents=True, exist_ok=True)
    opts_file = session_dir / f"opts{stage}.json"
    result_file = session_dir / f"result{stage}.json"
    opts_file.write_text(json.dumps(opts), encoding="utf-8")
    if result_file.exists():
        result_file.unlink()

    command = build_command(
        settings, stage, in_file, out_file, opts_file, result_file
    )

    t0 = time.time()
    timed_out = False
    try:
        proc = subprocess.run(
            [str(settings.matlab_exe), "-batch", command],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=settings.timeout_s,
            cwd=str(session_dir),
        )
        console = (proc.stdout or "") + (proc.stderr or "")
        returncode: int | None = proc.returncode
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        console = _text(exc.stdout) + _text(exc.stderr)
        returncode = None

    elapsed = round(time.time() - t0, 2)
    if log_path is not None:
        log_path.write_text(console, encoding="utf-8")

    # The result file is authoritative. webStage writes it inside its own
    # try/catch, so it exists even when the stage threw; only a MATLAB that
    # died before reaching the catch leaves it missing.
    if result_file.exists():
        try:
            result = json.loads(result_file.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            result = _synthetic(stage, [f"Result JSON was unreadable: {exc}"])
    elif timed_out:
        result = _synthetic(
            stage,
            [
                f"The stage exceeded the {settings.timeout_s} s limit and was "
                f"stopped. Raise WEB_TIMEOUT_S for large grids."
            ],
        )
    else:
        result = _synthetic(
            stage,
            [
                "MATLAB produced no result file. It probably failed to start "
                f"or crashed (exit code {returncode}). See the log.",
                *_error_lines(console),
            ],
        )

    result.setdefault("stage", stage)
    result["elapsedSeconds"] = elapsed
    result["timedOut"] = timed_out
    result["returncode"] = returncode
    result["log"] = console[-20000:]
    return result


def _synthetic(stage: int, errors: list[str]) -> dict[str, Any]:
    return {
        "ok": False,
        "stage": stage,
        "errors": errors,
        "warnings": [],
        "stats": {},
        "files": {},
    }


def _error_lines(console: str, limit: int = 8) -> list[str]:
    lines = [
        line.strip()
        for line in console.splitlines()
        if line.strip().startswith(("Error", "error:", "Undefined", "Unrecognized"))
    ]
    return lines[:limit]


def _text(value: Any) -> str:
    return value if isinstance(value, str) else ""
