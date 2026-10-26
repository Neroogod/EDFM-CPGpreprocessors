"""Where MATLAB, MRST and the working directory live.

Everything is overridable by environment variable and everything has a
default that works on this machine, so `python serve.py` needs no arguments.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path

# web/server/config.py -> web/server -> web -> repo root
REPO_ROOT = Path(__file__).resolve().parents[2]
WEB_ROOT = REPO_ROOT / "web"

_MATLAB_GLOBS = (
    r"C:\Program Files\MATLAB\*\bin\matlab.exe",
    r"C:\Program Files (x86)\MATLAB\*\bin\matlab.exe",
)


def find_matlab() -> Path | None:
    """Locate matlab.exe: env var, then PATH, then the usual install roots."""
    env = os.environ.get("MATLAB_EXE", "").strip()
    if env and Path(env).is_file():
        return Path(env)

    on_path = shutil.which("matlab")
    if on_path:
        return Path(on_path)

    for pattern in _MATLAB_GLOBS:
        root = Path(pattern).parts[0]
        rest = str(Path(*Path(pattern).parts[1:]))
        matches = sorted(Path(root).glob(rest))
        if matches:
            return matches[-1]          # newest release by name
    return None


@dataclass(frozen=True)
class Settings:
    mrst_root: Path
    web_root: Path
    work_dir: Path
    matlab_exe: Path | None
    timeout_s: int
    max_upload_bytes: int

    @classmethod
    def from_env(cls) -> "Settings":
        work = Path(
            os.environ.get("WEB_WORK_DIR", str(WEB_ROOT / "work"))
        ).resolve()
        work.mkdir(parents=True, exist_ok=True)
        return cls(
            mrst_root=Path(os.environ.get("MRST_ROOT", str(REPO_ROOT))).resolve(),
            web_root=WEB_ROOT,
            work_dir=work,
            matlab_exe=find_matlab(),
            timeout_s=int(os.environ.get("WEB_TIMEOUT_S", "1800")),
            max_upload_bytes=int(os.environ.get("WEB_MAX_UPLOAD", str(200 << 20))),
        )

    @property
    def matlab_dir(self) -> Path:
        return self.web_root / "matlab"

    @property
    def static_dir(self) -> Path:
        return self.web_root / "static"

    def describe(self) -> dict:
        return {
            "mrstRoot": str(self.mrst_root),
            "workDir": str(self.work_dir),
            "matlabExe": str(self.matlab_exe) if self.matlab_exe else None,
            "matlabAvailable": self.matlab_exe is not None,
            "timeoutSeconds": self.timeout_s,
        }
