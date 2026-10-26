"""Start the local server.

    python web/server/serve.py
    python web/server/serve.py --port 8123 --open

No dependencies beyond the standard library.
"""

from __future__ import annotations

import argparse
import sys
import webbrowser
from pathlib import Path

if __package__ in (None, ""):
    # Run directly as a path ("python web/server/serve.py"): put web/ on
    # sys.path so "server" resolves as a package and the relative imports
    # inside app.py keep working.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from server.app import make_server
    from server.config import Settings
else:
    from .app import make_server
    from .config import Settings


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="MRST preprocessor web app")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--open", action="store_true", help="open a browser")
    parser.add_argument(
        "--sweep-hours",
        type=float,
        default=24.0,
        help="delete session directories older than this at startup",
    )
    args = parser.parse_args(argv)

    settings = Settings.from_env()
    server = make_server(settings, args.host, args.port)
    removed = server.RequestHandlerClass.store.sweep(args.sweep_hours)

    url = f"http://{args.host}:{args.port}/"
    print("MRST -> CMG/Eclipse preprocessor")
    print(f"  serving   {settings.static_dir}")
    print(f"  work dir  {settings.work_dir}")
    if removed:
        print(f"  swept     {removed} old session(s)")
    if settings.matlab_exe:
        print(f"  MATLAB    {settings.matlab_exe}")
    else:
        print("  MATLAB    not found - browser-only pages still work; set")
        print("            MATLAB_EXE to enable the fracture stages")
    print(f"\n  {url}\n  Ctrl-C to stop\n")

    if args.open:
        webbrowser.open(url)

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nstopping")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
