"""Drive the page in a real headless browser and assert on the rendered DOM.

The API self-test proves the server works; this proves the page is wired to
it. It uses Chrome or Edge in headless mode with --dump-dom, which is enough
because app.js sets data-selftest="finished" on <body> when the scripted run
completes - so the dumped DOM is the finished state, not a snapshot of a page
still working.

    python web/server/browsertest.py             browser engine only (fast)
    python web/server/browsertest.py --matlab    also run the MATLAB pipeline
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from server.app import make_server
    from server.config import Settings
else:
    from .app import make_server
    from .config import Settings

CHROME_CANDIDATES = [
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
]

_results: list[tuple[str, bool, str]] = []


def find_browser() -> Path | None:
    for name in ("chrome", "msedge", "chromium"):
        found = shutil.which(name)
        if found:
            return Path(found)
    for candidate in CHROME_CANDIDATES:
        if Path(candidate).is_file():
            return Path(candidate)
    return None


def check(label: str, ok: bool, note: str = "") -> bool:
    _results.append((label, ok, note))
    print(f"  [{'PASS' if ok else 'FAIL'}] {label:<30} {note}")
    return ok


def render(browser: Path, url: str, budget_ms: int) -> str:
    with tempfile.TemporaryDirectory() as profile:
        proc = subprocess.run(
            [
                str(browser),
                "--headless",
                "--disable-gpu",
                "--no-sandbox",
                "--no-first-run",
                f"--user-data-dir={profile}",
                f"--virtual-time-budget={budget_ms}",
                "--dump-dom",
                url,
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=budget_ms / 1000 + 120,
        )
    return proc.stdout or ""


def attr(dom: str, element_id: str) -> str | None:
    """The full opening tag for an id, or None if it is not in the DOM."""
    m = re.search(r'id="%s"[^>]*>' % re.escape(element_id), dom)
    return m.group(0) if m else None


def text_of(dom: str, element_id: str) -> str:
    m = re.search(r'id="%s"[^>]*>(.*?)</' % re.escape(element_id), dom, re.S)
    return re.sub(r"<[^>]+>", "", m.group(1)).strip() if m else ""


def visible(dom: str, element_id: str) -> bool:
    tag = attr(dom, element_id)
    return tag is not None and "hidden" not in tag


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matlab", action="store_true",
                        help="also drive stages 1-7 through MATLAB")
    parser.add_argument("--keep", action="store_true",
                        help="print the dumped DOM on failure")
    args = parser.parse_args(argv)

    browser = find_browser()
    if browser is None:
        print("No Chrome or Edge found; skipping the browser test.")
        return 0

    settings = Settings.from_env()
    server = make_server(settings, "127.0.0.1", 0)
    port = server.server_address[1]
    base = f"http://127.0.0.1:{port}"
    threading.Thread(target=server.serve_forever, daemon=True).start()
    print(f"browser test\n  browser: {browser.name}\n  server:  {base}\n")

    try:
        # ---- browser engine: no server involvement beyond static files ---
        # engine=browser is explicit: the page now defaults to MATLAB mode
        # when the backend has it, so these checks must opt out.
        t0 = time.time()
        dom = render(browser, f"{base}/?grid=param&engine=browser&run=1", 15000)
        check("page renders", len(dom) > 2000, f"{len(dom)} bytes")
        check("modules loaded (engine card)", visible(dom, "card-engine"),
              "detect() found the backend")
        check("grid built in browser", visible(dom, "card-summary"),
              text_of(dom, "grid-status")[:70])
        check("export card shown", visible(dom, "card-export"))
        check("fracture card hidden in browser mode",
              not visible(dom, "card-frac"),
              "fractures need MATLAB")
        check("browser run is fast", time.time() - t0 < 60,
              f"{time.time() - t0:.1f}s")

        # ---- variable spacing -------------------------------------------
        # The lists are seeded from the uniform extent, so every block is the
        # same size and *CART is the right answer. Asserting *VARI here would
        # be asserting a bug: the writer is supposed to collapse a uniform
        # tensor grid rather than emit a redundant DI/DJ/DK block.
        dom = render(browser, f"{base}/?grid=param&engine=browser&spacing=variable&run=1", 15000)
        status = text_of(dom, "grid-status")
        check("variable spacing builds a grid", visible(dom, "card-summary"),
              status[:70])
        check("uniform sizes collapse to *CART", "CART" in status,
              "non-uniform lists are covered byte-for-byte in selftest.html")

        if not args.matlab:
            print("\n  (skipping the MATLAB pipeline; pass --matlab to include it)")
            return report()

        if not settings.matlab_exe:
            check("MATLAB available", False, "not found")
            return report()

        # ---- MATLAB engine: stages 1-7 end to end -----------------------
        print("\n  driving the MATLAB pipeline (a few minutes)…")
        t0 = time.time()
        url = f"{base}/?grid=param&engine=matlab&run=1&frac=1&export=1"
        dom = render(browser, url, 900000)
        elapsed = time.time() - t0

        finished = 'data-selftest="finished"' in dom
        check("scripted run finished", finished, f"{elapsed:.0f}s")
        check("MATLAB grid summary", visible(dom, "card-summary"),
              text_of(dom, "grid-status")[:70])
        check("fracture card shown", visible(dom, "card-frac"))

        frac = attr(dom, "frac-status") or ""
        check("fracture pipeline succeeded", "status ok" in frac,
              text_of(dom, "frac-status")[:80])

        stages_done = len(re.findall(r'<li class="ok">', dom))
        check("all five fracture stages green", stages_done >= 5,
              f"{stages_done} of 5")

        exp = attr(dom, "export-status") or ""
        check("deck written", "status ok" in exp,
              text_of(dom, "export-status")[:90])
        check("downloads offered", "api/session/" in dom and "Download" in dom)

        if args.keep and any(not ok for _, ok, _ in _results):
            print("\n----- DOM -----\n" + dom[:8000])

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
