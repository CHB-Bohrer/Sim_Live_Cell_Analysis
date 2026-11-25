"""Build docs/manual/SimLive_manual.pdf from the chapter files listed in manifest.yaml (takes seconds).

    python docs/manual/build.py

Needs pandoc and Playwright (Python package, or Node's `playwright` with its Chromium). Figures are NOT remade here:
run a script in figure_scripts/ only when the result it shows has changed (see the appendix of the manual)."""
import datetime as dt
import shutil
import subprocess
import tempfile
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
OUT = HERE / "SimLive_manual.pdf"


def html() -> Path:
    man = yaml.safe_load((HERE / "manifest.yaml").read_text(encoding="utf-8"))
    files = [str(HERE / c["file"]) for c in man["chapters"]]
    page = HERE / "_build.html"
    subprocess.run(["pandoc", *files, "-o", str(page), "--standalone", "--toc", "--toc-depth=1", "--number-sections",
                    "--css", "style.css", "--embed-resources", "--resource-path", str(HERE),
                    "--metadata", f"title={man['title']}", "--metadata", f"date={dt.date.today():%d %B %Y}"],
                   cwd=HERE, check=True)
    return page


def pdf(page: Path):
    footer = ('<div style="font-size:8px;width:100%;text-align:center;color:#6b7280">'
              '<span class="pageNumber"></span> / <span class="totalPages"></span></div>')
    try:
        from playwright.sync_api import sync_playwright
        with sync_playwright() as p:
            b = p.chromium.launch()
            pg = b.new_page()
            pg.goto(page.as_uri(), wait_until="load")
            pg.pdf(path=str(OUT), format="A4", print_background=True, display_header_footer=True,
                   header_template="<div></div>", footer_template=footer,
                   margin={"top": "18mm", "bottom": "18mm", "left": "16mm", "right": "16mm"})
            b.close()
    except ImportError:
        js = f"""const {{chromium}} = require('playwright');
(async () => {{ const b = await chromium.launch(); const pg = await b.newPage();
  await pg.goto({page.as_uri()!r}, {{waitUntil: 'load'}});
  await pg.pdf({{path: {str(OUT)!r}, format: 'A4', printBackground: true, displayHeaderFooter: true,
    headerTemplate: '<div></div>', footerTemplate: {footer!r},
    margin: {{top: '18mm', bottom: '18mm', left: '16mm', right: '16mm'}}}});
  await b.close(); }})();"""
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "print.js"
            f.write_text(js, encoding="utf-8")
            subprocess.run([shutil.which("node") or "node", str(f)], check=True)


if __name__ == "__main__":
    p = html()
    pdf(p)
    p.unlink()
    print("wrote", OUT)
