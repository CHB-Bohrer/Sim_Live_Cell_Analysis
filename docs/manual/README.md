# The manual (PDF)

`SimLive_manual.pdf` explains how the software works, chapter by chapter, with a figure or a number showing that each part works.

- Chapters: `sections/*.md` (one per stage; order and code map in `manifest.yaml`).
- Figures: `figures/` (made by `figure_scripts/fig_*.py`, CPU, under a minute each; run them in the simlive env, e.g.
  `scripts\run.cmd python docs\manual\figure_scripts\fig_spot_precision.py`). Chromatin and E-P figures come from saved results in `docs/img/`.
- Build: `python docs/manual/build.py` (needs pandoc + Playwright; seconds).

**Update rule:** run `python docs/manual/stale.py` after code changes. Edit only the chapters it lists. Rerun a figure script only when the
result it shows has changed. Rebuild the PDF once per PR that changed a chapter.
