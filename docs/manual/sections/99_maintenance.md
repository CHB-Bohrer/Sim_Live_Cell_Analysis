# Appendix: how this document is kept up to date {#maintenance}

The document is built from one Markdown file per chapter in `docs/manual/sections/`. Each figure has its own small script in `docs/manual/figure_scripts/`. `docs/manual/manifest.yaml` lists, for each chapter, the code it describes and the scripts that make its figures.

**Update rule (to keep the cost low):**

1. A chapter is edited only when a change on GitHub touches the code listed for it. `python docs/manual/stale.py` lists those chapters for the current branch.
2. A figure is remade only when the validation result behind it changes. Figures from saved results (chromatin, E–P) are copied in, never recomputed.
3. Rebuilding the PDF takes seconds: `python docs/manual/build.py`.

Each figure script runs on the CPU in under a minute. The polymer and E–P figures come from saved results and are never rerun just for this document.
