# The design document

`sim-live-paper.pdf` explains in plain language how the pipeline works, stage by stage, shows the evidence that each part works, lists
what is assumed rather than measured, and records why the main design choices were made. It uses the same system as the
chromatin-studio software paper: small section files, generated figures, one build script.

| Path | What it holds |
|---|---|
| `sections/NN-name.md` | One section each, in file-name order. The comment at the top names the code it describes (`Code:`) and the checks to re-run (`Check:`). |
| `analysis/checks.py` | Six short checks on the real code (probes, transcription, E-P curve recovery, blind calibration, spot precision, MS2 deconvolution samplers). Writes `data/checks.json` with a run manifest and `data/checks-table.md`. |
| `data/recorded-results.md` | Results quoted from larger runs whose data lives outside the repository, each with its source document. Edit by hand when a source changes. |
| `figures/make_figures.py` | Draws every figure from `data/`. Redraws only what changed. |
| `stale.py` | Lists the sections a code change made stale, from the `Code:` lines and the git difference since the last documented commit. |
| `build.py`, `style.css` | Assemble the sections and print the PDF. |

## Update after a code change (keep it cheap)

```
scripts\run.cmd python docs\paper\stale.py
scripts\run.cmd python docs\paper\analysis\checks.py <names it listed>
scripts\run.cmd python docs\paper\build.py --mark-done
```

Edit only the sections `stale.py` lists. The build needs pandoc, matplotlib and Playwright with Chromium (the tutorial's shot venv has
Playwright; or set `CHROMIUM_PATH` to Chrome or Edge). It never runs a simulation. Commit the section, data, figures and PDF with the code.

## Conventions

- Plain language, short sentences, a "Short version" box at the top of each section.
- **[Assumption]** marks a modelling choice or input without a measurement; **[Unverified]** marks a statement not checked against a run.
- Generated checks and recorded results are kept apart. A limit set before a check is not moved afterwards.
- A design choice gets an entry in `sections/12-decisions.md`: what, why, what was rejected. If the reason is not known, ask; never invent one.
