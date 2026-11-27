<!--
Section: Keeping this document current
Update when: the paper's build, checks or staleness tool change.
Code: docs/paper/build.py, docs/paper/stale.py, docs/paper/analysis/checks.py, docs/paper/figures/make_figures.py
Check: none
Sources: docs/paper/README.md.
-->

# Keeping this document current

::: {.box}
**Short version.** Each section names the code it describes. After a code change, one command lists the sections that are now stale and the checks to re-run. Only those are touched.
:::

1. `scripts\run.cmd python docs\paper\stale.py` lists stale sections and their checks.
2. Edit only the listed section files. Read the code difference, not whole files.
3. Re-run only the listed checks: `scripts\run.cmd python docs\paper\analysis\checks.py <name>`. Numbers in the validation table and figures come from the check's output, never from typing.
4. A design choice gets a row in Section 13. If its reason is not known, ask; do not invent one.
5. `scripts\run.cmd python docs\paper\build.py --mark-done` redraws stale figures, prints the PDF and records the commit as documented.

The build never runs a simulation. Recorded results change only when their source document does.
