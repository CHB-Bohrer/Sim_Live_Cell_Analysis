<!--
Section: Validation summary
Update when: any check is re-run, or a recorded result is superseded.
Code: docs/paper/analysis/checks.py, tests
Check: all
Sources: data/checks.json, data/recorded-results.md.
-->

# Validation summary

::: {.box}
**Short version.** Three checks are re-run from the code and write their own numbers. The larger results are quoted from documented runs. The unit tests (225 on CPU and 7 on GPU at 2026-11-10) cover the rest and are not repeated here.
:::

## Generated checks

Each row is produced by `docs/paper/analysis/checks.py`, with the limit fixed in the script before running. The run column is the manifest: when it ran and from which commit.

@include data/checks-table.md

## Recorded results

Runs whose data lives outside the repository. Each row names the document that holds the detail.

@include data/recorded-results.md

## What these do not show

- None of this is a comparison with a real movie. Agreement with theory and recovery of injected truth show the software is internally right, not that the simulation resembles the owner's microscope or cells.
- Several recorded results are single seeds and are marked so where quoted.
- The analysis stage is a placeholder for the owner's own analysis. Its scores describe this placeholder.
