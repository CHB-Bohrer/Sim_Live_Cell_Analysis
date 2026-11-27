<!--
Section: Summary
Update when: a stage changes state (built / stand-in / planned) or a headline number changes.
Code: status.json
Check: none
Sources: STATUS.md, docs/notes/status_detail.md.
-->

# Summary

::: {.box}
**Short version.** In a real microscope movie nobody knows the right answer. Here the movie is invented, so we do. The pipeline simulates chromatin, cells and a microscope, runs a complete analysis on the resulting images, and grades every step against the truth it started from.
:::

Sim Live Cell Analysis is research software for one question: **how much does each step of a live-cell chromatin analysis get wrong, and why?** The experiment it imitates has a promoter and an enhancer 100 kb apart, each labelled in its own colour, and MS2 bursts of transcription at the promoter, in nuclei that move, touch and divide.

The pipeline has seven stages (Figure 1). Stages 1 to 3 invent the data and keep the truth. Stages 4 to 6 are the analysis a real experiment would run. Stage 7 compares the two.

![The seven stages. Blue boxes invent the data; green boxes analyse and grade it.](figures/fig-pipeline.svg)

What this document gives you:

- how each stage works, with the equations where they matter and plain words everywhere else;
- for each part, the evidence that it works, with the limit that was set before looking;
- which numbers are measurements, which are assumptions, and which are placeholders waiting for real data (Section 12);
- why the main design choices were made (Section 13).

Two kinds of evidence appear, and they are kept apart. **Generated checks** are re-run from the code by `docs/paper/analysis/checks.py` and write their own numbers and figures. **Recorded results** come from larger runs whose data lives outside the repository; they are quoted with their date and the file that documents them.
