# Overview {#overview}

::: summary
**In one sentence.** We make fake microscope movies where the right answer is known. Then we run the analysis on them and measure how much each step gets wrong.
:::

## Why simulate?

In a real live-cell movie nobody knows exactly where a locus was, or when a promoter fired. So we can't tell whether an analysis is right. It might also look right for the wrong reason.

Here the movie is invented, so every true position, every true burst and every true cell identity is saved next to it. Each analysis step can be graded against that truth. We can also turn each microscope imperfection on or off, one at a time, and see what it does to the final biological answer. The biological question is how the enhancer–promoter (E–P) distance relates to transcription.

## The pipeline at a glance

```{=html}
<div class="flow">
  <div class="box built">1&nbsp; Chromatin<br><small>polymer + loop extrusion, probes, transcription</small></div>
  <div class="arrow">&rarr;</div>
  <div class="box built">2&nbsp; Cells<br><small>moving, dividing nuclei</small></div>
  <div class="arrow">&rarr;</div>
  <div class="box built">3&nbsp; Microscope<br><small>optics, camera, 28 error sources</small></div>
  <div class="arrow">&rarr;</div>
  <div class="box built">4&nbsp; Segment &amp; track<br><small>Cellpose, Trackastra</small></div>
  <div class="arrow">&rarr;</div>
  <div class="box partial">5&nbsp; Link<br><small>loci to cells</small></div>
  <div class="arrow">&rarr;</div>
  <div class="box partial">6&nbsp; Analysis<br><small>spots, MS2, Pol II</small></div>
  <div class="arrow">&rarr;</div>
  <div class="box built">7&nbsp; Grade<br><small>against the truth</small></div>
</div>
<p class="legend"><span class="key built"></span> built and tested &nbsp; <span class="key partial"></span> working, but a stand-in or partly built</p>
```

Each stage reads and writes plain files (tables, TIFF images), so any stage can be swapped out. All randomness comes from one seed per movie, so every movie can be remade exactly from its settings file and seed.

## How to read this document

Each chapter describes one part of the software. It covers what the part does, the main assumptions, **how we know it works** (a figure or a number checked against theory or a known result), and its limits. Numbers marked *placeholder* are invented values that wait for real measurements. You can skip the boxes marked *detail*.
