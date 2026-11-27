# Dashboard pages, performance, plots, walk-through and the 2026-10-20 fix round

Working notes moved verbatim from `CLAUDE.md` on 2026-11-27 (section order kept). Add new dated sections at the top; delete statements that are no longer true.

## Known problems fixed, tutorial link (2026-10-20; user: "Fix all known problems and make a link to the tutorial on the main page please")
- **Tutorial link on the dashboard home page** (both home versions): `style.tutorial_banner` (styled HTML, wraps on narrow windows) links to
  `app/static/tutorial/index.html` (opens in a new tab) and to the video release. `style.publish_tutorial()` copies `docs/tutorial` into
  `app/static/tutorial` (changed files only; gitignored): a page served over http may not link to `file://`, Streamlit refuses symlinks that
  leave `app/static`, and Streamlit 1.65 serves .html/.js with the right content types (checked). docs/tutorial stays the source.
- **Multi-copy linker:** `multicopy.untangle_by_place` + `repair` (docs/SPOT_FINDING.md "Copy swaps through wrong picks and gaps"): copy swaps made
  through wrong picks or gaps were invisible to the step costs (showcase, calibrated, tracked IDs: 7.9% on the other copy -> 0, rms 0.80 -> 0.15 um).
  25 sweep movies re-analysed, all three variants better (default: other copy at most 0.8%, was up to 2.5%; wrong picks 2.3 -> 0.6% at 1 stray/min),
  none worse by > 1 point. Dropped on the way (they made things worse): absolute territories in the nucleus frame (orientation flips by 180 deg on
  round nuclei), the pass on per-colour paths before pairing (stray paths at 4/min). 3D with two copies works (biplane run seed 72: 8/8 copies);
  3D scores now report stray / other-copy picks too.
- **E-P curve:** wrong-pick mixture (`with_outliers`, `assumed_outlier_rate`; corrected 4.11-4.13 vs reference 4.07 at 5-40% wrong picks when the
  share is supplied); bootstrap re-fits the prior per replicate; the "intervals 25% too narrow" claim WITHDRAWN (scan settings reused the same 20
  populations; 60 fresh ones give coverage 0.95). E-P page has the checkbox. docs/EP_CURVE.md "Update 2026-10-20".
- **Crash fixed:** cell placement drew nuclei one at a time and gave up when early ones blocked the rest (4 nuclei in the 40 um biplane field
  crashed; 8 never fit). `motion.place_cells` redraws the whole layout (100 tries; the first try uses the old random numbers, so every run that
  worked is unchanged, tested) and `placement_problem(cfg)` lets the sidebar refuse an impossible layout BEFORE starting, with the reason.
- **Transcription:** default coupling is now hill, half effect at 400 nm (transcription.DEFAULTS, loci_demo/loci_library/loci_3d configs, page):
  the user's E-P choice; the sharp 150 nm contact was reached in 0.75% of frames of the library distances (400 nm: 11.8%), so it barely coupled.
  BUG fixed: the sidebar hand-off fell back to the OLD TFF1-like gene when the Transcription tab had not been opened. Step 5 opens on the first burst.
- **Dashboard:** run summary + tracking scores only above pages that use the run (`NOT_ABOUT_RUN`: Scans, E-P curve, Library, Chromatin); one
  library selector on the Library page (and Chromatin opens on it); the live monitor checked during a real 2-simulation GPU library run, and a
  stalled generator no longer shows "now running" / time left; stage 7b table on the Metrics page (`propagation_view`); guide / glossary updated.
- **Code:** broad `except Exception` blocks narrowed to what they mean; `tests/test_dashboard_pages.py` (every page on the saved showcase runs + the
  empty home page; skipped where the runs do not exist; it caught a syntax slip in library_tab.py that this round introduced); `tests/test_segtrack_gpu.py`
  (Cellpose + Trackastra on a synthetic movie); new tests for placement, untangle_by_place, repair, wrong-pick correction.
- **Lesson (again):** do not edit `src/` while a pipeline run or a re-analysis is going: worker processes import modules from disk (a 3D run crashed
  with "name 'repair' is not defined" between two edits of one file). Do multi-step source edits in one write, or wait.
- **Needs the user (cannot be fixed here):** measured values for every placeholder (MS2 brightness R, instrument / camera, probe kinetics, stray-spot
  rate and brightness, wrong-pick share for the E-P fit), a measured MSD for the polymer calibration, a real movie; which meaning of nucleus
  "steadiness"; whether loops should live longer (more blocks per simulation) and whether a surrogate distance process should cover > 61 min movies;
  the plain-polymer library (73 h of GPU) needs a go-ahead; `gh auth refresh -s workflow` to enable docs/ci/syntax.yml. (The two stale
  `.git/worktrees` entries that made git print "Permission denied" are gone since the 2026-11-21 move: the copy had normal permissions.)
- **Planned work, not defects (unchanged):** stage 5 stitching of broken / switched tracks, 3D tracking and 3D cell isolation, Ultrack / TrackMate,
  density and frame-interval scans, choosing spots with promoter + enhancer + MS2 jointly, sister chromatids, E-P pool from the 100-simulation library.

## Dashboard walk-through findings (2026-10-09)
- Tested: all 24 coupling shape x rate combinations, frame intervals 5/20/120/600 s (warnings sensible), both sweeps, probe presets,
  Library/Chromatin tabs (Streamlit `AppTest` scripts + the browser), and 4 sidebar movies (toy, library+probes, library+MS2, toy+MS2) all
  finish with DET/CHOTA 1.0 (8 cells, 1 seed; the default MS2 spots, up to ~250 loops x 100 photons, did NOT disturb segmentation; more seeds
  and brighter spots not tested).
- BUG FIXED (important gotcha): every library-source movie died at the tracker with `OSError WinError 127 ... torch\lib\shm.dll`. Cause: the
  FIRST numpy BLAS call (any `@`/dot/`scipy Rotation.apply`) in the main process loads MKL's OpenMP runtimes (libiomp5md + libomp), after which
  `import torch` (Trackastra) fails. Code that runs before the tracker must avoid matmul; `library_loci.py` now uses a broadcast sum and its
  own `_random_rotation`. Diagnose with `psutil.Process().memory_maps()` filtered on "omp". Same class of problem as the Cellpose worker note.
- Fixed: sidebar probe/transcription checkboxes did nothing (never passed on); movie locus positions came from the config's fixed 5050/4950 kb
  (both loci on the last monomer in a 3 Mb library) -> now from the Transcription tab or the middle of the library region; Transcription tab default
  positions follow region size; run log decode crashed on non-ASCII (cp1252) -> UTF-8; figure fonts were tiny.
- `st.rerun(scope="fragment")` buttons work in the browser but raise under AppTest (test artefact only).
- All four "not done" items of this walk-through were fixed on 2026-10-20: the live monitor was checked during a real (test) library run; one library selector per
  page; the run's summary and tracking scores are no longer shown above pages that do not use the run (Scans, E-P curve, Library, Chromatin; `NOT_ABOUT_RUN`);
  Transcription step 5 opens on the cell's first Pol II; the coupling default is the Hill curve at 400 nm (the 150 nm contact was reached in 0.75% of frames).
- Literature checked against PMC full texts, see docs/PLAN_transcription.md: TFF1 ON 16 min (not 6), burst period 66-86 min, 1.5 transcripts/burst;
  Pol II 2.6 kb/min, 3' dwell 116 s (our default 30 s). Defaults not changed: needs the user's choice of gene.

- Scans page (fixed 2026-10-14): it assumed the tracking-sweep layout (`masks` column + CTC metrics) and crashed on opening because the newest sweep
  (`ms2_brightness`, from `scripts/sweep_ms2_brightness.py`) has another layout. Now generic: sweep CSVs are read as "varied settings, seed, [masks], run, results"
  (columns before `run` = settings, after = metrics, `seconds` ignored); tracking sweeps keep their named metrics and the masks switch. Not caused by the UI speed-up.
- Players (app/player.py): every in-browser movie now shows 'Loading frames… n / N' with a progress bar and a disabled button, then 'Ready: press Play'. Chromatin player: the grey circle is the simulation's confining sphere (the territory, NOT the nucleus); caption and in-image label say so. The territory is ~1-1.5 um across at 50 nm/unit vs a 12 um nucleus.

- Display scaling (2026-10-11): single-cell and all-cells players scale each locus colour once for the whole movie from the cell's own mask region (fastview.movie_limits); before, each frame was scaled to the crop maximum, so cosmic rays and fiducial beads in the crop corners made the real spots look missing (spot < 30% of crop max in 70-92% of frames in 3D runs). Measured locus D = MSD/(4t) at 10 s: library loci 0.0017-0.0019 um2/s (falls to 0.0011 at 100 s, nm_per_unit 50 and block_duration_s 10 are PLACEHOLDERS), toy loci 0.0029; the nucleus itself moves 0.43 um per frame.

## Dashboard background (2026-10-16)
- The page background is the Bohrer Lab website's faint chromatin field (user's choice among the site's backgrounds; source: the `.fieldbg` canvas of https://chb-bohrer.github.io, repo CHB-Bohrer.github.io).
  `app/bgfield.py` redraws the same recipe with PIL (seed 67, 1920x1080, transparent) into `app/static/bgfield.png` (committed, 450 KB); `style.CSS` puts it behind `.stApp` with `background-attachment: fixed`;
  `.streamlit/config.toml` has `[server] enableStaticServing = true` (serves app/static at /app/static/). Redraw: `scripts\run.cmd python scripts\make_bgfield.py [--seed --width --height]`. The dashboard colours (navy, teal)
  were NOT changed to the site's near-black/red palette. Restart the dashboard to see it (config change). Checked on port 8502: image served (200), no console errors, body text stays legible.

- Dashboard plot clean-up (2026-10-16; user: "anywhere else this could use improvement? ... fix all"). Method: render every page headlessly with Streamlit AppTest on a real run and save each figure by hooking
  `plt.close` (script kept out of the repo; save with `bbox_inches="tight"` like `style.fast_pyplot` does, otherwise titles/labels look cut off when they are not). Fixed: (1) `app/sweep_plot.py` replaces the one-line-per-combination
  sweep plot on 📈 Scans and the 🔗 E–P scans: each setting has ONE role (x axis / colour of lines / separate panels / held at a value that EXISTS in the sweep, sweeps are often not full grids), defaults by `default_roles`,
  widgets by `controls` (keys include the scan name), legend outside; all 60 (sweep x x-axis) combinations in $SIMLIVE_DATA/sweeps draw; (2) `style.fit_titles` wraps axes titles wider than their axes (called by `fast_pyplot` and
  `cached_figure`); (3) MS2 → Pol II legends above the panels, bigger heatmap numbers; (4) 🔬 Cells locus labels no longer overlap; (5) flat-field maps say "fit rejected" instead of an unexplained flat colour; (6) legends at readable
  size on Transcription/Probes. Tests: `tests/test_sweep_plot.py`. NOT looked at: sections behind toggles (Images build-up, players) and pages that draw no matplotlib figure (Metrics, Stage & fiducials, Library, Movie).
- Chromatin tab "Loop-extruding factors over time" (2026-10-16): the old plot (up to 60 factors x 2 legs x 8000 blocks as thin lines) was unreadable. Now `analysis.loop_coverage` (loops covering each 20 kb bin per
  time bin) as a heatmap + mean loop size (with a minutes axis from the assumed calibration and the boundary sites marked), and a zoom (`analysis.loop_traces`: genomic window x time window, at most 16 factors as shaded bands between the
  two legs, lines broken where a factor reloads; default 400 kb x 300 blocks, ~12-19 factors there). Tested in `tests/test_loop_view.py`; viewed on a real 8000-block library simulation. Streamlit does not
  hot-reload `app/chromatin_tab.py`: restart the dashboard to see edits.
