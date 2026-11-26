# Narrated tutorial, voice, playback, video releases

Working notes moved verbatim from `CLAUDE.md` on 2026-11-27 (section order kept). Add new dated sections at the top; delete statements that are no longer true.

## Narrated tutorial (2026-10-17; user: "in depth tutorial with animation and voice over ... the person should know how to use the UI")
- `docs/tutorial/index.html` (+ `engine.js`, `anims.js`, `scenes_a.js`, `scenes_b.js`, `tutorial.css`, `shots_data.js`, `img/`): self-contained, no libraries, works from `file://`.
  38 scenes / 246 sentences (~50 min; 2026-10-20). The dashboard home page links to it (served from a copy in app/static/tutorial, refreshed from docs/tutorial on every home-page load: edit docs/tutorial only). Narration = the browser's speech synthesis (Edge "Online (Natural)" voices are best; the picker is on the start screen). It is NOT
  Claude's own voice (not callable from a local page); the MP4 (below) uses a Microsoft neural voice instead.
- Screenshots are real dashboard pages: `scripts/tutorial/capture_shots.py` (Playwright driving the installed Edge, in a separate venv `C:\Users\cbohr\shotenv_tmp\v`, NOT the simlive env) captures TALL page images + the
  rectangle of each control; `scripts/tutorial/build_data.py <shots dir>` copies them to `docs/tutorial/img` and writes `shots_data.js`. Scenes refer to controls by name (`f("next", "Next error")`), so a
  layout change needs only a re-capture. Runs shown: `loci_demo_seed11_24bba2` (most pages), `tracking_demo_seed3_e15523` (tracking errors), `loci_demo_seed8_f5c88d` (Images part 2 measured effects), `loci_3d_seed23_2163f4` (Stage page),
  `loci_library_seed12_87a16d` (two copies + stray spots: shots `pg_cells_copies`, `pg_spots_copies`).
  The capture needs the test dashboard on port 8502 and takes ~10 min; it does not use the GPU. `python capture_shots.py <out> <shot> ...` re-takes single shots.
  Shots of 2026-10-20 in `C:\Users\cbohr\Dropbox\Live_Cell_Microscopy_Investigations\Sim_Live_Cell_Analysis\simulation_data\tutorial\shots_2026-10-20b`. Hand-set rectangles live in `build_data.py` `OVERRIDES` (used where the label text
  also appears higher up, e.g. a dataframe column called "simulation"); check them after a re-capture by drawing them onto the shot (they went stale once when a table was added).
  Check every sentence for script errors / missing rectangles with a Playwright loop over `T.scenes` (each beat: open `index.html#s=<id>&b=<n>`, collect console errors and
  the engine's "no rect" / "missing shot" warnings).
- Gotchas hit: (1) the in-app browser pane gives only 800x500 screenshots and shows `file://` pages as static snapshots (no scripts), so test the tutorial with Playwright (`tut.py`-style: open `index.html#s=<id>&b=<n>`
  and screenshot); (2) Windows long paths: the venv must be at a short path; (3) in Streamlit, `get_by_text` finds the intro boxes before the real widget (e.g. "Next error" appears in the intro text), so rectangles are
  found as the smallest visible container whose text contains the label; (4) a tall (9000 px) viewport stretches the fixed page background, so the capture injects `background-image:none`; (5) text in page JS strings
  written by a Python patch script loses `\n` escapes: use the Edit tool; (6) 2026-10-20: switching runs raced with the page capture (the old run's page was saved; now `select_run` waits for the new run),
  pages taller than the 9000 px viewport were cut (copies shots use 15000 px; `save_page` reads the real viewport height), and clicking an already-open page in the segmented control DESELECTS it
  (`goto_page` now skips when `aria-checked == "true"`); (7) Streamlit keeps the PREVIOUS page's elements, dimmed, until the new page replaces them, so a long page
  saved too early ends in the old page's leftovers (the first Spot finding shot ended in Trajectories plots, parts 3d / 4 were missing): `save_page` now waits until no
  `[data-stale="true"]` element is left. A scene with several pictures (`shot({ sets })`) stays on the last one chosen: switch back explicitly with `fs("main", ...)`.
- MP4s are GitHub releases, NEVER in the git tree (user 2026-10-17: "don't add another 59 MB to the history"; the first 59 MB version is still in history). CURRENT: release `tutorial-video-2026-11-26` (https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis/releases/tag/tutorial-video-2026-11-26;
  49 min, Kokoro af_heart clips at 1x, matches the tutorial at commit 57073b5; assets 128 MB + 98 MB 720p; published 2026-11-26 at the user's request; README, docs/tutorial/README, the dashboard banner
  (`app/style.py` TUTORIAL_VIDEO) and status.json point to it). Made from the recorded clips, not tts.py: `scripts/tutorial/video/beats_from_kokoro.py` (builds vid/audio + vid/beats.json,
  checks clip hashes), then record.py, mix.py, ffmpeg scale=1280:720 crf 27. Local working folder KEPT: `simulation_data\tutorial\video_2026-11-26\` (vid/, release/, stills/; ~1 GB, mostly the raw webm).
  Older: release `tutorial-video-2026-10-20`
  (https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis/releases/tag/tutorial-video-2026-10-20 (SUPERSEDED by tutorial-video-2026-10-25, then 2026-11-26); 52 min, with the copies scene; assets `sim_live_cell_analysis_tutorial.mp4` 110 MB and `_720p.mp4` 82 MB;
  published 2026-10-20 at the user's request; README links point to it). Previous: `tutorial-video-2026-10-17` (48 min, kept). Local copies (video_2026-10-20/-25/-25b folders) were DELETED 2026-11-21 to save space: the releases are the copies. Formerly in `simulation_data\tutorial\video_2026-10-20\`
  (+ `release/` with the asset names) and older ones in `C:\Users\cbohr\Dropbox\Live_Cell_Microscopy_Investigations\Sim_Live_Cell_Analysis\simulation_data\tutorial\`. Publishing a new release is outward-facing: ask the user first. `gh release create --target` needs a branch name or full SHA.
  Delete the audio cache before re-running `tts.py` (clips are keyed by scene:sentence, not by text). Found while checking the 2026-10-20 video: the recorder's frames lag the page's sentence log by 0.3-2.2 s, varying, so the voice LED the picture by ~0.7 s (also in the 2026-10-17
  video). `mix.py` now finds each caption change in the recording and places the audio there (voice now 0.17-0.63 s after its caption, median 0.22 s, 49 sentences checked).
  The player also pre-decodes the current and next scene's screenshots (`engine.js` `decodedImg` / `keepDecoded`) so switching to a 9000 px picture does not stall it.
  Voice = Microsoft `en-US-JennyNeural` via `edge-tts` (the user chose it; Claude's own voice cannot be produced). Scripts `scripts/tutorial/video/{tts,record,mix}.py`; the player accepts `#s=<id>&b=<n>&play=1`
  and `window.BEAT_DUR` (fixed sentence lengths) for recording. Regenerating takes ~50 min of real time. Playwright's recorder needed an ffmpeg copied to `%LOCALAPPDATA%\ms-playwright\ffmpeg-1013\ffmpeg-win64.exe`.
- When the dashboard changes: re-run the capture + build, then fix any scene whose control was renamed. Sentences must stay under ~250 characters (some online voices cut off long utterances). NOT checked: how the
  voices actually sound (no audio output was available to the author), and the tutorial on Chrome/Safari.

## Side views along z, natural tutorial voice, bug sweep (2026-10-25; user: "visualization of the spot along the X-Z axis in addition to X-Y", "update the tutorial in full", "default speed 1.5x", "still sounds like a robot")
- **Side views (Cells page, `app/cells_tab.side_views`)**: for z-stack runs, below the x-y figure of "One cell through time": x-z and y-z cuts of the chosen spot's own colour stack
  (`<ch>_zstack.tif`, frame of the slider; 3 pixel rows/columns averaged, zoomed 3 um sideways and 4 um up/down, each plane minus its own median background so the spot is not lost in
  the fixed-pattern noise), with found (white +, from `loci_positions_3d[_cal].csv`; radio when both exist), true (cyan o) and drawn (red x) positions, plus the height of every spot
  through time with the true heights as thick faint lines. Works for true and tracked IDs, with copies; plane heights from `imager_for_run(...).zk` (biplane too). 2D runs show a one-line notice.
  Tested on the 3D, biplane and two-copy 3D runs (`tests/test_dashboard_pages.py::test_cells_side_views`). NOT done: x-z view in the Spot finding page's per-step images and in the player.
- **Bugs found and fixed:** Cells page crashed on runs made before the `visible` column of `localization_errors.csv` (loci_demo_seed1_10849d); the tutorial said 26 image-error sources (there are 28; README said 22);
  `build_data.py` dropped every shot not in the folder it was given (it now merges into `shots_data.js`, so single shots can be re-taken).
- **Tutorial voice:** the browser voices sound robotic, so the player plays recorded clips (`docs/tutorial/audio/*.mp3` + `manifest.js`). FIRST version: Microsoft edge-tts Jenny (`make_clips.py`), judged still
  too synthetic by the user; CURRENT: Kokoro-82M, voice `af_heart` (user's pick from 4 samples), run locally on the CPU with NO text leaving the PC: `make_clips.py --export beats.json` (shot venv) then
  `make_clips_kokoro.py beats.json` (kokoro venv `C:/Users/cbohr/kokoro_tmp/v`, model files in `kokoro_tmp/model`, setup + SHA-256s in the script header; downloaded 2026-10-25 with the user's OK, 350 MB). 250 sentences, 22 MB.
  `PRON` in the script respells names the model misreads (Viterbi, Cellpose ...; spoken only, captions unchanged): extend it if a word sounds wrong. A clip is used only if its hash matches the sentence in `engine.js`
  (`hashText`), otherwise the browser voice speaks it, so after editing narration re-run `make_clips.py` (only changed sentences are synthesized). `docs/tutorial/audio/` is gitignored (like the MP4, not in
  the history): a fresh clone has the browser voice until the clips are made or published. Voice list has "Natural voice (recorded)" first. Default speed 1.5x then (2x since 2026-10-26, see "Tutorial playback" below).
  The MP4 was re-made the same day with this voice, the side-view sentences and the 28-source fix: release `tutorial-video-2026-10-25` (49 min, Kokoro af_heart voice; 126 MB + 88 MB 720p; first uploaded with the Jenny voice, assets replaced the same day; published at the user's request; the README, the dashboard home link and docs/tutorial/README pointed to it until 2026-11-26, when `tutorial-video-2026-11-26` replaced it as current; both earlier releases stay). Local files: `C:\Users\cbohr\Dropbox\Live_Cell_Microscopy_Investigations\Sim_Live_Cell_Analysis\simulation_data\tutorial\video_2026-10-25\`. The video narrates at normal speed (1x); only the web player defaults to 1.5x.
- New shot `pg_cells_3d` (capture_shots.py `cells_3d`, run `loci_3d_seed23_2163f4`); scene "cells" has 4 new sentences (250 sentences in 38 scenes now).

## Tutorial playback (2026-10-26; user: "the movie consistently lags ... the animations tend to go and then stop", "default speed ... 2x")
- **Cause found:** the PC itself. Windows' compositor (`dwm.exe`) kept the GPU's 3D engine at 81-97% with nothing else running (one 3440x1440 60 Hz screen;
  uptime 18 days). A Chrome trace showed the page's own work at ~1% per frame, but `DXGISwapChainImageBacking::Present` holding single frames 0.3-0.8 s
  (~18 frames/s). A blank page with one moving circle stalled the same way, so the tutorial code was not the cause. Check with
  `Get-Counter "\GPU Engine(*)\Utilization Percentage"` (dwm should be a few %); a sign-out or reboot usually clears it. We do NOT change system settings.
- **What the page did wrong:** its animation clock cut every frame's time step to 50 ms, so with late frames the canvas animations ran at ~40% of real
  time (measured: 1.6 s of animation per 4 s), fell behind the voice and looked frozen. Now `engine.js` passes the real elapsed time (capped at 1 s) and
  `Canvas` sub-steps its physics (`step()`) in pieces of <= 50 ms (`STEP`) and paints once: 3.8 s per 4 s at the same frame rate. Draw functions only use dt for
  tweens / accumulators, which are fine with a large dt.
- The next sentence's clip is now preloaded while the current one plays (`preloadAfter`, `takeClip`; the dashboard's static server takes ~0.25 s per file);
  measured readyState 4 (fully loaded) at every play.
- **Default speed 2x** (options now up to 2.5x; storage key `rate20`, so a 1.5x saved under `rate15` does not override the new default).
- **Scene "Distance, promoter, Pol II, MS2" (`A.transcription`) fixed** (user: "the graphs don't show anything"): its seed search found no example
  and silently fell back to a flat distance with no bursts. Now 40 min, distance = smooth random walk (90 s, 750 +- 280 nm, ~10% in contact), the first
  seed with 3-5 bursts (>= 2 starting in contact, >= 1 far, every burst >= 2 Pol II); faint full traces always drawn, contact shaded in every panel,
  burst / Pol II counters, axes in nm / loops / minutes; the sweep restarts when a panel appears and lasts 26 s / speed. Narration unchanged (clips still match).
- Test the tutorial with Playwright from the shot venv (`C:\Users\cbohr\shotenv_tmp\v`) over http (`python -m http.server` in docs/tutorial) with
  `--autoplay-policy=no-user-gesture-required`; `headless=False` uses the real GPU. To measure the animation clock: `T.onFrame(d => sum += d)` vs wall time.
