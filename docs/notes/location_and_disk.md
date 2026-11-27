# Where the project lives on disk; disk-space rules

Working notes moved verbatim from `CLAUDE.md` on 2026-11-27 (section order kept). Add new dated sections at the top; delete statements that are no longer true.

## Repo moved OUT of Dropbox (2026-11-27; user: "move only this code folder out of Dropbox to C:\Users\cbohr\Projects\Sim_Live_Cell_Analysis ... simulation_data stays in Dropbox")
- The git repo is now `C:\Users\cbohr\Projects\Sim_Live_Cell_Analysis` (why: Dropbox syncing `.git` can corrupt it; GitHub is the backup).
  Made with `Dropbox\project-kit\move_repo_out_of_dropbox.py` (`--check` first; it clones from GitHub, so everything must be pushed),
  then the clone `Projects\code` was renamed to `Sim_Live_Cell_Analysis`.
- The data did NOT move: `C:\Users\cbohr\Dropbox\Live_Cell_Microscopy_Investigations\Sim_Live_Cell_Analysis\simulation_data`. It is no longer
  next to the repo, so `scripts\run.cmd`, `simlive.io.runs.DATA_ROOT` and `scripts/tutorial/update_tutorial.py` search in this order:
  `$SIMLIVE_DATA`; `simulation_data` next to the repo folder; the Dropbox path above (via `%USERPROFILE%` / `Path.home()`); `~/SimLiveData`; `<repo>/data`.
- Fixed: the `.pth` and `direct_url.json` of both conda envs (`simlive`, `simlive-frozen-2026-10-05`) point to `Projects\Sim_Live_Cell_Analysis\src`;
  the Claude memory folder copied to `.claude\projects\C--Users-cbohr-Projects-Sim-Live-Cell-Analysis`; doc paths. Copied by hand (gitignored):
  `docs/tutorial/audio` (252 Kokoro clips), `app/static/tutorial`, and the empty `sheets2` folder. `.claude/launch.json` is tracked, so the clone has it.
- The old Dropbox `code\` folder was kept (the script's rename to "code (old repo - delete after checking)" failed with "Access denied":
  an open session held the folder). Its contents match commit 9d91ef2 plus the gitignored files above.
- The Dropbox `Live_Cell_Microscopy_Investigations\Sim_Live_Cell_Analysis\README.md` (the folder map, not in git) was updated to say so.

## Project location (moved 2026-11-21 into Dropbox, then 2026-11-25 one level down; the repo left Dropbox 2026-11-27, see above)
- 2026-11-25, user: "I want you to make a larger folder call Live_Cell_Microscopy investigations and store everything for this project in
  there". The whole `Dropbox\Sim_Live_Cell_Analysis` folder (code, simulation_data, README.md) was moved with one rename into
  `Dropbox\Live_Cell_Microscopy_Investigations\` (underscore, no space). Fixed afterwards: the SimLiveData junction, the `.pth` and
  `direct_url.json` of both conda envs, the doc paths, the Claude memory folder (copied). No code change was needed: `run.cmd` and
  `DATA_ROOT` find `simulation_data` relative to the repo. 229 CPU tests passed from the new place. Gotcha: the first move attempt failed
  with "Access denied" because Windows Explorer held a handle on `simulation_data\tutorial\voice_samples`; restarting Explorer fixed it.
- 2026-11-21, user: "I want all of this project stored within a folder in my dropbox".
- Until 2026-11-27 everything was in `C:\Users\cbohr\Dropbox\Live_Cell_Microscopy_Investigations\Sim_Live_Cell_Analysis\` (Dropbox syncs it; the user chose to include all ~115 GB of data):
  `code\` = the git repo (now in `C:\Users\cbohr\Projects`), `simulation_data\` (still there) = what used to be `C:\Users\cbohr\SimLiveData` (runs, sweeps, chromatin library, tutorial shots
  and videos), `README.md` there = a map of the folder (not in git). The old OneDrive copy (`OneDrive\Documents\Sim_Live_Cell_Analysis`) was
  sent to the Recycle Bin after the move was verified (225 CPU tests passed from the new place).
- `C:\Users\cbohr\SimLiveData` is now a JUNCTION to `simulation_data` so absolute paths saved in old run files keep working. Do not delete
  the junction's target through it.
- The `simlive` conda env (and `simlive-frozen-2026-10-05`) import the code through an editable install: their
  `site-packages\__editable__.simlive-0.1.0.pth` now holds `C:\Users\cbohr\Projects\Sim_Live_Cell_Analysis\src`. If the repo moves again, edit that
  line (or `pip install -e . --no-deps`), otherwise the env silently imports nothing / the old copy.
- NOT in Dropbox on purpose (path-bound, rebuildable): conda envs in `C:\Users\cbohr\miniforge3`, the tutorial venvs `C:\Users\cbohr\shotenv_tmp`
  and `C:\Users\cbohr\kokoro_tmp` (moving a venv breaks it).
- Claude Code: open sessions in `C:\Users\cbohr\Projects\Sim_Live_Cell_Analysis` (the auto-memory was copied to that project's memory folder on 2026-11-27).

## Disk space: compressed images, fewer copies (2026-11-21; user: "do you really need to save that many? ... cut it down ... use some logic")
- Measured before: 115 GB / 42,000 files, of which runs 96 GB (stage3 movies 40, stage5 per-cell crops 37, stage7 CTC copies 9.2 in 26,000
  files, stage4 masks 9.2). Most of it was copies of the same pixels or uncompressed masks.
- **Every TIFF the pipeline writes is now losslessly compressed** (`simlive.io.images.write_tiff`, zlib + predictor; bit-for-bit identical):
  label masks ~150x smaller, camera movies ~2.5x; one frame still reads in a few ms. Compressed files cannot be memory-mapped, so every
  former `tifffile.memmap` (isolate, bleach, flatfield, fiducials, biplane, localize3d) uses `io.images.open_stack`: memmap for old
  uncompressed files, else `PageStack` (decodes only the indexed frames / planes; picklable). Use `write_tiff` / `open_stack` in new code.
  Tests: `tests/test_images_io.py`. Old runs: `scripts\run.cmd python scripts\compress_tiffs.py [folder] [--dry-run]` (verifies each file).
- **CTC folders are no longer kept**: `tracking_demo` writes them to `stage7_validation/_ctc_scratch`, scores, and deletes them; only the track
  tables stay (`ctc_tracks_gt.txt`, `ctc_tracks_pred_<variant>.txt`). The scores were always in metrics.json.
- **Deleted with the user's OK (2026-11-21):** the run folders of 163 parameter-scan members (each scan's `results.csv` keeps every number; any run is
  remade from config + seed; runs named in docs / tutorial were kept), all CTC copies, `scratch_runs`, the tutorial video working folders
  (all three video versions are on GitHub Releases; only the 2026-10-17 MP4s stay locally). The 14 GB chromatin library was NOT touched.
- **Result:** `simulation_data` 115 GB / 42,415 files -> 28.1 GB / 9,178 files (runs 96 -> 13.5 GB, 53 runs kept; chromatin 14.4 GB unchanged).
  The 2,780 TIFFs of the kept runs went 55.2 -> 14.4 GB, every one verified identical after rewriting. Gotcha: Dropbox locks files while it
  uploads them, so `os.replace` can fail with WinError 5; `compress_tiffs.py` retries and keeps the original; just run it again.
- Scans whose member runs are gone: the 📈 Scans page reads only `sweeps/<name>/results.csv`, so it is unaffected; a run listed there can be remade
  with `scripts\sweep.py` (same grid + seeds -> same run ids). Exception: `scripts/sweep_copies_strays.py --reanalyse TAG` re-reads the scan's saved
  movies, which are gone: re-run the scan first (it re-makes them).
