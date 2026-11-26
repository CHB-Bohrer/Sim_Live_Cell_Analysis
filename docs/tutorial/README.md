# Narrated tutorial

`index.html` is a self-contained, animated, narrated tour of the whole program and its dashboard (about 50 minutes at 1x, about 25 minutes at the default 2x speed, in 38 short scenes;
every scene can be jumped to from the list on the left). **Open it from the dashboard's home page ("Open the narrated tutorial") or by
double-clicking `index.html`** (Microsoft Edge recommended). Nothing to install, nothing is sent anywhere. The dashboard serves a copy
(`app/static/tutorial`, refreshed from this folder whenever the home page loads; `app/style.publish_tutorial`), so edit the files HERE.

- **Voice:** the default is a recorded voice made with Kokoro (a free text-to-speech model that runs locally; voice `af_heart`): one MP3 per sentence in `audio/`, played by the page itself, so it works offline once the files
  exist. They are made by `scripts/tutorial/video/make_clips_kokoro.py` (see below) and are NOT in git (17 MB): without them the page falls back to your browser's speech synthesis (Web Speech API;
  Edge offers the most natural voices, names ending in "Online (Natural)"). Pick the voice and speed on the start screen or in the bottom bar; the default speed is 2x.
  A clip is used only while its sentence is unchanged (hash check), so after editing narration run the two steps below again (only edited sentences are re-recorded).
  This is NOT Claude's own voice: there is no way to call it from a local page.
- **Keys:** Space = play/pause, Left/Right = previous/next sentence, Shift + Left/Right = previous/next scene, M = voice on/off, C = captions.
- **Deep links:** `index.html#s=<scene id>&b=<sentence number>` opens a scene (ids are in `scenes_a.js` / `scenes_b.js`).
- **Your turn:** the last sentence of most scenes shows a short "Your turn" card with steps to try in the real dashboard.

## What it contains
| Part | Scenes |
|---|---|
| 0 The idea | what the program is for; the seven pipeline stages (animated) |
| 1 Starting the dashboard | launch (terminal), home screen, built-in map and glossary |
| 2 Making a simulation | the two sidebar forms (every setting), what happens when you press Run |
| 3 Cells and tracking | animation of tracking errors; Movie, Metrics, Trajectories, Scans pages |
| 4 The microscope | animation of a spot becoming an image; Images, Calibration, Stage and fiducials pages |
| 5 Chromatin | polymer + loop extrusion animation; Library, Chromatin pages |
| 6 Finding the loci | Cells page (incl. the x-z / y-z side views of 3D runs); animation of the Viterbi spot finder; animation of two promoter-enhancer copies, stray spots and the multi-copy linker; Spot finding page |
| 7 Transcription biology | probes, transcription, MS2 deconvolution and E-P curve animations; Probes, Transcription, MS2, E-P pages |
| 8 Configuration | Config page; command-line cheat sheet |
| 9 Putting it together | three workflows, cautions, cheat sheet |

## Updating it: one command
`scripts\tutorial\update.cmd` says which screenshots are out of date and checks every sentence; `scripts\tutorial\update.cmd all` re-takes the
stale screenshots, records the voice of edited sentences and checks again. The step-by-step guide for every kind of change is
**[UPDATING.md](UPDATING.md)**; [NUMBERS.md](NUMBERS.md) lists the sentences that state a number (re-read them when results change).

## How it is made (background; the command above does these steps)
- Screenshots are REAL dashboard pages captured headlessly with Edge + Playwright: `scripts/tutorial/capture_shots.py` (needs a dashboard on
  port 8502 and a Python with `playwright`, ideally a separate venv: `python -m venv C:\Users\<you>\shotenv_tmp\v` then `pip install playwright`;
  it drives the installed Edge, no browser download). It writes tall page images plus the on-screen rectangle of every control.
  Then `python scripts/tutorial/build_data.py <shots folder>` copies them to `img/` and writes `shots_data.js`. The scenes point at controls by name
  (e.g. `f("next", "Next error")`), so a layout change only needs a re-capture; a renamed or new control needs an edit of `scenes_*.js`.
- `engine.js`: player, speech, the screenshot "camera" (`Shot`), canvas helpers, terminal and card components. `anims.js`: the concept animations.
  `scenes_a.js`, `scenes_b.js`: the narration (one sentence per beat) and which visual cue goes with it. No libraries, works from `file://`.
- Keep each sentence under ~250 characters: some online voices stop reading long utterances.
- Run numbers in the screenshots come from `loci_demo_seed11_24bba2`, `tracking_demo_seed3_e15523`, `loci_demo_seed8_f5c88d`, `loci_3d_seed23_2163f4` and
  (two copies per cell, stray spots: shots `pg_cells_copies`, `pg_spots_copies`) `loci_library_seed12_87a16d`.

## Recording the voice clips (Kokoro, runs on this PC)
Two steps, two virtual environments (the setup of the second is in the header of `scripts/tutorial/video/make_clips_kokoro.py`; nothing is sent over the internet):
1. shot venv (playwright): `python scripts/tutorial/video/make_clips.py --export C:/Users/<you>/kokoro_tmp/beats.json`
2. kokoro venv: `python scripts/tutorial/video/make_clips_kokoro.py C:/Users/<you>/kokoro_tmp/beats.json af_heart` (other voices: af_bella, am_michael, bf_emma; optional third argument = speed). About 10 minutes for all 250 sentences, much less after an edit.
`make_clips.py` without `--export` is the older route with Microsoft's online voices (edge-tts, sends the text to Microsoft). To share the clips, attach `docs/tutorial/audio` to a release.

## The MP4 (49 min; 720p = 98 MB, full size 1600x900 = 128 MB)
Download it from the release: https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis/releases/tag/tutorial-video-2026-11-26
This version (2026-11-26) matches the interactive tutorial as of commit 57073b5: the recorded Kokoro af_heart clips at normal speed (1x; only the web player defaults to 2x),
the redrawn transcription scene (2026-10-26), the stray-spot wording fix and the 2026-11-10 screenshots of the re-run results. Earlier releases stay:
https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis/releases/tag/tutorial-video-2026-10-25 (49 min, Kokoro, older screenshots),
https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis/releases/tag/tutorial-video-2026-10-20 (52 min, Microsoft Jenny voice) and
https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis/releases/tag/tutorial-video-2026-10-17 (48 min, voice ~0.7 s ahead of the picture).
It is NOT in the git tree, to keep the repository small (an earlier 46-minute version, `docs/tutorial/tutorial_720p.mp4`, is still in the git history; it was removed from the tree on 2026-10-17).
The releases are the copies; the 2026-11-26 working folder (`simulation_data\tutorial\video_2026-11-26\`, with `vid/` and `release/`) is kept locally for now.
Voice = the same recorded Kokoro clips the player uses (`docs/tutorial/audio`, made on this PC, see above), NOT Claude's voice (there is none to call).
How it was made (scripts in `scripts/tutorial/video/`, run from the working folder with the screenshot venv plus `imageio-ffmpeg numpy scipy`):
1. Make sure every sentence has a matching clip (`scripts\tutorial\update.cmd` reports none missing). Then copy each `docs/tutorial/audio/<scene>_<n>.mp3` to `vid/audio/` and
   write `vid/beats.json` = `{"voice", "beats": [{id, text}], "dur": {"scene:n": seconds}}` with the sentence list read from the player (`T.scenes`), as `tts.py` does
   (`beats_from_kokoro.py` does exactly this and also checks each clip's hash against `hashText`). `tts.py` is the older route: it synthesizes
   the clips with Microsoft's online voice (edge-tts; sends the text to Microsoft) and caches them by scene and sentence number, not by text (delete `vid/audio` first).
2. `record.py`: plays `index.html#s=welcome&b=0&play=1` silently with every sentence lasting exactly its audio length (`window.BEAT_DUR`), screen-records it with Playwright and logs when each sentence started.
   Playwright's video recorder needs an ffmpeg at `%LOCALAPPDATA%\ms-playwright\ffmpeg-1013\ffmpeg-win64.exe`; the `imageio-ffmpeg` binary was copied there.
3. `mix.py <out.mp4>`: places each MP3 where that sentence's caption appears in the recording, then encodes H.264 + AAC. The recorder's frames arrive late after the page
   logs a sentence, by a varying amount, so audio placed at the logged time would lead the picture (the 2026-10-17 video has this); `--no-sync` restores that placement.
   Caption lag goes to `vid/caption_lag.json`; 2026-11-26: all 250 found, median 0.10 s, range 0.02-0.30 s.
4. 720p copy: `ffmpeg -i full.mp4 -vf scale=1280:720 -c:v libx264 -crf 27 -pix_fmt yuv420p -c:a copy -movflags +faststart <name>_720p.mp4`.
Takes about an hour in real time (the recording alone ~50 min). Re-do it whenever the narration or the dashboard changes.
