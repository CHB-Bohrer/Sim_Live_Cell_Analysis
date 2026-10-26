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

## How it is made (to update it when the dashboard changes)
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

## The MP4 (49 min; 720p = 88 MB, full size 1600x900 = 126 MB)
Download it from the release: https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis/releases/tag/tutorial-video-2026-10-25
This version (2026-10-25) uses the Kokoro af_heart clips (normal speed; the first 2026-10-25 upload used Microsoft Jenny and was replaced the same day), has the Cells side-view sentences and the 2026-10-25 screenshots. The 2026-10-20 one (52 min, browser-grade voice, no side views) stays at https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis/releases/tag/tutorial-video-2026-10-20. The previous one (48 min, without them, voice ~0.7 s ahead of the picture)
stays at https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis/releases/tag/tutorial-video-2026-10-17.
It is NOT in the git tree, to keep the repository small (an earlier 46-minute version, `docs/tutorial/tutorial_720p.mp4`, is still in the git history; it was removed from the tree on 2026-10-17).
Local copies: `C:\Users\cbohr\SimLiveData\tutorial\video_2026-10-25\` (and the older ones in `C:\Users\cbohr\SimLiveData\tutorial\`).
The same tutorial as a video with narration. Voice = Microsoft's online neural voice `en-US-JennyNeural` via the `edge-tts` package (NOT Claude's voice; there is none to call).
How it was made (scripts in `scripts/tutorial/video/`, run with the same separate venv as the screenshots plus `pip install edge-tts imageio-ffmpeg numpy scipy`):
1. `tts.py`: exports every sentence of the tutorial and synthesizes one MP3 each (sends the narration text, already public in this repo, to Microsoft's service); writes `beats.json` with durations.
2. `record.py`: plays `index.html#s=welcome&b=0&play=1` silently with every sentence lasting exactly its audio length (`window.BEAT_DUR`), screen-records it with Playwright and logs when each sentence started.
   Playwright's video recorder needs an ffmpeg at `%LOCALAPPDATA%\ms-playwright\ffmpeg-1013\ffmpeg-win64.exe`; the `imageio-ffmpeg` binary was copied there.
3. `mix.py`: places each MP3 where that sentence's caption appears in the recording, then encodes H.264 + AAC. The recorder's frames arrive 0.3-2.2 s
   (median 0.8 s) after the page logs a sentence, by a varying amount, so audio placed at the logged time led the picture (the 2026-10-17 video has this);
   `--no-sync` restores the old placement. Checked on the 2026-10-20 video at 49 sentences spread over the whole length: the voice starts 0.17-0.63 s after its caption
   (median 0.22 s, mostly the silence at the start of each voice clip).
Takes about 55 minutes in real time. Re-do it whenever the narration or the dashboard changes (delete `vid/audio` first: `tts.py` caches clips by scene and sentence number, not by text).
