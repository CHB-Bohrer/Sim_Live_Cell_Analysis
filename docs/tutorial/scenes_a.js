/* Scenes, part 1: the idea, starting the dashboard, making a simulation, cells and tracking, the microscope. */
"use strict";
const SC = (id, title, make, beats, tryit) => T.add({ id, title, make, beats: beats.map((b) => ({ s: b[0], c: b[1], pause: b[2] })), tryit });
const shot = (o) => (st) => new Shot(st, o);
const code = (s) => `<code>${s}</code>`;
/* focus helper: f("name", "label", opts) */
const f = (n, l, o) => (s) => s.focus(n, l, o);
const sec = (a, b, label, extra) => (s) => { const r = s.between(a, b); s.focus(r, label, Object.assign({ maxZoom: 2.2, padX: 1.35, padY: 1.1 }, extra)); };

/* ================================================================== 0 */
T.chapter("0 · The idea");
SC("welcome", "What this program is for", (st) => A.truthloop(st), [
  ["Welcome. This tutorial shows how the Sim Live Cell Analysis program works, and how to drive its dashboard. By the end you will be able to start it, make your own simulated movie, and read every page.", (a) => a.step(0)],
  ["The idea is simple. In a real microscope movie nobody knows the right answer. So we invent the movie. We simulate cells, chromatin, a microscope and a camera, and we keep an answer key with every true position.", (a) => a.step(1)],
  ["Then we pretend we do not know the answer. The analysis only sees the noisy movie. It finds the cells, follows them through time, and locates the promoter and the enhancer inside each one.", (a) => a.step(2)],
  ["It saves its estimates. These are what you would get from real data.", (a) => a.step(3)],
  ["Finally we compare the estimates with the key. The difference is a measured error, not a guess. And because every ingredient of the simulation can be switched on and off, we can learn exactly which microscope error costs how much accuracy.", (a) => a.step(4)],
  ["One honest warning. Many numbers in the program are placeholders until you supply measurements, for example a real mean squared displacement. Treat the results as a demonstration of the method, not as biology.", (a) => a.step(4)],
]);

SC("pipeline", "The pipeline in seven stages", (st) => A.pipeline(st), [
  ["The program is a pipeline of seven stages. Each stage is its own folder of code, and each one writes files that the next one reads, so you can inspect or replace any stage on its own.", (a) => a.stage(0)],
  ["Stage one is chromatin. A polymer simulation on the graphics card folds a ten megabase region and moves it, with loop extrusion and boundary elements. It is made once, saved, and reused. A simple stand-in motion exists for quick tests.", (a) => a.stage(1)],
  ["Stage two is the cells. Nuclei move around the field, bump into each other, change shape, and sometimes divide. Every cell keeps its own true ID number.", (a) => a.stage(2)],
  ["Stage three is the microscope. It turns cells and loci into camera images, with twenty-six switchable sources of image error, from drift and blur to hot pixels and bleaching. It also writes the answer key.", (a) => a.stage(3)],
  ["Stage four finds every nucleus in every frame, which is segmentation, and links them through time, which is tracking. It uses Cellpose and Trackastra, and it never sees the answer key.", (a) => a.stage(4)],
  ["Stage five cuts out one small movie per cell, so each cell can be analysed by itself, many in parallel.", (a) => a.stage(5)],
  ["Stage six finds the promoter and the enhancer in every frame, reads the MS2 spot, and turns it into Pol two start times. This is where your own analysis plugs in.", (a) => a.stage(6)],
  ["Stage seven compares everything with the answer key: tracking scores, locus errors, distance errors, and how much each microscope error costs.", (a) => a.stage(7)],
  ["The answer key travels straight to stage seven and skips the analysis. That is what makes the validation honest.", (a) => a.stage(0)],
]);

/* ================================================================== 1 */
T.chapter("1 · Starting the dashboard");
SC("launch", "Start it", (st) => {
  const o = { instant: false, mode: "term", st };
  o.term = new Term(st); o.term.instant = false;
  o.browser = null;
  o.destroy = () => { o.term.destroy(); if (o.browser) o.browser.destroy(); };
  o.run = async () => {
    o.term.instant = o.instant; o.term.clear();
    await o.term.type("scripts\\dashboard.cmd");
    await o.term.print([{ html: "Starting dashboard on http://127.0.0.1:8501 ... your browser will open in about 10 seconds.", cls: "dim" }, { html: "KEEP THIS WINDOW OPEN. Close it to stop the dashboard.", cls: "warn" }, "", { html: "  You can now view your Streamlit app in your browser.", cls: "ok" }, { html: "  Local URL: http://127.0.0.1:8501", cls: "ok" }]);
  };
  o.showBrowser = () => { o.term.el.style.display = "none"; if (!o.browser) { o.browser = new Shot(st, { shot: "home", label: "🏠 home" }); o.browser.instant = o.instant; o.browser.top(); } };
  return o;
}, [
  ["To start, open a terminal in the project folder and run scripts, backslash, dashboard dot c m d. You can also double-click that file.", (o) => { o.run(); }, 1800],
  ["Your browser opens at localhost, port eighty-five zero one. This address only exists on your own computer, nobody else can reach it.", (o) => { o.showBrowser(); }],
  ["Keep the terminal window open while you work. Closing it, or pressing Control C, stops the dashboard. It also prints any error messages there, so look at it if a page misbehaves.", (o) => { }],
  ["Everything else in the project runs through scripts, backslash, run dot c m d. It activates the correct Python environment for you. Never call the environment's python directly, it fails to load its libraries.", (o) => { }],
  ["The dashboard itself never runs the heavy simulation. When you press a run button it starts a separate process, so the web page stays responsive and never loads the deep-learning libraries.", (o) => { }],
  ["Your turn: start the dashboard now, and keep the page open next to this tutorial.", (o) => { }],
], [`In the project folder run ${code("scripts\\dashboard.cmd")}`, "Wait about ten seconds for your browser to open http://127.0.0.1:8501", "Keep the terminal window open", `To test a second copy on another port: ${code("scripts\\run.cmd python -m streamlit run app\\dashboard.py --server.port 8502")}`]);

SC("layout", "The home screen", shot({ shot: "home", label: "🏠 home" }), [
  ["This is the home screen. On the left is the sidebar. It is where you create new simulations, and where you choose which finished run you are looking at.", (s) => s.focus("sidebar", "Sidebar")],
  ["At the top, the title says which run is open, how many cells it has, how many frames, and the pixel size.", f("hero", "Which run is open")],
  ["Two folded boxes help newcomers: a map of the app, and a glossary of terms. Open them whenever a word is unfamiliar.", f(["guide", "glossary"], "Map and glossary", { multi: true })],
  ["What was simulated in this run lists the settings in plain words: how many nuclei, how long the movie is, the field of view, how many image errors are on, whether there are loci, probe binding, M S two transcription, and which segmentation and tracker were used. Hover a box for details.", f("summary", "Plain-language summary of the run")],
  ["Below are the headline tracking scores. The first row uses automatic segmentation, which is realistic. The second uses perfect masks, which isolates just the linking error. One point zero is perfect for both.", f(["scores_auto", "scores_gt"], "Headline scores", { multi: true })],
  ["Then comes the page bar. There are sixteen pages, and only the page you open is computed, which keeps the dashboard fast. Your choice survives reruns.", f("nav", "16 pages")],
  ["Run to view lists every finished movie, newest first. Pick another one and every page switches to it.", (s) => { s.focus("runsel", "Choose a run"); s.click("runsel"); }],
], ["Open the map and the glossary boxes", "Click through a few pages in the page bar", "Open Run to view and pick a different run; watch the title change"]);

SC("guide", "The built-in map and glossary", shot({ sets: { guide: "guide", glossary: "glossary" }, start: "guide", label: "🧭 New here?" }), [
  ["The map lists every part of the app as a card. The two sidebar forms make new runs. The tracking pages are Movie, Metrics, Trajectories and Scans.", (s) => s.focus([100, 60, 1100, 190], "Make runs · judge tracking")],
  ["Cells shows one small movie per nucleus. Spot finding explains how each locus and the M S two signal are found. Images, Stage and fiducials, and Calibration are about the microscope.", (s) => s.focus([100, 245, 1100, 195], "Cells · Spot finding · microscope pages")],
  ["Transcription and Probes build the models behind the M S two spot and the locus brightness. M S two to Pol two turns the trace into polymerase starts. Library and Chromatin cover the saved polymer simulations.", (s) => s.focus([100, 435, 1100, 130], "Models, deconvolution, chromatin")],
  ["At the bottom is the typical path: choose simulations in Library, tune the models, press Run in the sidebar, and read the results in Movie and Metrics. Every setting also has a small question-mark tooltip.", (s) => s.focus([100, 565, 1100, 50], "Typical path")],
  ["The glossary defines everything from ground truth and ID switch to point spread function, photon transfer curve, kymograph and deconvolution, in one or two lines each.", (s) => { s.use("glossary"); s.focus([95, 60, 1110, 330], "Glossary of terms"); }],
], ["Open the glossary and read three terms you did not know", "Find the page that shows what each image error does (Images)"]);

/* ================================================================== 2 */
T.chapter("2 · Making a simulation");
const NOTE = (h_, items) => `<h4>${h_}</h4><ul>${items.map((i) => "<li>" + i + "</li>").join("")}</ul>`;
SC("form_sim", "Form 1: New simulation", shot({ shot: "sb_sim", align: "left", minScale: 1.2, label: "sidebar · New simulation" }), [
  ["The first sidebar form is called New simulation. It runs the tracking-only experiment: moving nuclei, a basic microscope, segmentation and tracking. It takes about twenty seconds.", (s) => { s.focus(s.between("Random seed", "Acquisition"), "New simulation", { maxZoom: 2.2, padX: 1.3, padY: 1.2 }); s.note(NOTE("Tracking only", ["Defaults come from <code>configs/tracking_demo.yaml</code>", "About 20 seconds", "Everything is placed in one form: nothing runs until you press the button"])); }],
  ["The random seed comes first. The same settings and the same seed always give exactly the same movie.", (s) => { s.focus("Random seed", "Same settings + same seed = same movie"); s.note(null); }],
  ["Acquisition sets how many frames the movie has and the seconds between them. Longer intervals let cells move farther between images, which makes linking harder.", (s) => { sec("Acquisition", "Cells and motion", "Acquisition")(s); s.note(NOTE("Acquisition", ["Frames: how long the movie is", "Frame interval: seconds between images"])); }],
  ["Cells and motion sets the number of cells, the field of view, the nuclear radius, diffusion, directed speed and division. The speed scale is the most important knob for tracking errors: one is the default, one half is half as fast.", (s) => { sec("Cells and motion", "Nuclear shape", "Cells and motion")(s); s.note(NOTE("Cells and motion", ["More cells in the same field = more touching nuclei", "<b>Speed scale</b>: multiplies the cell speed", "Division probability per cell per frame"])); }],
  ["Nuclear shape and texture control how round and lumpy the nuclei are, how fast the shape changes, and the chromatin-like texture and dark nucleoli inside. In our tests texture broke the simple watershed segmenter, while Cellpose coped.", (s) => { sec("Nuclear shape", "Imaging", "Shape and texture")(s); s.note(NOTE("Shape and texture", ["Aspect ratio: 1 is round", "Shape persistence: high = stable shapes", "Texture contrast and nucleoli"])); }],
  ["Imaging has the brightness in photons per pixel per second, the camera's read noise, the pixel size and the numerical aperture. Lower brightness means a noisier picture.", (s) => { sec("Imaging", "Segmentation and tracker", "Imaging")(s); s.note(NOTE("Imaging", ["Photons per pixel per second", "Read noise (electrons)", "Pixel size (nm) and NA set the blur"])); }],
  ["Segmentation and tracker choose the software under test: Cellpose or a classic threshold and watershed, and one of Trackastra's modes: greedy, greedy without divisions, or I L P.", (s) => { sec("Segmentation and tracker", "Image errors", "Software under test")(s); s.note(NOTE("The software being tested", ["Segmentation: finds each nucleus", "Tracker: links nuclei through time"])); }],
  ["Image errors has one checkbox for the realistic set of error sources, and a list to switch individual sources off, so you can see what one error does. Then press Run simulation.", (s) => { s.focus(s.between("Image errors", null), "Then press Run", { maxZoom: 2.2, padX: 1.3, padY: 1.1 }); s.note(NOTE("Image errors", ["Tick: all sources with placeholder values", "Untick = the simple original model", "Switch single sources off to compare"])); }],
], ["Expand New simulation in the sidebar", "Set Speed scale to 1.5 and press Run simulation", "When it finishes, open Movie and Metrics", "Run again with 0.5 and compare the scores"]);

SC("form_loci", "Form 2: Single-cell loci run (basics)", shot({ shot: "sb_loci", align: "left", minScale: 1.2, label: "sidebar · New single-cell loci run" }), [
  ["The second form is the big one: New single-cell loci run. It adds two coloured loci to every nucleus, isolates each cell, and locates the loci in each one. It takes about two minutes, or up to ten with 3D.", (s) => { s.focus(s.between("Random seed", "Frame interval"), "Single-cell loci run", { maxZoom: 2.2, padX: 1.3, padY: 1.2 }); s.note(NOTE("What it adds", ["Two loci per nucleus: promoter and enhancer", "Cells are cut out one by one", "Loci found in every cell, in parallel"])); }],
  ["Cells in the field, frames, and frame interval set the size of the experiment. Make the movie long compared with the burst time if you study transcription: tens of minutes.", (s) => { s.focus(s.between("Cells in the field", "Cell speed scale"), "Size of the experiment", { maxZoom: 2.2, padX: 1.3, padY: 1.2 }); s.note(NOTE("Size", ["Cells in the field", "Frames", "Frame interval (s)"])); }],
  ["Cell speed scale of zero point four keeps tracking easy, so that locus analysis is not drowned by tracking errors. Photons per locus per second set how bright each spot is. Division probability lets cells split.", (s) => { s.focus(s.between("Cell speed scale", "Microscope"), "Speed, brightness, division", { maxZoom: 2.2, padX: 1.3, padY: 1.2 }); s.note(NOTE("Easy first run", ["Speed 0.4 keeps tracking easy", "Brightness is a placeholder until measured"])); }],
  ["Microscope switches between a generic placeholder and the Bruker Vutara V X L in its widefield mode, or in its biplane three-dimensional mode, with its objective, camera, flat illumination and colours taken one after another.", (s) => { s.focus(s.between("Microscope", "Field of view"), "Generic or Vutara VXL", { maxZoom: 2.2, padX: 1.3, padY: 1.2 }); s.note(NOTE("Microscope", ["Generic placeholder", "Vutara VXL: 60x / 1.3 NA, ~100 nm pixels, Orca Fusion BT camera", "Every value is sourced or marked as assumed in <code>docs/INSTRUMENTS.md</code>"])); }],
  ["Field of view sets the size of the imaged area. Room outside the field lets cells wander out of view and new ones walk in. And the checkbox for bleaching each cell by its own light dose is on by default.", (s) => { s.focus(s.between("Field of view", "3D z-stacks"), "Field, room outside, bleaching", { maxZoom: 2.2, padX: 1.3, padY: 1.1 }); s.note(NOTE("Field and bleaching", ["Small fields are much cheaper for z-stacks", "Cells outside the lit area do not bleach", "Daughters inherit the mother's dose"])); }],
], ["Expand New single-cell loci run", "Keep the defaults and press Run single-cell loci simulation (about 2 minutes)", "Open the Cells page when it finishes"]);

SC("form_loci2", "Form 2: 3D, calibration, models, spot finding", shot({ shot: "sb_loci", align: "left", minScale: 1.2, label: "sidebar · New single-cell loci run" }), [
  ["Further down is three-dimensional imaging. Record 3D z-stacks turns every frame into a stack of planes. Nuclei become ellipsoids, loci get a height, and the blur extends along z. The page warns you before a very large run.", (s) => { s.focus(s.between("3D z-stacks", "Calibration"), "3D z-stacks", { maxZoom: 2.2, padX: 1.3, padY: 1.05 }); s.note(NOTE("3D z-stacks", ["Planes per stack and plane spacing", "Fiducial beads on the slide", "Confocal instead of widefield", "Disk space = frames x planes x pixels x colours"])); }],
  ["Calibration simulates dark frames, a flat-field slide and a bead slide with the same instrument, estimates the instrument from them, and analyses the loci a second time with the correction applied.", (s) => { s.focus(s.between("Calibration", "Where do the locus movements come from?"), "Calibration", { maxZoom: 2.2, padX: 1.3, padY: 1.2 }); s.note(NOTE("Calibration", ["Dark frames, flat field, beads", "Blind estimate, then checked against the truth", "Locus analysis repeated with the correction"])); }],
  ["Where do the locus movements come from lets you choose invented stand-in motion, which is fast and needs no library, or the saved chromatin library, with your choice from the Library page.", (s) => { s.focus(s.between("Where do the locus movements come from?", "Model probe binding"), "Stand-in or library", { maxZoom: 2.2, padX: 1.3, padY: 1.2 }); s.note(NOTE("Source of locus motion", ["<b>Stand-in</b>: confined diffusion, invented", "<b>Library</b>: real polymer simulations (GPU, made once)"])); }],
  ["Two checkboxes add the biology. Model probe binding makes locus brightness follow a probe model. Model distance-dependent transcription adds M S two spots at the promoter, which switches on more often when the enhancer is close.", (s) => { s.focus(s.between("Model probe binding", "Spot finding (stage 6)"), "Probes and transcription", { maxZoom: 2.2, padX: 1.3, padY: 1.2 }); s.note(NOTE("Biology switches", ["Probe binding (settings from the Probes page)", "MS2 transcription (settings from the Transcription page)", "Works best with the library: real distances"])); }],
  ["Spot finding chooses the method: the Chen two thousand twenty-five difference of Gaussians with Viterbi, or the original brightest pixel. The candidate threshold is deliberately low, because the Viterbi path chooses among the candidates.", (s) => { s.focus(s.between("Spot finding (stage 6)", "Image errors"), "Stage 6 method", { maxZoom: 2.2, padX: 1.3, padY: 1.2 }); s.note(NOTE("Spot finding", ["Chen 2025 (DoG + Viterbi) or the old brightest pixel", "Candidate threshold in multiples of the noise"])); }],
  ["Finally the same image-errors controls as in the first form, and the big run button.", (s) => { s.focus(s.between("Image errors", null), "Run", { maxZoom: 2.2, padX: 1.3, padY: 1.1 }); s.note(NOTE("Run", ["Same image-error switches as form 1", "Press Run single-cell loci simulation"])); }],
], ["Tick Record 3D z-stacks with a 50 micron field and 15 planes, and run it", "Or tick probe binding and transcription with the library source", "Read the disk-size warning before running 3D"]);

SC("run", "What happens when you press Run", (st) => {
  const o = { instant: false }; o.term = new Term(st);
  o.destroy = () => o.term.destroy();
  o.lines = (from, to) => { const L = [
    ["Run folder: C:\\Users\\you\\SimLiveData\\runs\\loci_demo_seed11_24bba2", "dim"], ["Stage 2: simulating cell motion and shape", ""], ["Stage 1: simulating loci inside each nucleus", ""],
    ["Stage 1b: stochastic probe binding at the promoter and enhancer (sets each locus' brightness)", ""], ["Stage 1c: stochastic transcription driven by promoter-enhancer distance (the MS2 signal)", ""],
    ["Stage 3: rendering nucleus + locus0 + locus1 (MS2 spots are in the nuclear channel)", ""], ["Stage 4: segmenting", ""], ["Stage 4: tracking on gt_masks (trackastra)", ""], ["Stage 4: tracking on auto_masks (trackastra)", ""],
    ["Stage 5: isolating cells (identity source: truth)", ""], ["Stage 6: locating loci in 12 cells in parallel", ""], ["Stage 6: photobleaching correction (background removed first)", ""],
    ["Stage 6: deconvolving the MS2 traces into Pol II initiation events (Chen et al. 2025)", ""], ["Done", "ok"]].slice(from, to);
    o.term.print(L.map(([t, c]) => ({ html: t, cls: c, gap: 450 }))); };
  return o;
}, [
  ["When you press a run button, the dashboard starts a separate process and streams its progress lines into a status box. These lines are what it prints.", (o) => { o.term.clear(); o.term.instant = o.instant; o.lines(0, 1); }],
  ["Each line is a stage. Stage two simulates the cells. Stage one places the loci, and one B and one C add the probes and the transcription. Stage three renders the images.", (o) => o.lines(1, 6)],
  ["Stage four segments and tracks, once on perfect masks and once on automatic masks. Stage five isolates every cell. Stage six locates the loci, corrects bleaching and deconvolves the M S two traces.", (o) => o.lines(6, 13)],
  ["When it says Done, the new run is selected automatically and every page shows it. If the library source is chosen but no library is picked, the form asks you to choose one in the Library page first.", (o) => o.lines(13, 14)],
  ["Every run lives in a folder named after the configuration and the seed, in the Sim Live Data folder outside the repository. The same settings and seed always reproduce the same run.", (o) => { }],
  ["One practical rule: the graphics card is shared. Do not start a movie while a chromatin library is being generated, and run the polymer simulations and the movie pipeline one after the other.", (o) => { }],
], ["Run any form and watch the status box", "Find the new run in Run to view", `Open the Config page to see the exact settings and the run folder`]);

/* ================================================================== 3 */
T.chapter("3 · Cells and tracking");
SC("motion", "Cells move, the tracker follows", (st) => A.tracking(st), [
  ["Stage two moves the nuclei. Each has a persistent ID and an irregular shape that changes over time. They diffuse, drift, collide, and some divide. This is the ground truth.", (a) => a.panels(1)],
  ["Stage four then does what a real lab would. It segments each frame and links the nuclei through time with Trackastra, and gives them its own track numbers. Colours match when a cell is followed correctly.", (a) => a.panels(2)],
  ["When two nuclei touch, the tracker can swap their identities. That is an I D switch, marked by the red circle, and counted at the bottom.", (a) => { a.swaps(true).crowd(true); }],
  ["Faster cells make more mistakes. In our tests, errors came from crowding multiplied by speed.", (a) => a.speed(2.2)],
  ["Why does this matter for chromatin? A swapped cell carries its loci to the wrong track, so a locus trajectory would be built from two different cells. Stage seven measures how much this costs.", (a) => a.speed(1)],
], ["Open New simulation and run it at Speed scale 0.5, then at 1.5", "Compare the ID switches in Metrics"]);

SC("movie", "Page: Movie", shot({ sets: { main: "pg_movie", player: "pg_movie_player" }, label: "🎞 Movie" }), [
  ["The first page is Movie. At the top, two boxes say what the page shows and how to use it. Every page has them.", f("intro", "What it shows · how to use it")],
  ["Open Play the movie and switch on Prepare player for smooth playback in your browser. The frames are drawn once, then played without waiting. Space plays and pauses, the arrow keys step.", (s) => { s.use("player"); s.focus(["player", "prepare"], "Prepare player", { multi: true }); s.click("player"); }],
  ["Tracked on chooses which tracker result you see: automatic segmentation, which is realistic, or perfect ground-truth masks, which isolates the linking error.", (s) => { s.use("main"); s.focus("tracked", "Which tracker result"); }],
  ["Next error jumps to the next frame where something went wrong, and Previous error goes back. The small arrows step one frame.", (s) => { s.focus(["prev", "next"], "Jump between errors", { multi: true }); s.click("next"); }],
  ["Or drag the frame slider. Below it you see how many frames contain an error, and which ones.", (s) => s.focus(["slider", "errcount"], "Frame slider", { multi: true })],
  ["On the left are the true cells with their true ID numbers. On the right is what the tracker made, with its own track numbers. A track is coloured by the true cell it mostly follows, so matching colours mean it was followed correctly.", f("images", "True cells · tracked cells")],
  ["Red circles mark identity switches. Here the tracker swapped two cells that touch. Their numbers on the right no longer match their true colours.", f("images", "Red circle = identity switch", { maxZoom: 1.2 })],
  ["A red box spells it out: how many switches, missed cells and false-positive masks this frame has, and which track each cell had before and after.", f("errbox", "What went wrong in this frame")],
  ["The bar chart shows identity switches over time, and the folded box explains every outline colour on the pictures.", f(["howto", "bars"], "Switches over time", { multi: true })],
], ["Open the Movie page of any run", "Press Next error and read the red box", "Switch Tracked on to ground-truth masks: do the errors disappear?", "Tick Also show raw image panel to see what the segmenter sees"]);

SC("metrics", "Pages: Metrics and Trajectories", shot({ sets: { main: "pg_metrics", traj: "pg_traj" }, label: "📊 Metrics" }), [
  ["Metrics starts with a sentence in words: what fraction of nuclei were found, how many links were right, and how many identity switches occurred, for both automatic and perfect masks.", f("inwords", "In words")],
  ["Below is the table of every score for both variants. Compare the two columns: any gap between them is caused by the segmentation, not the tracker.", f("table", "Automatic vs perfect masks")],
  ["The folded box defines the numbers. DET is detection, LNK is linking, TRA is overall tracking, and CHOTA combines detection and association. Track purity and target effectiveness say how much of a track follows one cell.", f("meaning", "What do these numbers mean?")],
  ["The per-cell identity table shows which cells swapped. Sort it to find the troublemakers.", f("percell", "Per-cell identity table")],
  ["Trajectories draws the true path of every cell in colour and the tracker's paths as dashed lines, with red crosses at switches. A dashed line that jumps between coloured paths is a swap.", (s) => { s.use("traj"); s.focus("plot", "Coloured = truth · dashed = tracker · ✕ = switch"); }],
], ["Open Metrics and read the In words sentence", "Open Trajectories and find a red cross", "Compare both columns: which is larger, the segmentation or the linking error?"]);

SC("scans", "Page: Scans", shot({ shot: "pg_scans", label: "📈 Scans" }), [
  ["A single run is not evidence. The number of identity switches varies a lot from seed to seed, by five to ten, so compare conditions over at least five seeds. The Scans page shows such sweeps.", f("intro", "Many seeds, not one run")],
  ["Choose a scan first. Sweeps are made from a terminal with scripts, sweep dot p y, and each is stored in its own folder.", f(["scan", "caption"], "Pick a saved scan", { multi: true })],
  ["Then pick a metric, and which masks to score: automatic, ground truth, or both.", f(["metric", "masks"], "Metric and masks", { multi: true })],
  ["Every setting has exactly one role: the x axis, the colour of lines, a separate panel, or held at one value. Change the roles to look at the same sweep from different sides.", f("roles", "One role per setting")],
  ["Lines are means over seeds, bars are one standard deviation, and dots are single seeds. Differences smaller than the bars are not reliable. In this speed scan, scores stay high at slow speeds and fall as the cells speed up.", f("plot", "Mean ± 1 sd, dots = seeds")],
  ["A table of means and standard deviations is folded underneath, ready to copy.", f("table", "Table")],
], ["Run a sweep: scripts\\run.cmd python scripts\\sweep.py --name speed --seeds 1-5 --grid motion.speed_scale=1,0.5,0.25", "Open Scans and choose it", "Change the x axis and the colour of lines"]);

/* ================================================================== 4 */
T.chapter("4 · The microscope");
SC("forward", "From a spot to an image", (st) => A.microscope(st), [
  ["Stage three is the microscope. We start from an ideal point of light, exactly where the locus is. The cyan circle marks the truth.", (a) => a.step(0)],
  ["Diffraction blurs the point into the point spread function. The brightest pixel is still in the right place.", (a) => a.step(1)],
  ["Everything else in the nucleus that glows in this colour adds a haze under the spot.", (a) => a.step(2)],
  ["Light comes as photons, so there is shot noise. It is largest, relative to the signal, for dim spots.", (a) => a.step(3)],
  ["The camera adds its own electronic read noise to every pixel.", (a) => a.step(4)],
  ["Then the defects: hot pixels and a cosmic ray. Now the brightest pixel is red, because a naive detector would report a hot pixel instead of the locus.", (a) => a.step(5)],
  ["And even a perfect detector would be off, because stage drift and the chromatic shift between colours draw the spot in the wrong place. Some errors can be removed by the analysis, others only by calibration.", (a) => a.step(6)],
  ["The current version has twenty-six switchable sources. The Images page, next, explains each one and measures what it costs.", (a) => a.step(6)],
]);

SC("biplane", "Vutara biplane: two planes at once", (st) => A.biplane(st), [
  ["The Vutara can image in two ways. In its widefield mode it records one focal plane at a time. In its biplane mode it records two planes at the same moment.", (a) => a.phase(0)],
  ["A beam splitter sends half of the light to each half of one camera, so two planes six hundred nanometres apart in the sample appear side by side. The height of a spot comes from how bright it is in the two planes.", (a) => a.phase(1)],
  ["Two planes cover only a thin slice, and the method works best near or between them. A nucleus is about ten microns thick, so you still step through it, as in an ordinary z-stack.", (a) => a.phase(2)],
  ["Each stage position gives two planes, so a step of one point two microns covers the sample evenly every point six microns, with half as many stage moves. The price is that each plane gets half of the light.", (a) => a.phase(3)],
  ["The two camera halves are never registered perfectly. The second plane's image is shifted by a fraction of a pixel, which biases the positions found.", (a) => a.phase(4)],
  ["Beads in the image sit between the lowest pair of planes, so they appear in both. The difference of their two positions measures the shift in every frame, and the program undoes it. In a test it measured nought point three two and minus nought point two one pixels for a true nought point three five and minus nought point two five.", (a) => a.phase(5)],
  ["Two things are not measured this way: the true plane separation, and the z calibration of the point spread function, which the manufacturer does with a bead scan. They stay at their nominal values.", (a) => a.phase(5)],
  ["Your turn: in the single-cell loci form, choose Bruker Vutara VXL biplane 3D mode, keep the field at forty microns, and run a short movie. Then read the Spot finding and Stage pages.", (a) => a.phase(5)],
], ["In New single-cell loci run choose Microscope: Bruker Vutara VXL (biplane 3D mode)", "Use a 40 µm field and 30-60 frames (the page shows the disk size first)", "Open 🔎 Spot finding and 🖼 Images; open 📍 Stage &amp; fiducials", "Compare with a run in the widefield mode with Record 3D z-stacks ticked"]);

SC("images", "Page: Images", shot({ sets: { main: "pg_images", err: "pg_images_err", build: "pg_images_build" }, label: "🖼 Images" }), [
  ["The Images page has five parts. It puts every recorded image next to the ground truth behind it. Part one lets you pick a frame and see the camera images with the truth drawn on top.", f("p1", "Part 1: image + truth")],
  ["Choose the frame, the channels, which are the nucleus, the promoter and the enhancer, and the brightness ceiling. Log scale helps with dim things.", f(["frame", "channels"], "Frame and channels", { multi: true })],
  ["View switches between the whole field, one cell, and a zoom on a locus. Then tick what to draw on top: nucleus outlines, ID numbers, true locus positions, M S two spots, and beads.", f(["view", "draw"], "View and overlays", { multi: true })],
  ["Cyan circles are where a locus truly is. Red crosses are where that colour's spot is actually drawn, which differs by the chromatic shift. A diamond marks a promoter with an active M S two signal.", f("draw", "○ truth · ✕ drawn · ◇ MS2")],
  ["You can also show the error-free version, and the difference between recorded and error-free, which is exactly what the errors added.", f("errfree", "Error-free version")],
  ["The numbers below change from frame to frame: stage drift, focus error, the blur size, and the light level in each channel.", f("stats", "This frame's imaging conditions")],
  ["The ground-truth table lists the true locus positions and, for transcription, whether the promoter is on, the distance between promoter and enhancer, how many Pol two are on the gene, and the M S two signal.", f("truth", "Truth for one cell, this frame")],
  ["Part two is about one error source at a time. Choose a source, and the page explains what it is, what it does to the picture, and how it is corrected: by the analysis itself, with calibration images, or from the movie alone.", (s) => { s.use("err"); s.focus(["source", "what", "does"], "Pick one error source", { multi: true }); }],
  ["It also measures what the error does to the results, with and without each correction. In our showcase run the colour shift between channels dominated position errors, and a bead slide fixed it. Bleed-through of the bright M S two spot dominated intensities.", f("results", "Measured effect on the results")],
  ["Part three builds the image up one error at a time, from the ground truth to the recorded frame, with a table of what each source adds.", (s) => { s.use("build"); s.focus("p3", "Part 3: build-up", { maxZoom: 1.1 }); }],
  ["Part four ranks which error changes this frame most, and part five explains every source in full, with this run's settings.", (s) => { s.use("main"); s.focus(["p4", "p5"], "Parts 4 and 5", { multi: true }); }],
], ["Open Images, pick a frame, tick true locus positions", "Tick Also show the error-free version", "In part 2 choose Hot pixels, then Chromatic shift", "Press Build the image up in part 3"]);

SC("calibration", "Page: Calibration", shot({ shot: "pg_cal", label: "🎯 Calibration" }), [
  ["Real instruments are never perfect, so you calibrate them. This page simulates the calibration with the very instrument that recorded the movie, then checks every estimate against the truth.", f("intro", "Calibrating the microscope")],
  ["Four steps. Acquire pictures of things whose appearance is known. Estimate the instrument from those pictures alone. Check, because the simulation knows the truth. Apply the correction to the movie.", f("steps", "Acquire · Estimate · Check · Apply")],
  ["The table compares each estimate with the true value: gain, offset, read noise, dark current, hot pixels, illumination, blur, colour offsets and drift. In this run forty-six of forty-seven quantities were within tolerance.", f("table", "Estimated versus true")],
  ["Dark frames, shot with the shutter closed, show only the camera: its offset, fixed pattern, hot pixels, read noise and cosmic rays.", f("dark", "Dark frames")],
  ["The flat-field slide glows evenly. Several brightnesses give the photon transfer curve, whose slope is the gain, and the average shows the vignetting and pixel gain differences.", f("flat", "Flat-field slide")],
  ["The bead slide has tiny sources that appear in every colour. They reveal the blur, the offsets between colours, and the drift.", f("bead", "Bead slide")],
  ["And finally, does it help? The table compares the movie analysed without and with the calibration. In this run the median locus error fell from about forty-one to about twenty-seven nanometres.", f("help", "Does it help?")],
], ["Open Calibration on a run that has it (Simulate and apply a calibration was ticked)", "Look for the one warning in the table and read its note", "Compare the two rows of Does it help"]);

SC("stage", "Page: Stage and fiducials", shot({ shot: "pg_stage", label: "📍 Stage & fiducials" }), [
  ["This page needs a three-dimensional run. A z-stack is made by stepping a piezo through planes, and a real piezo is imperfect. The page explains how that matters.", f("intro", "z-piezo and beads")],
  ["Four cards tell the story. The piezo is asked for planes at nominal heights. But each step has a small random error, the response bows, differs going up and down, and lags after starting.", f("cards", "The four ideas")],
  ["Beads lie on the slide in every image and show how far the whole picture drifted. But they all sit at one height, so they only see the planes near the slide.", f("s2", "Beads on the slide")],
  ["Only the piezo's own position log corrects plane-to-plane errors. The last table shows what each correction does to locus distances and heights.", f("s3", "Does correcting help?")],
], ["Make a 3D run: tick Record 3D z-stacks and fiducial beads", "Open this page and compare the two corrections"]);
