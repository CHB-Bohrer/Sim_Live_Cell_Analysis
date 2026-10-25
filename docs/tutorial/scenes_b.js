/* Scenes, part 2: chromatin, finding loci, transcription biology, configuration, workflows. */
"use strict";

/* ================================================================== 5 */
T.chapter("5 · Chromatin");
SC("polymer", "Stage 1: the polymer", (st) => A.polymer(st), [
  ["Stage one is the chromatin. A polymer simulation of ten thousand monomers, one kilobase each, so ten megabases, run on the graphics card with polychrom and OpenMM, inside a confining sphere.", (a) => a.pe(false)],
  ["Two beads are tagged: the promoter and the enhancer, a hundred kilobases apart. How far apart they are, moment by moment, is what the transcription model needs.", (a) => a.pe(true)],
  ["Then loop extrusion. A factor lands on the chromatin, and both of its legs walk away from each other, growing a loop and pulling distant points together.", (a) => a.extrude(true)],
  ["Boundary elements, like C T C F sites, stop the factor. The loop is limited to a domain, and contacts across the boundary are reduced. Boundary positions are synthetic until you supply real ones.", (a) => a.boundaries(true)],
  ["We checked that the engine reproduces known behaviour, in a validation report. Eleven of twelve checks pass. The twelfth, insulation by the boundaries, is real but weaker than published, and the report says so honestly.", (a) => { }],
  ["A simulation only has polymer units. The calibration converts them: fifty-five nanometres per unit and nought point four six seconds per saved block. These come from two assumptions: a hundred kilobase pair is within four hundred nanometres ten percent of the time, and a single locus has a diffusion coefficient of two thousandths of a square micron per second. Replace them with a measured mean squared displacement.", (a) => a.scale(true), 900],
]);

SC("library", "Page: Library", shot({ shot: "pg_lib", label: "🗂 Library" }), [
  ["The Library page is about the saved polymer simulations that movies draw their loci from. A library is made once, over hours, then reused for many movies.", f("intro", "Saved chromatin simulations")],
  ["Four cards explain the path. The graphics card generates the simulations once. For a movie, each cell is given its own saved simulation. The two loci are read from it, converted to nanometres, rotated and placed in the nucleus. Then the frames are rendered and analysed.", f("how", "Generate · Choose · Place · Image")],
  ["One selector at the top chooses the library for the whole page. The live status refreshes every five seconds: simulations finished, the one running now and its time left, and the graphics card's load. If the computer sleeps it says Stalled.", f(["libsel", "live"], "One library, its live status", { multi: true })],
  ["Below, you choose how cells get their simulation: random from the whole library, random from the ones you tick, or pin particular cells to particular simulations.", f(["choose", "radio"], "Choose the simulations", { multi: true })],
  ["The preview table shows exactly what a movie with this seed will give each cell, if no cell divides.", f("preview", "Preview")],
  ["To use it, pick the saved chromatin library as the source in the loci form. The same choices work in config files with sim IDs and pin. Do not start a movie while the generator runs, the graphics card is shared.", f("use", "Use this in a movie")],
], ["Open Library and look at the status of the generator", "Choose Random from the ones I tick and tick three", "Run a loci simulation with the library source and check the preview matches"]);

SC("chromatin", "Page: Chromatin", shot({ shot: "pg_chrom", label: "🧬 Chromatin" }), [
  ["The Chromatin page looks inside one saved polymer simulation: the chain moving through time, its contact map, and the distance between any two loci.", f("intro", "Inside one simulation")],
  ["The folded validation report shows how the engine reproduces known behaviour from the Mirny lab, with every check and its result.", f("validation", "Validation report")],
  ["Pick a library and a simulation. The table lists each one with its speed, run time and mean loop size. The cards below show the monomers, the speed in steps per second, and the mean loop size.", f(["lib", "sim"], "Library and simulation", { multi: true })],
  ["Open the smooth player to watch the chromosome move in your browser, with two chosen loci and the loop-extruder lines. Or switch on the contact map, the contact probability against distance, and the three-dimensional shape.", f(["player", "contact"], "Watch it move", { multi: true })],
  ["Two loci: distance over time is the part that feeds the microscope step. Choose locus A's position and the separation to locus B in kilobases, and see how far apart they are in every saved block, with a histogram.", f(["two", "locusA", "sep"], "Two loci: distance over time", { multi: true })],
  ["Loop-extruding factors over time shows how many loops cover each position at each time, as a heatmap, with the mean loop size underneath. Bright bands are often inside a loop.", f("loops", "Loops over the whole simulation")],
  ["The zoom shows individual factors as bands between their two legs. A band widens as the factor extrudes, stops at a boundary, and breaks where the factor unbinds and reloads elsewhere.", f("zoom", "Zoom: individual factors")],
], ["Open Chromatin and pick the first simulation", "Choose locus A = 4000 kb and separation 100 kb", "Compare the distance trace with 500 kb", "Zoom into a 400 kb window in the loop view"]);

/* ================================================================== 6 */
T.chapter("6 · Finding the loci");
SC("cells", "Page: Cells", shot({ sets: { main: "pg_cells", copies: "pg_cells_copies", d3: "pg_cells_3d" }, label: "🔬 Cells" }), [
  ["Stage five cuts every cell out of the movie. Each cell becomes a fixed-size window that follows its centre, with the nucleus and one channel per locus. The Cells page shows them.", f("intro", "One small movie per cell")],
  ["First choose which identities to use: the true cell IDs from the simulation, or the tracker's. Running the same analysis on both is how tracking errors are carried into locus results.", f("ident", "True IDs or tracker IDs")],
  ["The summary table lists each cell's lifetime, gaps, border contact and nearest neighbour. Clean tracks are the ones to trust.", f("table", "Cell summary")],
  ["All cells at one frame shows every window side by side, and a player runs them together.", f(["playall", "allframe", "grid"], "All cells", { multi: true })],
  ["Then pick one cell to watch through time, with its located loci. White plus signs are what the analysis found, cyan circles are the truth, and red crosses are where each colour is really drawn.", f(["onecell", "cellsel"], "One cell through time", { multi: true })],
  ["The picture above is flat: it is the view from the top. In a three-dimensional run, a new part appears underneath it, called Side views, which shows where each spot sits along z, the height above the coverslip.", fs("d3", "side", "Side views: where the spots are along z")],
  ["Pick a spot in the list. The two left panels cut the recorded stack through it: x against z, then y against z, with the coverslip at the bottom. The blur is the axial point spread function, longer than in x and y.", fs("d3", ["spotsel", "sidefig"], "Cuts through the chosen spot", { multi: true })],
  ["White plus signs are the heights the analysis found. Cyan circles are the true heights, and red crosses are where each spot is really drawn. A plus far from its circle means a spot found at the wrong height.", fs("d3", "sidefig", "Found, true and drawn heights")],
  ["The right panel follows the height of every spot through time. Thick faint lines are the truth, and the grey band is the range of the stack. A spot near the edge of the band is easy to lose.", fs("d3", "sidefig", "Height through time")],
  ["Some cells now carry two promoter-enhancer pairs. Each found spot is labelled P for promoter or E for enhancer, plus its copy number, and pink squares mark stray spots, the decoys the finder has to ignore.", (s) => { s.use("copies"); s.focus("legend", "P0, E0, P1, E1 and stray spots"); }],
  ["The figure underneath is built to answer three questions. Top left: how far is every found spot from its true locus, frame by frame? A pink square there means the spot found was a stray, a diamond means it was the other copy.", (s) => { s.use("copies"); s.focus("fig", "Errors, offsets, distances"); }],
  ["Top right: every found position minus the true one. The width of the cloud is the precision, and its centre is the bias, for example an uncorrected colour shift. Below: the promoter-enhancer distance of each copy, found against true.", (s) => { s.use("copies"); s.focus("fig", "Precision, bias and each copy's distance"); }],
], ["Open Cells and choose True cell IDs", "Pick a cell with 2 promoter-enhancer copies", "On a 3D run (loci_3d), open Cells and look at Side views: choose a spot and compare the white plus with the cyan circle along z", "Find the frames with a pink square in the error plot, and look at them in the image above"]);

SC("trellis", "How a spot is found", (st) => A.trellis(st), [
  ["Stage six finds each locus in each cell, following Chen, Levo, Zoller and Gregor, two thousand twenty-five. In every frame we have a set of bright blobs, and one of them is the locus.", (a) => a.step(0)],
  ["First a median filter and a difference of Gaussians, with a deliberately low threshold. So there are several candidates per frame: bigger circles are brighter.", (a) => a.step(1)],
  ["Picking the brightest candidate in each frame jumps around, and it is wrong whenever a noisy blob outshines the real spot.", (a) => a.step(2)],
  ["The Viterbi algorithm instead finds the most likely whole trajectory. Jumping a long way is expensive, because a locus diffuses slowly, and bright candidates earn a reward.", (a) => a.step(3)],
  ["If a frame has no candidate at all, the path can bridge the gap, paying a miss cost, and even pick the spot up again after a long gap. Then a Gaussian fit refines the position.", (a) => a.step(4)],
]);

SC("copies", "Two copies and stray spots", (st) => A.copies(st), [
  ["Real cells are rarely as tidy as one promoter and one enhancer. A diploid cell can show both alleles, so the simulation can give a cell two promoter-enhancer pairs. Half the cells get two in the example settings.", (a) => a.phase(0)],
  ["Each copy has its own chromatin trajectory, from its own polymer simulation, and its own promoter, so each one bursts on its own, driven by its own promoter-enhancer distance. A copy can even get different rates, here twice as many bursts.", (a) => a.phase(1)],
  ["And every nucleus holds stray spots: probe or M S two aggregates that appear somewhere, wander, and vanish after about a minute. They are as bright as a locus, so a finder that takes the brightest spot will pick them.", (a) => a.phase(2)],
  ["The one-path finder follows one spot per colour, so it can never see two copies. The multi-copy linker extracts several paths per colour, each through the candidates the earlier ones left.", (a) => a.phase(3)],
  ["A path is kept only if it looks like a locus: there in most frames, about as bright as the best path, and as slow. Paths through noise and strays fail that test. Then each promoter path is paired with an enhancer path that stays within one and a half microns of it.", (a) => a.phase(4)],
  ["When two copies come close, two paths can swap identities. Exchanging their tails wherever that makes the steps shorter undoes it, and a pair prior anchors each promoter to its own enhancer. Over twenty-five movies the copy number was right in every cell without strays, and in nine of ten cells with four strays per minute.", (a) => a.phase(5)],
]);

SC("spots", "Page: Spot finding", shot({ sets: { main: "pg_spots", copies: "pg_spots_copies" }, label: "🔎 Spot finding" }), [
  ["The Spot finding page shows how stage six works. It has four parts: the method, every step on one cell, the results for this run, and a comparison of methods.", f("intro", "How loci are found")],
  ["Part one lays out the method next to the paper, step by step: one stack per nucleus, the median filter, the difference of Gaussians, the Viterbi trajectory, the position, and the intensity measured in an ellipsoid minus the ring around it.", f(["method", "table"], "Method and where it comes from", { multi: true })],
  ["Part two shows every intermediate image. Choose a cell, a colour, which is the promoter, the enhancer or the M S two spot, and a frame, or jump to a problem frame. Change the settings, press Apply, and watch what each step does.", f(["step2", "cellsel", "spotsel"], "Every step on one cell", { multi: true })],
  ["In a run with two copies, part two B shows what the linker did in the chosen cell. Every row is one path it tried, every tick a frame, coloured by what that spot really was: copy zero, copy one, a stray, or noise. A good result is one solid colour per kept path.", fs("copies", "s2b", "What the linker did")],
  ["Part three is the result for this run: the error of every frame against the truth, with detection rates, the share of picks on a stray spot or on the other copy, and the M S two trace at the tracked promoter.", fs("copies", "s3", "This run's result")],
  ["Three B corrects photobleaching, with the background removed first. Three C estimates uneven illumination and per-cell bleaching from the movie itself, with no flat-field image, and keeps the correction only if it predicts held-out cells better.", fs("main", ["s3b", "s3c"], "Bleaching and illumination", { multi: true })],
  ["Three D asks how bright an M S two spot must be before transcription is measurable at all, a sweep over many movies. The answer depends on the brightness of one polymerase compared with one nuclear pixel, so it must be measured on real data.", fs("main", "s3d", "How bright must the MS2 spot be?")],
  ["Three E is the sweep for copies and stray spots: the same movies analysed four ways, with more and more strays. One path per colour misses the second copy and jumps between copies; the multi-copy linker with the pair prior keeps them apart.", fs("copies", "s3e", "Copies and strays over many movies")],
  ["Part four compares methods on this whole run, for example the original brightest pixel against the new spot finder, so you can see what an analysis choice costs.", fs("main", "s4", "Compare methods")],
], ["Open Spot finding and read part 1", "In part 2 choose a cell with two copies, the promoter, and look at part 2b", "Read part 3e: which analysis keeps the copies apart?"]);

/* ================================================================== 7 */
T.chapter("7 · Transcription biology");
SC("probe_anim", "Probes bind a locus", (st) => A.probes(st), [
  ["A locus is labelled by many fluorescent probes. They float freely in the nucleus. Nothing sticks yet.", (a) => a.phase(0)],
  ["A free probe lands on the locus at rate k bind, and becomes bound.", (a) => a.phase(1)],
  ["A bound probe starts scanning along the locus at rate k scan, or lets go again at rate k unbind. A scanning probe is still fluorescent.", (a) => a.phase(2)],
  ["After about one over k off seconds a scanning probe falls off. Each probe cycles on its own and at random, so the number attached fluctuates, and the locus brightness is proportional to that number.", (a) => a.phase(3)],
  ["Faster turnover means a noisier trace. The promoter and the enhancer have separate settings. All rates here are placeholders until measured kinetics are supplied.", (a) => a.koff(0.5)],
]);

SC("probes", "Page: Probes", shot({ shot: "pg_probes", label: "🧪 Probes" }), [
  ["The Probes page lets you explore this model and put it into a movie.", f("intro", "Probe binding model")],
  ["The explainer repeats the cycle: free, bound, scanning, fall off. It is simulated exactly, event by event, with a Gillespie algorithm.", f("cycle", "How probes bind a locus")],
  ["Choose a preset or change the rates for the promoter. You see the average number attached, the noise, how long a probe stays, and how long the count fluctuates. Scrub through an example trace and watch the probes move through the states.", f("promoter", "Promoter (locus 0)")],
  ["The enhancer has its own panel, with its own rates.", f("enhancer", "Enhancer (locus 1)")],
  ["The parameter sweep runs the model over a range of one parameter, in seconds. Simulated points sit on the exact theory line when the model behaves.", f("sweep", "Parameter sweep")],
  ["To use these settings in a movie, tick Model probe binding in the loci form of the sidebar. The command-line version is shown in the last box.", f("use", "Use them in a movie")],
], ["Open Probes and pick a different preset for the promoter", "Drag Look at time and watch the states", "Run a sweep of k off"]);

SC("tx_anim", "Distance, promoter, Pol II, MS2", (st) => A.transcription(st), [
  ["Now the central question. The enhancer and promoter move around. The distance between them, from a saved polymer simulation, goes in at the top.", (a) => a.show(0)],
  ["The promoter switches on and off at random. While the enhancer is close, it switches on more often. Only one rate is raised, which one you choose: k on, the initiation rate, or k off.", (a) => a.show(1).shade(true)],
  ["While the promoter is on, Pol two molecules start and travel along the gene. Each copies the M S two cassette first, then the rest of the gene, waits at the end and is released.", (a) => a.show(2)],
  ["Each Pol two carries M S two loops, so the glow of the spot is the sum over all Pol two on the gene. That sum is what the microscope sees: a bright spot at the promoter, in the nuclear colour.", (a) => a.show(3)],
  ["The polymer does not feel transcription: the coupling goes one way. And all rates are per second of real time, so the frame interval can change freely.", (a) => a.show(3)],
]);

SC("tx", "Page: Transcription", shot({ shot: "pg_tx", label: "⚡ Transcription" }), [
  ["The Transcription page builds this chain in nine numbered steps. Work from top to bottom. Change any number and everything below updates.", f("intro", "Nine steps")],
  ["Step one sets how the movie samples time: frame interval, exposure and movie length. It warns you when there are too few frames per burst.", f("s1", "1 · Sampling")],
  ["Step two is the distance. Take it from a saved polymer simulation by choosing the library, the simulation, and the two genomic positions, or use an invented stand-in. The calibration numbers convert to nanometres and seconds.", f("s2", "2 · Distance")],
  ["Step three sets how the distance changes a rate: a sharp contact, a smooth Hill curve, a smooth exponential, or none as a control. The default is the Hill curve, half its effect at four hundred nanometres.", f("s3", "3 · How distance sets a rate")],
  ["Step four defines the promoter: mean off time, mean on time, how often a Pol two starts while on, the Pol two speed, the cassette and gene sizes, the number of loops and the wait at the end. A plot shows the three timescales against distance.", f("s4", "4 · Promoter switching")],
  ["Step five draws the Pol two molecules as lines on a kymograph, starting at this cell's first burst. Slope is speed, and the flat top is the wait before release.", f("s5", "5 · Pol II along the gene")],
  ["Step six shows the M S two signal over time with the exact expectation averaged over many cells.", f("s6", "6 · The MS2 signal")],
  ["Step seven shows what the microscope sees: tiles of the nucleus at successive times. A faint spot is invisible against the haze, which is where real transcription events get missed.", f("s7", "7 · What the microscope sees")],
  ["Step eight sweeps one setting over many simulated cells. Step nine gives the command, or the sidebar checkbox, to run it in a movie.", f(["s8", "s9"], "Sweep and hand-off", { multi: true })],
], ["Open Transcription and read the diagram at the top", "In step 3 switch the coupling shape between Sharp contact and Hill", "Press New random cell in step 4", "Tick the transcription checkbox in the loci form and run a library movie"]);

SC("deconv_anim", "From the MS2 trace to Pol II starts", (st) => A.deconv(st), [
  ["The recorded M S two trace is a blurred-in-time version of events. Every Pol two start, at a certain moment, is the real thing we want.", (a) => a.step(0)],
  ["One Pol two makes a step-like shape: the signal climbs as it copies the loops, stays flat, and drops at release. That is the kernel.", (a) => a.step(1)],
  ["The trace is the sum of all these shapes, plus noise. This sum is what the camera records.", (a) => a.step(2)],
  ["Deconvolution works backwards. Many possible histories are consistent with the trace, and where they agree the answer is certain. The result is a probability of a start in each second.", (a) => a.step(3)],
]);

SC("ms2", "Page: MS2 → Pol II", shot({ shot: "pg_ms2", label: "🧮 MS2 → Pol II" }), [
  ["The page M S two to Pol two shows this Bayesian deconvolution of Chen and colleagues, step by step, with the truth next to every step. Pick a run with M S two spots at the top.", f(["intro", "run"], "Pick a run with MS2 spots", { multi: true })],
  ["The yellow note under it matters. When a single Pol two is faint next to the nuclear background, the deconvolution mostly reflects noise, so pick a brighter run to see the method work.", f([90, 350, 1120, 70], "Read this note first")],
  ["Step one is what goes in: the spot intensity measured at the tracked promoter, with background and bleaching removed, divided by one fully loaded Pol two, so one means one polymerase. The unit check compares it with the truth.", f("s1", "1 · What goes in")],
  ["Step two shows the kernel, the shape one Pol two makes: ramp, plateau, wait, release. The assumptions are those of the paper. A random wait in the simulation versus a fixed wait in the model is one source of false events.", f("s2", "2 · The kernel")],
  ["Step three estimates the noise from frame-to-frame differences, since the signal changes slowly and the noise does not. Brighter spots are noisier.", f("s3", "3 · The noise model")],
  ["Step four shows one cell: the trace, the fit, the probability of a Pol two start in each second, and the bursts called by the rule. You can change the burst window and the number of starts needed, live.", f(["s4", "burst"], "4 · One cell", { multi: true })],
  ["Step five is the whole run against the truth: how many true starts were found, how many were false, and whether the rate agrees.", f("s5", "5 · The whole run")],
  ["Five B checks the sampler with four chains from different starts, and the usual convergence numbers. The default sampler is chosen per cell, and it switches to parallel tempering when the plain one fails.", f("s5b", "5b · Has the sampler converged?")],
  ["Step six maps the burst rule. A burst is a period of frequent starts. The paper's default misses bursts of a slow gene, so the maps show how well each window and threshold agrees with the truth.", f("s6", "6 · Choosing the burst rule")],
], ["Open MS2 → Pol II and read steps 1 to 3", "In step 4 pick a cell and change the burst window", "Compare the found starts with the true ones in step 5"]);

SC("ep_anim", "How well can the curve be recovered?", (st) => A.epcurve(st), [
  ["The final scientific question: how likely is the promoter to be in a burst, as a function of its distance to the enhancer? The true curve falls with distance.", (a) => a.n(10000).sigma(0).wrong(0).points(false).recovered(false)],
  ["An experiment gives only a limited number of cells, so the points scatter and the curve is noisy.", (a) => a.n(40).points(true).recovered(true)],
  ["With many more cells the recovered curve approaches the truth.", (a) => a.n(1200)],
  ["But measurement error flattens it. When distances are blurred by localization, colour shift and the projection from three to two dimensions, the effect looks smaller.", (a) => a.sigma(450)],
  ["Wrong spot picks lower the contrast further. A few percent of wrong picks matter much more than modest localization noise.", (a) => a.sigma(120).wrong(0.18)],
  ["The dashboard lets you set all of these and see what is recoverable, before you do the real experiment.", (a) => a.sigma(0).wrong(0)],
]);

SC("ep", "Page: E–P curve", shot({ shot: "pg_ep", label: "🔗 E–P curve" }), [
  ["The E P curve page does this for the polymer distances, seen through a two-dimensional, a z-stack, or a biplane microscope. It plots the true curve against the one an analysis recovers.", f("intro", "Probability of a burst vs distance")],
  ["Part one chooses the microscope, the number of cells, the delay between the two colours, the photons, the focus scatter, the wrong picks, and whether the single locus diffusion refers to a locus or to the distance. Then it simulates and fits, in about forty seconds.", f(["s1", "controls"], "Part 1: simulate and fit", { multi: true })],
  ["The result row compares the true contrast, the probability of a burst at one hundred fifty nanometres over that at a thousand, with the recovered one, and says whether the interval contains the truth. Compare the two before trusting an analysis.", f([90, 690, 1120, 95], "Truth vs recovered")],
  ["The promoter integrates the distance over about forty-five seconds, while the distance fluctuates within ten, so the curve is fitted from the distance history, not the distance now.", f("s1", "Why history matters")],
  ["Part one B is the earlier version with an invented slow distance process, kept for comparison.", f("s1b", "1b · Earlier version")],
  ["Part two shows saved scans: many populations per setting, with contrast, interval width, and whether the interval covers the reference. Choose the x axis, the colour of lines and panels, as in Scans.", f(["s2", "scan"], "Part 2: saved scans", { multi: true })],
], ["Open E–P curve and press Simulate and fit with the defaults", "Untick Enhancer acts, the null control: is a false effect found?", "Raise Wrong picks to 0.1 and simulate again"]);

/* ================================================================== 8 */
T.chapter("8 · Configuration and the command line");
SC("config", "Page: Config", shot({ shot: "pg_cfg", label: "⚙ Config" }), [
  ["The Config page shows the exact, fully resolved settings this run was made with, every default filled in, and where its files are.", f("intro", "Settings of this run")],
  ["Settings are grouped by section: acquisition, geometry, cells, motion, optics, the loci, image errors, spot finding, segmentation and tracking. Open one to read it as YAML.", f(["sections", "acq"], "Settings by section", { multi: true })],
  ["The loci section holds the source of motion, brightness, probes and transcription. Everything you can set in the sidebar is also a setting here.", f("loci", "loci")],
  ["At the bottom is the whole file. Copy it to a new file under configs to repeat or change the run from a terminal. The same config and seed always give the same run.", f(["whole", "files"], "The whole file + run folder", { multi: true })],
], ["Open Config and copy the whole file", "Save it as configs\\my_run.yaml and run it from a terminal"]);

SC("cli", "The command line", (st) => {
  const o = { instant: false }; o.term = new Term(st); o.destroy = () => o.term.destroy();
  const cmds = [
    ["scripts\\run.cmd python scripts\\run_tracking_demo.py --set motion.speed_scale=0.5 --set seed=2", "# one tracking run with changed settings"],
    ["scripts\\run.cmd python scripts\\run_tracking_demo.py --config configs\\loci_library.yaml --set loci.transcription.enabled=true", "# loci from the saved library, with MS2"],
    ["scripts\\run.cmd python scripts\\sweep.py --name speed --seeds 1-5 --grid motion.speed_scale=1,0.5,0.25", "# a grid x seeds, results in SimLiveData\\sweeps"],
    ["scripts\\run.cmd python scripts\\compare_localizers.py --run loci_demo_seed11_24bba2", "# redo only the locus analysis of a finished run"],
    ["scripts\\run.cmd python scripts\\deconvolve_ms2.py --run loci_demo_seed11_24bba2", "# MS2 trace -> Pol II starts"],
    ["scripts\\run.cmd pytest", "# the CPU tests; pytest -m gpu for the graphics-card tests"],
  ];
  o.add = async (i) => { o.term.instant = o.instant; for (const k of [].concat(i)) { await o.term.type(cmds[k][0], 10); o.term.line(cmds[k][1], "dim"); } };
  return o;
}, [
  ["Everything the dashboard does can also be run from a terminal, always through scripts, run dot c m d. Any setting can be overridden with dash dash set, key equals value.", (o) => { o.term.clear(); o.add(0); }],
  ["Pass a different config file to switch to the loci movie with the saved library, and turn on transcription.", (o) => o.add(1)],
  ["The sweep script runs a grid of settings across several seeds and saves a results file that the Scans page reads. Single runs are not evidence.", (o) => o.add(2)],
  ["You can redo only the locus analysis of a finished run, which takes a minute or two, and the M S two deconvolution.", (o) => o.add([3, 4])],
  ["And the tests. Plain pytest runs the CPU tests, and pytest with minus m gpu runs the graphics-card ones in a separate process.", (o) => o.add(5)],
], ["Copy one command above and run it", "Reload the dashboard and find the new run in Run to view"]);

/* ================================================================== 9 */
T.chapter("9 · Putting it together");
SC("recipes", "Three workflows", (st) => new Reveal(st, `<h3>Three typical questions, and the path through the dashboard</h3><div class="grid">
  <div class="card2"><h4>1 · Is my cell tracking good enough?</h4><p>New simulation with your cell speed and density → <b>Movie</b> (Next error) → <b>Metrics</b> (gap between masks = segmentation) → run a <b>sweep</b> over ≥ 5 seeds → <b>Scans</b>.</p></div>
  <div class="card2"><h4>2 · Which microscope error hurts my locus positions?</h4><p>New single-cell loci run (Vutara or generic, calibration on) → <b>Images</b> part 2: pick a source, measure it → <b>Calibration</b> → <b>Spot finding</b> (found, wrong picks, error in nm).</p></div>
  <div class="card2"><h4>3 · Can the E–P coupling be recovered from transcription?</h4><p><b>Library</b> (choose simulations) → <b>Transcription</b> (gene, coupling) → loci run with library + MS2 → <b>Spot finding</b> (MS2 trace) → <b>MS2 → Pol II</b> → <b>E–P curve</b>.</p></div></div>`), [
  ["Here are three typical questions and the path through the dashboard for each.", (r) => r.to(0)],
  ["First: is my cell tracking good enough? Make a tracking simulation with your own cell speed and density, step through the errors in Movie, read Metrics, then sweep over at least five seeds and read Scans.", (r) => r.to(1)],
  ["Second: which microscope error hurts my locus positions? Make a loci run with calibration on, pick error sources in Images and measure them, check Calibration, and read found, wrong picks and the error in nanometres on Spot finding.", (r) => r.to(2)],
  ["Third: can the enhancer promoter coupling be recovered from transcription? Choose simulations in Library, set the gene in Transcription, run a library movie with M S two, read the trace, deconvolve it, and test the curve in E P curve.", (r) => r.to(3)],
], ["Pick the question closest to yours and follow the path once, end to end"]);

SC("cautions", "Things to remember", (st) => new Reveal(st, `<h3>Read the results with these in mind</h3><div class="grid" style="grid-template-columns:1fr 1fr">
  <div class="card2 warn"><h4>Placeholders</h4><p>Brightness, kinetics, error sizes and the polymer calibration are placeholders or assumptions until measured. The numbers show how the method behaves, not biology.</p></div>
  <div class="card2 warn"><h4>One run is not evidence</h4><p>Compare conditions over several seeds. Differences smaller than the bars in Scans are not reliable.</p></div>
  <div class="card2 warn"><h4>One graphics card</h4><p>Chromatin simulations and the movie pipeline share it: run them one after the other. A library takes hours.</p></div>
  <div class="card2 warn"><h4>Restart after code edits</h4><p>Streamlit does not reliably reload helper modules. Restart the dashboard after editing files in <code>app/</code>.</p></div>
  <div class="card2"><h4>Reproducible by design</h4><p>Same config and seed → same run folder. Everything is in <code>SimLiveData</code>, outside the repository.</p></div>
  <div class="card2"><h4>Where to read more</h4><p><code>README.md</code>, <code>docs/ARCHITECTURE.md</code>, <code>docs/SPOT_FINDING.md</code>, <code>docs/IMAGE_ERRORS.md</code>, <code>docs/EP_CURVE.md</code>.</p></div></div>`), [
  ["Before you start, a few things to keep in mind.", (r) => r.to(0)],
  ["Many values are placeholders or assumptions until you supply measurements. The results show how the method behaves, not biology.", (r) => r.only([0])],
  ["A single run is not evidence. Compare conditions over several seeds.", (r) => r.only([0, 1])],
  ["The graphics card is shared, so run chromatin simulations and movies one after the other. And if you edit files in the app folder, restart the dashboard.", (r) => r.only([0, 1, 2, 3])],
  ["Everything is reproducible from a config and a seed, and the documentation folder has the details of each method.", (r) => r.to(6)],
]);

SC("cheat", "Cheat sheet: which page answers what", (st) => new Reveal(st, `<h3>Which page answers which question</h3><table class="q">
  <tr><td>Did the tracker follow the cells? Where did it fail?</td><td>Movie · Metrics · Trajectories</td></tr>
  <tr><td>Is a result reliable over many seeds?</td><td>Scans</td></tr>
  <tr><td>What does each cell look like by itself, with its loci?</td><td>Cells</td></tr>
  <tr><td>How is each locus found, and how often is it wrong?</td><td>Spot finding</td></tr>
  <tr><td>What does the image error X do, and how is it fixed?</td><td>Images</td></tr>
  <tr><td>How good is the instrument calibration?</td><td>Calibration</td></tr>
  <tr><td>Do z-piezo errors matter in 3D?</td><td>Stage &amp; fiducials</td></tr>
  <tr><td>How does distance control the MS2 spot?</td><td>Transcription</td></tr>
  <tr><td>How do probes set the locus brightness?</td><td>Probes</td></tr>
  <tr><td>When did Pol II start? Where are the bursts?</td><td>MS2 → Pol II</td></tr>
  <tr><td>Can the E–P curve be recovered, and how well?</td><td>E–P curve</td></tr>
  <tr><td>Which saved simulations feed a movie?</td><td>Library · Chromatin</td></tr>
  <tr><td>What exactly was this run?</td><td>Config · the summary on top</td></tr></table>`), [
  ["Here is a cheat sheet: the question, and the page that answers it.", (r) => r.to(0)],
  ["For tracking questions go to Movie, Metrics and Trajectories, and to Scans for reliability.", (r) => r.to(2)],
  ["For loci go to Cells and Spot finding. For microscope errors go to Images, Calibration and Stage and fiducials.", (r) => r.to(7)],
  ["For biology go to Transcription, Probes, M S two to Pol two, and E P curve. For the saved polymer simulations, Library and Chromatin. For what exactly a run was, Config.", (r) => r.to(14)],
  ["That is the whole program. Start the dashboard, press Run, and work through the pages with the Your Turn steps. Welcome aboard.", (r) => r.to(14)],
]);
