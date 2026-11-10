# Sentences of the tutorial that state a number

Made by `scripts\tutorial\update.cmd numbers`. When results or defaults change, re-read these against the dashboard and the docs;
edit the sentence in `scenes_a.js` / `scenes_b.js` (search for its first words), then run `scripts\tutorial\update.cmd all`.


## pipeline
- `pipeline:0` The program is a pipeline of seven stages. Each stage is its own folder of code, and each one writes files that the next one reads, so you can inspect or replace any stage on its own.
- `pipeline:1` Stage one is chromatin. A polymer simulation on the graphics card folds a ten megabase region and moves it, with loop extrusion and boundary elements. It is made once, saved, and reused. A simple stand-in motion exists for quick tests.

## layout
- `layout:5` Below are the headline tracking scores. The first row uses automatic segmentation, which is realistic. The second uses perfect masks, which isolates just the linking error. One point zero is perfect for both.

## form_sim
- `form_sim:0` The first sidebar form is called New simulation. It runs the tracking-only experiment: moving nuclei, a basic microscope, segmentation and tracking. It takes about twenty seconds.

## form_loci
- `form_loci:0` The second form is the big one: New single-cell loci run. It adds two coloured loci to every nucleus, isolates each cell, and locates the loci in each one. It takes about two minutes, or up to ten with 3D.
- `form_loci:2` Cell speed scale of zero point four keeps tracking easy, so that locus analysis is not drowned by tracking errors. Photons per locus per second set how bright each spot is. Division probability lets cells split.

## movie
- `movie:6` Red circles mark identity switches. Here the tracker swapped two cells that touch. Their numbers on the right no longer match their true colours.

## biplane
- `biplane:0` The Vutara can image in two ways. In its widefield mode it records one focal plane at a time. In its biplane mode it records two planes at the same moment.
- `biplane:1` A beam splitter sends half of the light to each half of one camera, so two planes six hundred nanometres apart in the sample appear side by side. The height of a spot comes from how bright it is in the two planes.
- `biplane:2` Two planes cover only a thin slice, and the method works best near or between them. A nucleus is about ten microns thick, so you still step through it, as in an ordinary z-stack.
- `biplane:3` Each stage position gives two planes, so a step of one point two microns covers the sample evenly every point six microns, with half as many stage moves. The price is that each plane gets half of the light.
- `biplane:5` Beads in the image sit between the lowest pair of planes, so they appear in both. The difference of their two positions measures the shift in every frame, and the program undoes it. In a test it measured nought point three two and minus nought point two one pixels for a true nought point three five and minus nought point two five.
- `biplane:7` Your turn: in the single-cell loci form, choose Bruker Vutara VXL biplane 3D mode, keep the field at forty microns, and run a short movie. Then read the Spot finding and Stage pages.

## calibration
- `calibration:6` And finally, does it help? The table compares the movie analysed without and with the calibration. In this run the median locus error fell from about forty-one to about twenty-seven nanometres.

## polymer
- `polymer:0` Stage one is the chromatin. A polymer simulation of ten thousand monomers, one kilobase each, so ten megabases, run on the graphics card with polychrom and OpenMM, inside a confining sphere.
- `polymer:1` Two beads are tagged: the promoter and the enhancer, a hundred kilobases apart. How far apart they are, moment by moment, is what the transcription model needs.
- `polymer:4` We checked that the engine reproduces known behaviour, in a validation report. Eleven of twelve checks pass. The twelfth, insulation by the boundaries, is real but weaker than published, and the report says so honestly.
- `polymer:5` A simulation only has polymer units. The calibration converts them: fifty-five nanometres per unit and nought point four six seconds per saved block. These come from two assumptions: a hundred kilobase pair is within four hundred nanometres ten percent of the time, and a single locus has a diffusion coefficient of two thousandths of a square micron per second. Replace them with a measured mean squared displacement.

## library
- `library:2` One selector at the top chooses the library for the whole page. The live status refreshes every five seconds: simulations finished, the one running now and its time left, and the graphics card's load. If the computer sleeps it says Stalled.

## copies
- `copies:4` A path is kept only if it looks like a locus: there in most frames, about as bright as the best path, and as slow. Paths through noise and strays fail that test. Then each promoter path is paired with an enhancer path that stays within one and a half microns of it.
- `copies:5` When two copies come close, two paths can swap identities. Exchanging their tails wherever that makes the steps shorter undoes it, and a pair prior anchors each promoter to its own enhancer. Over twenty-five movies the copy number was right in every cell without strays, and in nine of ten cells with four strays per minute.

## tx
- `tx:3` Step three sets how the distance changes a rate: a sharp contact, a smooth Hill curve, a smooth exponential, or none as a control. The default is the Hill curve, half its effect at four hundred nanometres.

## ms2
- `ms2:7` Five B checks the sampler with four chains from different starts, and the usual convergence numbers. The default sampler is chosen per cell, and it switches to parallel tempering when the plain one fails.

## ep_anim
- `ep_anim:4` Wrong spot picks lower the contrast further. A few percent of wrong picks matter much more than modest localization noise.

## ep
- `ep:1` Part one chooses the microscope, the number of cells, the delay between the two colours, the photons, the focus scatter, the wrong picks, and whether the single locus diffusion refers to a locus or to the distance. Then it simulates and fits, in about forty seconds.
- `ep:2` The result row compares the true contrast, the probability of a burst at one hundred fifty nanometres over that at a thousand, with the recovered one, and says whether the interval contains the truth. Compare the two before trusting an analysis.
- `ep:3` The promoter integrates the distance over about forty-five seconds, while the distance fluctuates within ten, so the curve is fitted from the distance history, not the distance now.
