# Plan (not built): use all three signals to choose spots, and tell apart loci on different chromosomes

Asked by the user 2026-10-11 ("use all three loci to select for good spots, and separate the loci from different chromosomes");
written down for later. Nothing here is implemented yet.

## Where things stand
- Each colour is found on its own: one Viterbi path per cell for the promoter (locus 0), one for the enhancer (locus 1)
  (`stage6_analysis/dog_viterbi.py`).
- The MS2 signal is read at the tracked promoter (`ms2.py`, `ms2_position: promoter`). That is already one way the signals are used
  together: the promoter's own colour places the MS2 measurement.
- The simulation has ONE copy of each locus per cell. Diploid cells have two (four after replication: sister chromatids), so the
  question "which promoter belongs with which enhancer" never arises in the simulated data yet.

## A. Choose spots jointly (one chromosome copy)
1. **A pair term in the path.** Run one Viterbi over candidate PAIRS (promoter candidate i, enhancer candidate j), at most 5 x 5 = 25
   states per frame, with each spot's own reward and link cost plus a pair cost -log p(d_ij). Here d_ij is the promoter-enhancer
   distance (corrected for the chromatic shift) and p(d) is the distance distribution.
   - p(d) can be learnt from frames where both spots are unambiguous (one strong candidate each), so it is not assumed.
   - The pair term rejects a ghost / hot pixel / cosmic ray in one colour if it sits microns away from the other colour's spot.
   - Caution: this biases the very distance the experiment measures. Score it on simulations whose true distance distribution differs
     from the prior (e.g. coupling on / off), and report results with and without the pair term.
2. **MS2 as a third witness.** When the MS2 spot is bright (ON frames), its candidate must coincide with the promoter's within the
   registration error (~50-100 nm after calibration). Add a reward for a promoter candidate that has an MS2 candidate on top of it in
   the nuclear channel. That confirms the promoter during bursts, and it is independent evidence because it is a different colour and
   fluorophore.
3. **Quality flags per frame:**
   - only one colour found;
   - the pair distance is in the tail of p(d);
   - MS2 is ON but there is no promoter candidate under it;
   - the cell's path re-acquired the spot after a long gap.

   Downstream analyses can drop flagged frames, and the simulation can measure how much each flag lowers the error rate.
4. **Validation:** wrong-pick rate, median error and the bias of the measured d(t) distribution, against the truth, with each piece
   on/off. Use more than 5 seeds and the image-error sources that make single-colour picks fail: hot pixels, cosmic rays, cross-talk.

## B. Two (or four) copies: separate homologs and sisters
1. **Simulation first.**
   - Two polymer copies per cell: two library simulations in separate territories. The distance between homologous territories is
     typically several microns; take the value from the user's cells.
   - Promoter and enhancer of the same copy keep their polymer (cis) relationship.
   - MS2 is on one allele or on both: the user's construct decides.
   - Optionally S/G2: sister chromatids, cis-like pairs that are close (< ~1 um) and cohesed.
   - Truth: copy ID of every spot.
2. **Analysis.**
   - Per colour, track TWO spots per cell: two disjoint Viterbi paths (k-best disjoint paths / min-cost flow) or one Viterbi over
     candidate pairs with an exclusion rule.
   - Then pair promoter tracks with enhancer tracks: a Hungarian assignment on the time-averaged cis distance. Cis pairs stay within
     ~1 um; trans pairs are as far apart as the territories.
   - Flag the cell as ambiguous when the two copies come closer than ~2 x the typical cis distance.
   - Allele-specific MS2: assign each MS2 spot to the promoter it sits on.
3. **Validation:**
   - pairing error rate and identity swaps between copies, as a function of homolog separation and locus mobility;
   - how often a cell must be flagged;
   - how much a wrong pairing biases d(t).

## Needed from the user
- Copy number of the labelled loci in their cells (diploid? aneuploid line?), and whether both alleles carry the labels and the MS2 cassette.
- Typical distance between the two homologous copies, and whether cells are imaged in S/G2 (sister chromatids).
- Whether the downstream question uses both alleles or one per cell.
