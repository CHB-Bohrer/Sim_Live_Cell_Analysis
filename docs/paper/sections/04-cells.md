<!--
Section: Cells (stage 2)
Update when: the motion model, collisions, division or the nuclear shape model change.
Code: src/simlive/stage2_cells
Check: none
Sources: src/simlive/stage2_cells/motion.py, docs/notes/status_detail.md.
-->

# Cells: nuclei that move, touch and divide

::: {.box}
**Short version.** Each nucleus does a persistent random walk, bounces off its neighbours and the walls, slowly changes shape, and may divide. Chromatin is added afterwards in the nucleus's own frame, so cell motion never changes locus dynamics.
:::

Per frame, a cell takes a diffusive step (diffusion coefficient $D$, so the step has standard deviation $\sqrt{2 D\,\Delta t}$ per axis) plus a directed drift whose heading itself diffuses. Cells collide as hard circles and reflect at the walls. A margin outside the field of view lets cells leave and enter.

A nucleus is an ellipse whose outline is perturbed by a few low-order harmonics,

$$\rho(\theta) = 1 + \sum_k a_k \cos(k\theta + \phi_k),$$

and the aspect ratio, orientation, amplitudes $a_k$ and phases $\phi_k$ drift slowly, so nuclei rotate and deform. In 3D runs each cell also has a height and a vertical half-axis, drawn from their own random stream so the lateral simulation is unchanged.

Division follows the Cell Tracking Challenge convention: the parent track ends, and two daughters with new IDs and a recorded parent start in the next frame.

`speed_scale` multiplies all displacements. It is the main knob for how hard tracking is (Section 7). The default is deliberately hard; the real cell speed and frame interval have not been supplied **[Assumption]**.
