/* Concept animations (canvas, logical size 1280 x 720). Each factory returns a controller whose methods are the "cues" used by scenes.js. */
"use strict";
const A = {};
const COLS = ["#60a5fa", "#f472b6", "#fbbf24", "#4ade80", "#a78bfa", "#fb923c", "#2dd4bf", "#f87171", "#c084fc", "#facc15"];

function panel(c, x, y, w, h_, title, col) {
  rrect(c, x, y, w, h_, 16, "#12203b", col || PAL.line, 2);
  if (title) txt(c, title, x + 18, y + 34, { size: 19, weight: 700, color: col && col !== PAL.line ? col : PAL.muted });
}
function show(v) { return clamp(v, 0, 1); }

/* ================================================================ 1. the validation idea */
A.truthloop = (stage) => {
  const S = { step: new Tween(0) };
  const r = rng(5);
  const path = (t, k) => [Math.sin(t * 0.7 + k) * 40 + Math.sin(t * 1.9 + k * 2) * 14, Math.cos(t * 0.6 + k * 1.3) * 28 + Math.sin(t * 2.3 + k) * 10];
  const cv = new Canvas(stage, (c, t) => {
    const st = S.step.step(0.016, 4), a = (n) => show(st - n + 0.4);
    txt(c, "Simulate a movie where the answer is known, then test the analysis on it", 640, 54, { size: 30, weight: 700, align: "center" });
    // left: simulation
    panel(c, 40, 100, 380, 520, "1 · SIMULATION", PAL.teal);
    c.save(); c.globalAlpha = a(0.5) || 0.0; c.beginPath(); c.ellipse(230, 330, 130, 100, 0, 0, 7); c.fillStyle = "#1c2f52"; c.fill(); c.strokeStyle = PAL.blue; c.lineWidth = 3; c.stroke();
    const p0 = path(t, 0), p1 = path(t, 1.7);
    const P = [230 + p0[0], 330 + p0[1]], E = [230 + p1[0] * 0.8 + 16, 330 + p1[1] * 0.8];
    circ(c, P[0], P[1], 9, PAL.purple); circ(c, E[0], E[1], 9, PAL.amber);
    c.restore();
    txt(c, "nucleus with a promoter (purple)\nand an enhancer (orange)", 230, 480, { size: 17, align: "center", color: PAL.muted, alpha: a(0.5) });
    rrect(c, 70, 520, 320, 70, 12, "#0e2a2a", PAL.teal); txt(c, "ANSWER KEY\ntrue positions, saved", 230, 548, { size: 18, align: "center", weight: 700, color: PAL.teal, alpha: a(0.5) });
    // middle: movie + analysis
    c.globalAlpha = a(1.4);
    arrow(c, 430, 330, 470, 330, PAL.muted);
    panel(c, 480, 100, 380, 520, "2 · ANALYSIS (never sees the key)", PAL.blue);
    rrect(c, 520, 160, 300, 220, 10, "#05080f");
    for (let i = 0; i < 900; i++) { const x = 520 + r() * 300, y = 160 + r() * 220, g = r(); c.fillStyle = `rgba(160,180,210,${g * 0.35})`; c.fillRect(x, y, 3, 3); }
    c.beginPath(); c.ellipse(670, 270, 90, 70, 0, 0, 7); c.fillStyle = "rgba(120,150,200,.18)"; c.fill();
    circ(c, 670 + (P[0] - 230) * 0.7, 270 + (P[1] - 330) * 0.7, 5, "rgba(255,255,255,.85)"); circ(c, 670 + (E[0] - 230) * 0.7, 270 + (E[1] - 330) * 0.7, 5, "rgba(255,255,255,.7)");
    txt(c, "noisy recorded movie", 670, 405, { size: 17, align: "center", color: PAL.muted });
    txt(c, "find cells → follow them →\nlocate the loci", 670, 460, { size: 20, align: "center", weight: 600 });
    c.globalAlpha = 1;
    // estimates
    c.globalAlpha = a(2.4);
    const jx = (k) => Math.sin(t * 5 + k) * 6 + Math.sin(t * 11 + k * 3) * 3, jy = (k) => Math.cos(t * 4.3 + k) * 6;
    const eP = [P[0] + jx(1) + 3, P[1] + jy(1)], eE = [E[0] + jx(2) - 4, E[1] + jy(2)];
    rrect(c, 520, 520, 300, 70, 12, "#0c1d3a", PAL.blue); txt(c, "RESULT\nestimated positions", 670, 548, { size: 18, align: "center", weight: 700, color: PAL.blue });
    c.globalAlpha = 1;
    c.globalAlpha = a(3.0);
    arrow(c, 870, 330, 910, 330, PAL.muted);
    panel(c, 920, 100, 320, 520, "3 · COMPARE", PAL.amber);
    c.beginPath(); c.ellipse(1080, 300, 120, 90, 0, 0, 7); c.fillStyle = "#1c2f52"; c.fill();
    const q = (p) => [1080 + (p[0] - 230) * 0.9, 300 + (p[1] - 330) * 0.9];
    for (const [tr, es, col] of [[P, eP, PAL.purple], [E, eE, PAL.amber]]) {
      const a1 = q(tr), b1 = q(es); circ(c, a1[0], a1[1], 8, null, "#67e8f9", 3); line(c, b1[0] - 8, b1[1] - 8, b1[0] + 8, b1[1] + 8, PAL.red, 3); line(c, b1[0] - 8, b1[1] + 8, b1[0] + 8, b1[1] - 8, PAL.red, 3); line(c, a1[0], a1[1], b1[0], b1[1], col, 2);
    }
    txt(c, "○ truth   ✕ analysis", 1080, 430, { size: 18, align: "center", color: PAL.muted });
    txt(c, "error = distance between them\nmeasured, never guessed", 1080, 480, { size: 20, align: "center", weight: 600 });
    rrect(c, 950, 530, 260, 60, 12, "#2b2410", PAL.amber); txt(c, "e.g. median error 26 nm", 1080, 568, { size: 20, align: "center", weight: 700, color: PAL.amber });
    c.globalAlpha = 1;
  });
  return { instant: false, step(n) { S.step.set(n, this.instant); return this; }, destroy() { cv.destroy(); } };
};

/* ================================================================ 2. the pipeline */
A.pipeline = (stage) => {
  const S = { k: new Tween(0), cur: 0 };
  const boxes = [
    { n: "1", t: "Chromatin", s: "polymer simulations\n(GPU library)", x: 50, y: 110, tag: "built · toy stand-in too" },
    { n: "2", t: "Cells", s: "motion, collisions,\nshape, division", x: 360, y: 110, tag: "built" },
    { n: "3", t: "Microscope", s: "recorded images,\n24 error sources", x: 670, y: 110, tag: "built" },
    { n: "4", t: "Segment, track", s: "Cellpose, Trackastra", x: 980, y: 110, tag: "built" },
    { n: "5", t: "Cut out cells", s: "one small movie\nper cell", x: 205, y: 390, tag: "built" },
    { n: "6", t: "Locate loci", s: "spot finder, MS2,\nPol II deconvolution", x: 515, y: 390, tag: "built · yours plugs in" },
    { n: "7", t: "Validate", s: "compare with the\nanswer key", x: 825, y: 390, tag: "built" },
  ];
  const detail = [
    "", "Stage 1 simulates how the chromatin folds and moves. Each simulation is saved once, then reused. Output: positions of every monomer over time.",
    "Stage 2 moves the nuclei around the field, bumps them, changes their shape and lets some divide. Every cell keeps a true ID number.",
    "Stage 3 turns the cells and loci into camera images: blur, noise, drift, haze and the rest. Output: the movie, plus the answer key.",
    "Stage 4 finds every nucleus in every frame (segmentation) and links them through time (tracking). It never sees the answer key.",
    "Stage 5 cuts one small movie per cell, so each cell can be analysed by itself, in parallel.",
    "Stage 6 finds the promoter and the enhancer in every frame, reads the MS2 spot, and turns it into Pol II starts. This is where your own analysis plugs in.",
    "Stage 7 compares everything with the answer key: tracking scores, locus error, distance error, and how much each microscope error costs.",
  ];
  const cv = new Canvas(stage, (c, t) => {
    const k = S.k.step(0.016, 5);
    txt(c, "The pipeline: seven stages, each writes files the next one reads", 640, 52, { size: 30, weight: 700, align: "center" });
    const cen = (b) => [b.x + 125, b.y + 80];
    const links = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6]];
    for (const [a, b] of links) {
      const A_ = boxes[a], B_ = boxes[b]; let x1, y1, x2, y2;
      if (a === 3) { x1 = A_.x + 125; y1 = A_.y + 160; x2 = B_.x + 125; y2 = B_.y; } else { x1 = A_.x + 250; y1 = A_.y + 80; x2 = B_.x; y2 = B_.y + 80; }
      if (a === 3) { line(c, x1, y1, x1, 300, PAL.grey, 3); line(c, x1, 300, x2, 300, PAL.grey, 3); arrow(c, x2, 300, x2, y2, PAL.grey, 3); }
      else arrow(c, x1 + 4, y1, x2 - 6, y2, PAL.grey, 3);
    }
    // answer key
    c.setLineDash([8, 6]); c.strokeStyle = PAL.teal; c.lineWidth = 2.5; c.beginPath(); c.moveTo(670, 270); c.lineTo(670, 340); c.lineTo(950, 340); c.lineTo(950, 390); c.stroke(); c.setLineDash([]);
    txt(c, "answer key (true IDs and positions)", 690, 330, { size: 16, color: PAL.teal });
    boxes.forEach((b, i) => {
      const on = Math.abs(k - (i + 1)) < 0.5, act = show(1 - Math.abs(k - (i + 1)));
      rrect(c, b.x, b.y, 250, 160, 18, on ? "#14355a" : "#12203b", on ? PAL.teal : PAL.line, on ? 4 : 2);
      circ(c, b.x + 34, b.y + 36, 18, on ? PAL.teal : "#27406b"); txt(c, b.n, b.x + 34, b.y + 44, { size: 22, weight: 800, align: "center", color: on ? "#052" : PAL.text });
      txt(c, b.t, b.x + 64, b.y + 44, { size: 22, weight: 700 });
      txt(c, b.s, b.x + 22, b.y + 96, { size: 19, color: PAL.muted });
      txt(c, b.tag, b.x + 22, b.y + 146, { size: 14, color: b.tag.includes("stand-in") || b.tag.includes("yours") ? PAL.amber : PAL.green });
      if (act > 0.05) { c.globalAlpha = act * 0.35; rrect(c, b.x - 6, b.y - 6, 262, 172, 22, null, PAL.teal, 6); c.globalAlpha = 1; }
    });
    // moving packet
    const seq = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6]]; const ph = (t * 0.5) % seq.length; const [a, b] = seq[Math.floor(ph)], f = ph % 1;
    const pa = cen(boxes[a]), pb = cen(boxes[b]); circ(c, lerp(pa[0], pb[0], f), lerp(pa[1], pb[1], f), 7, PAL.amber);
    if (S.cur > 0) {
      rrect(c, 40, 590, 1200, 100, 14, "#0c1830", PAL.teal);
      txt(c, detail[S.cur], 66, 626, { size: 21, color: PAL.text, lh: 1.3, w: 1130 });
    }
  });
  return { instant: false, stage(n) { S.cur = n; S.k.set(n, this.instant); return this; }, destroy() { cv.destroy(); } };
};

/* ================================================================ 3. cells moving, tracker swapping IDs */
A.tracking = (stage) => {
  const r = rng(11), N = 8;
  const cells = Array.from({ length: N }, (_, i) => ({ x: 0.12 + (i % 4) * 0.25 + r() * 0.06, y: 0.2 + Math.floor(i / 4) * 0.5 + r() * 0.1, vx: 0, vy: 0, lab: i, cool: 0, ang: r() * 6, sh: 0.5 + r() * 0.3 }));
  const S = { panels: 1, speed: 1, swaps: false, n: 0, flash: [], crowd: 0 };
  const cv = new Canvas(stage, (c, t, dt) => {
    txt(c, "Stage 2 moves the cells. Stage 4 has to follow them.", 640, 48, { size: 28, weight: 700, align: "center" });
    const drawP = (ox, title, tracker) => {
      panel(c, ox, 80, 580, 560, title, tracker ? PAL.blue : PAL.teal);
      for (let i = 0; i < N; i++) {
        const q = cells[i], x = ox + 30 + q.x * 520, y = 120 + q.y * 480, col = COLS[(tracker ? q.lab : i) % COLS.length];
        c.save(); c.translate(x, y); c.rotate(q.ang + t * 0.05);
        c.beginPath(); c.ellipse(0, 0, 46, 46 * (0.7 + 0.15 * Math.sin(t * 0.7 + i)), 0, 0, 7); c.fillStyle = col + "cc"; c.fill(); c.strokeStyle = col; c.lineWidth = 3; c.stroke(); c.restore();
        txt(c, String(tracker ? q.lab : i), x, y + 8, { size: 24, weight: 800, align: "center", color: "#08111f" });
      }
      for (const f of S.flash) if (tracker) { const a = f.life; c.globalAlpha = clamp(a, 0, 1); circ(c, ox + 30 + f.x * 520, 120 + f.y * 480, 54 - 20 * a, null, PAL.red, 5); c.globalAlpha = 1; }
    };
    drawP(40, "TRUE cells (ground truth IDs)", false);
    c.globalAlpha = show(S.panels - 1 + 0.001) > 0 ? 1 : 0; if (S.panels > 1) drawP(660, "TRACKED cells (tracker's IDs)", true); c.globalAlpha = 1;
    if (S.panels > 1) { rrect(c, 840, 650, 220, 52, 12, "#2a1414", PAL.red); txt(c, "ID switches: " + S.n, 950, 685, { size: 22, weight: 800, align: "center", color: PAL.red }); }
    txt(c, "speed ×" + S.speed.toFixed(1), 70, 685, { size: 20, color: PAL.muted });
  });
  cv.step = (dt) => {
    const sp = S.speed;
    for (let i = 0; i < N; i++) {
      const q = cells[i]; q.vx += (r() - 0.5) * 0.35 * dt * sp; q.vy += (r() - 0.5) * 0.35 * dt * sp; q.vx *= 0.97; q.vy *= 0.97;
      const m = Math.hypot(q.vx, q.vy), lim = 0.06 * sp; if (m > lim) { q.vx *= lim / m; q.vy *= lim / m; }
      q.cool -= dt;
    }
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
      const a = cells[i], b = cells[j], dx = b.x - a.x, dy = (b.y - a.y) * 0.86, d = Math.hypot(dx, dy) || 1e-3, min = 0.19;
      if (d < min) { const f = (min - d) * 0.5; a.vx -= dx / d * f * dt * 8; a.vy -= dy / d * f * dt * 8; b.vx += dx / d * f * dt * 8; b.vy += dy / d * f * dt * 8; }
      if (S.swaps && d < min * 0.97 && a.cool <= 0 && b.cool <= 0 && S.panels > 1) {
        [a.lab, b.lab] = [b.lab, a.lab]; a.cool = b.cool = 3; S.n++; S.flash.push({ x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, life: 1.5 });
      }
    }
    for (const q of cells) { q.x += q.vx * dt; q.y += q.vy * dt; if (q.x < 0.07 || q.x > 0.93) q.vx *= -1; if (q.y < 0.1 || q.y > 0.92) q.vy *= -1; q.x = clamp(q.x, 0.05, 0.95); q.y = clamp(q.y, 0.08, 0.94); }
    S.flash.forEach((f) => (f.life -= dt * 0.5)); S.flash = S.flash.filter((f) => f.life > 0);
  };
  return {
    instant: false, panels(n) { S.panels = n; return this; }, speed(v) { S.speed = v; return this; },
    swaps(on) { S.swaps = on; return this; }, crowd(v) { if (v) cells.forEach((q, i) => { q.x = 0.3 + (i % 4) * 0.13; q.y = 0.3 + Math.floor(i / 4) * 0.2; }); return this; },
    destroy() { cv.destroy(); },
  };
};

/* ================================================================ 4. forward model: a spot becomes an image */
A.microscope = (stage) => {
  const N = 41, S = { step: 0, accum: 0 }, r0 = rng(3);
  const off = document.createElement("canvas"); off.width = off.height = N; const octx = off.getContext("2d"), img = octx.createImageData(N, N);
  const tx = 20.3, ty = 19.6;
  const hot = [[8, 31], [33, 10]], cosmic = [[30, 33], [31, 33], [32, 34]];
  const titles = ["Ideal: a point of light, where the locus really is", "Diffraction: the optics blur the point (PSF)", "Plus nuclear haze from the other fluorescent protein",
    "Plus shot noise: light arrives as photons", "Plus camera read noise", "Plus hot pixels and a cosmic ray", "Plus stage drift and colour shift: drawn in the wrong place"];
  let rr = rng(9), lastNoise = -1; let frame = new Float32Array(N * N), brightest = [0, 0];
  const build = (st) => {
    const sh = st >= 6 ? [1.6, -1.2] : [0, 0], sg = 1.7;
    for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
      let v = 0;
      if (st === 0) { v = (Math.abs(x - Math.round(tx + sh[0])) < 1 && Math.abs(y - Math.round(ty + sh[1])) < 1) ? 255 : 0; }
      else { const d2 = (x - tx - sh[0]) ** 2 + (y - ty - sh[1]) ** 2; v = 255 * Math.exp(-d2 / (2 * sg * sg)); }
      if (st >= 2) v += 55 * Math.exp(-(((x - 20) ** 2 + (y - 20) ** 2) / (2 * 14 * 14))) + 22;
      if (st >= 3) v += gauss(rr) * Math.sqrt(Math.max(v, 1)) * 1.1;
      if (st >= 4) v += gauss(rr) * 9;
      frame[y * N + x] = v;
    }
    if (st >= 5) { for (const [x, y] of hot) frame[y * N + x] = 330; for (const [x, y] of cosmic) frame[y * N + x] = 420; }
    let b = -1; for (let i = 0; i < N * N; i++) if (frame[i] > b) { b = frame[i]; brightest = [i % N, Math.floor(i / N)]; }
  };
  const cv = new Canvas(stage, (c, t, dt) => {
    S.accum += dt; if (S.step >= 3 && S.accum > 0.35) { S.accum = 0; rr = rng((Math.random() * 1e9) | 0); build(S.step); } else if (lastNoise !== S.step) { lastNoise = S.step; rr = rng(9); build(S.step); }
    for (let i = 0; i < N * N; i++) { const v = clamp(frame[i] / 330, 0, 1), k = i * 4; const g = Math.pow(v, 0.8); img.data[k] = 30 + 225 * g; img.data[k + 1] = 40 + 215 * g; img.data[k + 2] = 60 + 195 * g; img.data[k + 3] = 255; }
    octx.putImageData(img, 0, 0);
    txt(c, titles[S.step], 640, 50, { size: 28, weight: 700, align: "center" });
    const X0 = 70, Y0 = 90, Z = 14;
    c.imageSmoothingEnabled = false; c.drawImage(off, X0, Y0, N * Z, N * Z);
    rrect(c, X0 - 3, Y0 - 3, N * Z + 6, N * Z + 6, 4, null, PAL.line);
    const px = (x) => X0 + (x + 0.5) * Z, py = (y) => Y0 + (y + 0.5) * Z;
    const sh = S.step >= 6 ? [1.6, -1.2] : [0, 0];
    circ(c, px(tx), py(ty), 22, null, "#67e8f9", 3);
    txt(c, "○ true locus position", X0, Y0 + N * Z + 34, { size: 19, color: "#67e8f9" });
    if (S.step >= 1) { const b = brightest; const col = Math.hypot(b[0] - tx, b[1] - ty) > 5 ? PAL.red : PAL.green; line(c, px(b[0]) - 14, py(b[1]) - 14, px(b[0]) + 14, py(b[1]) + 14, col, 4); line(c, px(b[0]) - 14, py(b[1]) + 14, px(b[0]) + 14, py(b[1]) - 14, col, 4);
      txt(c, "✕ brightest pixel (naive detector)", X0 + 290, Y0 + N * Z + 34, { size: 19, color: col }); }
    // profile through row ty
    const PX = 740, PY = 140, PW = 480, PH = 300;
    panel(c, PX - 30, PY - 50, PW + 60, PH + 190, "Brightness along the row through the locus", PAL.line);
    line(c, PX, PY + PH, PX + PW, PY + PH, PAL.grey, 2); line(c, PX, PY, PX, PY + PH, PAL.grey, 2);
    c.beginPath(); const row = Math.round(ty);
    for (let x = 0; x < N; x++) { const v = clamp(frame[row * N + x] / 440, 0, 1); const X = PX + (x / (N - 1)) * PW, Y = PY + PH - v * PH; x ? c.lineTo(X, Y) : c.moveTo(X, Y); }
    c.strokeStyle = PAL.teal; c.lineWidth = 3; c.stroke();
    const notes = ["", "The blur is a property of the lens: wider for longer wavelengths and lower NA.", "Everything the cell glows with in this colour adds a haze under the spot.", "Dim spots are the noisiest relative to their brightness.", "Every pixel gets a little electronic noise as well.", "A hot pixel can be brighter than the real spot, so a naive detector locks onto it.", "The spot you analyse is not where the locus is: chromatic shift and drift are errors too."];
    if (S.step > 0) txt(c, notes[S.step], PX, PY + PH + 50, { size: 19, color: PAL.muted, lh: 1.35, w: 470 });
    txt(c, `${S.step + 1} of 7`, PX + PW, PY + PH + 140, { size: 18, color: PAL.grey, align: "right" });
  });
  return { instant: false, step(n) { S.step = n; return this; }, destroy() { cv.destroy(); } };
};

/* ================================================================ 5. polymer + loop extrusion */
A.polymer = (stage) => {
  const N = 72, R = 250, CX = 330, CY = 400, r = rng(21), L = 19, BR = 8;
  const P = Array.from({ length: N }, (_, i) => { const a = i / N * 6.283 * 2.2, rad = 30 + i * 2.4; return [CX + Math.cos(a) * rad * 0.9, CY + Math.sin(a) * rad * 0.9]; });
  const S = { pe: false, ext: false, bnd: false, scale: false, a: 36, b: 37, tick: 0, trace: [], dtr: 0 };
  const ip = 21, ie = 52, bL = 17, bR = 56;
  const step = (dt) => {
    const h_ = 0.004, steps = Math.max(1, Math.round(dt / 0.016 * 4));
    for (let s = 0; s < steps; s++) {
      const F = P.map(() => [0, 0]);
      for (let i = 0; i < N - 1; i++) { const dx = P[i + 1][0] - P[i][0], dy = P[i + 1][1] - P[i][1], d = Math.hypot(dx, dy) || 1, f = 60 * (d - L); F[i][0] += dx / d * f; F[i][1] += dy / d * f; F[i + 1][0] -= dx / d * f; F[i + 1][1] -= dy / d * f; }
      for (let i = 0; i < N - 2; i++) { const dx = P[i + 2][0] - P[i][0], dy = P[i + 2][1] - P[i][1], d = Math.hypot(dx, dy) || 1, f = 6 * (d - 2 * L * 0.93); F[i][0] += dx / d * f; F[i][1] += dy / d * f; F[i + 2][0] -= dx / d * f; F[i + 2][1] -= dy / d * f; }
      if (S.ext) { const a = P[S.a], b = P[S.b], dx = b[0] - a[0], dy = b[1] - a[1], d = Math.hypot(dx, dy) || 1, f = 90 * (d - 12); F[S.a][0] += dx / d * f; F[S.a][1] += dy / d * f; F[S.b][0] -= dx / d * f; F[S.b][1] -= dy / d * f; }
      for (let i = 0; i < N; i++) for (let j = i + 2; j < N; j++) { const dx = P[j][0] - P[i][0], dy = P[j][1] - P[i][1], d2 = dx * dx + dy * dy; if (d2 < 27 * 27 && d2 > 0.01) { const d = Math.sqrt(d2), f = 40 * (27 - d); F[i][0] -= dx / d * f; F[i][1] -= dy / d * f; F[j][0] += dx / d * f; F[j][1] += dy / d * f; } }
      for (let i = 0; i < N; i++) {
        const [x, y] = P[i], dx = x - CX, dy = y - CY, d = Math.hypot(dx, dy); if (d > R - 12) { const f = 120 * (d - (R - 12)); F[i][0] -= dx / d * f; F[i][1] -= dy / d * f; }
        P[i][0] += F[i][0] * h_ + gauss(r) * 2.1; P[i][1] += F[i][1] * h_ + gauss(r) * 2.1;
      }
    }
  };
  let acc = 0, ext_t = 0;
  const cv = new Canvas(stage, (c, t) => {
    txt(c, "Stage 1: a chromatin fibre folding and moving inside a confining sphere", 640, 46, { size: 26, weight: 700, align: "center" });
    circ(c, CX, CY, R, "#0c162b", PAL.line, 3); txt(c, "simulation sphere (a territory, not the nucleus)", CX, CY + R + 28, { size: 15, color: PAL.grey, align: "center" });
    if (S.ext && S.bnd) { /* boundary markers */ }
    c.lineWidth = 7; c.lineCap = "round"; c.strokeStyle = "#5d86c9"; c.beginPath(); P.forEach((p, i) => (i ? c.lineTo(p[0], p[1]) : c.moveTo(p[0], p[1]))); c.stroke(); c.lineCap = "butt";
    if (S.ext) { const a = P[S.a], b = P[S.b]; c.lineWidth = 4; c.strokeStyle = PAL.green; c.beginPath(); for (let i = S.a; i <= S.b; i++) i === S.a ? c.moveTo(P[i][0], P[i][1]) : c.lineTo(P[i][0], P[i][1]); c.stroke();
      circ(c, a[0], a[1], 9, PAL.green, "#052e16", 2); circ(c, b[0], b[1], 9, PAL.green, "#052e16", 2); line(c, a[0], a[1], b[0], b[1], PAL.green, 2, [3, 3]); }
    if (S.bnd) for (const i of [bL, bR]) { const p = P[i]; c.beginPath(); c.moveTo(p[0], p[1] - 17); c.lineTo(p[0] - 11, p[1] + 5); c.lineTo(p[0] + 11, p[1] + 5); c.closePath(); c.fillStyle = PAL.red; c.fill(); }
    if (S.pe) { circ(c, P[ip][0], P[ip][1], 11, PAL.purple, "#fff", 2); circ(c, P[ie][0], P[ie][1], 11, PAL.amber, "#fff", 2); txt(c, "promoter", P[ip][0] - 16, P[ip][1] - 16, { size: 16, color: PAL.purple, align: "right", weight: 700 }); txt(c, "enhancer", P[ie][0] + 16, P[ie][1] - 14, { size: 16, color: PAL.amber, weight: 700 }); }
    if (S.ext) txt(c, "● loop-extruding factor (green): both legs move apart, growing a loop", 40, 100, { size: 16, color: PAL.green });
    if (S.bnd) txt(c, "▲ boundary element: the factor stops here", 40, 124, { size: 16, color: PAL.red });
    // trace
    const TX = 700, TY = 120, TW = 520, TH = 250;
    if (S.pe) {
      panel(c, TX - 30, TY - 50, TW + 60, TH + 100, "Promoter–enhancer distance over time", PAL.line);
      const dmax = 360; line(c, TX, TY + TH, TX + TW, TY + TH, PAL.grey, 2); line(c, TX, TY, TX, TY + TH, PAL.grey, 2);
      const thr = TY + TH - (95 / dmax) * TH; line(c, TX, thr, TX + TW, thr, PAL.amber, 2, [6, 5]); txt(c, "contact (400 nm)", TX + TW, thr - 8, { size: 15, color: PAL.amber, align: "right" });
      c.beginPath(); S.trace.forEach((v, i) => { const X = TX + (i / 300) * TW, Y = TY + TH - clamp(v / dmax, 0, 1) * TH; i ? c.lineTo(X, Y) : c.moveTo(X, Y); }); c.strokeStyle = PAL.teal; c.lineWidth = 3; c.stroke();
      txt(c, "time →", TX + TW, TY + TH + 28, { size: 15, color: PAL.grey, align: "right" });
      const d = Math.hypot(P[ie][0] - P[ip][0], P[ie][1] - P[ip][1]); txt(c, "now: " + Math.round(d * 1.2 * 1) + " units", TX, TY + TH + 28, { size: 16, color: PAL.muted });
    }
    if (S.scale) {
      rrect(c, 680, 520, 540, 150, 14, "#0c1830", PAL.teal);
      txt(c, "Calibration turns simulation units into real ones:", 700, 556, { size: 19, weight: 700, color: PAL.teal });
      txt(c, "1 monomer = 1 kb of chromatin\n55 nm per length unit, 0.46 s per saved block\n(assumed: 100 kb pair within 400 nm 10% of the time,\n locus D = 0.002 µm²/s)", 700, 586, { size: 17, color: PAL.text, lh: 1.3 });
    }
  });
  cv.step = (dt) => {
    acc += dt; step(dt);
    if (S.ext) { ext_t += dt; if (ext_t > 0.22) { ext_t = 0; const lo = S.bnd ? bL : 2, hi = S.bnd ? bR : N - 3; if (S.a > lo) S.a--; if (S.b < hi) S.b++; if (S.a <= lo && S.b >= hi && !S.bnd) { S.a = 36; S.b = 37; } } }
    if (acc > 0.1) { acc = 0; const d = Math.hypot(P[ie][0] - P[ip][0], P[ie][1] - P[ip][1]); S.trace.push(d); if (S.trace.length > 300) S.trace.shift(); }
  };
  return {
    instant: false, pe(on) { S.pe = on; return this; }, extrude(on) { S.ext = on; if (on) { S.a = 36; S.b = 37; } return this; }, boundaries(on) { S.bnd = on; return this; },
    scale(on) { S.scale = on; return this; }, destroy() { cv.destroy(); },
  };
};

/* ================================================================ 6. spot finding: candidates + Viterbi */
A.trellis = (stage) => {
  const T_ = 13, r = rng(41), S = { step: 0, prog: new Tween(0) };
  const X0 = 130, X1 = 1160, Y0 = 130, Y1 = 560;
  const truth = []; let y = 0.5; for (let i = 0; i < T_; i++) { y = clamp(y + (r() - 0.5) * 0.12, 0.3, 0.7); truth.push(y); }
  const cand = [];
  for (let i = 0; i < T_; i++) {
    const cs = [];
    if (i !== 7) cs.push({ y: truth[i], b: i === 3 || i === 9 ? 0.45 : 0.85 + r() * 0.1, real: true });
    const nf = 2 + Math.floor(r() * 2);
    for (let k = 0; k < nf; k++) cs.push({ y: clamp(0.1 + r() * 0.8, 0.08, 0.92), b: (i === 3 || i === 9) && k === 0 ? 0.8 : 0.25 + r() * 0.35, real: false });
    cand.push(cs);
  }
  const px = (i) => X0 + (i / (T_ - 1)) * (X1 - X0), py = (v) => Y0 + v * (Y1 - Y0);
  const bright = cand.map((cs) => cs.reduce((a, b) => (b.b > a.b ? b : a)));
  const cv = new Canvas(stage, (c, t) => {
    const p = S.prog.step(0.016, 3);
    txt(c, ["Every nucleus becomes a column of candidates per frame", "A low threshold keeps several candidates per frame", "Picking the brightest candidate each frame jumps around", "The Viterbi path asks for the most likely whole trajectory", "A frame with no spot is bridged, at a cost"][S.step] || "", 640, 48, { size: 27, weight: 700, align: "center" });
    panel(c, 60, 80, 1160, 520, null);
    c.save(); c.translate(86, 460); c.rotate(-Math.PI / 2); txt(c, "position in the cell →", 0, 0, { size: 16, color: PAL.grey }); c.restore();
    for (let i = 0; i < T_; i++) { line(c, px(i), Y0 - 20, px(i), Y1 + 20, "#1b2c4d", 1); txt(c, "t" + i, px(i), Y1 + 44, { size: 14, color: PAL.grey, align: "middle" === "x" ? "center" : "center" }); }
    for (let i = 0; i < T_; i++) cand[i].forEach((k) => { const rad = 5 + k.b * 15, a = S.step >= 1 ? 1 : 0.0; c.globalAlpha = a; circ(c, px(i), py(k.y), rad, S.step >= 1 ? "rgba(147,164,195,.35)" : null, "#93a4c3", 2); c.globalAlpha = 1; });
    if (S.step === 0) { for (let i = 0; i < T_; i++) txt(c, "frame", px(i), 118, { size: 13, color: PAL.grey, align: "center" }); }
    if (S.step >= 2 && S.step < 3) { const pts = bright.map((k, i) => [px(i), py(k.y)]); c.beginPath(); pts.forEach(([x, y_], i) => (i ? c.lineTo(x, y_) : c.moveTo(x, y_))); c.strokeStyle = PAL.red; c.lineWidth = 3.5; c.stroke(); pts.forEach(([x, y_], i) => circ(c, x, y_, 7, bright[i].real ? PAL.green : PAL.red)); txt(c, "wrong picks when a noisy bright spot beats the real one", 640, 640, { size: 20, color: PAL.red, align: "center" }); }
    if (S.step >= 3) {
      const pts = []; truth.forEach((v, i) => { if (i !== 7) pts.push([px(i), py(v), i]); });
      c.beginPath(); pts.forEach(([x, y_], i) => (i ? c.lineTo(x, y_) : c.moveTo(x, y_))); c.strokeStyle = PAL.green; c.lineWidth = 5; c.stroke();
      pts.forEach(([x, y_]) => circ(c, x, y_, 7, PAL.green));
      txt(c, "link cost = −log of the diffusion propagator (big jumps cost a lot)   ·   reward = log SNR, minus the clutter of false candidates", 640, 640, { size: 18, color: PAL.green, align: "center" });
    }
    if (S.step >= 4) { const a = [px(6), py(truth[6])], b = [px(8), py(truth[8])]; c.setLineDash([8, 6]); line(c, a[0], a[1], b[0], b[1], PAL.amber, 4); c.setLineDash([]); circ(c, px(7), py((truth[6] + truth[8]) / 2), 16, null, PAL.amber, 3); txt(c, "no candidate in frame 7:\nthe path pays a miss cost\nand carries on", px(7) + 26, py((truth[6] + truth[8]) / 2) - 70, { size: 17, color: PAL.amber, lh: 1.25 }); }
  });
  return { instant: false, step(n) { S.step = n; return this; }, destroy() { cv.destroy(); } };
};

/* ================================================================ 7. probes binding a locus */
A.probes = (stage) => {
  const r = rng(77), N = 26, S = { phase: 0, koff: 0.18, trace: [], acc: 0 };
  const ps = Array.from({ length: N }, () => ({ st: 0, x: 90 + r() * 560, y: 150 + r() * 380, vx: 0, vy: 0, s: r() }));
  const LX0 = 260, LX1 = 480, LY = 340;
  const cv = new Canvas(stage, (c, t) => {
    txt(c, "Probes land on a locus, slide along it, and fall off", 640, 46, { size: 27, weight: 700, align: "center" });
    c.beginPath(); c.ellipse(370, 340, 300, 220, 0, 0, 7); c.fillStyle = "#0c162b"; c.fill(); c.strokeStyle = PAL.line; c.lineWidth = 3; c.stroke();
    c.beginPath(); c.moveTo(LX0, LY); c.lineTo(LX1, LY); c.strokeStyle = PAL.purple; c.lineWidth = 12; c.lineCap = "round"; c.stroke(); c.lineCap = "butt"; txt(c, "locus", (LX0 + LX1) / 2, LY + 34, { size: 17, color: PAL.purple, align: "center", weight: 700 });
    const col = ["#6b7fa8", PAL.amber, PAL.green];
    for (const p of ps) { circ(c, p.x, p.y, p.st ? 8 : 6, col[p.st]); }
    // counts
    const cnt = [0, 0, 0]; ps.forEach((p) => cnt[p.st]++);
    const BX = 740; txt(c, "attached probes", BX, 120, { size: 20, weight: 700 });
    [["FREE", "#6b7fa8", 0], ["BOUND", PAL.amber, 1], ["SCANNING", PAL.green, 2]].forEach(([n, cl, k], i) => { rrect(c, BX, 140 + i * 50, 28 + cnt[k] * 16, 32, 6, cl); txt(c, `${n}  ${cnt[k]}`, BX + 40 + cnt[k] * 16, 163 + i * 50, { size: 18, color: PAL.text }); });
    if (S.phase >= 3) {
      const TX = 740, TY = 370, TW = 460, TH = 200; panel(c, TX - 20, TY - 50, TW + 40, TH + 90, "Locus brightness ∝ probes attached", PAL.line);
      line(c, TX, TY + TH, TX + TW, TY + TH, PAL.grey, 2);
      c.beginPath(); S.trace.forEach((v, i) => { const X = TX + (i / 240) * TW, Y = TY + TH - (v / 26) * TH; i ? c.lineTo(X, Y) : c.moveTo(X, Y); }); c.strokeStyle = PAL.teal; c.lineWidth = 3; c.stroke();
      txt(c, S.koff > 0.3 ? "fast turnover: noisier trace" : "slow turnover: smoother trace", TX + TW, TY + TH + 30, { size: 16, color: PAL.muted, align: "right" });
    }
    const leg = ["", "FREE → BOUND at rate k_bind", "BOUND → SCANNING at k_scan (or back to FREE at k_unbind)", "SCANNING → FREE at k_off", ""][Math.min(S.phase, 3)];
    if (leg) txt(c, leg, 370, 640, { size: 20, color: PAL.amber, align: "center" });
  });
  cv.step = (dt) => {
    for (const p of ps) {
      if (p.st === 0) {
        p.vx += (r() - 0.5) * 900 * dt; p.vy += (r() - 0.5) * 900 * dt; p.vx *= 0.9; p.vy *= 0.9; p.x += p.vx * dt; p.y += p.vy * dt;
        const dx = (p.x - 370) / 290, dy = (p.y - 340) / 210; if (dx * dx + dy * dy > 1) { p.x -= p.vx * dt * 2; p.y -= p.vy * dt * 2; p.vx *= -1; p.vy *= -1; }
        if (S.phase >= 1 && r() < 0.12 * dt) { p.st = 1; p.x = LX0 + r() * (LX1 - LX0); p.y = LY - 18; }
      } else if (p.st === 1) { if (r() < 0.9 * dt) { p.st = 2; } else if (r() < 0.25 * dt) p.st = 0; p.y = LY - 18; }
      else { p.x += (r() - 0.5) * 120 * dt; p.x = clamp(p.x, LX0, LX1); p.y = LY - 8 + Math.sin(p.s * 9 + cv.t * 2) * 1; if (S.phase >= 2 && r() < S.koff * dt) { p.st = 0; p.vx = (r() - 0.5) * 80; p.vy = -60; } }
      if (S.phase < 1 && p.st) p.st = 0;
    }
    S.acc += dt; if (S.acc > 0.05) { S.acc = 0; S.trace.push(ps.filter((p) => p.st).length); if (S.trace.length > 240) S.trace.shift(); }
  };
  return { instant: false, phase(n) { if (n === 3 && S.phase !== 3) S.trace = []; S.phase = n; return this; }, koff(v) { S.koff = v; S.trace = []; return this; }, destroy() { cv.destroy(); } };
};

/* ================================================================ 8. transcription: distance -> promoter -> Pol II -> MS2 */
A.transcription = (stage) => {
  // pre-simulate 60 minutes at 1 s
  const NS = 3600, F = 14, k0 = 1 / 700, koff = 1 / 55, kinit = 0.12, thr = 400;
  let best = null;
  for (let seed = 1; seed < 400 && !best; seed++) {
    const r = rng(seed * 7919), d = new Float32Array(NS), on = new Uint8Array(NS), starts = [];
    let x = 600, st = 0;
    for (let i = 0; i < NS; i++) { x += (700 - x) * 0.08 + gauss(r) * 130; x = Math.max(60, x); d[i] = x; const near = x < thr ? 1 : 0; if (!st) { if (r() < k0 * (1 + (F - 1) * near)) st = 1; } else { if (r() < koff) st = 0; else if (r() < kinit) starts.push(i); } on[i] = st; }
    const bursts = on.reduce((a, v, i) => a + (v && !on[i - 1] ? 1 : 0), 0);
    if (bursts >= 3 && bursts <= 4 && starts.length > 12) best = { d, on, starts, r };
  }
  if (!best) { const d = new Float32Array(NS).fill(700), on = new Uint8Array(NS); best = { d, on, starts: [], r: rng(1) }; }
  const rel = best.starts.map((s) => ({ s, hold: 90 + best.r() * 40 }));
  const ms2 = new Float32Array(NS); rel.forEach(({ s, hold }) => { for (let i = s; i < NS; i++) { const a = i - s; if (a > hold + 31) break; ms2[i] += a < 31 ? 24 * a / 31 : (a < hold + 31 ? 24 : 0); } });
  const S = { show: 0, tw: new Tween(0), play: 0, shade: false };
  const X0 = 90, X1 = 1220, W = X1 - X0, px = (i) => X0 + i / NS * W;
  const cv = new Canvas(stage, (c, t, dt) => {
    S.play = (S.play + dt / 20) % 1.1; const pi = Math.min(NS - 1, Math.floor(Math.min(S.play, 1) * NS));
    txt(c, "Stage 1c: distance → promoter state → Pol II → MS2 spot", 640, 44, { size: 26, weight: 700, align: "center" });
    // distance
    panel(c, 50, 70, 1180, 170, null); txt(c, "enhancer–promoter distance", 70, 100, { size: 17, color: PAL.teal, weight: 700 });
    line(c, X0, 175, X1, 175, "#22365c", 1); const yd = (v) => 225 - clamp(v / 1400, 0, 1) * 130;
    line(c, X0, yd(thr), X1, yd(thr), PAL.amber, 1.5, [6, 5]); txt(c, "contact", X1 - 4, yd(thr) - 6, { size: 14, color: PAL.amber, align: "right" });
    if (S.shade) for (let i = 0; i < pi; i++) if (best.d[i] < thr) { c.fillStyle = "rgba(45,212,191,.22)"; c.fillRect(px(i), 85, W / NS + 0.6, 150); }
    c.beginPath(); for (let i = 0; i <= pi; i += 4) { i ? c.lineTo(px(i), yd(best.d[i])) : c.moveTo(px(i), yd(best.d[i])); } c.strokeStyle = PAL.blue; c.lineWidth = 2; c.stroke();
    // promoter
    c.globalAlpha = show(S.tw.step(dt, 3)); // promoter
    if (S.show >= 1) {
      panel(c, 50, 250, 1180, 100, null); txt(c, "promoter ON / OFF", 70, 280, { size: 17, color: PAL.purple, weight: 700 });
      for (let i = 0; i < pi; i += 2) if (best.on[i]) { c.fillStyle = PAL.purple; c.fillRect(px(i), 292, W / NS * 2 + 0.6, 40); }
      line(c, X0, 332, X1, 332, "#22365c", 1); txt(c, "k_on is raised while the enhancer is close", X1 - 8, 280, { size: 15, color: PAL.muted, align: "right" });
    }
    c.globalAlpha = 1;
    if (S.show >= 2) {
      panel(c, 50, 360, 1180, 130, null); txt(c, "Pol II on the gene (each carries MS2 loops)", 70, 390, { size: 17, color: PAL.amber, weight: 700 });
      line(c, 140, 450, 1160, 450, PAL.grey, 4); rrect(c, 140, 440, 330, 20, 4, "rgba(74,222,128,.35)"); txt(c, "MS2 cassette", 150, 432, { size: 13, color: PAL.green });
      txt(c, "TSS", 128, 456, { size: 13, color: PAL.muted, align: "right" });
      rel.forEach(({ s, hold }) => { const a = pi - s; if (a >= 0 && a < hold + 31) { const x = 140 + clamp(a / 120, 0, 1) * 1020; circ(c, x, 450, 9, PAL.amber, "#7c5a05", 2); } });
    }
    if (S.show >= 3) {
      panel(c, 50, 500, 1180, 200, null); txt(c, "MS2 signal = loops carried by all Pol II on the gene (the bright spot)", 70, 530, { size: 17, color: PAL.green, weight: 700 });
      const yy = (v) => 680 - clamp(v / 150, 0, 1) * 130; line(c, X0, 680, X1, 680, "#22365c", 1);
      c.beginPath(); for (let i = 0; i <= pi; i += 2) { i ? c.lineTo(px(i), yy(ms2[i])) : c.moveTo(px(i), yy(ms2[i])); } c.strokeStyle = PAL.green; c.lineWidth = 2.5; c.stroke();
    }
    c.fillStyle = "rgba(255,255,255,.5)"; c.fillRect(px(pi), 80, 2, 620);
    txt(c, "0", X0, 715, { size: 13, color: PAL.grey, align: "center" }); txt(c, "60 min", X1, 715, { size: 13, color: PAL.grey, align: "right" });
  });
  return { instant: false, show(n) { S.show = n; S.tw.set(1, this.instant); return this; }, shade(on) { S.shade = on; return this; }, destroy() { cv.destroy(); } };
};

/* ================================================================ 9. from the MS2 trace to Pol II starts */
A.deconv = (stage) => {
  const r = rng(55), NT = 900, S = { step: 0, tw: new Tween(0) };
  const starts = [120, 150, 175, 205, 560, 585, 612, 640];
  const ker = (a) => (a < 0 ? 0 : a < 31 ? a / 31 : a < 120 ? 1 : 0);
  const clean = new Float32Array(NT), noisy = new Float32Array(NT);
  for (let i = 0; i < NT; i++) { let v = 0; starts.forEach((s) => (v += ker(i - s))); clean[i] = v; noisy[i] = v + gauss(r) * 0.55; }
  const post = new Float32Array(NT);
  const sm = (a, w) => a.map((_, i) => { let s = 0, n = 0; for (let k = -w; k <= w; k++) { const j = i + k; if (j >= 0 && j < NT) { s += a[j]; n++; } } return s / n; });
  const d1 = Array.from(sm(noisy, 8)).map((v, i, a) => (i > 15 ? (a[i] - a[i - 16]) / 16 : 0)); for (let i = 0; i < NT; i++) post[i] = Math.max(0, d1[i + 16] || 0);
  const pmax = Math.max(...post); for (let i = 0; i < NT; i++) post[i] = post[i] / pmax * 0.85;
  const X0 = 90, X1 = 1220, px = (i) => X0 + i / NT * (X1 - X0);
  const cv = new Canvas(stage, (c, t, dt) => {
    const a = (n) => show(S.tw.step(dt, 3) - n + 1);
    txt(c, "Reading the MS2 trace backwards: when did Pol II start?", 640, 46, { size: 27, weight: 700, align: "center" });
    panel(c, 50, 80, 1180, 150, null); txt(c, "true Pol II starts (the simulation knows)", 70, 112, { size: 17, color: PAL.amber, weight: 700 });
    starts.forEach((s) => { line(c, px(s), 130, px(s), 210, PAL.amber, 3); circ(c, px(s), 130, 5, PAL.amber); });
    // kernel inset
    if (S.step >= 1) { c.globalAlpha = a(1); txt(c, "one Pol II → ramp, plateau, release", 760, 112, { size: 16, color: PAL.muted });
      c.beginPath(); for (let i = 0; i <= 150; i++) { const X = 760 + i * 2.2, Y = 205 - ker(i) * 60; i ? c.lineTo(X, Y) : c.moveTo(X, Y); } c.strokeStyle = PAL.green; c.lineWidth = 3; c.stroke(); c.globalAlpha = 1; }
    panel(c, 50, 245, 1180, 210, null); txt(c, S.step >= 2 ? "what the microscope records: sum of all Pol II + noise" : "sum of all the ramps (no noise)", 70, 277, { size: 17, color: PAL.green, weight: 700 });
    if (S.step >= 1) { c.globalAlpha = a(1); const yy = (v) => 440 - clamp((v + 1) / 9, 0, 1) * 150; line(c, X0, yy(0), X1, yy(0), "#22365c", 1);
      c.beginPath(); for (let i = 0; i < NT; i += 2) { const v = S.step >= 2 ? noisy[i] : clean[i]; i ? c.lineTo(px(i), yy(v)) : c.moveTo(px(i), yy(v)); } c.strokeStyle = S.step >= 2 ? "#86efac" : PAL.green; c.lineWidth = 2; c.stroke(); c.globalAlpha = 1; }
    if (S.step >= 3) { c.globalAlpha = a(3); panel(c, 50, 470, 1180, 200, null); txt(c, "inferred: probability of a start in each second (illustration; the dashboard draws the real posterior)", 70, 502, { size: 17, color: PAL.pink, weight: 700 });
      for (let i = 0; i < NT; i += 3) { const v = post[i]; if (v > 0.03) { c.fillStyle = "rgba(244,114,182,.75)"; c.fillRect(px(i), 650 - v * 110, 4, v * 110); } }
      starts.forEach((s) => line(c, px(s), 610, px(s), 655, "#fff", 2)); txt(c, "white ticks = truth", X1 - 10, 700, { size: 15, color: PAL.muted, align: "right" }); c.globalAlpha = 1; }
  });
  return { instant: false, step(n) { S.step = n; S.tw.set(n, this.instant); return this; }, destroy() { cv.destroy(); } };
};

/* ================================================================ 10. E-P curve: how well can it be recovered? */
A.epcurve = (stage) => {
  const S = { N: new Tween(0), sig: new Tween(0), wrong: new Tween(0), pts: new Tween(0), show: new Tween(0) };
  const r = rng(91), noise = Array.from({ length: 14 }, () => gauss(r));
  const X0 = 140, X1 = 1180, Y0 = 110, Y1 = 560, F = 12, dh = 400;
  const truthP = (d) => 0.045 * (1 + (F - 1) / (1 + Math.pow(d / dh, 4)));
  const px = (d) => X0 + d / 1500 * (X1 - X0), py = (p) => Y1 - p / 0.6 * (Y1 - Y0);
  const cv = new Canvas(stage, (c, t, dt) => {
    S.N.step(dt, 2); S.sig.step(dt, 2); S.wrong.step(dt, 2); S.pts.step(dt, 3); S.show.step(dt, 3);
    txt(c, "How well can the whole P(burst ON) curve be recovered?", 640, 46, { size: 27, weight: 700, align: "center" });
    panel(c, 60, 80, 1160, 530, null);
    line(c, X0, Y1, X1, Y1, PAL.grey, 2); line(c, X0, Y0, X0, Y1, PAL.grey, 2);
    for (let d = 0; d <= 1500; d += 250) { txt(c, String(d), px(d), Y1 + 26, { size: 14, color: PAL.grey, align: "center" }); }
    txt(c, "promoter–enhancer distance (nm)", (X0 + X1) / 2, Y1 + 52, { size: 17, color: PAL.muted, align: "center" }); c.save(); c.translate(86, 340); c.rotate(-Math.PI / 2); txt(c, "P(promoter in a burst)", 0, 0, { size: 17, color: PAL.muted, align: "center" }); c.restore();
    const flat = clamp(S.sig.v / 600, 0, 0.9), mean = 0.13, low = 1 - S.wrong.v * 1.1;
    const shown = (d) => (lerp(truthP(d), mean, flat)) * (1 - S.wrong.v * 0.9 * 0.6);
    // true curve
    c.beginPath(); for (let d = 0; d <= 1500; d += 10) { d ? c.lineTo(px(d), py(truthP(d))) : c.moveTo(px(d), py(truthP(d))); } c.strokeStyle = PAL.teal; c.lineWidth = 4; c.stroke();
    txt(c, "true curve", px(80), py(truthP(80)) - 14, { size: 18, color: PAL.teal, weight: 700 });
    // points
    const nb = 14;
    if (S.pts.v > 0.02) { c.globalAlpha = clamp(S.pts.v, 0, 1); for (let i = 0; i < nb; i++) { const d = 60 + i * 105, p = shown(d), n = Math.max(2, S.N.v) / nb, sd = Math.sqrt(p * (1 - p) / n); const v = clamp(p + noise[i] * sd * 0.8, 0.003, 0.58);
      line(c, px(d), py(clamp(v + sd, 0, 0.6)), px(d), py(clamp(v - sd, 0, 0.6)), "#fde68a", 2); circ(c, px(d), py(v), 6, PAL.amber); } c.globalAlpha = 1; }
    if (S.show.v > 0.02) { c.globalAlpha = S.show.v; c.beginPath(); for (let d = 0; d <= 1500; d += 10) { d ? c.lineTo(px(d), py(shown(d))) : c.moveTo(px(d), py(shown(d))); } c.strokeStyle = PAL.pink; c.lineWidth = 4; c.setLineDash([10, 6]); c.stroke(); c.setLineDash([]);
      txt(c, "recovered curve", px(900), py(shown(900)) - 22, { size: 18, color: PAL.pink, weight: 700 }); c.globalAlpha = 1; }
    const contrast = shown(150) / shown(1000);
    rrect(c, 690, 118, 490, 130, 14, "#0c1830", PAL.line);
    txt(c, "N cells: " + (S.N.v > 5000 ? "very many" : Math.round(S.N.v)) + "     localization error: " + Math.round(S.sig.v) + " nm\nwrong picks: " + Math.round(S.wrong.v * 100) + "%", 710, 150, { size: 17, color: PAL.muted, lh: 1.4 });
    if (S.show.v > 0.5) txt(c, "effect = P(150 nm) / P(1000 nm) = " + contrast.toFixed(1) + "   (truth: " + (truthP(150) / truthP(1000)).toFixed(1) + ", cartoon curve)", 710, 226, { size: 17, color: PAL.pink, weight: 700 });
  });
  return { instant: false, n(v) { S.N.set(v, this.instant); return this; }, sigma(v) { S.sig.set(v, this.instant); return this; }, wrong(v) { S.wrong.set(v, this.instant); return this; },
    points(on) { S.pts.set(on ? 1 : 0, this.instant); return this; }, recovered(on) { S.show.set(on ? 1 : 0, this.instant); return this; }, destroy() { cv.destroy(); } };
};
