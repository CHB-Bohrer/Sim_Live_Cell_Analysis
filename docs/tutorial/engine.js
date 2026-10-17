/* Narrated tutorial engine: scenes made of beats (a sentence of narration + a visual cue), spoken with the browser's
   speech synthesis (Microsoft Edge has the most natural voices). No libraries, works from file://. */
"use strict";
const T = {
  scenes: [], chapters: [], si: -1, bi: 0, playing: false, frozen: false, started: false,
  rate: 1, voice: null, voices: [], mute: false, caps: true, ctl: null, stage: null, doneSet: new Set(),
  token: 0, timers: [],
};
const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => Array.from(el.querySelectorAll(s));
function h(tag, attrs, ...kids) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (k === "class") e.className = v; else if (k === "html") e.innerHTML = v; else if (k === "style") e.style.cssText = v;
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), v); else e.setAttribute(k, v);
  }
  for (const c of kids.flat()) if (c != null) e.append(c.nodeType ? c : document.createTextNode(c));
  return e;
}
const sleep = (ms) => new Promise((r) => T.timers.push(setTimeout(r, ms)));
function store(k, v) { try { if (v === undefined) return JSON.parse(localStorage.getItem("simlive_tut_" + k)); localStorage.setItem("simlive_tut_" + k, JSON.stringify(v)); } catch (e) { return null; } }

/* ---------------------------------------------------------------- registration */
T.chapter = (title) => { T.chapters.push({ title, scenes: [] }); };
T.add = (scene) => {
  scene.ch = T.chapters.length - 1; scene.idx = T.scenes.length;
  T.chapters[scene.ch].scenes.push(scene.idx); T.scenes.push(scene);
};

/* ---------------------------------------------------------------- voices */
const PREFER = ["aria", "jenny", "sonia", "libby", "guy", "ryan", "emma", "andrew", "brian", "ava", "michelle", "steffan"];
function rankVoice(v) {
  const n = v.name.toLowerCase();
  if (!/^en[-_]/i.test(v.lang)) return -1;
  let s = 1;
  if (n.includes("natural")) s += 100;
  if (n.includes("online")) s += 40;
  if (n.includes("neural")) s += 60;
  if (n.includes("google")) s += 25;
  if (n.includes("samantha") || n.includes("ava") || n.includes("allison")) s += 20;
  PREFER.forEach((p, i) => { if (n.includes(p)) s += 12 - i * 0.5; });
  if (/en-us/i.test(v.lang)) s += 5;
  if (n.includes("david") || n.includes("zira") || n.includes("mark")) s -= 20;     // old robotic desktop voices
  return s;
}
function loadVoices() {
  const vs = (window.speechSynthesis ? speechSynthesis.getVoices() : []).filter((v) => /^en[-_]/i.test(v.lang));
  vs.sort((a, b) => rankVoice(b) - rankVoice(a));
  T.voices = vs;
  const saved = store("voice");
  T.voice = vs.find((v) => v.name === saved) || vs[0] || null;
  for (const sel of $$("select.voice")) {
    sel.innerHTML = "";
    if (!vs.length) sel.append(h("option", { value: "" }, "No voices found (captions only)"));
    vs.forEach((v) => sel.append(h("option", { value: v.name }, v.name.replace(/^Microsoft /, "").replace(/ - English.*/, "") + (v.localService ? "" : "  (online)"))));
    if (T.voice) sel.value = T.voice.name;
  }
}
function wordsOf(text) { return text.split(/(\s+)/); }

/* ---------------------------------------------------------------- speech */
let utter = null;
function stopSpeech() { try { if (window.speechSynthesis) speechSynthesis.cancel(); } catch (e) { } utter = null; }
function speak(text, spans, token) {
  return new Promise((resolve) => {
    const nWords = text.split(/\s+/).length;
    const est = (nWords / 2.55 / T.rate) * 1000 + 700;            // read-aloud time, also the fallback when there is no voice
    let done = false;
    const finish = () => { if (!done) { done = true; clearTimeout(wd); resolve(); } };
    const wd = setTimeout(finish, T.mute || !T.voice ? est : est * 1.8 + 4000);
    T.timers.push(wd);
    if (T.mute || !T.voice || !window.speechSynthesis) {            // no audio: reveal the words at reading speed
      const t0 = performance.now();
      const tick = () => {
        if (done || token !== T.token) return;
        const k = Math.floor(((performance.now() - t0) / est) * spans.length);
        spans.forEach((s, i) => s.classList.toggle("said", i < k));
        requestAnimationFrame(tick);
      };
      tick();
      return;
    }
    stopSpeech();
    const u = new SpeechSynthesisUtterance(text);
    u.voice = T.voice; u.lang = T.voice.lang; u.rate = T.rate; u.pitch = 1;
    u.onboundary = (e) => {
      if (token !== T.token) return;
      let pos = 0, k = 0; const parts = wordsOf(text);
      for (const p of parts) { if (pos >= e.charIndex) break; pos += p.length; if (p.trim()) k++; }
      spans.forEach((s, i) => s.classList.toggle("said", i < k));
    };
    u.onend = () => { spans.forEach((s) => s.classList.add("said")); finish(); };
    u.onerror = () => finish();
    utter = u;
    try { T.timers.push(setTimeout(() => { if (!done && token === T.token) speechSynthesis.speak(u); }, 80)); } catch (e) { finish(); }   // speak() right after cancel() is dropped by some browsers
  });
}

/* ---------------------------------------------------------------- playback */
function setCaption(text) {
  const cap = $("#cap-text"); cap.innerHTML = "";
  const spans = [];
  for (const p of wordsOf(text)) {
    if (!p.trim()) { cap.append(p); continue; }
    const s = h("span", { class: "w" }, p); cap.append(s); spans.push(s);
  }
  return spans;
}
function updateUI() {
  const sc = T.scenes[T.si]; if (!sc) return;
  $("#crumb").textContent = T.chapters[sc.ch].title + "  ›  " + sc.title;
  $$("#rail button").forEach((b) => { b.classList.toggle("on", +b.dataset.i === T.si); b.classList.toggle("done", T.doneSet.has(+b.dataset.i)); });
  const total = T.scenes.reduce((a, s) => a + s.beats.length, 0);
  const before = T.scenes.slice(0, T.si).reduce((a, s) => a + s.beats.length, 0) + T.bi + 1;
  $(".prog i").style.width = (100 * before / total) + "%";
  $("#info").textContent = `scene ${T.si + 1}/${T.scenes.length} · beat ${T.bi + 1}/${sc.beats.length}`;
  $("#play").textContent = T.playing ? "⏸ Pause" : "▶ Play";
  const act = $(`#rail button[data-i="${T.si}"]`); if (act && act.scrollIntoView) act.scrollIntoView({ block: "nearest" });
  try { history.replaceState(null, "", `#s=${sc.id}&b=${T.bi}`); } catch (e) { }
}
function showTry(on) {
  const sc = T.scenes[T.si], el = $("#try");
  if (on && sc && sc.tryit) {
    el.innerHTML = "<h5>Your turn</h5><ol>" + sc.tryit.map((t) => "<li>" + t + "</li>").join("") + "</ol>";
    el.classList.add("show");
  } else el.classList.remove("show");
}

function mountScene(i, upToBeat) {
  T.token++; stopSpeech(); T.timers.forEach(clearTimeout); T.timers = [];
  if (T.ctl && T.ctl.destroy) T.ctl.destroy();
  const stage = $("#stage"); stage.innerHTML = ""; showTry(false);
  T.si = i; T.bi = 0;
  const sc = T.scenes[i];
  T.ctl = sc.make(stage);
  if (T.ctl) T.ctl.instant = true;
  for (let b = 0; b < upToBeat; b++) { try { sc.beats[b].c && sc.beats[b].c(T.ctl); } catch (e) { console.error(e); } }
  if (T.ctl) T.ctl.instant = false;
}

async function runBeat(i, b, token) {
  const sc = T.scenes[i], beat = sc.beats[b];
  T.si = i; T.bi = b; updateUI();
  const spans = setCaption(beat.s);
  try { beat.c && beat.c(T.ctl); } catch (e) { console.error(e); }
  if (b === sc.beats.length - 1) showTry(true);
  if (!T.playing) { spans.forEach((s) => s.classList.add("said")); return; }
  await speak(beat.s, spans, token);
  if (token !== T.token || !T.playing) return;
  await sleep(beat.pause != null ? beat.pause : 450);
  if (token !== T.token || !T.playing) return;
  if (b < sc.beats.length - 1) return runBeat(i, b + 1, token);
  T.doneSet.add(i); updateUI();
  await sleep(sc.tryit ? 2600 : 900);
  if (token !== T.token || !T.playing) return;
  if (i < T.scenes.length - 1) { T.goto(i + 1, 0, true); } else { T.pause(); setCaption("That is the end of the tutorial. Everything shown here is in the dashboard, so open it and try the Your Turn steps."); $("#cap-text").querySelectorAll(".w").forEach((s) => s.classList.add("said")); }
}

T.goto = (i, b = 0, keepPlaying = false) => {
  i = Math.max(0, Math.min(T.scenes.length - 1, i));
  const was = keepPlaying || T.playing;
  T.playing = false; T.frozen = false;
  mountScene(i, b);
  T.playing = was;
  const token = ++T.token;
  runBeat(i, b, token);
};
T.play = () => { if (T.playing) return; T.playing = true; T.frozen = false; updateUI(); T.goto(T.si < 0 ? 0 : T.si, T.si < 0 ? 0 : T.bi, true); };
T.pause = () => { T.playing = false; T.frozen = true; T.token++; stopSpeech(); T.timers.forEach(clearTimeout); T.timers = []; updateUI(); };
T.toggle = () => (T.playing ? T.pause() : T.play());
T.nextBeat = () => { const sc = T.scenes[T.si]; if (T.bi < sc.beats.length - 1) T.goBeat(T.bi + 1); else if (T.si < T.scenes.length - 1) T.goto(T.si + 1, 0); };
T.prevBeat = () => { if (T.bi > 0) T.goBeat(T.bi - 1); else if (T.si > 0) T.goto(T.si - 1, 0); };
T.goBeat = (b) => { T.goto(T.si, b); };

/* ---------------------------------------------------------------- global animation clock */
const frameFns = new Set(); let lastT = 0;
T.onFrame = (fn) => { frameFns.add(fn); return () => frameFns.delete(fn); };
function loop(now) {
  const dt = Math.min(0.05, (now - lastT) / 1000 || 0); lastT = now;
  const d = T.frozen ? 0 : dt;
  for (const fn of Array.from(frameFns)) { try { fn(d, now / 1000); } catch (e) { console.error(e); } }
  requestAnimationFrame(loop);
}
requestAnimationFrame(loop);

/* ---------------------------------------------------------------- screenshot camera */
class Shot {
  /* sets: {key: shotName}; the browser frame shows the picture at its natural size in a moving camera. */
  constructor(stage, opts) {
    this.stage = stage; this.opts = opts; this.sets = opts.sets || { main: opts.shot }; this.key = opts.start || Object.keys(this.sets)[0];
    this.root = h("div", { class: "browser" });
    this.url = h("span", { class: "url" }, opts.url || "localhost:8501"); this.pgl = h("span", { class: "pg" }, opts.label || "");
    this.root.append(h("div", { class: "chrome" }, h("i"), h("i"), h("i"), this.url, this.pgl));
    this.view = h("div", { class: "view" }); this.world = h("div", { class: "world" });
    this.hls = []; this.tip = h("div", { class: "tip" });
    this.cur = h("div", { class: "cur", html: '<svg viewBox="0 0 24 24" width="26" height="26"><path d="M3 2l7 18 2.6-7.4L20 10z" fill="#fff" stroke="#111" stroke-width="1.4" stroke-linejoin="round"/></svg>' });
    this.view.append(this.world); this.root.append(this.view); stage.append(this.root);
    this.use(this.key, true);
  }
  get data() { return window.SHOTS[this.sets[this.key]]; }
  use(key, first) {
    this.key = key; const d = this.data; if (!d) { console.warn("missing shot", key); return; }
    this.world.innerHTML = "";
    this.img = h("img", { src: d.src, style: `width:${d.w}px;height:${d.h}px` });
    this.world.append(this.img);
    this.hls = []; this.world.append(this.tip, this.cur);
    this.tip.classList.remove("show"); this.cur.classList.remove("show");
    this.dpr = d.dpr || 1;
    this.camera(null, { top: true }, true);
  }
  size() { const b = this.view.getBoundingClientRect(); return [b.width, b.height]; }
  rectOf(spec) {
    const d = this.data, R = d.rects;
    if (typeof spec === "string") return R[spec] || null;
    if (Array.isArray(spec) && typeof spec[0] === "number") return spec;
    if (Array.isArray(spec)) {
      const rs = spec.map((n) => R[n]).filter(Boolean); if (!rs.length) return null;
      const x0 = Math.min(...rs.map((r) => r[0])), y0 = Math.min(...rs.map((r) => r[1]));
      const x1 = Math.max(...rs.map((r) => r[0] + r[2])), y1 = Math.max(...rs.map((r) => r[1] + r[3]));
      return [x0, y0, x1 - x0, y1 - y0];
    }
    return null;
  }
  camera(rect, o = {}, instant) {
    const d = this.data, [vw, vh] = this.size();
    const fitW = vw / d.w, fitH = vh / d.h;
    const s0 = this.opts.minScale || fitW;
    let s, tx, ty;
    if (!rect || o.top) { s = s0; tx = (vw - d.w * s) / 2; ty = 0; }
    else {
      const padX = o.padX || 1.3, padY = o.padY || 1.4;
      s = Math.min(vw / (rect[2] * padX), vh / (rect[3] * padY), (o.maxZoom || 1.3) * (this.opts.minScale ? 1 : 1));
      s = Math.max(s, s0);
      tx = vw / 2 - (rect[0] + rect[2] / 2) * s;
      ty = (rect[3] * s > vh * 0.82 || o.top || (rect[3] < 110 && rect[2] > 1000)) ? vh * 0.15 - rect[1] * s : vh / 2 - (rect[1] + rect[3] / 2) * s;
    }
    if (d.w * s > vw) tx = Math.min(0, Math.max(vw - d.w * s, tx)); else tx = this.opts.align === "left" ? vw * 0.05 : (vw - d.w * s) / 2;
    if (d.h * s > vh) ty = Math.min(0, Math.max(vh - d.h * s, ty)); else ty = 0;
    this.s = s;
    this.root.classList.toggle("instant", !!(instant || this.instant));
    this.world.style.setProperty("--inv", String(1 / s));
    this.world.style.transform = `translate(${tx}px,${ty}px) scale(${s})`;
    this.tx = tx; this.ty = ty;
  }
  clear() {
    this.hls.forEach((e) => e.remove()); this.hls = [];
    this.tip.classList.remove("show");
  }
  /* focus: move the camera to a named rectangle (or several), outline it and show a label. */
  focus(spec, label, o = {}) {
    const ok = this.instant;
    const rect = this.rectOf(spec); if (!rect) { console.warn("no rect", spec); return this; }
    this.clear();
    const names = Array.isArray(spec) && typeof spec[0] === "string" && o.multi ? spec : null;
    this.camera(o.view ? this.rectOf(o.view) : rect, o);
    const place = () => {
      const inv = 1 / this.s, pad = 5 * inv;
      const rs = names ? names.map((n) => this.rectOf(n)).filter(Boolean) : [rect];
      for (const r of rs) {
        const e = h("div", { class: "hl" + (names ? " multi" : "") });
        e.style.cssText = `left:${r[0] - pad}px;top:${r[1] - pad}px;width:${r[2] + 2 * pad}px;height:${r[3] + 2 * pad}px`;
        this.world.insertBefore(e, this.tip); this.hls.push(e);
        requestAnimationFrame(() => e.classList.add("show"));
      }
      if (label) {
        this.tip.textContent = label;
        const sc = this.s, [vw, vh] = this.size();
        const rx = rect[0] * sc + this.tx, ry = rect[1] * sc + this.ty;
        const tipH = 38, tipW = Math.min(360, label.length * 9 + 28);
        const vx = Math.max(10, Math.min(rx, vw - tipW - 10));
        let vy = ry - tipH - 10;
        if (vy < 8) vy = Math.max(8, Math.min(ry + 10, vh - tipH - 8));
        this.tip.style.transform = `scale(${inv})`;
        this.tip.style.left = `${(vx - this.tx) / sc}px`;
        this.tip.style.top = `${(vy - this.ty) / sc}px`;
        requestAnimationFrame(() => this.tip.classList.add("show"));
      }
    };
    if (ok) place(); else this._t = setTimeout(place, 750);
    return this;
  }
  between(a, b, pad = 10) {
    const A_ = this.data.rects[a], B_ = b ? this.data.rects[b] : null; if (!A_) return null;
    const y0 = A_[1] - pad, y1 = B_ ? B_[1] - pad : this.data.h;
    return [0, y0, this.data.w, Math.max(40, y1 - y0)];
  }
  note(html) {
    if (!this.noteEl) { this.noteEl = h("div", { class: "note" }); this.view.append(this.noteEl); }
    this.noteEl.innerHTML = html || ""; this.noteEl.classList.toggle("show", !!html); return this;
  }
  top() { this.clear(); this.camera(null, { top: true }); return this; }
  /* the mouse pointer moves to a rectangle and clicks. */
  click(spec, o = {}) {
    const r = this.rectOf(spec); if (!r) return this;
    const x = r[0] + (o.fx != null ? o.fx : 0.5) * r[2], y = r[1] + (o.fy != null ? o.fy : 0.5) * r[3];
    const place = () => {
      this.cur.classList.add("show");
      this.cur.style.transform = `scale(${1 / this.s})`;
      this.cur.style.left = x + "px"; this.cur.style.top = y + "px";
      if (!this.instant && !o.nodown) setTimeout(() => {
        const rp = h("div", { class: "ripple" }); rp.style.left = x + "px"; rp.style.top = y + "px"; this.world.append(rp); setTimeout(() => rp.remove(), 900);
      }, 900);
    };
    if (this.instant) place(); else setTimeout(place, o.delay != null ? o.delay : 150);
    return this;
  }
  hideCursor() { this.cur.classList.remove("show"); return this; }
  destroy() { clearTimeout(this._t); }
}

/* ---------------------------------------------------------------- canvas scenes (logical 1280 x 720) */
const PAL = { bg: "#0f1a31", panel: "#16233f", line: "#2a3c64", text: "#e8edf6", muted: "#93a4c3", teal: "#2dd4bf", blue: "#60a5fa", amber: "#fbbf24",
  red: "#f87171", green: "#4ade80", purple: "#a78bfa", pink: "#f472b6", grey: "#64748b" };
function rng(seed) { let s = seed >>> 0 || 1; return () => { s ^= s << 13; s >>>= 0; s ^= s >>> 17; s ^= s << 5; s >>>= 0; return s / 4294967296; }; }
function gauss(r) { let u = 0, v = 0; while (!u) u = r(); v = r(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
const lerp = (a, b, t) => a + (b - a) * t;
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
const ease = (t) => t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2;
class Canvas {
  constructor(stage, drawFn, opts = {}) {
    this.stage = stage; this.W = opts.W || 1280; this.H = opts.H || 720;
    this.cv = h("canvas"); stage.append(this.cv); this.ctx = this.cv.getContext("2d");
    this.t = 0; this.drawFn = drawFn; this.instant = false;
    this.ro = new ResizeObserver(() => this.resize()); this.ro.observe(stage); this.resize();
    this.off = T.onFrame((dt) => { this.t += dt; this.step && this.step(dt); this.paint(dt); });
  }
  resize() {
    const r = this.stage.getBoundingClientRect(), dpr = Math.min(2, window.devicePixelRatio || 1);
    this.cv.width = Math.max(2, Math.round(r.width * dpr)); this.cv.height = Math.max(2, Math.round(r.height * dpr));
    this.k = this.cv.width / this.W;
  }
  paint(dt) {
    const c = this.ctx; c.setTransform(this.k, 0, 0, this.k, 0, 0);
    c.fillStyle = PAL.bg; c.fillRect(0, 0, this.W, this.H);
    this.drawFn(c, this.t, dt);
  }
  destroy() { this.off(); this.ro.disconnect(); }
}
function txt(c, s, x, y, o = {}) {
  const pa = c.globalAlpha;
  c.font = `${o.weight || 500} ${o.size || 22}px Inter, "Segoe UI", system-ui, sans-serif`;
  c.fillStyle = o.color || PAL.text; c.textAlign = o.align || "left"; c.textBaseline = o.base || "alphabetic";
  if (o.alpha != null) c.globalAlpha = pa * o.alpha;
  let lines = String(s).split("\n");
  if (o.w) {                                   // word-wrap to a maximum width
    const out = [];
    for (const ln of lines) { let cur = ""; for (const wd of ln.split(" ")) { const t = cur ? cur + " " + wd : wd; if (c.measureText(t).width > o.w && cur) { out.push(cur); cur = wd; } else cur = t; } out.push(cur); }
    lines = out;
  }
  const lh = (o.size || 22) * (o.lh || 1.3);
  lines.forEach((l, i) => c.fillText(l, x, y + i * lh));
  c.globalAlpha = pa;
}
function rrect(c, x, y, w, h_, r, fill, stroke, lw = 2) {
  c.beginPath(); c.roundRect ? c.roundRect(x, y, w, h_, r) : c.rect(x, y, w, h_);
  if (fill) { c.fillStyle = fill; c.fill(); } if (stroke) { c.strokeStyle = stroke; c.lineWidth = lw; c.stroke(); }
}
function circ(c, x, y, r, fill, stroke, lw = 2) {
  c.beginPath(); c.arc(x, y, r, 0, Math.PI * 2);
  if (fill) { c.fillStyle = fill; c.fill(); } if (stroke) { c.strokeStyle = stroke; c.lineWidth = lw; c.stroke(); }
}
function line(c, x1, y1, x2, y2, col, lw = 2, dash) {
  c.beginPath(); c.moveTo(x1, y1); c.lineTo(x2, y2); c.strokeStyle = col; c.lineWidth = lw; c.setLineDash(dash || []); c.stroke(); c.setLineDash([]);
}
function arrow(c, x1, y1, x2, y2, col, lw = 3) {
  line(c, x1, y1, x2, y2, col, lw);
  const a = Math.atan2(y2 - y1, x2 - x1);
  c.beginPath(); c.moveTo(x2, y2); c.lineTo(x2 - 14 * Math.cos(a - 0.45), y2 - 14 * Math.sin(a - 0.45)); c.lineTo(x2 - 14 * Math.cos(a + 0.45), y2 - 14 * Math.sin(a + 0.45)); c.closePath(); c.fillStyle = col; c.fill();
}
/* smooth value that moves toward a target (for cue-driven animation) */
class Tween {
  constructor(v = 0) { this.v = v; this.to = v; }
  set(x, instant) { this.to = x; if (instant) this.v = x; return this; }
  step(dt, speed = 3) { this.v += (this.to - this.v) * Math.min(1, dt * speed); return this.v; }
}

/* ---------------------------------------------------------------- terminal */
class Term {
  constructor(stage, parent) {
    this.el = h("div", { class: "term" }); (parent || stage).append(this.el); this.alive = true; this.instant = false;
  }
  line(html, cls = "") { const d = h("div", { class: "ln " + cls, html }); this.el.append(d); this.el.scrollTop = 1e6; return d; }
  clear() { this.el.innerHTML = ""; }
  async type(cmd, speed = 38) {
    const d = this.line('<span class="p">PS C:\\Sim_Live_Cell_Analysis&gt;</span> <span class="t"></span><span class="caret"></span>');
    const t = $(".t", d);
    if (this.instant) { t.textContent = cmd; $(".caret", d).remove(); return; }
    for (const ch of cmd) { if (!this.alive) return; t.textContent += ch; await new Promise((r) => setTimeout(r, speed + Math.random() * 25)); }
    $(".caret", d).remove();
  }
  async print(lines, gap = 220) {
    for (const l of lines) { if (!this.alive) return; this.line(l.html || l, l.cls || ""); if (!this.instant) await new Promise((r) => setTimeout(r, l.gap != null ? l.gap : gap)); }
  }
  destroy() { this.alive = false; }
}

/* ---------------------------------------------------------------- reveal cards / tables */
class Reveal {
  constructor(stage, html) { this.el = h("div", { class: "cards", html }); stage.append(this.el); this.n = 0; this.items = () => $$(".card2, tr", this.el); }
  to(n) { this.items().forEach((e, i) => e.classList.toggle("show", i < n)); this.n = n; return this; }
  more() { return this.to(this.n + 1); }
  only(idx) { this.items().forEach((e, i) => e.classList.toggle("show", idx.includes(i))); return this; }
}

/* ---------------------------------------------------------------- boot */
T.boot = () => {
  const rail = $("#rail");
  T.chapters.forEach((ch, ci) => {
    rail.append(h("h4", {}, ch.title));
    ch.scenes.forEach((si) => rail.append(h("button", { "data-i": si, onclick: () => { T.started = true; T.goto(si, 0); $("#rail").classList.remove("open"); } }, T.scenes[si].title)));
  });
  $("#play").onclick = () => T.toggle();
  $("#prev").onclick = () => T.prevBeat(); $("#next").onclick = () => T.nextBeat();
  $("#prevs").onclick = () => T.goto(T.si - 1, 0); $("#nexts").onclick = () => T.goto(T.si + 1, 0);
  $("#replay").onclick = () => T.goto(T.si, 0);
  $("#rail-toggle").onclick = () => $("#rail").classList.toggle("open");
  const sp = store("rate"); if (sp) T.rate = sp;
  for (const sel of $$("select.rate")) { sel.value = String(T.rate); sel.onchange = (e) => { T.rate = +e.target.value; store("rate", T.rate); $$("select.rate").forEach((s) => (s.value = String(T.rate))); }; }
  $("#mute").onclick = () => { T.mute = !T.mute; $("#mute").classList.toggle("on", T.mute); $("#mute").textContent = T.mute ? "🔇 Silent" : "🔊 Voice"; if (T.mute) stopSpeech(); };
  $("#capbtn").onclick = () => { T.caps = !T.caps; $("#caption").classList.toggle("off", !T.caps); $("#capbtn").classList.toggle("on", T.caps); };
  $("#capbtn").classList.add("on");
  for (const sel of $$("select.voice")) sel.onchange = (e) => { T.voice = T.voices.find((v) => v.name === e.target.value) || null; store("voice", e.target.value); $$("select.voice").forEach((s) => (s.value = e.target.value)); };
  $("#testvoice").onclick = () => { stopSpeech(); if (T.voice) { const u = new SpeechSynthesisUtterance("Hello. This is how the narration will sound."); u.voice = T.voice; u.rate = T.rate; speechSynthesis.speak(u); } };
  if (window.speechSynthesis) { loadVoices(); speechSynthesis.onvoiceschanged = loadVoices; setTimeout(loadVoices, 600); setTimeout(loadVoices, 2000); }
  document.addEventListener("keydown", (e) => {
    if (/select|input|textarea/i.test(e.target.tagName)) return;
    if (e.code === "Space") { e.preventDefault(); if (!T.started) return; T.toggle(); }
    else if (e.code === "ArrowRight") e.shiftKey ? T.goto(T.si + 1, 0) : T.nextBeat();
    else if (e.code === "ArrowLeft") e.shiftKey ? T.goto(T.si - 1, 0) : T.prevBeat();
    else if (e.key === "m") $("#mute").click();
    else if (e.key === "c") $("#capbtn").click();
  });
  const m = /s=([\w-]+)(?:&b=(\d+))?/.exec(location.hash);
  if (m) {
    const i = T.scenes.findIndex((s) => s.id === m[1]);
    if (i >= 0) { $("#start").remove(); T.started = true; T.frozen = false; T.goto(i, +(m[2] || 0)); return; }
  }
  $("#startbtn").onclick = () => { $("#start").remove(); T.started = true; T.goto(0, 0, true); };
  $("#browsebtn").onclick = () => { $("#start").remove(); T.started = true; T.goto(0, 0); };
  mountScene(0, 0); T.si = 0; updateUI();
};
