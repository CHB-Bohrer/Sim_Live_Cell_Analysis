"""One look for the whole dashboard: page CSS, matplotlib defaults, and small layout helpers.

Colours are defined once here (PALETTE). Matplotlib figures use the same palette so plots match the page.
"""
import html

import matplotlib as mpl

# dark theme (the user asked for a darker interface, 2026-10-12); .streamlit/config.toml holds the same colours for Streamlit's widgets
PALETTE = {"ink": "#e5e7eb", "muted": "#94a3b8", "line": "#26324a", "card": "#111a2c", "bg": "#0b1220",
           "accent": "#14b8a6", "accent_soft": "#123b3a", "good": "#4ade80", "warn": "#f59e0b", "bad": "#f87171",
           "blue": "#60a5fa", "violet": "#a78bfa", "fig": "#0f1829", "grid": "#1d2940"}
CYCLE = ["#2dd4bf", "#60a5fa", "#f59e0b", "#a78bfa", "#f87171", "#94a3b8", "#22d3ee", "#facc15"]

CSS = f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
:root {{ --ink:{PALETTE['ink']}; --muted:{PALETTE['muted']}; --line:{PALETTE['line']}; --accent:{PALETTE['accent']}; --card:{PALETTE['card']}; }}
.block-container {{ padding-top: 2.2rem; padding-bottom: 4rem; max-width: 1280px; }}
h1 {{ font-weight: 700; letter-spacing: -0.02em; font-size: 2rem !important; margin-bottom: .1rem; }}
h2, h3 {{ font-weight: 650; letter-spacing: -0.01em; }}
h3 {{ margin-top: 2rem; }}
[data-testid="stCaptionContainer"] {{ color: var(--muted); }}

/* tabs as quiet pills */
[data-baseweb="tab-list"] {{ gap: .25rem; border-bottom: 1px solid var(--line); }}
[data-baseweb="tab"] {{ height: 2.6rem; padding: 0 .9rem; border-radius: .6rem .6rem 0 0; font-weight: 550; }}
[data-baseweb="tab"][aria-selected="true"] {{ color: var(--accent); }}
[data-baseweb="tab-highlight"] {{ background-color: var(--accent) !important; height: 3px; }}

/* metrics as cards */
[data-testid="stMetric"] {{ background:var(--card); border:1px solid var(--line); border-radius:.8rem; padding:.8rem 1rem;
  box-shadow: 0 1px 2px rgba(0,0,0,.35); }}
[data-testid="stMetricLabel"] {{ color: var(--muted); }}
[data-testid="stMetricValue"] {{ font-weight: 700; letter-spacing:-0.02em; }}

/* expanders, forms, dataframes, images get the same card edge */
[data-testid="stExpander"] {{ border:1px solid var(--line); border-radius:.8rem; background:var(--card); }}
[data-testid="stForm"] {{ border:1px solid var(--line); border-radius:.8rem; background:var(--card); }}
[data-testid="stDataFrame"], [data-testid="stImage"] img, [data-testid="stPyplot"] img {{ border-radius:.6rem; }}
[data-testid="stSidebar"] {{ background:#0e1626; border-right:1px solid var(--line); }}
.stButton > button, .stFormSubmitButton > button {{ border-radius:.6rem; font-weight:600; }}

/* helper components */
.sl-hero {{ margin: 0 0 1.1rem 0; }}
.sl-hero p {{ color: var(--muted); margin:.1rem 0 0 0; font-size:1.02rem; }}
.sl-card {{ background:var(--card); border:1px solid var(--line); border-radius:.9rem; padding:1rem 1.1rem; min-height:12.5rem;
  box-shadow: 0 1px 2px rgba(0,0,0,.35); }}
.sl-card h4 {{ margin:.1rem 0 .35rem 0; font-size:1rem; }}
.sl-card p {{ margin:0; color:var(--muted); font-size:.92rem; line-height:1.45; }}
.sl-step {{ display:inline-flex; width:1.7rem; height:1.7rem; border-radius:50%; background:{PALETTE['accent_soft']};
  color:var(--accent); font-weight:700; align-items:center; justify-content:center; margin-right:.5rem; }}
.sl-badge {{ display:inline-block; padding:.15rem .65rem; border-radius:999px; font-size:.82rem; font-weight:650; }}
.sl-badge.good {{ background:#14532d; color:#bbf7d0; }} .sl-badge.warn {{ background:#78350f; color:#fde68a; }}
.sl-badge.bad {{ background:#7f1d1d; color:#fecaca; }} .sl-badge.idle {{ background:#1f2937; color:#cbd5e1; }}
.sl-badge.run {{ background:{PALETTE['accent_soft']}; color:#99f6e4; }}
.sl-stats {{ display:grid; grid-template-columns: repeat(auto-fill, minmax(200px, 1fr)); gap:.6rem; margin:.5rem 0 .2rem 0; }}
.sl-stat {{ background:var(--card); border:1px solid var(--line); border-radius:.7rem; padding:.55rem .8rem; }}
.sl-stat b {{ display:block; font-size:1.35rem; letter-spacing:-0.02em; line-height:1.3; }}
.sl-stat span {{ color:var(--muted); font-size:.8rem; }}
.sl-mono {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size:.85rem; }}

/* ---- polish layer ---- */
.stApp, .stMarkdown, p, label, input, textarea, button, [data-baseweb="tab"] {{ font-family: 'Inter', system-ui, -apple-system, 'Segoe UI', sans-serif; }}
[data-testid="stIconMaterial"], [class*="material-symbols"], .material-icons {{ font-family: 'Material Symbols Rounded', 'Material Icons' !important; }}
.stApp {{ background: radial-gradient(1200px 400px at 20% -10%, #12324a 0%, rgba(11,18,32,0) 60%), #0b1220; }}
.sl-hero {{ border-left: 5px solid var(--accent); padding: .35rem 0 .35rem 1rem; margin-bottom: 1.4rem; }}
.sl-hero h1 {{ background: linear-gradient(90deg, #2dd4bf, #60a5fa); -webkit-background-clip: text; background-clip: text;
  -webkit-text-fill-color: transparent; }}
[data-baseweb="tab"] {{ transition: background .15s, color .15s; }}
[data-baseweb="tab"]:hover {{ background: {PALETTE['accent_soft']}; color: var(--accent); }}
[data-testid="stMetric"], [data-testid="stExpander"], [data-testid="stForm"], .sl-card, .sl-stat {{
  transition: box-shadow .15s, transform .15s; }}
[data-testid="stMetric"]:hover, .sl-card:hover {{ box-shadow: 0 6px 18px rgba(0,0,0,.45); transform: translateY(-1px); }}
[data-testid="stExpander"] summary {{ font-weight: 600; }}
.stButton > button, .stFormSubmitButton > button {{ transition: box-shadow .15s, transform .15s; }}
.stButton > button:hover, .stFormSubmitButton > button:hover {{ box-shadow: 0 4px 12px rgba(20,184,166,.30); transform: translateY(-1px); }}
.stButton > button[kind="primary"], .stFormSubmitButton > button[kind="primary"] {{
  background: linear-gradient(135deg, #0f766e, #0e7490); border: 0; color: #f0fdfa; }}
[data-testid="stAlert"] {{ border-radius: .8rem; }}
/* compact metric cards (the headline scores were too large) */
[data-testid="stMetric"] {{ padding: .35rem .7rem; border-radius: .6rem; }}
[data-testid="stMetricLabel"] p {{ font-size: .78rem; }}
[data-testid="stMetricValue"] {{ font-size: 1.35rem !important; line-height: 1.25; }}
[data-testid="stMetric"]:hover {{ transform: none; }}
.stMarkdown strong {{ font-weight: 600; }}
[data-baseweb="input"], [data-baseweb="select"] > div {{ border-radius: .55rem; }}
[data-testid="stSidebar"] h2 {{ font-size: 1.25rem; letter-spacing: -0.01em; }}
[data-testid="stSidebar"] [data-testid="stExpander"] {{ background: var(--card); }}
.sl-intro {{ display:grid; grid-template-columns: 1fr 1fr; gap: .8rem; margin: .2rem 0 1.2rem 0; }}
.sl-intro > div {{ border-radius: .85rem; padding: .85rem 1.05rem; border: 1px solid var(--line); background:var(--card); }}
.sl-intro .what {{ border-left: 4px solid var(--accent); }}
.sl-intro .how {{ border-left: 4px solid {PALETTE['blue']}; }}
.sl-intro b.t {{ display:block; font-size:.74rem; letter-spacing:.08em; text-transform:uppercase; margin-bottom:.25rem; }}
.sl-intro .what b.t {{ color: var(--accent); }} .sl-intro .how b.t {{ color: {PALETTE['blue']}; }}
.sl-intro p {{ margin:0; color: var(--ink); font-size:.93rem; line-height:1.5; }}
@media (max-width: 900px) {{ .sl-intro {{ grid-template-columns: 1fr; }} }}
.sl-guide {{ display:grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap:.6rem; }}
.sl-guide div {{ background:var(--card); border:1px solid var(--line); border-radius:.75rem; padding:.6rem .8rem; font-size:.88rem;
  line-height:1.4; color: var(--muted); }}
.sl-guide b {{ color: var(--ink); display:block; font-size:.95rem; }}
</style>
"""


PYPLOT_DPI = 110      # st.pyplot's own default is 200 dpi: a 17-inch figure became a 3400 px PNG (rasterised + compressed on every rerun)


def fast_pyplot(st) -> None:
    """Replace st.pyplot by a version that saves at PYPLOT_DPI with light PNG compression and shows the bytes with st.image. Same call
    signature as st.pyplot(fig, width=...). The page is ~1250 px wide, so 110 dpi (1870 px for the widest figures) loses nothing visible
    while rendering ~3x fewer pixels; measured page reruns got ~2x faster."""
    import io

    if getattr(st.pyplot, "_fast", False):
        return

    def pyplot(fig=None, clear_figure=None, width="content", **kw):
        import matplotlib.pyplot as plt

        fig = fig or plt.gcf()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=kw.pop("dpi", PYPLOT_DPI), bbox_inches="tight", pil_kwargs={"compress_level": 1}, **kw)
        st.image(buf.getvalue(), width=width if width in ("stretch", "content") or isinstance(width, int) else "content")
        if clear_figure:
            fig.clf()

    pyplot._fast = True
    st.pyplot = pyplot


def cached_figure(st, key: tuple, draw, width="stretch", keep: int = 30) -> None:
    """Show a matplotlib figure that depends only on `key` (e.g. run path + file time + settings): drawn once per session, then re-shown as
    PNG bytes. `draw()` returns the figure. Keeps the last `keep` figures per browser session."""
    import io

    import matplotlib.pyplot as plt

    cache = st.session_state.setdefault("_fig_cache", {})
    png = cache.pop(key, None)
    if png is None:
        fig = draw()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=PYPLOT_DPI, bbox_inches="tight", pil_kwargs={"compress_level": 1})
        plt.close(fig)
        png = buf.getvalue()
    cache[key] = png                                  # most recently used last
    while len(cache) > keep:
        cache.pop(next(iter(cache)))
    st.image(png, width=width)


def apply(st) -> None:
    """Inject the page CSS and set matplotlib defaults (call once, right after st.set_page_config)."""
    st.markdown(CSS, unsafe_allow_html=True)
    fast_pyplot(st)
    mpl.rcParams.update({
        "figure.facecolor": PALETTE["fig"], "axes.facecolor": PALETTE["fig"], "savefig.facecolor": PALETTE["fig"],
        "axes.edgecolor": "#3b4a63", "axes.labelcolor": PALETTE["ink"], "text.color": PALETTE["ink"],
        "xtick.color": PALETTE["muted"], "ytick.color": PALETTE["muted"], "legend.labelcolor": PALETTE["ink"],
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": PALETTE["grid"],
        "grid.linewidth": 0.8, "axes.axisbelow": True, "axes.titlesize": 17, "axes.titleweight": "bold",
        "axes.titlelocation": "left", "axes.labelsize": 14, "legend.frameon": False, "font.size": 14,
        "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 12,   # figures are ~14 in wide shown at ~70 dpi: small fonts were unreadable
        "axes.prop_cycle": mpl.cycler(color=CYCLE), "lines.linewidth": 2.0,
    })


def hero(st, title: str, subtitle: str = "") -> None:
    st.markdown(f'<div class="sl-hero"><h1>{html.escape(title)}</h1><p>{html.escape(subtitle)}</p></div>',
                unsafe_allow_html=True)


def intro(st, what: str, how: str) -> None:
    """Two short boxes at the top of a tab: WHAT it shows and HOW to use it (plain text; no HTML needed)."""
    st.markdown(f'<div class="sl-intro"><div class="what"><b class="t">What this shows</b><p>{html.escape(what)}</p></div>'
                f'<div class="how"><b class="t">How to use it</b><p>{html.escape(how)}</p></div></div>', unsafe_allow_html=True)


GUIDE = [("➕ New simulation (sidebar)", "Cells move, are segmented, tracked and scored. About 20 s."),
         ("🔬 New single-cell loci run (sidebar)", "Adds two coloured loci per nucleus (optional probes, MS2, 3D z-stacks and a calibration). 2 to 10 min."),
         ("🎞 Movie · 📊 Metrics · 🧭 Trajectories", "How well the tracker did, frame by frame, with the errors marked."),
         ("🔬 Cells · 📈 Scans", "Single-cell movies with their loci; sweeps over many seeds and settings."),
         ("🔎 Spot finding", "How each locus and the MS2 signal are found (Chen et al. 2025), every intermediate image, the error per frame, bleaching and illumination corrections, and a comparison of methods."),
         ("🖼 Images", "Every recorded image next to the ground truth behind it; pick ONE error source to see what it is, how it is corrected and what it does to the results; the image built up one error at a time."),
         ("📍 Stage & fiducials", "z-piezo plane errors, the beads on the slide that track drift and focus, and what correcting with them (or the piezo's log) does."),
         ("🎯 Calibration", "Simulated dark frames, flat-field slide and bead slide: what the instrument is, estimated vs true, and how much correcting for it helps."),
         ("⚡ Transcription · 🧪 Probes", "Build and explore the models behind the MS2 spot and the locus brightness."),
         ("🧮 MS2 → Pol II", "Turn each cell's MS2 trace into the moments Pol II started and into bursts (Chen et al. 2025 deconvolution), against the truth."),
         ("🗂 Library · 🧬 Chromatin", "The saved polymer simulations: progress, choosing them, and looking inside one.")]


def guide(st) -> None:
    """Collapsed 'start here' map of the whole app."""
    with st.expander("🧭 New here? A map of this app", expanded=False):
        st.markdown('<div class="sl-guide">' + "".join(f"<div><b>{html.escape(a)}</b>{html.escape(b)}</div>" for a, b in GUIDE)
                    + "</div>", unsafe_allow_html=True)
        st.caption("Typical path: (1) choose which saved simulations to use in 🗂 Library, (2) tune the models in ⚡ / 🧪, "
                   "(3) press Run in the sidebar, (4) read the results in the movie and metrics tabs. Every setting has a ⓘ tooltip.")


GLOSSARY = {
    "Ground truth": "What the simulation knows to be true (true cell masks and IDs, true locus positions). Real data never has it.",
    "Segmentation": "Finding which pixels belong to each nucleus in one image.",
    "Tracking": "Linking the nuclei found in each frame into tracks, so the same cell keeps the same ID over time.",
    "ID switch": "A true cell that changes track ID between two frames (a swap or a broken track).",
    "CHOTA / LNK / DET / TRA": "Cell Tracking Challenge scores: DET = detection, LNK = linking, TRA = overall tracking, CHOTA = detection and association "
                               "combined. 1.0 is perfect.",
    "Locus": "A labelled spot on the chromosome (here a promoter and an enhancer, each in its own colour).",
    "MS2 spot": "A bright spot at an active gene: many fluorescent MS2 loops sit on the RNAs being made. It is in the nuclear channel's colour.",
    "PSF": "Point spread function: the blur a microscope turns a point of light into. Sets the size of every spot.",
    "Photon / electron / ADU": "Light arrives as photons; the sensor turns some into electrons (quantum efficiency); the camera reports counts (ADU).",
    "Shot noise": "Random variation in the number of photons detected, unavoidable even for a steady light; grows as the square root of the signal.",
    "Read noise": "Electronic noise added every time a pixel is read; independent of the signal.",
    "Chromatic shift": "Different colours are imaged at slightly different positions, so two loci of different colours appear displaced.",
    "Flat-field": "How even the illumination is across the field of view.",
    "Library": "A collection of saved polymer (chromosome) simulations; each cell in a movie gets its own.",
    "Polymer unit": "The length unit of the polymer simulation until it is calibrated to nanometres and seconds (still a placeholder).",
    "z-stack": "One 3D picture made of several planes recorded at different heights, one after the other.",
    "Maximum projection": "A 2D picture made by keeping the brightest value along z at every position (what segmentation and tracking use in 3D runs).",
    "Axial": "Along z, the direction of the microscope's optical axis (up and down); lateral means in the image plane.",
    "Dark frame": "A picture taken with the shutter closed: shows only the camera's own offset, noise, hot pixels and cosmic rays.",
    "Flat field": "A picture of a uniformly glowing slide: shows the unevenness of the illumination and of the pixels' gains.",
    "Bead slide": "A slide of tiny fluorescent beads (smaller than the blur) used to measure the PSF, the colour offsets and the drift.",
    "Photon transfer curve": "Noise variance plotted against signal; its slope is the camera gain (counts per electron).",
    "DSNU / PRNU": "Fixed pixel-to-pixel differences: in the dark offset (DSNU) and in the gain (PRNU).",
    "Piezo": "The small motor that moves the focus up and down in steps to record the planes of a z-stack. Its real positions differ a little from the requested ones.",
    "Hysteresis": "A piezo ends up at a slightly different place when it arrives from below than when it arrives from above.",
    "Fiducial bead": "A bright tiny bead lying on the slide, visible in every image, used to measure how far the whole picture has drifted.",
    "Position log": "The piezo's own sensor reading of where each plane really was (a few nm of noise).",
    "Kymograph": "A picture with time along one axis and position along a line on the other.",
    "Pol II": "RNA polymerase II, the enzyme that copies a gene into RNA. Each Pol II on the gene carries the MS2 loops it has made, so it glows.",
    "Initiation (Pol II start)": "The moment one Pol II starts copying the gene. Bursts are periods in which they start often.",
    "Deconvolution": "Working backwards from a blurred-in-time signal (the MS2 trace: the sum of every Pol II's glow) to the events that made it.",
    "Posterior samples": "Many possible answers, each consistent with the data and its noise; where they agree the answer is certain.",
    "Kernel": "The signal one single event produces over time (here: one Pol II's glow from start to release).",
    "Self-calibration": "Estimating an instrument property (background, hot pixels, illumination, colour shift) from the movie itself instead of from calibration images.",
    "Error-free reference": "The same movie rendered by a perfect microscope that still has the unavoidable blur and photon / read noise.",
}


def glossary(st) -> None:
    """Collapsed list of the terms used in the dashboard."""
    with st.expander("📖 Glossary of terms", expanded=False):
        st.markdown("\n".join(f"- **{k}**: {v}" for k, v in GLOSSARY.items()))


def card(st, title: str, body: str, step: int | None = None) -> None:
    """A small bordered box with a heading and a short paragraph (body may contain simple HTML)."""
    num = f'<span class="sl-step">{step}</span>' if step else ""
    st.markdown(f'<div class="sl-card"><h4>{num}{html.escape(title)}</h4><p>{body}</p></div>', unsafe_allow_html=True)


def badge(st, text: str, kind: str = "idle") -> None:
    st.markdown(f'<span class="sl-badge {kind}">{html.escape(text)}</span>', unsafe_allow_html=True)


def stat_grid(st, items: list[tuple[str, str, str]]) -> None:
    """Compact two-column grid of (label, value, tooltip) numbers; unlike st.metric it never truncates text."""
    cells = "".join(f'<div class="sl-stat" title="{html.escape(tip)}"><b>{html.escape(val)}</b><span>{html.escape(lab)}</span></div>'
                    for lab, val, tip in items)
    st.markdown(f'<div class="sl-stats">{cells}</div>', unsafe_allow_html=True)
