"""One look for the whole dashboard: page CSS, matplotlib defaults, and small layout helpers.

Colours are defined once here (PALETTE). Matplotlib figures use the same palette so plots match the page.
"""
import html

import matplotlib as mpl

PALETTE = {"ink": "#1f2933", "muted": "#64748b", "line": "#dfe4e1", "card": "#ffffff", "bg": "#fbfbf9",
           "accent": "#0f766e", "accent_soft": "#d9efec", "good": "#15803d", "warn": "#b45309", "bad": "#b91c1c",
           "blue": "#2563eb", "violet": "#7c3aed"}
CYCLE = ["#0f766e", "#2563eb", "#b45309", "#7c3aed", "#b91c1c", "#475569", "#0891b2", "#a16207"]

CSS = f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
:root {{ --ink:{PALETTE['ink']}; --muted:{PALETTE['muted']}; --line:{PALETTE['line']}; --accent:{PALETTE['accent']}; }}
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
[data-testid="stMetric"] {{ background:#fff; border:1px solid var(--line); border-radius:.8rem; padding:.8rem 1rem;
  box-shadow: 0 1px 2px rgba(15,23,42,.04); }}
[data-testid="stMetricLabel"] {{ color: var(--muted); }}
[data-testid="stMetricValue"] {{ font-weight: 700; letter-spacing:-0.02em; }}

/* expanders, forms, dataframes, images get the same card edge */
[data-testid="stExpander"] {{ border:1px solid var(--line); border-radius:.8rem; background:#fff; }}
[data-testid="stForm"] {{ border:1px solid var(--line); border-radius:.8rem; background:#fff; }}
[data-testid="stDataFrame"], [data-testid="stImage"] img, [data-testid="stPyplot"] img {{ border-radius:.6rem; }}
[data-testid="stSidebar"] {{ background:#f1f3f1; border-right:1px solid var(--line); }}
.stButton > button, .stFormSubmitButton > button {{ border-radius:.6rem; font-weight:600; }}

/* helper components */
.sl-hero {{ margin: 0 0 1.1rem 0; }}
.sl-hero p {{ color: var(--muted); margin:.1rem 0 0 0; font-size:1.02rem; }}
.sl-card {{ background:#fff; border:1px solid var(--line); border-radius:.9rem; padding:1rem 1.1rem; min-height:12.5rem;
  box-shadow: 0 1px 2px rgba(15,23,42,.04); }}
.sl-card h4 {{ margin:.1rem 0 .35rem 0; font-size:1rem; }}
.sl-card p {{ margin:0; color:var(--muted); font-size:.92rem; line-height:1.45; }}
.sl-step {{ display:inline-flex; width:1.7rem; height:1.7rem; border-radius:50%; background:{PALETTE['accent_soft']};
  color:var(--accent); font-weight:700; align-items:center; justify-content:center; margin-right:.5rem; }}
.sl-badge {{ display:inline-block; padding:.15rem .65rem; border-radius:999px; font-size:.82rem; font-weight:650; }}
.sl-badge.good {{ background:#dcfce7; color:#166534; }} .sl-badge.warn {{ background:#fef3c7; color:#92400e; }}
.sl-badge.bad {{ background:#fee2e2; color:#991b1b; }} .sl-badge.idle {{ background:#e5e7eb; color:#374151; }}
.sl-badge.run {{ background:{PALETTE['accent_soft']}; color:#115e59; }}
.sl-stats {{ display:grid; grid-template-columns: repeat(2, 1fr); gap:.6rem; margin:.5rem 0 .2rem 0; }}
.sl-stat {{ background:#fff; border:1px solid var(--line); border-radius:.7rem; padding:.55rem .8rem; }}
.sl-stat b {{ display:block; font-size:1.35rem; letter-spacing:-0.02em; line-height:1.3; }}
.sl-stat span {{ color:var(--muted); font-size:.8rem; }}
.sl-mono {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size:.85rem; }}

/* ---- polish layer ---- */
.stApp, .stMarkdown, p, label, input, textarea, button, [data-baseweb="tab"] {{ font-family: 'Inter', system-ui, -apple-system, 'Segoe UI', sans-serif; }}
[data-testid="stIconMaterial"], [class*="material-symbols"], .material-icons {{ font-family: 'Material Symbols Rounded', 'Material Icons' !important; }}
.stApp {{ background: linear-gradient(180deg, #f7faf9 0%, #fbfbf9 260px); }}
.sl-hero {{ border-left: 5px solid var(--accent); padding: .35rem 0 .35rem 1rem; margin-bottom: 1.4rem; }}
.sl-hero h1 {{ background: linear-gradient(90deg, #0f766e, #2563eb); -webkit-background-clip: text; background-clip: text;
  -webkit-text-fill-color: transparent; }}
[data-baseweb="tab"] {{ transition: background .15s, color .15s; }}
[data-baseweb="tab"]:hover {{ background: {PALETTE['accent_soft']}; color: var(--accent); }}
[data-testid="stMetric"], [data-testid="stExpander"], [data-testid="stForm"], .sl-card, .sl-stat {{
  transition: box-shadow .15s, transform .15s; }}
[data-testid="stMetric"]:hover, .sl-card:hover {{ box-shadow: 0 6px 18px rgba(15,23,42,.08); transform: translateY(-1px); }}
[data-testid="stExpander"] summary {{ font-weight: 600; }}
.stButton > button, .stFormSubmitButton > button {{ transition: box-shadow .15s, transform .15s; }}
.stButton > button:hover, .stFormSubmitButton > button:hover {{ box-shadow: 0 4px 12px rgba(15,118,110,.25); transform: translateY(-1px); }}
.stButton > button[kind="primary"], .stFormSubmitButton > button[kind="primary"] {{
  background: linear-gradient(135deg, #0f766e, #0e7490); border: 0; }}
[data-testid="stAlert"] {{ border-radius: .8rem; }}
/* compact metric cards (the headline scores were too large) */
[data-testid="stMetric"] {{ padding: .35rem .7rem; border-radius: .6rem; }}
[data-testid="stMetricLabel"] p {{ font-size: .78rem; }}
[data-testid="stMetricValue"] {{ font-size: 1.35rem !important; line-height: 1.25; }}
[data-testid="stMetric"]:hover {{ transform: none; }}
.stMarkdown strong {{ font-weight: 600; }}
[data-baseweb="input"], [data-baseweb="select"] > div {{ border-radius: .55rem; }}
[data-testid="stSidebar"] h2 {{ font-size: 1.25rem; letter-spacing: -0.01em; }}
[data-testid="stSidebar"] [data-testid="stExpander"] {{ background: #fff; }}
.sl-intro {{ display:grid; grid-template-columns: 1fr 1fr; gap: .8rem; margin: .2rem 0 1.2rem 0; }}
.sl-intro > div {{ border-radius: .85rem; padding: .85rem 1.05rem; border: 1px solid var(--line); background:#fff; }}
.sl-intro .what {{ border-left: 4px solid var(--accent); }}
.sl-intro .how {{ border-left: 4px solid {PALETTE['blue']}; }}
.sl-intro b.t {{ display:block; font-size:.74rem; letter-spacing:.08em; text-transform:uppercase; margin-bottom:.25rem; }}
.sl-intro .what b.t {{ color: var(--accent); }} .sl-intro .how b.t {{ color: {PALETTE['blue']}; }}
.sl-intro p {{ margin:0; color: var(--ink); font-size:.93rem; line-height:1.5; }}
@media (max-width: 900px) {{ .sl-intro {{ grid-template-columns: 1fr; }} }}
.sl-guide {{ display:grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap:.6rem; }}
.sl-guide div {{ background:#fff; border:1px solid var(--line); border-radius:.75rem; padding:.6rem .8rem; font-size:.88rem;
  line-height:1.4; color: var(--muted); }}
.sl-guide b {{ color: var(--ink); display:block; font-size:.95rem; }}
</style>
"""


def apply(st) -> None:
    """Inject the page CSS and set matplotlib defaults (call once, right after st.set_page_config)."""
    st.markdown(CSS, unsafe_allow_html=True)
    mpl.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "axes.edgecolor": PALETTE["line"], "axes.labelcolor": PALETTE["ink"], "text.color": PALETTE["ink"],
        "xtick.color": PALETTE["muted"], "ytick.color": PALETTE["muted"],
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e8ece9",
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
         ("🔬 New single-cell loci run (sidebar)", "Adds two coloured loci per nucleus (and optional probes / MS2). About 2 min."),
         ("🎞 Movie · 📊 Metrics · 🧭 Trajectories", "How well the tracker did, frame by frame, with the errors marked."),
         ("🔬 Cells · 📈 Scans", "Single-cell movies with their loci; sweeps over many seeds and settings."),
         ("⚡ Transcription · 🧪 Probes", "Build and explore the models behind the MS2 spot and the locus brightness."),
         ("🗂 Library · 🧬 Chromatin", "The saved polymer simulations: progress, choosing them, and looking inside one.")]


def guide(st) -> None:
    """Collapsed 'start here' map of the whole app."""
    with st.expander("🧭 New here? A map of this app", expanded=False):
        st.markdown('<div class="sl-guide">' + "".join(f"<div><b>{html.escape(a)}</b>{html.escape(b)}</div>" for a, b in GUIDE)
                    + "</div>", unsafe_allow_html=True)
        st.caption("Typical path: (1) choose which saved simulations to use in 🗂 Library, (2) tune the models in ⚡ / 🧪, "
                   "(3) press Run in the sidebar, (4) read the results in the movie and metrics tabs. Every setting has a ⓘ tooltip.")


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
