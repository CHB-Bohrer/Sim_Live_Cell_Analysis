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
        "grid.linewidth": 0.8, "axes.axisbelow": True, "axes.titlesize": 13, "axes.titleweight": "bold",
        "axes.titlelocation": "left", "axes.labelsize": 11, "legend.frameon": False, "font.size": 10.5,
        "axes.prop_cycle": mpl.cycler(color=CYCLE), "lines.linewidth": 2.0,
    })


def hero(st, title: str, subtitle: str = "") -> None:
    st.markdown(f'<div class="sl-hero"><h1>{html.escape(title)}</h1><p>{html.escape(subtitle)}</p></div>',
                unsafe_allow_html=True)


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
