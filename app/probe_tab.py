"""The dashboard's 🧪 Probes tab: the stochastic probe-binding model at the promoter and the enhancer.

Pick parameters (with a live example), see the model drawn for each locus with its probe counts, sweep a parameter, and
send the chosen parameters to a movie run. Reads/writes nothing on disk except the optional run view
(<run>/stage1_chromatin/probe_occupancy.csv). Imports only numpy / pandas / numba (via simlive), never torch.
"""
import html
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

import style
from simlive.stage1_chromatin import probes as P

KEYS = ("n_probes", "k_bind", "k_scan", "k_unbind", "k_off")
LOCUS_LABEL = {"promoter": "Promoter (locus 0)", "enhancer": "Enhancer (locus 1)"}
PRESETS = {
    "Placeholder default (stays ~10 s)": {},
    "Short-lived binder (stays ~2 s)": {"k_scan": 1.0, "k_off": 0.6, "k_unbind": 0.3},
    "Long-lived binder (stays ~60 s)": {"k_scan": 0.3, "k_off": 0.02, "k_unbind": 0.05},
    "Few probes (noisy signal)": {"n_probes": 8},
    "Many probes (smooth signal)": {"n_probes": 200, "k_bind": 0.01},
    "Saturating locus (20 sites)": {"capacity": 20, "k_bind": 0.4, "n_probes": 60},
}
COL = {"free": "#94a3b8", "bound": "#2563eb", "scan": "#0f766e", "bad": "#b45309"}


def _k(locus: str, name: str) -> str:
    return f"pr_{locus}_{name}"


def _apply_preset(locus: str) -> None:
    """Callback: fill the locus' widgets from the chosen preset (on top of the placeholder defaults)."""
    vals = {**P.DEFAULTS[locus], **PRESETS[st.session_state[_k(locus, "preset")]]}
    for n in KEYS:
        st.session_state[_k(locus, n)] = vals[n]
    st.session_state[_k(locus, "cap_on")] = vals["capacity"] is not None
    st.session_state[_k(locus, "cap")] = int(vals["capacity"] or 20)


def _init(locus: str) -> None:
    d = P.DEFAULTS[locus]
    for n in KEYS:
        st.session_state.setdefault(_k(locus, n), d[n])
    st.session_state.setdefault(_k(locus, "cap_on"), False)
    st.session_state.setdefault(_k(locus, "cap"), 20)


def current(locus: str) -> dict:
    """The parameters currently set in the tab for one locus (dict, defaults if the tab was never opened)."""
    d = P.DEFAULTS[locus]
    out = {n: st.session_state.get(_k(locus, n), d[n]) for n in KEYS}
    out["capacity"] = int(st.session_state.get(_k(locus, "cap"), 20)) if st.session_state.get(_k(locus, "cap_on")) else None
    return out


def overrides_for_movie() -> list[str]:
    """`--set` overrides that switch the probe model on in a movie run with the settings from this tab."""
    ov = ["loci.probes.enabled=true", f"loci.probes.photons_per_probe_s={st.session_state.get('pr_photons', 300)}"]
    for loc in P.LOCUS_NAMES:
        for k, v in current(loc).items():
            ov.append(f"loci.probes.{loc}.{k}={'null' if v is None else v}")
    return ov


# ----------------------------------------------------------------------------- the picture of the model
def model_svg(p: P.ProbeParams, counts: tuple[int, int, int], mean: tuple[float, float, float], when: str) -> str:
    """The three states as boxes with the probes drawn as dots, the rates on the arrows, counts and averages inside."""
    f, b, s = counts
    boxes = [("FREE", "floating", f, mean[0], COL["free"], 20),
             ("BOUND", "just landed", b, mean[1], COL["bound"], 290), ("SCANNING", "sliding", s, mean[2],
                                                                      COL["scan"], 560)]
    out = ['<svg viewBox="0 0 780 290" xmlns="http://www.w3.org/2000/svg" style="width:100%;font-family:system-ui,sans-serif">',
           '<defs><marker id="ar" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">'
           '<path d="M0,0 L10,5 L0,10 z" fill="#475569"/></marker></defs>']
    for name, sub, n, m, color, x in boxes:
        out.append(f'<rect x="{x}" y="40" width="200" height="150" rx="14" fill="{color}" fill-opacity="0.10" stroke="{color}" '
                   f'stroke-width="2"/>')
        out.append(f'<text x="{x + 12}" y="62" font-size="14" font-weight="700" fill="{color}">{name}</text>')
        out.append(f'<text x="{x + 12}" y="79" font-size="11" fill="#64748b">{sub}</text>')
        out.append(f'<text x="{x + 188}" y="66" font-size="30" font-weight="700" text-anchor="end" fill="#1f2933">{n}</text>')
        out.append(f'<text x="{x + 188}" y="82" font-size="11" text-anchor="end" fill="#64748b">average {m:.1f}</text>')
        if p.n_probes <= 80:
            for i in range(n):
                out.append(f'<circle cx="{x + 18 + (i % 12) * 14.5}" cy="{104 + (i // 12) * 14.5}" r="5" fill="{color}"/>')
    arrows = [("M224,100 L286,100", 255, 92, f"k_bind {p.k_bind:g}"), ("M494,100 L556,100", 525, 92, f"k_scan {p.k_scan:g}")]
    for d, x, y, t in arrows:
        out.append(f'<path d="{d}" stroke="#475569" stroke-width="2.2" fill="none" marker-end="url(#ar)"/>')
        out.append(f'<text x="{x}" y="{y}" font-size="11" text-anchor="middle" fill="#334155">{t}</text>')
    out.append('<path d="M660,192 C660,262 120,262 120,194" stroke="#475569" stroke-width="2.2" fill="none" marker-end="url(#ar)"/>')
    out.append(f'<text x="390" y="274" font-size="12" text-anchor="middle" fill="#334155">falls off: k_off {p.k_off:g} per s '
               f'(scans {1 / p.k_off:.0f} s on average)</text>')
    out.append('<path d="M390,192 C390,222 150,222 150,194" stroke="#b45309" stroke-width="1.6" fill="none" stroke-dasharray="5 4" '
               'marker-end="url(#ar)"/>')
    out.append(f'<text x="390" y="236" font-size="11" text-anchor="middle" fill="#b45309">lets go again: k_unbind {p.k_unbind:g}</text>')
    cap = f" · {p.capacity} binding sites" if p.capacity else ""
    out.append(f'<text x="20" y="22" font-size="12" fill="#64748b">{p.n_probes} probes in the pool{cap} · counts at {html.escape(when)}</text>')
    out.append("</svg>")
    return "".join(out)


# ----------------------------------------------------------------------------- one locus panel
def _locus_panel(locus: str) -> P.ProbeParams | None:
    _init(locus)
    st.markdown(f"#### {LOCUS_LABEL[locus]}")
    left, right = st.columns([5, 7], gap="large")
    left.selectbox("Start from a preset", list(PRESETS), key=_k(locus, "preset"), on_change=_apply_preset, args=(locus,),
                 help="Fills the boxes below; change any number afterwards.")
    c = left.columns(2)
    c[0].number_input("Probes in the pool", 1, 2000, key=_k(locus, "n_probes"), help=P.PARAM_HELP["n_probes"])
    c[1].number_input("k_bind (1/s)", 0.0001, 100.0, key=_k(locus, "k_bind"), step=0.01, format="%.4g", help=P.PARAM_HELP["k_bind"])
    c[0].number_input("k_scan (1/s)", 0.0001, 100.0, key=_k(locus, "k_scan"), step=0.05, format="%.4g", help=P.PARAM_HELP["k_scan"])
    c[1].number_input("k_unbind (1/s)", 0.0, 100.0, key=_k(locus, "k_unbind"), step=0.05, format="%.4g", help=P.PARAM_HELP["k_unbind"])
    c[0].number_input("k_off (1/s)", 0.0001, 100.0, key=_k(locus, "k_off"), step=0.01, format="%.4g", help=P.PARAM_HELP["k_off"])
    with c[1]:
        st.checkbox("Limit the binding sites", key=_k(locus, "cap_on"), help=P.PARAM_HELP["capacity"])
        st.number_input("Binding sites", 1, 5000, key=_k(locus, "cap"), disabled=not st.session_state[_k(locus, "cap_on")],
                        label_visibility="collapsed")
    try:
        p = P.ProbeParams.from_dict(current(locus), locus)
    except ValueError as e:
        st.error(str(e))
        return None
    sm = P.summary(p)
    with left:
        style.stat_grid(st, [("probes attached (average)", f"{sm['mean_attached']:.1f}", P.PARAM_HELP["n_probes"]),
                             ("signal noise (sd / mean)", f"{100 * sm['cv_attached']:.0f}%", "relative fluctuation of the attached count"),
                             ("a probe stays attached", f"{sm['residence_s']:.1f} s", "average time per visit"),
                             ("count fluctuates over", f"{sm['correlation_time_s']:.1f} s", "how long before the count forgets its value")])
    if p.capacity and sm["mean_attached"] > 0.7 * p.capacity:
        left.caption("⚠ The site limit is close to the average count, so the exact numbers above overestimate it; the example "
                   "trace below includes the limit.")

    tr_key = f"pr_{locus}_seed"
    seed = st.session_state.setdefault(tr_key, 1)
    window = st.session_state.get("pr_window", 300)
    tr = P.simulate_trace(p, window, seed)
    ts = np.arange(0, window + 1e-9, 1.0)
    d = P.sample_trace(tr, ts, 1e-3)
    when = right.slider("Look at time (s) in the example trace", 0, int(window), 0, key=f"pr_{locus}_when")
    i = min(int(when), len(d) - 1)
    mean = (sm["mean_free"], sm["mean_bound"], sm["mean_scanning"])
    right.markdown(model_svg(p, (int(d.n_free[i]), int(d.n_bound[i]), int(d.n_scanning[i])), mean, f"{when} s"), unsafe_allow_html=True)
    if right.button("New random example", key=f"pr_{locus}_new"):
        st.session_state[tr_key] = seed + 1
        st.rerun(scope="fragment")

    fig, ax = plt.subplots(1, 2, figsize=(14, 3.1), gridspec_kw={"width_ratios": [4, 1]})
    ax[0].fill_between(ts, 0, d.n_scanning, step="post", color=COL["scan"], alpha=0.35, label="scanning", lw=0)
    ax[0].fill_between(ts, d.n_scanning, d.n_attached, step="post", color=COL["bound"], alpha=0.35, label="bound", lw=0)
    ax[0].step(ts, d.n_attached, where="post", color="#1f2933", lw=1.3, label="attached (what you see)")
    ax[0].axhline(sm["mean_attached"], color=COL["bad"], ls="--", lw=1.2, label="theory average")
    ax[0].axvline(when, color="#64748b", lw=1)
    ax[0].set_xlabel("time (s)"); ax[0].set_ylabel("probes at the locus"); ax[0].set_title("Probes attached through time")
    ax[0].legend(ncol=4, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.28))
    long = P.sample_trace(P.simulate_trace(p, 2000, seed + 7), np.arange(0, 2000), 1e-3).n_attached
    ax[1].hist(long, bins=np.arange(long.min(), long.max() + 2) - 0.5, density=True, color=COL["bound"], alpha=0.6)
    ax[1].set_xlabel("probes attached"); ax[1].set_title("How often each count occurs", fontsize=11)
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    return p


# ----------------------------------------------------------------------------- sweep
def _sweep_section() -> None:
    st.markdown("### Parameter sweep")
    st.caption("Run the model over a range of one parameter (fast, seconds) and see how the signal changes. Simulated points sit "
               "on the exact theory line when the model is behaving; the same sweep can be run from a terminal.")
    c = st.columns(4)
    locus = c[0].selectbox("Locus", P.LOCUS_NAMES, format_func=LOCUS_LABEL.get, key="pr_sw_locus")
    param = c[1].selectbox("Parameter to vary", list(KEYS), index=4, key="pr_sw_param")
    cur = float(current(locus)[param]) or 0.1                      # the range defaults to 0.2x .. 5x the current value
    kk = f"{locus}_{param}"
    lo = c[2].number_input("From", 0.0, 100000.0, max(cur / 5, 1.0) if param == "n_probes" else cur / 5, format="%.4g", key=f"pr_sw_lo_{kk}")
    hi = c[3].number_input("To", 0.0, 100000.0, cur * 5, format="%.4g", key=f"pr_sw_hi_{kk}")
    c = st.columns(4)
    n = int(c[0].number_input("Number of values", 2, 25, 6, key="pr_sw_n"))
    log = c[1].checkbox("Space them logarithmically", True, key="pr_sw_log")
    reps = int(c[2].number_input("Repeats per value", 1, 20, 3, key="pr_sw_reps"))
    dur = float(c[3].number_input("Seconds simulated per repeat", 200.0, 100000.0, 3000.0, 200.0, key="pr_sw_dur"))
    st.caption(f"The other settings stay at what is set above for the {locus} (currently {param} = {cur:g}).")
    go = st.button("▶ Run sweep", type="primary", key="pr_sw_go")
    if go:
        if lo <= 0 and log:
            st.error("A logarithmic range must start above 0.")
        elif hi <= lo:
            st.error("'To' must be larger than 'From'.")
        else:
            vals = np.geomspace(lo, hi, n) if log else np.linspace(lo, hi, n)
            vals = np.unique(np.maximum(np.round(vals), 1).astype(int)) if param == "n_probes" else vals
            with st.spinner("Simulating…"):
                try:
                    st.session_state["pr_sweep_res"] = (locus, param, P.sweep(current(locus), {param: [float(v) if param != "n_probes" else int(v) for v in vals]},
                                                                              locus, reps, dur, seed=1))
                except ValueError as e:
                    st.error(str(e))
    res = st.session_state.get("pr_sweep_res")
    if not res:
        return
    locus, param, df = res
    g = df.groupby(param).agg(**{c: (c, "mean") for c in df.columns if c.startswith(("th_", "sim_"))}).reset_index()
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.6))
    for a, (th, sim, title) in zip(ax, [("th_mean_attached", "sim_mean_attached", "Probes attached (average)"),
                                         ("th_cv_attached", "sim_cv", "Signal noise (CV = sd / mean)"),
                                         ("th_p_none_attached", "sim_p_none", "Fraction of time with no probe attached")]):
        a.plot(g[param], g[th], "-", color=COL["bad"], label="exact theory")
        a.scatter(df[param], df[sim], s=18, color=COL["bound"], alpha=0.6, label="simulation (each repeat)")
        a.set_xlabel(param); a.set_title(title, fontsize=11)
        if g[param].max() / max(g[param].min(), 1e-12) > 20:
            a.set_xscale("log")
    ax[0].legend(fontsize=8)
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    with st.expander("Table and download"):
        st.dataframe(g.round(4), hide_index=True, width="stretch")
        st.download_button("Download all repeats (CSV)", df.to_csv(index=False), f"probe_sweep_{locus}_{param}.csv", "text/csv")
    vals = ",".join(f"{v:.4g}" for v in g[param])
    st.markdown("**Run the same sweep from a terminal**, or push it through the whole movie pipeline to see how binding noise "
                "changes the locus analysis:")
    st.code(f"scripts\\run.cmd python scripts\\sweep_probes.py --name my_sweep --locus {locus} --grid {param}={vals}\n"
            f"scripts\\run.cmd python scripts\\sweep.py --name probes_{param} --config configs\\loci_demo.yaml "
            f"--set loci.probes.enabled=true --grid loci.probes.{locus}.{param}={vals} --seeds 1-5", language="bash")


# ----------------------------------------------------------------------------- tab
@st.fragment
def render(run_path: Path | None = None) -> None:
    style.intro(st, "A stochastic model of fluorescent probes landing on, sliding along and leaving a locus. The number attached sets the locus brightness, so binding noise appears in the movie.", "Pick a preset or change the rates for the promoter and enhancer, scrub the example trace, then run a sweep to see how one rate changes the signal. Tick the probe option in the sidebar loci form to use these settings in a movie.")
    st.markdown("### How probes bind a locus")
    c = st.columns(4)
    with c[0]:
        style.card(st, "Free", "Fluorescent probes float around in the nucleus. A free probe lands on the locus at rate "
                   "<b>k_bind</b>.", 1)
    with c[1]:
        style.card(st, "Bind", "A probe that has just landed is <b>bound</b>. It either starts scanning (<b>k_scan</b>) or "
                   "lets go again (<b>k_unbind</b>).", 2)
    with c[2]:
        style.card(st, "Scan", "A scanning probe slides along the locus. It is still fluorescent, so it counts towards the "
                   "brightness.", 3)
    with c[3]:
        style.card(st, "Fall off", "After about 1 / <b>k_off</b> seconds the probe falls off and becomes free again. "
                   "Then the cycle repeats.", 4)
    st.write("")
    st.caption("Each probe cycles independently and at random (exact Gillespie simulation), so the number attached "
               "fluctuates. In a movie the locus brightness is proportional to the number attached, which is what adds "
               "binding noise to the locus signal. The promoter and the enhancer have separate settings. All rates are "
               "placeholders until measured kinetics are supplied.")

    sc = st.columns([2, 2, 3])
    sc[0].number_input("Photons per attached probe per s", 1, 100000, 300, 50, key="pr_photons",
                       help="Brightness of one probe. Used only in movies (the locus gets this x probes attached).")
    sc[1].number_input("Example trace length (s)", 60, 3000, 300, 60, key="pr_window")
    for locus in P.LOCUS_NAMES:
        with st.container(border=True):
            _locus_panel(locus)
    _sweep_section()

    st.markdown("### Use these settings in a movie")
    st.markdown("In the sidebar, open **New single-cell loci run** and tick **Model probe binding at the two loci**. "
                "Or from a terminal:")
    st.code("scripts\\run.cmd python scripts\\run_tracking_demo.py --config configs\\loci_demo.yaml " +
            " ".join(f'--set "{o}"' for o in overrides_for_movie()), language="bash")

    occ_f = run_path / "stage1_chromatin" / "probe_occupancy.csv" if run_path else None
    if occ_f is not None and occ_f.exists():
        st.markdown("### In the selected run")
        occ = pd.read_csv(occ_f)
        cid = st.selectbox("Cell", sorted(occ.cell_id.unique()), key="pr_cell")
        fig, ax = plt.subplots(figsize=(12, 3.2))
        for loc, color in zip(P.LOCUS_NAMES, (COL["scan"], COL["bound"])):
            g = occ[(occ.cell_id == cid) & (occ.locus_name == loc)]
            if g.empty:
                continue
            ax.step(g.t, g.mean_attached, where="mid", color=color, label=f"{loc} (as imaged)")
            ax.axhline(P.summary(P.ProbeParams.from_dict(current(loc), loc))["mean_attached"], color=color, ls=":", lw=1)
        ax.set_xlabel("frame"); ax.set_ylabel("probes attached"); ax.legend()
        ax.set_title("Probes attached per frame in this run (dotted = theory average for the settings above)")
        fig.tight_layout()
        st.pyplot(fig, width="stretch")
        plt.close(fig)
