"""The dashboard's ⚡ Transcription tab: promoter-enhancer distance -> promoter switching -> Pol II -> MS2 signal -> the picture.

Every step of stage 1c is shown with its own plot, in the order the model works. Distances come from a saved polymer simulation
(or an invented stand-in when no library exists yet). All rates are per second of real time, so the frame interval can be changed
freely; the model is only sampled at the frame times. Imports numpy / pandas / numba / scipy through simlive, never torch.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

import library_tab
import style
from simlive.calibration.polymer import default_scale
from simlive.stage1_chromatin import library as L
from simlive.stage1_chromatin import transcription as T
from simlive.stage3_microscopy.render import _spot_photons

C = {"dist": "#60a5fa", "contact": "#14b8a6", "on": "#f59e0b", "ms2": "#4ade80", "theory": "#e5e7eb", "grey": "#94a3b8"}
KIND = {"contact": "Sharp contact (on / off)", "hill": "Smooth (Hill curve)", "exp": "Smooth (exponential)", "none": "None (control)"}
RATE = {"k_on": "turning ON (k_on)", "k_init": "initiating Pol II while ON (k_init)", "k_off": "turning OFF (k_off): contact keeps it ON longer"}


def _set(key, default):
    return st.session_state.setdefault(key, default)


# ----------------------------------------------------------------------------- pictures
def overview_svg() -> str:
    dots = "".join(f'<circle cx="{512 + i * 15}" cy="118" r="5" fill="#4ade80" opacity="0.85"/>' for i in range(8))
    return f"""<svg viewBox="0 0 800 270" xmlns="http://www.w3.org/2000/svg" style="width:100%;font-family:system-ui,sans-serif">
<defs><marker id="a2" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="#94a3b8"/></marker></defs>
<path d="M30,160 L770,160" stroke="#94a3b8" stroke-width="3"/>
<rect x="90" y="140" width="46" height="40" rx="8" fill="#a78bfa" fill-opacity=".85"/><text x="113" y="203" font-size="12" text-anchor="middle" fill="#cbd5e1">enhancer</text>
<text x="113" y="218" font-size="11" text-anchor="middle" fill="#94a3b8">(locus 1)</text>
<rect x="430" y="140" width="46" height="40" rx="8" fill="#f59e0b" fill-opacity=".9"/><text x="453" y="203" font-size="12" text-anchor="middle" fill="#cbd5e1">promoter</text>
<text x="453" y="218" font-size="11" text-anchor="middle" fill="#94a3b8">(locus 0)</text>
<path d="M113,138 C113,78 453,78 453,138" stroke="#60a5fa" stroke-width="2.4" fill="none" stroke-dasharray="6 5"/>
<text x="283" y="34" font-size="13" text-anchor="middle" fill="#60a5fa" font-weight="700">distance d  (from the polymer simulation)</text>
<text x="283" y="54" font-size="11" text-anchor="middle" fill="#94a3b8">small d  →  promoter switches ON more often</text>
<path d="M476,160 L760,160" stroke="#4ade80" stroke-width="5"/>
<text x="560" y="148" font-size="11" text-anchor="middle" fill="#4ade80">MS2 cassette (loops)</text>
{dots}
<text x="660" y="148" font-size="11" text-anchor="middle" fill="#94a3b8">gene body</text>
<g fill="#94a3b8"><circle cx="520" cy="180" r="7"/><circle cx="600" cy="180" r="7"/><circle cx="700" cy="180" r="7"/></g>
<text x="650" y="212" font-size="11" text-anchor="middle" fill="#94a3b8">Pol II travel along the gene</text>
<text x="400" y="256" font-size="12" text-anchor="middle" fill="#cbd5e1">MS2 brightness = loops carried by all Pol II on the gene → a bright spot at the promoter</text>
</svg>"""


def shade_contact(ax, t, contact, color=C["contact"]):
    """Shade the time spans where `contact` (bool per sample) is True."""
    if not len(t) or not contact.any():
        return
    edges = np.flatnonzero(np.diff(np.concatenate([[0], contact.astype(int), [0]])))
    dt = np.median(np.diff(t)) if len(t) > 1 else 1.0
    for a, b in zip(edges[::2], edges[1::2]):
        ax.axvspan(t[a], t[min(b, len(t) - 1)] + dt * (b >= len(t)), color=color, alpha=0.14, lw=0)


# ----------------------------------------------------------------------------- sections
def _step(n: int, title: str, text: str) -> None:
    st.markdown(f"### {n}. {title}")
    st.caption(text)


def _distance_source(dt: float, minutes: float):
    """Step 2: returns (d_nm array, dt_block_s, label, loops_warning)."""
    libs = library_tab.libraries()
    opts = (["library"] if libs else []) + ["standin"]
    src = st.radio("Where does the distance come from?", opts, horizontal=True, key="tx_src",
                   format_func={"library": "A saved polymer simulation", "standin": "Invented stand-in (not from the polymer simulations)"}.get)
    if src == "library":
        c = st.columns(4)
        lib = c[0].selectbox("Library", libs, key="tx_lib")
        df = L.list_library(lib)
        if df.empty:
            st.warning("This library has no finished simulations yet.")
            return None
        sim = c[1].selectbox("Simulation", df.sort_values("seed").sim_id.tolist(), key="tx_sim")
        n_kb = int(df[df.sim_id == sim].n_monomers.iloc[0])      # 1 monomer = 1 kb; defaults sit 100 kb apart in the middle of THIS region
        mid = (n_kb // 2) // 10 * 10
        prom = int(c[2].number_input(f"Promoter position (kb, 0-{n_kb - 1})", 0, n_kb - 1, mid + 50, 10, key=f"tx_prom_{n_kb}"))
        enh = int(c[3].number_input(f"Enhancer position (kb, 0-{n_kb - 1})", 0, n_kb - 1, mid - 50, 10, key=f"tx_enh_{n_kb}"))
        c = st.columns(3)
        nm0, blk0, scale_src = default_scale()
        nm = float(c[0].number_input("nm per polymer unit", 1.0, 1000.0, round(nm0, 1), 5.0, key="tx_nm",
                                     help=f"Default from {scale_src}: assumes a 100 kb pair is within 400 nm 10% of the time."))
        blk = float(c[1].number_input("Seconds per saved block", 0.01, 3600.0, round(blk0, 3), 0.1, key="tx_blk", format="%.3f",
                                      help=f"Default from {scale_src}: assumes a single locus has D = 0.002 um^2/s at a 1 min lag."))
        c[2].metric("Genomic separation", f"{abs(prom - enh)} kb")
        st.session_state["tx_movie_loci"] = {"library": lib, "positions_kb": [prom, enh], "nm_per_unit": nm, "block_duration_s": blk}   # read by the sidebar
        traj, beads = L.locus_trajectories(df[df.sim_id == sim].path.iloc[0], [prom, enh])
        d = np.linalg.norm(traj[:, 0] - traj[:, 1], axis=1) * nm
        st.caption(f"Monomers actually used: promoter {beads[0]} kb, enhancer {beads[1]} kb (positions snap to the stored 10 kb grid). "
                   f"The saved simulation covers {len(d) * blk / 60:.0f} min at {blk:g} s per block.")
        if abs(int(beads[0]) - int(beads[1])) == 0:
            st.error("Promoter and enhancer snap to the same monomer; move them apart.")
        return d, blk, f"{sim}", len(d) * blk < minutes * 60
    c = st.columns(3)
    frac = c[0].slider("Fraction of time in contact", 0.0, 0.6, 0.1, 0.01, key="tx_sfrac")
    loop = c[1].number_input("Mean time per contact (s)", 10.0, 3600.0, 240.0, 30.0, key="tx_sloop")
    seed = int(c[2].number_input("Random seed", 0, 10000, 1, key="tx_sseed"))
    blk = 5.0
    d = T.standin_distance(int(minutes * 60 / blk) + 10, blk, np.random.default_rng(seed), frac, loop)
    st.caption("This distance is INVENTED (a looped/open switch with confined wiggling); it is not from polymer simulations or data.")
    return d, blk, "stand-in", False


def _plot_distance(d, blk, p, minutes):
    t = np.arange(len(d)) * blk / 60.0
    m = t <= minutes
    fig, ax = plt.subplots(1, 2, figsize=(14, 3.0), gridspec_kw={"width_ratios": [4, 1]})
    ax[0].plot(t[m], d[m], color=C["dist"], lw=1)
    if p is not None and p.coupling == "contact":
        ax[0].axhline(p.d_contact_nm, color=C["contact"], ls="--", lw=1.3, label=f"contact threshold {p.d_contact_nm:g} nm")
        shade_contact(ax[0], t[m], d[m] <= p.d_contact_nm)
        ax[0].legend(loc="upper right", fontsize=11)
    ax[0].set_xlabel("time (min)"); ax[0].set_ylabel("distance (nm)"); ax[0].set_title("Promoter–enhancer distance through time")
    ax[1].hist(d, bins=40, color=C["dist"], alpha=0.6, density=True)
    if p is not None and p.coupling == "contact":
        ax[1].axvline(p.d_contact_nm, color=C["contact"], ls="--")
    ax[1].set_xlabel("distance (nm)"); ax[1].set_title("How often", fontsize=11)
    fig.tight_layout(); st.pyplot(fig, width="stretch"); plt.close(fig)


def _coupling_inputs() -> dict:
    c = st.columns([3, 3, 2])
    kind = c[0].radio("How does distance act?", list(KIND), format_func=KIND.get, key="tx_kind", horizontal=False)
    rate = c[1].radio("Which rate depends on distance?", list(RATE), format_func=RATE.get, key="tx_rate")
    fold = c[2].number_input("Fold change at full contact", 1.0, 1000.0, 12.0, 1.0, key="tx_fold", help=T.HELP["fold"])
    out = {"coupling": kind, "coupled_rate": rate, "fold": float(fold)}
    c = st.columns(4)
    if kind == "contact":
        out["d_contact_nm"] = float(c[0].number_input("Contact distance (nm)", 10.0, 2000.0, 150.0, 10.0, key="tx_dc", help=T.HELP["d_contact_nm"]))
    elif kind == "hill":
        out["d_half_nm"] = float(c[0].number_input("Distance of half effect (nm)", 10.0, 2000.0, 200.0, 10.0, key="tx_dh"))
        out["hill_n"] = float(c[1].number_input("Steepness (Hill n)", 0.5, 20.0, 4.0, 0.5, key="tx_hn"))
    elif kind == "exp":
        out["d_decay_nm"] = float(c[0].number_input("Decay length (nm)", 10.0, 2000.0, 150.0, 10.0, key="tx_dd"))
    else:
        out["f_const"] = float(c[0].slider("Constant coupling", 0.0, 1.0, 0.1, 0.05, key="tx_fc",
                                           help="Control: the rate is raised by the same factor whatever the distance."))
    return out


def _kinetics_inputs() -> dict:
    c = st.columns(3)
    off = c[0].number_input("Mean OFF time when the enhancer is far (min)", 0.1, 10000.0, 20.0, 5.0, key="tx_off")
    on = c[1].number_input("Mean ON time (min)", 0.1, 1000.0, 0.75, 0.25, key="tx_on")
    ini = c[2].number_input("While ON, a Pol II starts every (s)", 1.0, 3600.0, 7.0, 1.0, key="tx_ini")
    c = st.columns(5)
    out = {"k_on": 1 / (off * 60.0), "k_off": 1 / (on * 60.0), "k_init": 1 / ini}
    out["elongation_kb_min"] = float(c[0].number_input("Pol II speed (kb/min)", 0.2, 20.0, 2.5, 0.5, key="tx_v", help=T.HELP["elongation_kb_min"]))
    out["cassette_kb"] = float(c[1].number_input("MS2 cassette (kb)", 0.0, 20.0, 1.3, 0.1, key="tx_cas", help=T.HELP["cassette_kb"]))
    out["gene_kb"] = float(c[2].number_input("Gene after cassette (kb)", 0.0, 200.0, 2.0, 0.5, key="tx_gene", help=T.HELP["gene_kb"]))
    out["n_loops"] = int(c[3].number_input("MS2 loops", 1, 200, 24, 1, key="tx_loops", help=T.HELP["n_loops"]))
    out["dwell_s"] = float(c[4].number_input("Wait at gene end (s)", 1.0, 600.0, 20.0, 5.0, key="tx_dw", help=T.HELP["dwell_s"]))
    return out


def current_params() -> dict:
    """Parameters as currently set in the tab (defaults if never opened), for the hand-off to a movie run."""
    g = st.session_state.get
    p = dict(T.DEFAULTS)
    p.update(coupling=g("tx_kind", "contact"), coupled_rate=g("tx_rate", "k_on"), fold=float(g("tx_fold", 12.0)),
             d_contact_nm=float(g("tx_dc", 150.0)), d_half_nm=float(g("tx_dh", 200.0)), hill_n=float(g("tx_hn", 4.0)),
             d_decay_nm=float(g("tx_dd", 150.0)), f_const=float(g("tx_fc", 0.1)),
             k_on=1 / (g("tx_off", 50.0) * 60.0), k_off=1 / (g("tx_on", 16.0) * 60.0), k_init=1 / g("tx_ini", 600.0),
             elongation_kb_min=float(g("tx_v", 2.5)), cassette_kb=float(g("tx_cas", 1.3)), gene_kb=float(g("tx_gene", 5.0)),
             n_loops=int(g("tx_loops", 24)), dwell_s=float(g("tx_dw", 116.0)))
    return p


def overrides_for_movie() -> list[str]:
    ov = ["loci.transcription.enabled=true", f"loci.transcription.photons_per_loop_s={st.session_state.get('tx_phot', 100)}"]
    ov += [f"loci.transcription.params.{k}={v:.6g}" if isinstance(v, float) else f"loci.transcription.params.{k}={v}"
           for k, v in current_params().items()]
    return ov


# ----------------------------------------------------------------------------- the tab
@st.fragment
def render(run_path: Path | None = None) -> None:
    style.intro(st, "How the distance between a promoter and its enhancer, taken from a saved polymer simulation, makes the promoter switch ON, "
                    "fire polymerases, and produce a bright MS2 spot. Every step has its own plot.",
                "Work through steps 1 to 9 from top to bottom. Change any number and everything below updates. Step 8 sweeps one setting; "
                "step 9 gives the command (or sidebar checkbox) to run it in a movie. Defaults mimic the TFF1 gene; rates are per second of real time.")
    st.markdown("### From chromosome distance to the MS2 spot")
    st.markdown(overview_svg(), unsafe_allow_html=True)
    st.caption("The polymer simulation tells us how far the enhancer is from the promoter at every moment. That distance sets how "
               "readily the promoter switches ON. While ON it fires Pol II, each of which carries MS2 loops that light up a spot at the "
               "promoter. Nine steps below follow this chain; every number and curve updates as you change a setting. The polymer "
               "does not feel transcription (one-way coupling).")

    # 1 sampling
    _step(1, "How the movie samples time", "All rates below are per SECOND of real time. Changing the frame interval never needs "
                                             "a rate to be re-entered: the model runs in continuous time and is only read at the frame times.")
    c = st.columns(3)
    dt = float(c[0].number_input("Frame interval (s)", 0.5, 3600.0, 20.0, 5.0, key="tx_dt"))
    expo = float(c[1].number_input("Exposure (s)", 0.001, 60.0, 0.1, 0.05, key="tx_expo", format="%.3f"))
    minutes = float(c[2].number_input("Movie length (min)", 5.0, 600.0, 90.0, 10.0, key="tx_min"))
    n_frames = int(minutes * 60 / dt)
    slot = st.container()

    # 2 distance
    _step(2, "Distance between promoter and enhancer", "Taken from a saved polymer simulation at the genomic positions below "
                                                        "(100 kb apart by default), converted to nanometres.")
    src = _distance_source(dt, minutes)
    if src is None:
        return
    d, blk, label, too_short = src
    if too_short:
        st.warning("The saved trajectory is shorter than the movie, so it is repeated back and forth (no jump at the seams). "
                   "Correlations longer than the trajectory are then artificial.")

    _plot_distance(d, blk, None, minutes)

    # 3 coupling
    _step(3, "How the distance sets a rate", "Choose the shape and which rate responds. The curve shows how strongly the enhancer "
                                              "acts (0 to 1) at each distance, drawn over how often each distance occurs.")
    cp = _coupling_inputs()
    p0 = T.TranscriptionParams.from_dict(cp)
    xs = np.linspace(0, max(600.0, np.percentile(d, 99.5)), 400)
    fig, ax = plt.subplots(1, 2, figsize=(14, 2.9), gridspec_kw={"width_ratios": [2, 3]})
    ax[0].plot(xs, T.contact_fraction(p0, xs), color=C["contact"], lw=2.2)
    a2 = ax[0].twinx()
    a2.hist(d, bins=50, color=C["dist"], alpha=0.25, density=True); a2.set_yticks([]); a2.grid(False)
    ax[0].set_xlabel("distance (nm)"); ax[0].set_ylabel("coupling: 0 far, 1 full contact"); ax[0].set_ylim(-0.03, 1.05)
    ax[0].set_title("How strongly the enhancer acts (blue = how often)", fontsize=11)
    tt = np.arange(len(d)) * blk / 60.0
    ax[1].plot(tt, d, color=C["dist"], lw=0.8)
    inc = (d <= p0.d_contact_nm) if p0.coupling == "contact" else (T.contact_fraction(p0, d) > 0.5)
    thr = p0.d_contact_nm if p0.coupling == "contact" else (float(d[inc].max()) if inc.any() else float(np.percentile(d, 5)))
    if p0.coupling != "none":
        ax[1].axhline(thr, color=C["contact"], ls="--", lw=1.2)
        shade_contact(ax[1], tt, inc)
        ax[1].set_ylim(0, max(4 * thr, float(np.percentile(d, 25))))      # zoom on the low end: contact is where the curve dips
    ax[1].set_xlabel("time (min)"); ax[1].set_ylabel("distance (nm), low end")
    ax[1].set_title(f"In contact {100 * inc.mean():.1f}% of the time (teal; dashed = threshold); zoomed on short distances"
                    if p0.coupling != "none" else "No distance dependence: the rate is raised by a constant", fontsize=11)
    fig.tight_layout(); st.pyplot(fig, width="stretch"); plt.close(fig)
    # 4 promoter
    _step(4, "Promoter switching (OFF ⇄ ON)", "Mean times, not rates, so they are easy to relate to published bursts. The enhancer-far "
                                               "values are the baseline; contact multiplies the chosen rate.")
    kin = _kinetics_inputs()
    try:
        p = T.TranscriptionParams.from_dict({**cp, **kin})
    except ValueError as e:
        st.error(str(e))
        return
    with slot:
        for lvl, msg in T.check_sampling(p, dt, expo):
            {"good": st.success, "warn": st.warning, "bad": st.error}[lvl](msg)
        st.caption(f"{n_frames} frames of {dt:g} s = {minutes:g} min. Exposure {expo:g} s.")

    far, near = T.summary(p, 1e9), T.summary(p, 1.0)
    kon, koff, kinit = T.rates_from_distance(p, xs)
    fig, ax = plt.subplots(figsize=(14, 2.6))
    for arr, nm_, col in ((1 / kon / 60, "mean OFF time (min)", C["on"]), (1 / koff / 60, "mean ON time (min)", C["contact"]),
                          (1 / kinit / 60, "Pol II interval (min)", "#a78bfa")):
        ax.plot(xs, arr, color=col, label=nm_, lw=2)
    ax.set_yscale("log"); ax.set_xlabel("distance between promoter and enhancer (nm)"); ax.legend(fontsize=11, ncol=3)
    ax.set_title("The three timescales of the promoter as the enhancer gets closer (a flat line = does not depend on distance)", fontsize=11)
    fig.tight_layout(); st.pyplot(fig, width="stretch"); plt.close(fig)
    style.stat_grid(st, [("ON fraction, enhancer far", f"{100 * far['on_fraction']:.0f}%", "exact"),
                         ("ON fraction, full contact", f"{100 * near['on_fraction']:.0f}%", "exact"),
                         ("bursts per hour, far → contact", f"{far['bursts_per_hour']:.1f} → {near['bursts_per_hour']:.1f}", "exact"),
                         ("Pol II per burst", f"{far['pol_per_burst']:.0f}", "exact: initiation rate / ON→OFF rate")])

    # run the model once
    seed = int(st.session_state.get("tx_seed", 1))
    ft = np.arange(0, n_frames) * dt
    cell = T.simulate_cell(p, d, blk, ft, expo, seed)
    t_min = np.arange(len(cell.d_nm)) * blk / 60.0
    horizon = minutes
    c = st.columns([1, 5])
    if c[0].button("New random cell", key="tx_new"):
        st.session_state["tx_seed"] = seed + 1
        st.rerun(scope="fragment")
    c[1].caption(f"Random seed {seed}: one simulated cell. Every cell gets its own random history; the distance trace above is shared.")

    # promoter state
    fig, ax = plt.subplots(2, 1, figsize=(14, 3.9), sharex=True, gridspec_kw={"height_ratios": [1, 1]})
    contact = (T.contact_fraction(p, cell.d_nm) > 0.5) if p.coupling != "none" else np.zeros(len(cell.d_nm), bool)
    ax[0].plot(t_min, cell.d_nm, color=C["dist"], lw=0.9); ax[0].set_ylabel("distance (nm)")
    shade_contact(ax[0], t_min, contact)
    st_t = np.concatenate([cell.switch_t, [horizon * 60]]) / 60
    ax[1].step(st_t, np.concatenate([cell.switch_on, [cell.switch_on[-1]]]), where="post", color=C["on"], lw=1.6, label="promoter (one cell)")
    p_start, pmean = T.expected_on(p, cell.d_nm, blk)
    ax[1].plot(np.arange(len(p_start)) * blk / 60, p_start, color=C["theory"], ls="--", lw=1.2, label="probability ON (exact, averaged over many cells)")
    shade_contact(ax[1], t_min, contact)
    ax[1].set_ylabel("promoter ON"); ax[1].set_xlabel("time (min)"); ax[1].set_xlim(0, horizon); ax[1].legend(fontsize=11, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2)
    blk_contact = contact[: len(pmean)]
    if blk_contact.any() and (~blk_contact).any():
        ax[0].set_title(f"Step 4: P(ON) {100 * pmean[blk_contact].mean():.1f}% while in contact (teal) vs {100 * pmean[~blk_contact].mean():.1f}% far "
                        f"(exact, this distance history; contact {100 * blk_contact.mean():.1f}% of the time)")
    else:
        ax[0].set_title("Step 4: the promoter switches ON more often while the enhancer is close (teal shading; no contact in this trace)")
    fig.tight_layout(); st.pyplot(fig, width="stretch"); plt.close(fig)

    # 5 polymerases
    _step(5, "Polymerases travel along the gene", "Each initiation sends one Pol II along the MS2 cassette and the gene body at the "
                                                   "speed set above. The picture is a kymograph: time across, position along the gene up.")
    span = float(st.slider("Show a window of (min)", 5, int(max(10, minutes)), int(min(30, minutes)), 5, key="tx_win"))
    start = float(st.slider("Window starts at (min)", 0, int(max(0, minutes - span)), 0, 5, key="tx_start"))
    ev = cell.events
    L_kb = p.cassette_kb + p.gene_kb
    fig, ax = plt.subplots(2, 1, figsize=(14, 4.4), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
    ax[0].axhspan(0, p.cassette_kb, color=C["ms2"], alpha=0.12)
    ax[0].text(start + span * 0.005, p.cassette_kb / 2, "MS2 cassette", va="center", fontsize=8, color=C["ms2"])
    for r in ev.itertuples():
        t0, t1 = r.init_s / 60, (r.init_s + p.t_gene_s) / 60
        if t1 < start or t0 > start + span:
            continue
        ax[0].plot([t0, t1, r.release_s / 60], [0, L_kb, L_kb], color=C["theory"], lw=1.1, alpha=0.8)
    ax[0].set_ylim(0, L_kb * 1.04); ax[0].set_ylabel("kb along gene")
    ax[0].set_title(f"Step 5: each line is one Pol II ({len(ev)} in this cell). Slope = speed; flat top = waiting before release")
    ax[1].step(cell.frames.t_s / 60, cell.frames.n_polII, where="post", color=C["theory"], lw=1.2)
    ax[1].set_ylabel("Pol II"); ax[1].set_xlabel("time (min)"); ax[1].set_xlim(start, start + span)
    fig.align_ylabels(ax)
    fig.tight_layout(); st.pyplot(fig, width="stretch"); plt.close(fig)
    style.stat_grid(st, [("Pol II transit time", f"{p.t_gene_s / 60:.1f} min", "from initiation to the end of the gene"),
                         ("MS2 ramp-up time", f"{p.t_cassette_s:.0f} s", "while the cassette is being transcribed"),
                         ("Pol II on the gene (average)", f"{far['mean_polII']:.1f} far → {near['mean_polII']:.1f} contact", "exact"),
                         ("One Pol II's full signal", f"{p.n_loops} loops", "")])

    # 6 MS2 signal
    _step(6, "The MS2 signal", "Total loops carried by all Pol II on the gene. Grey = the true signal at every moment; dots = what each "
                                "frame records (averaged over the exposure); dashed = the exact expectation over many cells.")
    fine = np.arange(0, minutes * 60, max(1.0, dt / 8))
    Te = (ev.release_s - ev.init_s).to_numpy()
    order = np.argsort(ev.init_s.to_numpy())
    _, inst = T._sample_polymerases(fine, 1e-3, ev.init_s.to_numpy()[order], Te[order], p.t_cassette_s, float(p.n_loops))
    exp_t = np.arange(0, minutes * 60, max(dt, blk))
    expected = T.expected_ms2(p, cell.d_nm, blk, exp_t)
    warm = p.t_gene_s + p.dwell_s
    fig, ax = plt.subplots(figsize=(14, 3.4))
    ax.plot(fine / 60, inst, color=C["grey"], lw=0.9, label="true signal")
    ax.plot(cell.frames.t_s / 60, cell.frames.ms2_loops, "o", color=C["ms2"], ms=3, label=f"recorded every {dt:g} s")
    ax.plot(exp_t / 60, expected, color=C["theory"], ls="--", lw=1.2, label="exact expectation (over cells)")
    ax.axvspan(0, warm / 60, color="#e2e8f0", alpha=0.6, lw=0)
    shade_contact(ax, t_min, contact)
    ax.set_xlim(0, minutes); ax.set_xlabel("time (min)"); ax.set_ylabel("MS2 signal (loops)"); ax.legend(fontsize=11, loc="upper right")
    ax.set_title("Step 6: bursts of initiation appear as a ramp-up, a plateau and a decay (the grey band is the model's warm-up)")
    fig.tight_layout(); st.pyplot(fig, width="stretch"); plt.close(fig)
    fr = cell.frames
    corr = np.corrcoef(fr.ms2_loops, fr.in_contact.astype(float))[0, 1] if fr.in_contact.any() and not fr.in_contact.all() and fr.ms2_loops.std() > 0 else float("nan")
    style.stat_grid(st, [("MS2 mean in this cell", f"{fr.ms2_loops.mean():.1f} loops", ""),
                         ("MS2 expected mean (any cell)", f"{np.mean(expected[exp_t > warm]):.1f} loops" if (exp_t > warm).any() else "—",
                          "exact"),
                         ("Correlation with being in contact", "n/a" if np.isnan(corr) else f"{corr:.2f}",
                          "weakened by the delay: contact → ON → Pol II → signal"),
                         ("Frames brighter than half a Pol II", f"{100 * (fr.ms2_loops > p.n_loops * 0.5).mean():.0f}%", "share of frames with more than n_loops / 2 loops of signal")])

    # 7 microscope
    _step(7, "What the microscope sees", "MS2-GFP carries a nuclear localisation signal, so it also glows through the whole nucleus (same colour as "
                                         "the nuclear label). The spot is the extra light at the promoter; here it sits on that haze.")
    c = st.columns(4)
    ph = float(c[0].number_input("Photons per MS2 loop per s", 1.0, 100000.0, 100.0, 10.0, key="tx_phot", help="PLACEHOLDER brightness."))
    bg = float(c[1].number_input("Nuclear haze (photons/pixel/s)", 0.0, 100000.0, 3000.0, 100.0, key="tx_bg"))
    pxn = float(c[2].number_input("Pixel size (nm)", 50.0, 1000.0, 130.0, 10.0, key="tx_px"))
    sigma = float(c[3].number_input("PSF sigma (pixels)", 0.5, 5.0, 1.4, 0.1, key="tx_sig"))
    st.caption(f"Relative brightness R = one Pol II with all {p.n_loops} loops / one nuclear pixel = {p.n_loops * ph / max(bg, 1e-9):.2f} "
               "(placeholder 0.8). The 🔎 Spot finding page (part 3d) shows how well transcription can be measured as a function of R; "
               "measure R on a real movie as the background-subtracted spot photons per Pol II over the nuclear pixel value (offset removed).")
    nt = 8
    on_frames = np.nonzero(fr.ms2_loops.to_numpy() > 0.5)[0]
    f0 = int(min(max(0, (on_frames[0] if len(on_frames) else 0) - 1), max(0, n_frames - nt)))
    first = int(st.slider("First frame shown", 0, max(0, n_frames - nt), f0, key="tx_f0",
                          help="Starts just before this cell's first burst."))
    idx = np.arange(first, min(first + nt, n_frames))
    rng = np.random.default_rng(seed)
    fig, axes = plt.subplots(1, nt, figsize=(14, 2.3))
    S = 25
    tiles = []
    for i in idx:
        spot = pd.DataFrame({"y_px": [S / 2 - 0.5], "x_px": [S / 2 - 0.5], "photons": [ph * expo * fr.ms2_loops.iloc[i]]})
        exp_img = _spot_photons((S, S), spot, sigma) + bg * expo
        tiles.append(rng.poisson(exp_img).astype(float))
    vmax = max(np.max(t) for t in tiles) if tiles else 1
    for ax, i, im in zip(np.atleast_1d(axes), idx, tiles):
        ax.imshow(im, cmap="gray", vmin=bg * expo * 0.6, vmax=max(vmax, bg * expo * 1.5)); ax.axis("off")
        ax.set_title(f"{fr.t_s.iloc[i] / 60:.1f} min\n{fr.ms2_loops.iloc[i]:.0f} loops", fontsize=9)
    for ax in np.atleast_1d(axes)[len(idx):]:
        ax.axis("off")
    fig.tight_layout(); st.pyplot(fig, width="stretch"); plt.close(fig)
    st.caption("Same grey scale in every tile. Few loops = invisible against the haze; this is where transcription events get missed.")

    # 8 sweep
    _sweep(p, d, blk, dt, expo, minutes)

    # 9 hand-off
    _step(9, "Use these settings in a movie", "In the sidebar, open New single-cell loci run and tick the transcription option. "
                                              "The movie then uses its own frame interval with the same rates.")
    st.code("scripts\\run.cmd python scripts\\run_tracking_demo.py --config configs\\loci_library.yaml " +
            " ".join(f'--set "{o}"' for o in overrides_for_movie()) + f" --set acquisition.frame_interval_s={dt:g}", language="bash")
    _run_view(run_path)


SWEEP = {"fold": ("Fold change at full contact", 12.0), "d_contact_nm": ("Contact distance (nm)", 150.0),
         "off_min": ("Mean OFF time when far (min)", 60.0), "on_min": ("Mean ON time (min)", 6.0), "init_s": ("Pol II interval (s)", 30.0),
         "elongation_kb_min": ("Pol II speed (kb/min)", 2.5)}


def _sweep(p, d, blk, dt, expo, minutes):
    _step(8, "Parameter sweep", "Vary one setting, keep the rest as above, simulate many cells on the same distance trace, and see how the "
                                 "activity and its link to contact change. Add this over a whole movie with scripts/sweep.py.")
    c = st.columns(5)
    key = c[0].selectbox("Setting to vary", list(SWEEP), format_func=lambda k: SWEEP[k][0], key="tx_sw_key")
    base = {"fold": p.fold, "d_contact_nm": p.d_contact_nm, "off_min": 1 / p.k_on / 60, "on_min": 1 / p.k_off / 60,
            "init_s": 1 / p.k_init, "elongation_kb_min": p.elongation_kb_min}[key]
    lo = c[1].number_input("From", 0.0, 1e6, float(base) / 4, format="%.4g", key=f"tx_sw_lo_{key}")
    hi = c[2].number_input("To", 0.0, 1e6, float(base) * 4, format="%.4g", key=f"tx_sw_hi_{key}")
    n = int(c[3].number_input("Values", 2, 15, 6, key="tx_sw_n"))
    reps = int(c[4].number_input("Cells per value", 2, 100, 20, key="tx_sw_reps"))
    if p.coupling != "contact" and key == "d_contact_nm":
        st.info("The contact distance only matters with the sharp-contact shape.")
    if st.button("▶ Run sweep", type="primary", key="tx_sw_go"):
        if lo <= 0 or hi <= lo:
            st.error("Need 0 < From < To.")
        else:
            vals = np.geomspace(lo, hi, n)
            rows, ft = [], np.arange(0, int(minutes * 60 / dt)) * dt
            prog = st.progress(0.0, text="Simulating…")
            for i, v in enumerate(vals):
                q = {"fold": dict(fold=max(1.0, v)), "d_contact_nm": dict(d_contact_nm=v), "off_min": dict(k_on=1 / (v * 60)),
                     "on_min": dict(k_off=1 / (v * 60)), "init_s": dict(k_init=1 / v), "elongation_kb_min": dict(elongation_kb_min=v)}[key]
                pp = T.TranscriptionParams(**{**p.__dict__, **q})
                for r in range(reps):
                    c_ = T.simulate_cell(pp, d, blk, ft, expo, 1000 * i + r)
                    f_ = c_.frames
                    sw = c_.switch_t[1:][c_.switch_on[1:] == 1]
                    sw = sw[(sw >= 0) & (sw <= minutes * 60)]
                    warm = pp.t_gene_s + pp.dwell_s
                    m = f_.t_s > warm
                    ct = f_.in_contact.astype(float)
                    rows.append({"value": v, "ON fraction": f_.promoter_on.mean(), "bursts per hour": len(sw) / (minutes / 60),
                                 "MS2 mean (loops)": f_.ms2_loops[m].mean() if m.any() else np.nan,
                                 "corr(MS2, contact)": np.corrcoef(f_.ms2_loops, ct)[0, 1] if ct.std() > 0 and f_.ms2_loops.std() > 0 else np.nan})
                prog.progress((i + 1) / n, text=f"{i + 1} of {n} values")
            prog.empty()
            st.session_state["tx_sweep"] = (key, pd.DataFrame(rows))
    res = st.session_state.get("tx_sweep")
    if res and res[0] == key:
        df = res[1]
        cols = [c for c in df.columns if c != "value"]
        fig, ax = plt.subplots(1, 4, figsize=(15, 3.2))
        for a, col in zip(ax, cols):
            g = df.groupby("value")[col].agg(["mean", "std"]).reset_index()
            a.errorbar(g.value, g["mean"], yerr=g["std"].fillna(0), marker="o", capsize=3, color=C["dist"])
            a.set_xscale("log"); a.set_xlabel(SWEEP[key][0]); a.set_title(col, fontsize=11)
        fig.tight_layout(); st.pyplot(fig, width="stretch"); plt.close(fig)
        st.caption("Dots with bars = mean ± 1 sd over cells. The correlation is weak when bursts are much longer than the time in contact "
                   "or when the frame interval is coarse.")
        with st.expander("Table and download"):
            g = df.groupby("value").mean().round(4).reset_index()
            st.dataframe(g, hide_index=True, width="stretch")
            st.download_button("Download all cells (CSV)", df.to_csv(index=False), f"transcription_sweep_{key}.csv", "text/csv")


def _run_view(run_path):
    f = run_path / "stage1_chromatin" / "transcription.csv" if run_path else None
    if f is None or not f.exists():
        return
    st.markdown("### In the selected run")
    tx = pd.read_csv(f)
    tx = tx if "copy" in tx else tx.assign(copy=0)
    n_cp = tx.groupby("cell_id")["copy"].nunique()
    cells = sorted(tx.cell_id.unique())
    two = [c for c in cells if n_cp[c] > 1]
    busy = tx.groupby("cell_id").ms2_loops.apply(lambda v: float((v > 0.5).mean()))
    pick = max(cells, key=lambda c: (n_cp[c] > 1, busy.get(c, 0.0)))          # several copies and the most bursting first
    cid = st.selectbox("Cell", cells, index=cells.index(pick), key="tx_cell",
                       format_func=lambda c: f"cell {c}" + (f" ({n_cp[c]} promoter-enhancer copies)" if n_cp[c] > 1 else ""))
    g = tx[tx.cell_id == cid]
    cols = ["#e5e7eb", "#2dd4bf", "#f59e0b", "#a78bfa"]
    fig, ax = plt.subplots(3, 1, figsize=(14, 5.6), sharex=True)
    for c, q in g.groupby("copy"):
        q = q.sort_values("t_s")
        lab = f"copy {c}" + ("" if q.get("ms2_labelled", pd.Series([True])).iloc[0] else " (no MS2 cassette)")
        ax[0].plot(q.t_s / 60, q.d_nm, color=cols[int(c) % 4], lw=1.2, label=lab)
        ax[1].step(q.t_s / 60, q.promoter_on + 0.06 * c, where="mid", color=cols[int(c) % 4], lw=1.4)
        ax[2].plot(q.t_s / 60, q.ms2_loops, color=cols[int(c) % 4], marker="o", ms=2.5, lw=1.0)
    if g["copy"].nunique() == 1:
        shade_contact(ax[0], (g.t_s / 60).to_numpy(), g.in_contact.to_numpy().astype(bool))
    ax[0].set_ylabel("distance (nm)"); ax[1].set_ylabel("promoter ON"); ax[1].set_yticks([0, 1]); ax[1].set_yticklabels(["OFF", "ON"])
    ax[2].set_ylabel("MS2 (loops)"); ax[2].set_xlabel("time (min)")
    ax[0].legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=4, fontsize=10)
    fig.tight_layout(); st.pyplot(fig, width="stretch"); plt.close(fig)
    if g["copy"].nunique() > 1:
        st.caption("Each copy has its own chromatin trajectory, so its own promoter-enhancer distance, and its own promoter: their bursts are "
                   "independent (and follow their own rates if loci.copies.transcription gives per-copy overrides). The ON lines are offset "
                   "slightly so they do not hide each other.")
