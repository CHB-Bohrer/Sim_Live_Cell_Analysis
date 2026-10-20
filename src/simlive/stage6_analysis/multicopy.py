"""Stage 6: several promoter-enhancer copies per cell, and stray spots that are not loci.

The one-copy spot finder (dog_viterbi.py) follows ONE path per colour through the candidates. A cell with two copies of the pair (homologous
alleles, `loci.copies`) shows two promoter spots and two enhancer spots, and any nucleus can hold short-lived stray spots of either colour
(probe aggregates; image error source `stray_spots`). `track_copies` turns the candidates of all colours of one cell into copies:

  1  paths: per colour, up to max_copies + extra_tracks Viterbi paths, one after the other, each through the candidates the earlier paths
     left unused (same link costs and rewards as for one path; D re-estimated for every path);
  2  which paths are loci: beyond the first path of each colour a path counts only if it looks like that first path: coverage (share of
     the cell's frames with a spot) >= min_track_coverage, median SNR >= min_snr_ratio x the first path's, and diffusion D <= max_D_ratio x
     the first path's. A labelled locus is there in nearly every frame, as bright as the other copy and slow; once the loci are taken, the
     leftover candidates are noise peaks and stray spots (each ~1 min), and a path through them is faint and hops (measured on the polymer
     library movies: loci median SNR 7.6-10.8, D 2.7-5.6 px^2 per frame; leftover paths SNR 3.3-6.7, D 5-117, yet coverage 0.65-0.98 because
     there is always some candidate to hop to, so coverage alone does not separate them);
  2b untangling: the paths are found one after another, so the first one is free to jump from one copy to the other where that buys a
     brighter stretch, and the second then takes the complement: the two copies swap identity (measured: whole cells swapped for the last
     30% of the movie, copies 1.3 um apart). Wherever crossing two paths' identities makes their frame-to-frame links cheaper (sum of squared
     steps / gap), their tails are exchanged (`untangle`). Both paths keep their spots, only who-is-who changes;
  3  pairing: the promoter and enhancer of one copy stay within about a micrometre of each other, two copies (or a locus and a stray) do not.
     Cost = median distance between a promoter path and an enhancer path over the frames where both have a spot; the Hungarian assignment
     keeps pairs with cost <= pair_max_um, at most max_copies. Paths without a partner are dropped; if no pair qualifies at all, the first
     path of each colour is kept (the one-copy result). More than two colours: every further colour joins the copy whose promoter is nearest;
  4  pair prior (optional): each copy's paths are found again, copy by copy, the promoter anchored on the copy's enhancer path and then the
     enhancer on the new promoter path, the candidates of the copies already finished excluded:
       pair_prior: gate   every candidate farther than pair_gate_um from the anchor is penalised by ((d - gate) / softness)^2 / 2 nats. It
                          rejects a stray spot that a path picked up after losing its locus. Flat inside, so distances shorter than the gate
                          are not pulled towards smaller values (a pair that really separates by more than the gate is);
       pair_prior: gauss  -d^2 / (2 s^2) with s from the cell's own pairs (median pair distance / 1.18 in 2D, / 1.54 in 3D) or pair_sigma_um.
                          It also decides which copy a candidate belongs to when two copies come close (a swap of two promoters costs the
                          extra distance to the wrong enhancer), but it pulls every measured distance towards smaller values (biased; report
                          with and without it).
  5  the copies are untangled once more, both colours together (a swap of a whole copy, promoter and enhancer at once, is undone by the
     summed steps of both).
  5 also runs `untangle_by_place` on the paired copies (2026-10-20): step costs cannot see a swap made through a wrong pick or a gap (in one cell both
     copies' paths ended up on each other's locus from frame ~72 on, 2.6 um apart). There every spot goes to the track whose RECENT place
     (median of its last few spots, relative to the nucleus centroid) it is nearest, where that is clearly better.
Copies are numbered by the rank of their promoter path (copy 0 = the best path). The number says nothing about WHICH true copy (allele) it
is; scoring matches them (localize.match_copies).
"""
from __future__ import annotations

import numpy as np

from simlive.stage6_analysis import dog_viterbi as DV

PLACE_MIN_SEP_UM = 0.8     # untangle_by_place acts only while the tracks' recent places are farther apart than this
SETTINGS = dict(max_copies=1, extra_tracks=2, min_track_coverage=0.5, min_snr_ratio=0.6, max_D_ratio=3.0, pair_max_um=1.5, pair_prior="gate",
                pair_gate_um=None, gate_softness_px=2.0, pair_sigma_um=None)
HELP = {
    "max_copies": "Most promoter-enhancer copies a cell can have (2 = diploid, 4 = after replication). 1 = the one-copy finder.",
    "extra_tracks": "Paths tried per colour beyond max_copies (to find the loci among stray spots).",
    "min_track_coverage": "A second or later path counts as a locus only if it has a spot in at least this share of the frames.",
    "min_snr_ratio": "... and only if its median SNR is at least this fraction of the colour's best path (a second copy is about as bright).",
    "max_D_ratio": "... and only if it moves at most this many times as fast as the best path (a path through leftover noise hops).",
    "pair_max_um": "A promoter path and an enhancer path belong to one copy only if their median distance is at most this (um).",
    "pair_prior": "How the other colour of a copy is used to choose spots: gate (flat window of pair_gate_um), gauss (smooth, biased), none.",
    "pair_gate_um": "If set: candidates farther than this from the other colour's spot of the same copy are penalised (um). Empty = off.",
    "pair_sigma_um": "Scale of the gauss pair prior (um); empty = from the cell's own promoter-enhancer pairs.",
    "gate_softness_px": "How quickly the penalty grows beyond the gate (px): ((d - gate) / softness)^2 / 2 nats.",
}


def _prior_kind(pair_prior, pair_gate_um) -> str | None:
    if pair_prior == "gauss":
        return "gauss"
    return "gate" if pair_prior in ("gate", None) and pair_gate_um else None


def multicopy_on(kw: dict) -> bool:
    """Do these locus_analysis settings ask for the multi-copy linker (more than one copy, or a pair prior)?"""
    return int(kw.get("max_copies") or 1) > 1 or _prior_kind(kw.get("pair_prior", "gate"), kw.get("pair_gate_um")) is not None


def path_quality(choice: list, valid: np.ndarray | None, max_gap: int, snr: list | None = None) -> dict:
    """Coverage (share of valid frames with a spot), number of connected stretches (separated by more than max_gap missed frames), and the
    median SNR of its spots (if `snr` is given)."""
    n = len(choice)
    valid = np.ones(n, bool) if valid is None else np.asarray(valid, bool)
    hit = np.array([c is not None for c in choice]) & valid
    f = np.nonzero(hit)[0]
    seg = int(1 + (np.diff(f) > max_gap + 1).sum()) if len(f) else 0
    out = {"coverage": float(hit.sum() / max(valid.sum(), 1)), "n_frames": int(hit.sum()), "segments": seg}
    if snr is not None:
        v = [float(np.asarray(snr[i])[c]) for i, c in enumerate(choice) if c is not None]
        out["snr_median"] = float(np.median(v)) if v else 0.0
    return out


def rel_positions(peaks: list, choice: list, centre: np.ndarray, z_scale: float = 1.0) -> np.ndarray:
    """(n, ndim) positions of a path relative to the nucleus centroid (z scaled to pixels), NaN where the path has no spot."""
    nd = next((np.asarray(p).shape[-1] for p in peaks if len(p)), 2)
    out = np.full((len(peaks), nd), np.nan)
    for i, j in enumerate(choice):
        if j is None:
            continue
        q = np.asarray(peaks[i][j], float).copy()
        q[-2:] -= centre[i]
        if nd == 3:
            q[0] *= z_scale
        out[i] = q
    return out


def fill_gaps(a: np.ndarray) -> np.ndarray:
    """Linear interpolation over the frames a path has no spot (held at the ends); all NaN stays NaN."""
    a = a.copy()
    ok = np.isfinite(a[:, 0])
    if ok.sum() == 0:
        return a
    t = np.arange(len(a))
    for k in range(a.shape[1]):
        a[:, k] = np.interp(t, t[ok], a[ok, k])
    return a


def pair_cost(a: np.ndarray, b: np.ndarray, min_common: int) -> float:
    """Median distance (px) between two paths over the frames where both have a spot; inf if fewer than min_common such frames."""
    both = np.isfinite(a[:, 0]) & np.isfinite(b[:, 0])
    if both.sum() < min_common:
        return np.inf
    return float(np.median(np.linalg.norm(a[both] - b[both], axis=1)))


def gate_offsets(peaks: list, centre: np.ndarray, partner: np.ndarray, gate_px: float, softness_px: float, z_scale: float = 1.0,
                 sigma_px: float | None = None) -> list:
    """Reward offset (nats, <= 0) of every candidate, d = distance to the partner path: -((d - gate) / softness)^2 / 2 beyond the gate, or
    with sigma_px (the gauss prior) -d^2 / (2 sigma^2) everywhere."""
    out = []
    for i, p in enumerate(peaks):
        p = np.asarray(p, float)
        if not len(p) or not np.isfinite(partner[i, 0]):
            out.append(np.zeros(len(p)))
            continue
        q = p.copy()
        q[:, -2:] -= centre[i]
        if q.shape[1] == 3:
            q[:, 0] *= z_scale
        d = np.linalg.norm(q - partner[i], axis=1)
        if sigma_px:
            out.append(-0.5 * (d / sigma_px) ** 2)
        else:
            out.append(-0.5 * (np.maximum(d - gate_px, 0.0) / max(softness_px, 1e-6)) ** 2)
    return out


def extract_paths(peaks: list, snr: list, centre: np.ndarray, n_paths: int, valid: np.ndarray | None = None, used: list | None = None,
                  offsets: list | None = None, **vkw) -> list[dict]:
    """Up to n_paths Viterbi paths, each through the candidates the earlier ones (and `used`) did not take. Each: {choice (index into the
    frame's candidates or None), D_px2, coverage, n_frames, segments}."""
    used = [np.zeros(len(p), bool) for p in peaks] if used is None else [np.asarray(u, bool).copy() for u in used]
    out = []
    for _ in range(int(n_paths)):
        idx = [np.nonzero(~u)[0] for u in used]
        pk = [np.asarray(p)[i] if len(p) else np.asarray(p) for p, i in zip(peaks, idx)]
        sn = [np.asarray(s)[i] if len(s) else np.asarray(s) for s, i in zip(snr, idx)]
        if not any(len(x) for x in pk):
            break
        off = None if offsets is None else [np.asarray(o)[i] if len(o) else np.zeros(0) for o, i in zip(offsets, idx)]
        ch, D = DV.viterbi_track(pk, sn, centre, reward_offset=off, **vkw)
        choice = [None if j is None else int(idx[i][j]) for i, j in enumerate(ch)]
        if all(c is None for c in choice):
            break
        for i, c in enumerate(choice):
            if c is not None:
                used[i][c] = True
        out.append({"choice": choice, "D_px2": D, **path_quality(choice, valid, int(vkw.get("viterbi_max_gap", 5)), snr)})
    return out


def untangle(tracks: list, peaks: list, centre: np.ndarray, z_scale: float = 1.0) -> list:
    """Exchange the tails of tracks wherever that makes the frame-to-frame links cheaper. `tracks[a]` = {colour: choice list} (one colour for a
    colour's own paths, all colours of a copy for whole copies); `peaks[colour][i]` = candidates of frame i. At every frame, every track that
    has a spot there (in any of its colours) is linked to its last spot (per colour); the cost of continuing track a with track b's spots
    is the sum over colours of |step|^2 / gap. If the cheapest assignment (Hungarian) is not the identity, the tracks' tails from that frame
    on are permuted. Returns new {colour: choice} dicts (the same spots, only their grouping changes)."""
    from scipy.optimize import linear_sum_assignment

    K = len(tracks)
    if K < 2:
        return tracks
    tr = [{k: list(v) for k, v in t.items()} for t in tracks]
    cols = sorted({k for t in tr for k in t})
    n = len(centre)
    pos = [{k: rel_positions(peaks[k], t[k], centre, z_scale) for k in t} for t in tr]
    last = [{k: None for k in t} for t in tr]
    for t in range(n):
        here = [a for a in range(K) if any(k in pos[a] and np.isfinite(pos[a][k][t, 0]) for k in cols)]
        ready = [a for a in here if any(last[a].get(k) is not None for k in cols)]
        if len(ready) >= 2:
            C = np.zeros((len(ready), len(ready)))
            for i, a in enumerate(ready):
                for j, b in enumerate(ready):
                    cost, used = 0.0, 0
                    for k in cols:
                        lk = last[a].get(k)
                        if lk is None or k not in pos[b] or not np.isfinite(pos[b][k][t, 0]):
                            continue
                        cost += float(np.sum((pos[b][k][t] - lk[1]) ** 2)) / max(t - lk[0], 1)
                        used += 1
                    C[i, j] = cost if used else 1e6
            rr, cc = linear_sum_assignment(C)
            if C[rr, cc].sum() < np.trace(C) - 1e-9:
                perm = {ready[r]: ready[c] for r, c in zip(rr, cc)}
                new_tr = [{k: list(v) for k, v in x.items()} for x in tr]
                new_pos = [{k: v.copy() for k, v in x.items()} for x in pos]
                for a, b in perm.items():
                    for k in cols:
                        if k in tr[a] and k in tr[b]:
                            new_tr[a][k][t:] = tr[b][k][t:]
                            new_pos[a][k][t:] = pos[b][k][t:]
                tr, pos = new_tr, new_pos
        for a in range(K):
            for k in pos[a]:
                if np.isfinite(pos[a][k][t, 0]):
                    last[a][k] = (t, pos[a][k][t].copy())
    return tr


def untangle_by_place(tracks: list, peaks: list, centre: np.ndarray, z_scale: float = 1.0, memory: int = 6, min_sep_px: float = 6.0,
                      margin: float = 0.25) -> list:
    """Give every track the spots that lie where IT has recently been. `untangle` links each spot to the track's last spot (squared step /
    gap), so a wrong pick moves the reference and a swap made through a wrong pick or a gap goes unseen (found 2026-10-20: in one cell the
    enhancer paths changed places at frames 70/72 and the promoter paths at 72/76, copies 2.6 um apart, with wrong picks in between; the
    copies stayed swapped for the rest of the movie). Here each track's reference is the median of its last `memory` spots (relative to the
    nucleus centroid, so a moving nucleus does not matter, and over a few frames a turning one hardly does): one wrong pick cannot move it,
    and it waits unchanged through any gap. Frame by frame, per colour, the spots present are assigned to the tracks (Hungarian on the
    distance to the references, tracks without a spot included); the assignment replaces the current one only if it is shorter by more than
    `margin` x the distance between the references, and only while the references are more than `min_sep_px` apart (copies that overlap are
    left to the step costs). An absolute 'territory' version was tried first and dropped: on near-round nuclei the outline's orientation can
    flip by 180 degrees, which exchanges two copies that sit on opposite sides of the centre. Same input and output as `untangle`."""
    from scipy.optimize import linear_sum_assignment

    K = len(tracks)
    if K < 2:
        return tracks
    tr = [{k: list(v) for k, v in t.items()} for t in tracks]
    n = len(centre)
    for k in sorted({q for t in tr for q in t}):
        if not all(k in t for t in tr):
            continue
        P = [rel_positions(peaks[k], t[k], centre, z_scale) for t in tr]
        choice = [list(t[k]) for t in tr]
        hist = [[] for _ in range(K)]                          # recent positions of every track (after any reassignment)
        for i in range(n):
            here = [a for a in range(K) if np.isfinite(P[a][i, 0])]
            if not here:
                continue
            if all(len(hist[a]) >= memory for a in range(K)):      # a place known from a few spots only is no reference (e.g. one wrong pick)
                ref = [np.median(np.array(h[-memory:]), axis=0) for h in hist]
                sep = min(float(np.linalg.norm(ref[x] - ref[y])) for x in range(K) for y in range(x + 1, K))
                if sep > min_sep_px:
                    D = np.array([[np.linalg.norm(P[a][i] - ref[b]) for b in range(K)] for a in here])
                    r2, c2 = linear_sum_assignment(D)
                    dest = {here[r]: int(c) for r, c in zip(r2, c2)}
                    own = sum(D[j, a] for j, a in enumerate(here))
                    if any(dest[a] != a for a in here) and own - D[r2, c2].sum() > margin * sep:
                        spots = {a: (choice[a][i], P[a][i].copy()) for a in here}
                        for a in here:
                            choice[a][i] = None
                            P[a][i] = np.nan
                        for a in here:
                            choice[dest[a]][i], P[dest[a]][i] = spots[a]
            for a in range(K):
                if np.isfinite(P[a][i, 0]):
                    hist[a].append(P[a][i].copy())
        for a in range(K):
            tr[a][k] = choice[a]
    return tr


def repair(copies: list, peaks: list, centre: np.ndarray, z_scale: float, min_common: int, margin: float = 0.8) -> list:
    """Pair every further colour with the promoters again after `untangle_by_place`. That pass repairs each colour's tracks on its own, so
    a promoter and an enhancer that were paired while one of them was swapped for part of the movie end up on DIFFERENT true copies for the
    whole movie (seen 2026-10-20 without a pair prior: one cell, 19% -> 50% of a copy's spots on the other copy). Hungarian on the median
    distance between the promoter tracks and the colour's tracks; the new pairing replaces the old one only if its total distance is below
    `margin` x the old one."""
    from scipy.optimize import linear_sum_assignment

    K = len(copies)
    if K < 2 or any(0 not in c for c in copies):
        return copies
    out = [dict(c) for c in copies]
    P = [rel_positions(peaks[0], c[0], centre, z_scale) for c in copies]
    for k in sorted({q for c in copies for q in c} - {0}):
        if not all(k in c for c in copies):
            continue
        Q = [rel_positions(peaks[k], c[k], centre, z_scale) for c in copies]
        C = np.array([[pair_cost(P[a], Q[b], min_common) for b in range(K)] for a in range(K)])
        C = np.where(np.isfinite(C), C, 1e6)
        rr, cc = linear_sum_assignment(C)
        if C[rr, cc].sum() < margin * np.trace(C):
            for a, b in zip(rr, cc):
                out[a][k] = copies[b][k]
    return out


def assign_pairs(C: np.ndarray, limit: float) -> list[tuple[int, int]]:
    """Rows paired with columns (cost matrix C, e.g. distances) so that as many pairs as possible have cost <= limit, and among those the
    total cost is smallest. Every pair beyond the limit costs the same large amount: with the raw costs the assignment would break a good
    pair to shorten an invalid one. Returns the valid (row, column) pairs."""
    from scipy.optimize import linear_sum_assignment

    if not C.size:
        return []
    rr, cc = linear_sum_assignment(np.where(C <= limit, C, 1e6))
    return [(int(a), int(b)) for a, b in zip(rr, cc) if C[a, b] <= limit]


def _used_by(paths: list, n_frames: int, n_cand: list) -> list:
    used = [np.zeros(k, bool) for k in n_cand]
    for p in paths:
        for i, c in enumerate(p["choice"]):
            if c is not None:
                used[i][c] = True
    return used


def track_copies(peaks: list, snr: list, centre: np.ndarray, valid: np.ndarray | None = None, *, px_um: float = 0.13, z_scale: float = 1.0,
                 max_copies: int = 1, extra_tracks: int = 2, min_track_coverage: float = 0.5, min_snr_ratio: float = 0.6,
                 max_D_ratio: float = 3.0, pair_max_um: float = 1.5,
                 pair_prior: str | None = "gate", pair_gate_um: float | None = None, gate_softness_px: float = 2.0,
                 pair_sigma_um: float | None = None, **vkw) -> dict:
    """All colours of one cell -> copies. `peaks[k][i]` / `snr[k][i]`: candidates of colour k in frame i (crop px; (z, y, x) for stacks),
    `centre[i]` the nucleus centroid (crop px), `valid[i]` whether frame i has a cell region. Returns {"copies": [{colour: choice list}],
    "paths": {colour: [path dicts]}, "pairs": [(promoter path, enhancer path, median distance um)], "fallback": True if no pair qualified,
    "gated": whether the pair gate re-ran the paths}. Other keywords: the Viterbi settings (dog_viterbi.viterbi_track)."""
    L, n = len(peaks), len(centre)
    max_copies = max(int(max_copies or 1), 1)
    vkw = {**vkw, "z_scale": z_scale}
    paths = {k: extract_paths(peaks[k], snr[k], centre, max_copies + int(extra_tracks), valid, **vkw) for k in range(L)}

    def looks_like_a_locus(k: int, r: int) -> bool:
        if r == 0:
            return True
        p, best = paths[k][r], paths[k][0]
        return (p["coverage"] >= min_track_coverage and p.get("snr_median", 0) >= min_snr_ratio * best.get("snr_median", 0)
                and p["D_px2"] <= max_D_ratio * max(best["D_px2"], 0.5))

    ok = {k: [looks_like_a_locus(k, r) for r in range(len(paths[k]))] for k in range(L)}
    for k in range(L):                     # undo swaps between this colour's locus paths (see 2b above); noise paths are left alone
        idx = [r for r in range(len(paths[k])) if ok[k][r]]
        if len(idx) >= 2:
            fixed = untangle([{k: paths[k][r]["choice"]} for r in idx], peaks, centre, z_scale)
            for r, f in zip(idx, fixed):
                paths[k][r] = {**paths[k][r], "choice": f[k], **path_quality(f[k], valid, int(vkw.get("viterbi_max_gap", 5)), snr[k])}
    pos = {k: [rel_positions(peaks[k], p["choice"], centre, z_scale) for p in paths[k]] for k in range(L)}
    n_valid = int(np.sum(valid)) if valid is not None else n
    min_common = max(5, int(0.1 * n_valid))
    pair_px = float(pair_max_um) / px_um
    copies, pairs, fallback = [], [], False
    if L == 1:
        copies = [{0: p["choice"]} for r, p in enumerate(paths[0]) if ok[0][r]][:max_copies]
    elif L >= 2:
        P, E = [r for r in range(len(paths[0])) if ok[0][r]], [r for r in range(len(paths[1])) if ok[1][r]]
        C = np.full((len(P), len(E)), np.inf)
        for a, i in enumerate(P):
            for b, j in enumerate(E):
                C[a, b] = pair_cost(pos[0][i], pos[1][j], min_common)
        if C.size:
            good = sorted((P[a], E[b], C[a, b]) for a, b in assign_pairs(C, pair_px))
            pairs = [(i, j, float(c) * px_um) for i, j, c in good[:max_copies]]
        if not pairs and paths[0] and paths[1]:                  # nothing pairs up: keep the first path of each colour (one-copy result)
            fallback = True
            pairs = [(0, 0, float(pair_cost(pos[0][0], pos[1][0], 1)) * px_um)]
        copies = [{0: paths[0][i]["choice"], 1: paths[1][j]["choice"]} for i, j, _ in pairs]
        if not copies and (paths[0] or paths[1]):                 # one colour has no path at all
            copies = [{k: paths[k][0]["choice"] for k in range(2) if paths[k]}]
        for k in range(2, L):                                    # further colours join the copy with the nearest promoter
            cand = [r for r in range(len(paths[k])) if ok[k][r]]
            if not cand or not copies:
                continue
            C = np.array([[pair_cost(pos[0][pairs[c][0]], pos[k][r], min_common) if c < len(pairs) else np.inf for r in cand]
                          for c in range(len(copies))])
            for a, b in assign_pairs(C, pair_px):
                copies[a][k] = paths[k][cand[b]]["choice"]
    kind = _prior_kind(pair_prior, pair_gate_um)
    sigma_px = None
    if kind == "gauss":
        nd = next((np.asarray(p).shape[-1] for p in peaks[0] if len(p)), 2)
        med = [c / px_um for _, _, c in pairs if np.isfinite(c)]
        sigma_px = float(pair_sigma_um) / px_um if pair_sigma_um else (float(np.median(med)) / (1.18 if nd == 2 else 1.54) if med else pair_px / 2)
        sigma_px = max(sigma_px, 1.0)
    if kind is not None and copies and L >= 2:
        n_cand = {k: [len(p) for p in peaks[k]] for k in range(L)}
        gate_px = float(pair_gate_um or pair_max_um) / px_um
        for c, cp in enumerate(copies):
            for k in [q for q in (0, 1) if q in cp] + [q for q in cp if q > 1]:     # promoter first (anchored on the enhancer), then the rest
                partner_k = 1 if k == 0 else 0
                if partner_k not in cp:
                    continue
                partner = fill_gaps(rel_positions(peaks[partner_k], cp[partner_k], centre, z_scale))
                used = _used_by([{"choice": copies[o][k]} for o in range(c) if k in copies[o]], n, n_cand[k])   # finished copies only
                off = gate_offsets(peaks[k], centre, partner, gate_px, gate_softness_px, z_scale, sigma_px)
                re = extract_paths(peaks[k], snr[k], centre, 1, valid, used=used, offsets=off, **vkw)
                if re:
                    cp[k] = re[0]["choice"]
    if len(copies) >= 2:                   # a whole copy (promoter and enhancer together) can still be swapped: untangle the copies jointly
        copies = untangle(copies, peaks, centre, z_scale)
        copies = untangle_by_place(copies, peaks, centre, z_scale, min_sep_px=PLACE_MIN_SEP_UM / px_um)
        copies = repair(copies, peaks, centre, z_scale, min_common)
    return {"copies": copies, "paths": paths, "pairs": pairs, "fallback": fallback, "gated": kind is not None, "prior": kind,
            "sigma_px": sigma_px}
