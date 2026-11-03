"""Numerical post-analysis used by the report: grid-convergence (Richardson / GCI),
interpolation of the slot-length design chart, sensitivities."""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


def read_csv(path: Path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k, v in r.items():
            try:
                r[k] = float(v)
            except (TypeError, ValueError):
                pass
    return rows


def richardson(h, f):
    """Celik et al. (2008) / Roache GCI for three grids.

    h: cell sizes fine -> coarse (h1 < h2 < h3); f: matching values.
    Returns dict with order p, extrapolated value, GCI of the finest grid (absolute,
    same units as f, safety factor 1.25) and the convergence type.
    """
    h1, h2, h3 = h
    f1, f2, f3 = f
    r21, r32 = h2 / h1, h3 / h2
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0 or e32 == 0:
        return {"p": float("nan"), "f_ext": f1, "gci": 0.0, "type": "converged"}
    ratio = e32 / e21
    s = 1.0 if ratio > 0 else -1.0
    kind = "monotonic" if ratio > 0 else "oscillatory"
    p = 1.0
    for _ in range(200):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p_new = abs(math.log(abs(ratio)) + q) / math.log(r21)
        if abs(p_new - p) < 1e-10:
            break
        p = p_new
    p_used = min(max(p, 0.5), 4.0)  # guard against unphysical orders
    f_ext = (r21 ** p_used * f1 - f2) / (r21 ** p_used - 1.0)
    gci = 1.25 * abs(f1 - f2) / (r21 ** p_used - 1.0)
    return {"p": p, "p_used": p_used, "f_ext": f_ext, "gci": gci, "type": kind,
            "r21": r21, "r32": r32}


def interp_L_for(target_db, Ls, Cs):
    """Slot length that gives target coupling, by monotone cubic-free linear interpolation in L
    of the measured C(L) curve (C decreases with L)."""
    order = np.argsort(Cs)
    c = np.asarray(Cs)[order]
    if not (c[0] <= target_db <= c[-1]):
        return None   # outside the simulated range: do not extrapolate
    return float(np.interp(target_db, c, np.asarray(Ls)[order]))


def slope(xs, ys, x0):
    """Central finite-difference slope dy/dx at x0 from the tabulated points around it."""
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    i = int(np.argmin(np.abs(xs - x0)))
    lo, hi = max(i - 1, 0), min(i + 1, len(xs) - 1)
    return float((ys[hi] - ys[lo]) / (xs[hi] - xs[lo]))


def band_where(freq, values, threshold, above=True):
    """Contiguous band around the centre where values stay above (or below) threshold."""
    freq, values = np.asarray(freq), np.asarray(values)
    ok = values >= threshold if above else values <= threshold
    if ok.all():
        return float(freq[0]), float(freq[-1]), True
    c = len(freq) // 2
    if not ok[c]:
        return float("nan"), float("nan"), False
    lo = c
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    hi = c
    while hi < len(freq) - 1 and ok[hi + 1]:
        hi += 1
    return float(freq[lo]), float(freq[hi]), False


def fit_power_law(h, f):
    """Least-squares fit f = f0 + k h^p over all points (p searched on a fine grid)."""
    h, f = np.asarray(h, float), np.asarray(f, float)
    best = None
    for p in np.linspace(0.1, 3.0, 581):
        A = np.vstack([np.ones_like(h), h ** p]).T
        coef, res, *_ = np.linalg.lstsq(A, f, rcond=None)
        rss = float(np.sum((A @ coef - f) ** 2))
        if best is None or rss < best[0]:
            best = (rss, p, coef)
    rss, p, (f0, k) = best
    return {"f0": float(f0), "k": float(k), "p": float(p), "rms": math.sqrt(rss / len(h))}


def joint_fit(series, p_grid=None):
    """Fit several mesh series that share one limit C0 and one order p:
    C_i(h) = C0 + k_i h^p. series = [[(h, C), ...], ...]."""
    p_grid = np.linspace(0.2, 3.0, 561) if p_grid is None else p_grid
    best = None
    n_s = len(series)
    for p in p_grid:
        A, y = [], []
        for i, ser in enumerate(series):
            for h, c in ser:
                row = [1.0] + [0.0] * n_s
                row[1 + i] = h ** p
                A.append(row)
                y.append(c)
        A, y = np.array(A), np.array(y)
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        rss = float(np.sum((A @ coef - y) ** 2))
        if best is None or rss < best[0]:
            best = (rss, p, coef, len(y))
    rss, p, coef, n = best
    return {"C0": float(coef[0]), "p": float(p), "k": [float(k) for k in coef[1:]], "rms": math.sqrt(rss / n)}


def joint_fit_with_spread(series_full, min_points=(3, 3)):
    """Joint fit plus a robustness spread: leave-one-out over all points, and dropping
    the coarsest points of the first series down to min_points."""
    base = joint_fit(series_full)
    alts = []
    for i, ser in enumerate(series_full):
        for j in range(len(ser)):
            sub = [s if k != i else s[:j] + s[j + 1:] for k, s in enumerate(series_full)]
            if all(len(s) >= 2 for s in sub):
                alts.append(joint_fit(sub)["C0"])
    first = sorted(series_full[0], key=lambda q: q[0])
    for m in range(min_points[0], len(first) + 1):
        alts.append(joint_fit([first[:m]] + series_full[1:])["C0"])
    spread = max(abs(a - base["C0"]) for a in alts)
    return base, alts, spread
