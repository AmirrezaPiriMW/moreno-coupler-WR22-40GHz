"""Figures for the technical report. All data come from design/design.json
and the openEMS result folders; nothing is typed in by hand."""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402

INK = "#17202B"
MUTED = "#5B6472"
BLUE = "#2A5788"
RUST = "#A8461A"
GREEN = "#3C7A4A"
GREY = "#8A8F98"
FILL1 = "#E3E8EE"
FILL2 = "#F3EFE6"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 9.5,
    "axes.labelsize": 9, "axes.edgecolor": INK, "axes.linewidth": 0.8,
    "xtick.color": INK, "ytick.color": INK, "legend.fontsize": 8,
    "legend.frameon": False, "savefig.dpi": 220, "axes.grid": False,
})


def _dim(ax, p0, p1, text, off=(0, 0), color=INK, fs=8, rot=0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="<|-|>", mutation_scale=7, lw=0.7, color=color))
    mx, my = (p0[0] + p1[0]) / 2 + off[0], (p0[1] + p1[1]) / 2 + off[1]
    ax.text(mx, my, text, ha="center", va="center", fontsize=fs, color=color, rotation=rot,
            bbox=dict(fc="white", ec="none", pad=0.6))


def geometry(r: dict, path: Path):
    a, b, t = r["a_mm"], r["b_mm"], r["wall_mm"]
    L, w, d, sp = r["slot_length_mm"], r["slot_width_mm"], r["d_mm"], r["spacing_mm"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.9), gridspec_kw={"width_ratios": [1.25, 1]})

    # --- top view of the junction -----------------------------------------
    ext = 1.25 * a
    ax.add_patch(Rectangle((0, -ext), a, a + 2 * ext, fc=FILL1, ec=INK, lw=1))          # guide 1 (z)
    ax.add_patch(Rectangle((-ext, 0), a + 2 * ext, a, fc=FILL2, ec=INK, lw=1, alpha=0.85))  # guide 2 (x)
    ax.add_patch(Rectangle((0, 0), a, a, fc="none", ec=RUST, lw=1, ls="--"))
    ax.plot([0, a], [0, a], ls=":", lw=0.8, c=MUTED)
    for cx, cz in r["slot_centres_mm"]:
        ax.add_patch(Rectangle((cx - L / 2, cz - w / 2), L, w, fc=INK, ec="none"))
        ax.add_patch(Rectangle((cx - w / 2, cz - L / 2), w, L, fc=INK, ec="none"))
    kw = dict(fontsize=8.5, color=INK, ha="center", va="center", fontweight="bold")
    ax.text(a / 2, -ext + 0.9, "Port 1\ninput", **kw)
    ax.text(a / 2, a + ext - 0.9, "Port 2\nthrough", **kw)
    ax.text(-ext / 2, a / 2, "Port 3\ncoupled", **{**kw, "color": RUST})
    ax.text(a + ext / 2, a / 2, "Port 4\nisolated", **{**kw, "color": MUTED})
    ax.annotate("", xy=(a * 0.82, -ext + 2.6), xytext=(a * 0.82, -ext + 0.4),
                arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.2))
    _dim(ax, (0, -0.9), (a, -0.9), "a = %.3f" % a)
    _dim(ax, (-0.55, 0), (-0.55, d), "d", off=(-0.35, 0))
    _dim(ax, (d, a + 0.75), (a - d, a + 0.75), "")
    ax.text(a / 2, a + 1.0, "Δ = %.3f" % sp, ha="center", va="bottom", fontsize=8)
    ax.set_xlim(-ext, a + ext)
    ax.set_ylim(-ext, a + ext)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title("(a) Top view: guide 1 along z (below), guide 2 along x (above)")

    # --- zoom on the shared wall ------------------------------------------
    ax2.add_patch(Rectangle((0, 0), a, a, fc=FILL1, ec=INK, lw=1))
    ax2.plot([0, a], [0, a], ls=":", lw=0.8, c=MUTED)
    for cx, cz in r["slot_centres_mm"]:
        ax2.add_patch(Rectangle((cx - L / 2, cz - w / 2), L, w, fc=INK, ec="none"))
        ax2.add_patch(Rectangle((cx - w / 2, cz - L / 2), w, L, fc=INK, ec="none"))
        ax2.plot(cx, cz, "o", ms=2.5, c="white")
    c1, c2 = r["slot_centres_mm"][0][0], r["slot_centres_mm"][1][0]
    # L of slot 1's x-arm, drawn below the cross (its z-arm ends at c1 - L/2)
    yL = c1 - L / 2 - 0.28
    for xe in (c1 - L / 2, c1 + L / 2):
        ax2.plot([xe, xe], [c1 - w / 2, yL - 0.08], c=INK, lw=0.4)
    _dim(ax2, (c1 - L / 2, yL), (c1 + L / 2, yL), "", fs=7.5)
    ax2.text(c1 + L / 2 + 0.12, yL, "L = %.2f" % L, fontsize=7.5, ha="left", va="center")
    # w across slot 2's z-arm, just above its top end
    zw = c2 + L / 2 + 0.22
    for xe in (c2 - w / 2, c2 + w / 2):
        ax2.plot([xe, xe], [c2 + L / 2, zw + 0.08], c=INK, lw=0.4)
    ax2.add_patch(FancyArrowPatch((c2 - w / 2 - 0.3, zw), (c2 - w / 2, zw), arrowstyle="-|>", mutation_scale=6,
                                  lw=0.6, color=INK))
    ax2.add_patch(FancyArrowPatch((c2 + w / 2 + 0.3, zw), (c2 + w / 2, zw), arrowstyle="-|>", mutation_scale=6,
                                  lw=0.6, color=INK))
    ax2.text(c2 + w / 2 + 0.36, zw, "w = %.2f" % w, fontsize=7.5, ha="left", va="center",
             bbox=dict(fc=FILL1, ec="none", pad=0.5))
    # arm end to wall
    _dim(ax2, (0, c1), (c1 - L / 2, c1), "", fs=7)
    ax2.text((c1 - L / 2) / 2, c1 + 0.12, "%.2f" % (c1 - L / 2), fontsize=7, ha="center", va="bottom")
    # Delta between the centres, above the wall
    yD = a + 0.3
    for xe in (c1, c2):
        ax2.plot([xe, xe], [c1 if xe == c1 else c2, yD + 0.08], c=MUTED, lw=0.4, ls=":")
    _dim(ax2, (c1, yD), (c2, yD), "", fs=7)
    ax2.text((c1 + c2) / 2, yD + 0.1, "Δ = %.3f" % sp, fontsize=7.5, ha="center", va="bottom")
    # centre coordinates, leaders kept clear of the slots
    ax2.annotate("(%.3f, %.3f)" % (c1, c1), xy=(c1, c1), xytext=(c1 + 0.95, c1 - 0.62), fontsize=7.5,
                 color=BLUE, bbox=dict(fc=FILL1, ec="none", pad=0.5), arrowprops=dict(arrowstyle="-", lw=0.6, color=BLUE, shrinkA=0, shrinkB=2))
    ax2.annotate("(%.3f, %.3f)" % (c2, c2), xy=(c2, c2), xytext=(c2 - 2.9, c2 + 0.55), fontsize=7.5,
                 color=BLUE, bbox=dict(fc=FILL1, ec="none", pad=0.5), arrowprops=dict(arrowstyle="-", lw=0.6, color=BLUE, shrinkA=0, shrinkB=2))
    ax2.set_xlim(-0.3, a + 0.75)
    ax2.set_ylim(-0.3, a + 0.85)
    ax2.set_aspect("equal")
    ax2.set_xlabel("x (mm)")
    ax2.set_ylabel("z (mm)")
    ax2.set_title("(b) Common wall, dimensions in mm")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def cross_section(r: dict, path: Path):
    a, b, t = r["a_mm"], r["b_mm"], r["wall_mm"]
    fig, ax = plt.subplots(figsize=(7.2, 1.9))
    wall = 0.6
    ax.add_patch(Rectangle((-wall, -wall), a + 2 * wall, 2 * b + t + 2 * wall, fc="#C9CDD3", ec=INK, lw=0.8))
    ax.add_patch(Rectangle((0, 0), a, b, fc=FILL1, ec=INK, lw=0.8))
    ax.add_patch(Rectangle((0, b + t), a, b, fc=FILL2, ec=INK, lw=0.8))
    c1 = r["slot_centres_mm"][0][0]
    L = r["slot_length_mm"]
    # a cut at z = d passes only through slot 1's x-arm
    ax.add_patch(Rectangle((c1 - L / 2, b), L, t, fc="white", ec=INK, lw=0.5))
    ax.text(a / 2, b / 2, "guide 1\n(TE$_{10}$ along z)", ha="center", va="center", fontsize=6.5)
    ax.text(a / 2, b + t + b / 2, "guide 2\n(TE$_{10}$ along x)", ha="center", va="center", fontsize=6.5)
    xd = a + wall + 0.35
    _dim(ax, (xd, 0), (xd, b), "", fs=7)
    ax.text(xd + 0.25, b / 2, "b = %.3f" % b, fontsize=7, va="center")
    _dim(ax, (xd, b + t), (xd, 2 * b + t), "", fs=7)
    ax.text(xd + 0.25, b + t + b / 2, "b = %.3f" % b, fontsize=7, va="center")
    ax.annotate("t = %.2f (shared wall)" % t, xy=(a - 0.1, b + t / 2), xytext=(xd + 0.25, b + t / 2),
                fontsize=7, va="center", arrowprops=dict(arrowstyle="-", lw=0.6))
    ax.annotate("slot arm, L = %.2f" % L, xy=(c1, b + t / 2), xytext=(-1.0, b + t / 2), fontsize=7, ha="right",
                va="center", arrowprops=dict(arrowstyle="-", lw=0.6))
    ax.set_xlim(-4.8, a + 4.6)
    ax.set_ylim(-1.0, 2 * b + t + 1.0)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(path, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def coupling_law(r: dict, path: Path):
    a = r["a_mm"]
    x = np.linspace(0, a, 241)
    X, Z = np.meshgrid(x, x)
    to3 = np.abs(np.sin(np.pi * (X + Z) / a))
    to4 = np.abs(np.sin(np.pi * (X - Z) / a))
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.1))
    for ax, data, title in ((axs[0], to3, "Towards port 3 (coupled)\n|sin π(x+z)/a|"),
                            (axs[1], to4, "Towards port 4 (isolated)\n|sin π(x−z)/a|")):
        im = ax.imshow(data, origin="lower", extent=[0, a, 0, a], cmap="Blues", vmin=0, vmax=1)
        ax.plot([0, a], [0, a], c=RUST, lw=0.8, ls="--")
        for cx, cz in r["slot_centres_mm"]:
            ax.plot(cx, cz, "+", ms=11, mew=1.8, c=RUST)
        ax.set_title(title)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("z (mm)")
    fig.colorbar(im, ax=axs, shrink=0.85, label="relative magnetic coupling")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def sparams(metrics_csv: Path, f0: float, path: Path, s4p: Path | None = None):
    import csv

    rows = list(csv.DictReader(open(metrics_csv)))
    f = np.array([float(r_["freq_ghz"]) for r_ in rows])
    g = lambda k: np.array([float(r_[k]) for r_ in rows])  # noqa: E731
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 5.0), sharex=True, gridspec_kw={"height_ratios": [1.6, 1]})
    ax1.plot(f, -g("insertion_loss_db"), c=INK, lw=1.4, label="S21 through")
    ax1.plot(f, -g("coupling_db"), c=BLUE, lw=1.8, label="S31 coupled")
    ax1.plot(f, -g("return_loss_db"), c=GREEN, lw=1.2, label="S11 input match")
    ax1.plot(f, -g("isolation_db"), c=RUST, lw=1.2, label="S41 isolated")
    ax1.set_ylabel("|S| (dB)")
    ax1.set_ylim(-90, 5)
    ax1.grid(alpha=0.25)
    ax1.legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    ax2.plot(f, g("directivity_db"), c=RUST, lw=1.6)
    ax2.set_ylabel("Directivity (dB)")
    ax2.set_xlabel("Frequency (GHz)")
    ax2.grid(alpha=0.25)
    for ax in (ax1, ax2):
        ax.axvline(f0, ls="--", c=GREY, lw=0.8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def bar_line(xs, ys, path, xlabel, ylabel, title, extra=None, invert=False, xlog=False):
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    for label, y, c, m in ys:
        ax.plot(xs, y, marker=m, c=c, lw=1.4, ms=4.5, label=label)
    if extra:
        extra(ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.25)
    if invert:
        ax.invert_yaxis()
    if xlog:
        ax.set_xscale("log")
    if len(ys) > 1:
        ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figures for the verification and design-chart sections
# ---------------------------------------------------------------------------
def convergence(conv, stair, legacy, rich, path):
    """(a) coupling vs h, corrected model; (b) corrected vs superseded model; (c) directivity and return loss."""
    fig, (ax1, ax0, ax2) = plt.subplots(1, 3, figsize=(7.2, 2.8), gridspec_kw={"width_ratios": [1.25, 1, 1]})
    h = [c["slot_cell_mm"] for c in conv]
    C = [c["coupling_db_at_f0"] for c in conv]
    if rich:
        ax1.axhspan(rich["f_ext"] - rich["gci"], rich["f_ext"] + rich["gci"], color=BLUE, alpha=0.12, lw=0)
        ax1.axhline(rich["f_ext"], c=BLUE, lw=0.6, ls=":")
        ax1.plot([0], [rich["f_ext"]], "*", c=RUST, ms=10, label="converged (bracket)", zorder=5)
    ax1.plot(h, C, "o-", c=BLUE, lw=1.6, ms=4.5, label="1/3–2/3 edge rule")
    if stair:
        st_ = sorted(stair, key=lambda s_: s_["slot_cell_mm"])
        ax1.plot([s_["slot_cell_mm"] for s_ in st_], [s_["coupling_db_at_f0"] for s_ in st_], "s-", c=GREEN,
                 lw=1.4, ms=4.5, label="lines on slot edges")
    ax1.set_xlim(-0.003, max(h) * 1.08)
    ax1.set_xlabel("slot / wall cell h (mm)")
    ax1.set_ylabel("coupling at 40 GHz (dB)")
    ax1.invert_yaxis()
    ax1.grid(alpha=0.25)
    ax1.legend(fontsize=6.2, loc="upper left", handlelength=1.5, borderaxespad=0.3)
    ax1.set_title("(a) Metal-only model", fontsize=8.5)
    if legacy:
        ax0.plot([s_["slot_cell_mm"] for s_ in legacy], [s_["coupling_db_at_f0"] for s_ in legacy], "^--", c=GREY,
                 lw=1, ms=4, label="superseded")
    ax0.plot(h, C, "o-", c=BLUE, lw=1.4, ms=3.5, label="metal-only")
    ax0.set_xlabel("slot / wall cell h (mm)")
    ax0.set_ylabel("coupling at 40 GHz (dB)")
    ax0.invert_yaxis()
    ax0.grid(alpha=0.25)
    ax0.legend(fontsize=6.6, loc="upper left")
    ax0.set_title("(b) Versus superseded model", fontsize=8.5)
    ax2.plot(h, [c["directivity_db_at_f0"] for c in conv], "o-", c=RUST, lw=1.4, ms=4, label="directivity")
    ax2.plot(h, [c["return_loss_db_at_f0"] for c in conv], "o-", c=GREEN, lw=1.2, ms=4, label="return loss")
    ax2.set_xlim(left=0)
    ax2.set_xlabel("slot / wall cell h (mm)")
    ax2.set_ylabel("dB at 40 GHz")
    ax2.grid(alpha=0.25)
    ax2.legend(fontsize=6.6)
    ax2.set_title("(c) Directivity, return loss", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def design_chart(Ls, C_sim, C_corr, est, targets, path):
    fig, ax = plt.subplots(figsize=(7.2, 3.3))
    ax.plot(Ls, C_sim, "o", c=GREY, ms=4.5, label="openEMS, h = 0.035 mm")
    xs = np.linspace(min(Ls), max(Ls), 200)
    ax.plot(xs, C_corr(xs), "-", c=BLUE, lw=1.8, label="mesh-corrected (h → 0)")
    ax.plot(Ls, est, "--", c=MUTED, lw=1.1, label="closed-form estimate (Sec. 3.5)")
    for C, L in targets:
        if L is None:
            continue
        ax.plot([L, L], [C, ax.get_ylim()[0] if False else 40], ":", c=RUST, lw=0.8)
        ax.plot(L, C, "D", c=RUST, ms=5)
        ax.annotate("%g dB → L = %.3f mm" % (C, L), xy=(L, C), xytext=(L + 0.05, C + 2.2), fontsize=7.5,
                    color=RUST, bbox=dict(fc="white", ec="none", pad=1.0),
                    arrowprops=dict(arrowstyle="-", lw=0.6, color=RUST))
    ax.set_xlabel("slot arm length L (mm)")
    ax.set_ylabel("coupling at 40 GHz (dB)")
    ax.set_ylim(max(max(C_sim), max(est)) + 1.5, min(min(C_sim), min(est)) - 1.5)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def freq_family(curves, path, key="coupling_db", ylabel="coupling (dB)", invert=True, f0=40.0):
    """curves: list of (label, freq_ghz array, values array)."""
    fig, ax = plt.subplots(figsize=(3.5, 3.2))
    cmap = plt.get_cmap("viridis")
    for i, (label, f, v) in enumerate(curves):
        ax.plot(f, v, c=cmap(i / max(1, len(curves) - 1) * 0.85), lw=1.4, label=label)
    ax.axvline(f0, ls="--", c=GREY, lw=0.8)
    ax.set_xlabel("frequency (GHz)")
    ax.set_ylabel(ylabel)
    if invert:
        ax.invert_yaxis()
    ax.grid(alpha=0.25)
    ax.legend(fontsize=6.5, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.22), handlelength=1.4,
              columnspacing=0.9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def field_maps(npz_path, r, path, npml=8):
    """|Ey| at f0 on the mid-planes of both guides, dB relative to the incident wave peak.

    Both guides are drawn along their own axis of propagation; the 8-cell PML at each
    feed end is left out.
    """
    D = np.load(npz_path)
    a = float(D["a_mm"])
    L, w = r["slot_length_mm"], r["slot_width_mm"]
    ref = np.max(np.abs(D["Ef_guide1_Ey"]))
    fig, axs = plt.subplots(2, 1, figsize=(7.2, 3.9), sharex=True)
    sl = slice(npml, -npml)
    # (a) guide 1: horizontal axis z (propagation), vertical axis x
    x, z, E = D["Ef_guide1_x_mm"], D["Ef_guide1_z_mm"], np.abs(D["Ef_guide1_Ey"])
    ix = (x >= 0) & (x <= a)
    x, z, E = x[ix], z[sl], E[ix][:, sl]
    im = axs[0].pcolormesh(z, x, 20 * np.log10(np.maximum(E / ref, 1e-6)), shading="auto", cmap="magma",
                           vmin=-70, vmax=0, rasterized=True)
    # (b) guide 2: horizontal axis x (propagation), vertical axis z
    x2, z2, E2 = D["Ef_guide2_x_mm"], D["Ef_guide2_z_mm"], np.abs(D["Ef_guide2_Ey"])
    iz = (z2 >= 0) & (z2 <= a)
    x2, z2, E2 = x2[sl], z2[iz], E2[sl][:, iz]
    axs[1].pcolormesh(x2, z2, 20 * np.log10(np.maximum(E2 / ref, 1e-6)).T, shading="auto", cmap="magma",
                      vmin=-70, vmax=0, rasterized=True)
    for k, ax in enumerate(axs):
        for cx, cz in r["slot_centres_mm"]:
            # guide 1 is drawn with (z, x) axes, guide 2 with (x, z)
            u, v = (cz, cx) if k == 0 else (cx, cz)
            for (du, dv) in ((L, w), (w, L)):
                ax.add_patch(Rectangle((u - du / 2, v - dv / 2), du, dv, fill=False, ec="#9fd3ff", lw=0.6))
        ax.set_aspect("equal")
    axs[0].set_ylabel("x (mm)")
    axs[0].set_title("(a) Guide 1, y = b/2: port 1 (left) → port 2 (right)",
                     fontsize=8, loc="left")
    axs[1].set_ylabel("z (mm)")
    axs[1].set_xlabel("position along the guide (mm): z for guide 1, x for guide 2")
    axs[1].set_title("(b) Guide 2, y = b + t + b/2: port 3 coupled (left), port 4 isolated (right)",
                     fontsize=8, loc="left")
    cb = fig.colorbar(im, ax=axs, shrink=0.9, pad=0.02)
    cb.set_label("|E$_y$| (dB rel. incident peak)")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
