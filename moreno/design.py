"""Moreno cross-guide directional coupler: design equations.

Two rectangular waveguides cross at right angles and share one broad wall.
Two "+"-shaped crossed slots cut through that wall couple energy from the
lower guide (guide 1) into the upper guide (guide 2).

Coordinates (all lengths in mm, frequencies in GHz)
---------------------------------------------------
* Guide 1 runs along +z:  x in [0, a], y in [0, b].
* Shared wall:            y in [b, b + t].
* Guide 2 runs along x:   z in [0, a], y in [b + t, 2b + t].
* The shared square is x, z in [0, a].

Port map: 1 = input (-z end of guide 1), 2 = through (+z),
          3 = coupled (-x end of guide 2), 4 = isolated (+x).

Physics (small-aperture, magnetic-dipole model)
-----------------------------------------------
For TE10 in both guides, a slot whose magnetic polarizability is the same
in every in-plane direction, centred at (x0, z0), excites guide 2 with

    towards port 3:  A3 ~ sin(pi (x0 + z0) / a) * exp(-j beta (x0 + z0))
    towards port 4:  A4 ~ sin(pi (x0 - z0) / a) * exp(-j beta (z0 - x0))

so any slot on the diagonal x0 = z0 is ideally directional.  The two slots
sit at (d, d) and (a - d, a - d).  On their own they would couple with
opposite sign; spacing them by Delta = a - 2d = lambda_g / 4 adds the 180
degree path difference that makes them reinforce at port 3.

Absolute coupling depends on the slot polarizability, which this model
holds constant.  Set it by sweeping the slot length in a full-wave solver.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

C0_MM_GHZ = 299.792458  # speed of light in mm * GHz

# name: (a_mm, b_mm, band_low_GHz, band_high_GHz)
WAVEGUIDES = {
    "WR-137": (34.849, 15.799, 5.85, 8.20),
    "WR-90": (22.860, 10.160, 8.20, 12.40),
    "WR-62": (15.799, 7.899, 12.40, 18.00),
    "WR-42": (10.668, 4.318, 18.00, 26.50),
    "WR-28": (7.112, 3.556, 26.50, 40.00),
    "WR-22": (5.690, 2.845, 33.00, 50.00),
    "WR-15": (3.759, 1.880, 50.00, 75.00),
}


# --------------------------------------------------------------------------
# Basic waveguide relations
# --------------------------------------------------------------------------
def te10_cutoff_ghz(a_mm: float) -> float:
    return C0_MM_GHZ / (2.0 * a_mm)


def free_space_wavelength_mm(f_ghz: float) -> float:
    return C0_MM_GHZ / f_ghz


def guide_wavelength_mm(f_ghz: float, a_mm: float) -> float:
    """TE10 guide wavelength; NaN below cut-off."""
    ratio = te10_cutoff_ghz(a_mm) / f_ghz
    if ratio >= 1.0:
        return float("nan")
    return free_space_wavelength_mm(f_ghz) / math.sqrt(1.0 - ratio * ratio)


def slot_wall_loss_db(f_ghz: float, slot_length_mm: float, wall_mm: float) -> float:
    """Loss through a thick wall, treating the slot as a below-cut-off guide.

    The slot's TE10 cut-off wavelength is 2L.  Above resonance (L > lambda0/2)
    the slot propagates and this first-order loss is zero.
    """
    k0 = 2.0 * math.pi / free_space_wavelength_mm(f_ghz)
    kc = math.pi / slot_length_mm
    if kc <= k0:
        return 0.0
    return 8.686 * math.sqrt(kc * kc - k0 * k0) * wall_mm


def slot_resonance_factor(f_ghz: float, slot_length_mm: float) -> float:
    """First-order growth of a slot's magnetic polarizability towards resonance.

    alpha_m(f) = alpha_m(0) / (1 - (2 L / lambda0)^2): the slot's half-wave
    resonance is at L = lambda0 / 2 (Cohn's resonant-aperture correction).
    """
    x = 2.0 * slot_length_mm / free_space_wavelength_mm(f_ghz)
    return float("inf") if x >= 1.0 else 1.0 / (1.0 - x * x)


def conductor_loss_db_per_m(f_ghz: float, a_mm: float, b_mm: float, sigma_s_per_m: float = 5.8e7) -> float:
    """TE10 conductor attenuation of a rectangular guide (Pozar, Microwave Engineering, Sec. 3.3).

    alpha_c = Rs / (a^3 b beta k eta) * (2 b pi^2 + a^3 k^2)   [Np/m], smooth walls.
    """
    a, b = a_mm * 1e-3, b_mm * 1e-3
    f = f_ghz * 1e9
    mu0, eta = 4e-7 * math.pi, 376.730313
    k = 2.0 * math.pi * f / (C0_MM_GHZ * 1e6)
    kc = math.pi / a
    if k <= kc:
        return float("inf")
    beta = math.sqrt(k * k - kc * kc)
    rs = math.sqrt(math.pi * f * mu0 / sigma_s_per_m)
    alpha = rs / (a ** 3 * b * beta * k * eta) * (2.0 * b * math.pi ** 2 + a ** 3 * k ** 2)
    return 8.685889638 * alpha


def _rect_gap(r1, r2) -> float:
    """Shortest distance between two axis-aligned rectangles (x0, x1, z0, z1); 0 if they touch."""
    dx = max(r2[0] - r1[1], r1[0] - r2[1], 0.0)
    dz = max(r2[2] - r1[3], r1[2] - r2[3], 0.0)
    return math.hypot(dx, dz)


def slot_arms(a_mm: float, d_mm: float, L: float, w: float):
    """The four arms of the two crosses as (x0, x1, z0, z1), slot 1 first."""
    arms = []
    for c in (d_mm, a_mm - d_mm):
        arms.append((c - L / 2, c + L / 2, c - w / 2, c + w / 2))
        arms.append((c - w / 2, c + w / 2, c - L / 2, c + L / 2))
    return arms


# --------------------------------------------------------------------------
# Design
# --------------------------------------------------------------------------
@dataclass
class MorenoInputs:
    waveguide: str = "WR-22"
    f0_ghz: float = 40.0
    slot_length_mm: float = 2.0
    slot_width_mm: float = 0.35
    wall_mm: float = 0.25
    spacing_mm: Optional[float] = None  # None -> lambda_g / 4
    feed_lengths_lg: float = 1.0        # feed length beyond the shared square
    sweep_fraction: float = 0.125       # simulated band = f0 * (1 +/- this)

    def guide(self):
        if self.waveguide not in WAVEGUIDES:
            raise ValueError(
                "Unknown waveguide %r; choose one of %s" % (self.waveguide, ", ".join(WAVEGUIDES))
            )
        return WAVEGUIDES[self.waveguide]


def pair_factor(f_ghz: float, a_mm: float, d_mm: float, spacing_mm: float) -> float:
    """|A3| of the slot pair relative to one slot at the best spot (d = a/4)."""
    beta = 2.0 * math.pi / guide_wavelength_mm(f_ghz, a_mm)
    return 2.0 * abs(math.sin(2.0 * math.pi * d_mm / a_mm)) * abs(math.sin(beta * spacing_mm))


def compute(inp: MorenoInputs) -> dict:
    a, b, band_lo, band_hi = inp.guide()
    f0 = inp.f0_ghz
    L, w, t = inp.slot_length_mm, inp.slot_width_mm, inp.wall_mm

    fc10 = te10_cutoff_ghz(a)
    lam0 = free_space_wavelength_mm(f0)
    lamg = guide_wavelength_mm(f0, a)
    if math.isnan(lamg):
        raise ValueError("f0 = %.3f GHz is below the TE10 cut-off (%.3f GHz)" % (f0, fc10))
    beta = 2.0 * math.pi / lamg

    spacing = inp.spacing_mm if inp.spacing_mm else lamg / 4.0
    d = (a - spacing) / 2.0
    pair = pair_factor(f0, a, d, spacing)
    pair_db = 20.0 * math.log10(max(pair, 1e-12))
    wall_db = slot_wall_loss_db(f0, L, t)
    l_res = lam0 / 2.0
    # The arm end face (x = d - L/2, or z = d - L/2 for the other arm) is the
    # feature closest to a side wall; by symmetry slot 2 has the same margin.
    clearance = d - L / 2.0
    # Narrowest metal web between the two crosses: shortest distance between any
    # arm of slot 1 and any arm of slot 2 (0 when they touch or overlap).
    arms = slot_arms(a, d, L, w)
    slot_gap = min(_rect_gap(p, q) for p in arms[:2] for q in arms[2:])
    min_clearance = max(0.1, 0.02 * a)
    spacing_deg = math.degrees(beta * spacing)
    feed = round(inp.feed_lengths_lg * lamg, 3)

    sweep_lo = max(band_lo, f0 * (1.0 - inp.sweep_fraction), fc10 * 1.05)
    sweep_hi = min(band_hi, f0 * (1.0 + inp.sweep_fraction), 2.0 * fc10 * 0.98)

    checks = [
        {"name": "single_mode", "ok": 1.25 * fc10 <= f0 <= 0.95 * 2.0 * fc10,
         "detail": "f0 within 1.25 fc .. 0.95 * 2 fc (%.2f .. %.2f GHz)" % (1.25 * fc10, 1.9 * fc10)},
        {"name": "wall_clearance", "ok": clearance >= min_clearance,
         "detail": ("slot arms clear side walls by %.3f mm (want >= %.3f)" % (clearance, min_clearance)
                    if clearance >= 0 else
                    "slot arms run %.3f mm past the side walls (want >= %.3f clear)" % (-clearance, min_clearance))},
        {"name": "below_resonance", "ok": L <= 0.9 * l_res,
         "detail": "L = %.0f%% of resonant length lambda0/2 = %.3f mm" % (100 * L / l_res, l_res)},
        {"name": "electrical_spacing", "ok": 70.0 <= spacing_deg <= 110.0,
         "detail": "beta * Delta = %.1f deg (aim for 90)" % spacing_deg},
        {"name": "slots_separate", "ok": slot_gap >= min_clearance,
         "detail": "narrowest metal web between the two crosses %.3f mm (want >= %.3f)" % (slot_gap, min_clearance)},
        {"name": "narrow_slot", "ok": w <= L / 4.0,
         "detail": "arm aspect ratio L/w = %.1f (want >= 4 for a narrow-slot cross)" % (L / w)},
    ]

    return {
        "inputs": asdict(inp),
        "a_mm": a, "b_mm": b, "band_ghz": [band_lo, band_hi],
        "f0_ghz": f0,
        "fc_te10_ghz": fc10, "fc_te20_ghz": 2.0 * fc10,
        "lambda0_mm": lam0, "lambda_g_mm": lamg, "beta_rad_per_mm": beta,
        "spacing_mm": spacing, "spacing_diagonal_mm": spacing * math.sqrt(2.0),
        "d_mm": d,
        "slot_centres_mm": [[d, d], [a - d, a - d]],
        "slot_length_mm": L, "slot_width_mm": w, "wall_mm": t,
        "slot_resonant_length_mm": l_res,
        "pair_vs_single_slot_db": pair_db,
        "copper_loss_db_per_m": conductor_loss_db_per_m(f0, a, b),
        "wall_loss_db": wall_db,
        "clearance_mm": clearance,
        "slot_gap_mm": slot_gap,
        "electrical_spacing_deg": spacing_deg,
        "feed_mm": feed,
        "sweep_ghz": [round(sweep_lo, 4), round(sweep_hi, 4)],
        "checks": checks,
        "all_checks_ok": all(c["ok"] for c in checks),
    }


def ripple(result: dict, n: int = 101):
    """Coupling variation (dB, relative to f0) from the closed-form model across the band.

    For crossed guides the magnetic coupling of a slot does not depend on frequency
    by itself (the guide impedance cancels; see docs/report), so the variation comes
    from the slot polarizability growing towards resonance, the array factor of the
    pair and the attenuation through the wall.
    """
    a, d, sp = result["a_mm"], result["d_mm"], result["spacing_mm"]
    L, t = result["slot_length_mm"], result["wall_mm"]
    lo, hi = result["band_ghz"]
    f0 = result["f0_ghz"]

    def level(f):
        return (20 * math.log10(max(pair_factor(f, a, d, sp), 1e-12) * slot_resonance_factor(f, L))
                - slot_wall_loss_db(f, L, t))

    ref = level(f0)
    out = []
    for i in range(n):
        f = lo + (hi - lo) * i / (n - 1)
        out.append((f, level(f) - ref))
    return out


# --------------------------------------------------------------------------
# Output helpers
# --------------------------------------------------------------------------
def summary_table(r: dict) -> str:
    rows = [
        ("Waveguide", "%s  a = %.3f mm, b = %.3f mm" % (r["inputs"]["waveguide"], r["a_mm"], r["b_mm"])),
        ("f0", "%.3f GHz" % r["f0_ghz"]),
        ("TE10 / TE20 cut-off", "%.3f / %.3f GHz" % (r["fc_te10_ghz"], r["fc_te20_ghz"])),
        ("lambda0 / lambda_g", "%.3f / %.3f mm" % (r["lambda0_mm"], r["lambda_g_mm"])),
        ("Slot spacing Delta", "%.3f mm per axis (%.3f mm on diagonal)" % (r["spacing_mm"], r["spacing_diagonal_mm"])),
        ("Slot offset d", "%.3f mm" % r["d_mm"]),
        ("Slot centres (x, z)", "(%.3f, %.3f), (%.3f, %.3f) mm" % tuple(sum(r["slot_centres_mm"], []))),
        ("Slot arm L x w", "%.3f x %.3f mm, wall t = %.3f mm" % (r["slot_length_mm"], r["slot_width_mm"], r["wall_mm"])),
        ("Pair vs one slot at a/4", "%+.2f dB" % r["pair_vs_single_slot_db"]),
        ("Thick-wall slot loss", "%.2f dB" % r["wall_loss_db"]),
        ("Feed length", "%.3f mm" % r["feed_mm"]),
        ("Simulation sweep", "%.3f - %.3f GHz" % tuple(r["sweep_ghz"])),
    ]
    width = max(len(k) for k, _ in rows)
    lines = ["%s  %s" % (k.ljust(width), v) for k, v in rows]
    lines.append("")
    for c in r["checks"]:
        lines.append("[%s] %s" % ("ok " if c["ok"] else "!! ", c["detail"]))
    return "\n".join(lines)


def plot_layout(r: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    a, L, w = r["a_mm"], r["slot_length_mm"], r["slot_width_mm"]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    ax.add_patch(Rectangle((0, 0), a, a, fill=True, fc="#EEF2F6", ec="#17202B", lw=2))
    ax.plot([0, a], [0, a], ls="--", lw=1, c="#4A5361")
    for cx, cz in r["slot_centres_mm"]:
        ax.add_patch(Rectangle((cx - L / 2, cz - w / 2), L, w, fc="#17202B"))
        ax.add_patch(Rectangle((cx - w / 2, cz - L / 2), w, L, fc="#17202B"))
    ax.set_xlim(-0.08 * a, 1.08 * a)
    ax.set_ylim(1.08 * a, -0.08 * a)  # z grows downward, input at top
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)  ->  port 4 isolated")
    ax.set_ylabel("z (mm)  ->  port 2 through")
    ax.set_title("Shared wall, top view (port 1 input at top, port 3 coupled at left)", fontsize=9)

    data = ripple(r)
    ax2.plot([f for f, _ in data], [v for _, v in data], c="#2A5788", lw=2)
    ax2.axvline(r["f0_ghz"], ls="--", c="#A8461A", lw=1)
    ax2.set_xlabel("Frequency (GHz)")
    ax2.set_ylabel("|S31| relative to f0 (dB)")
    ax2.set_title("Closed-form trend of |S31| across the band (up = stronger coupling)", fontsize=9)
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Moreno cross-guide coupler calculator")
    p.add_argument("--wg", default="WR-22", choices=sorted(WAVEGUIDES))
    p.add_argument("--f0", type=float, default=40.0, help="centre frequency, GHz")
    p.add_argument("--L", type=float, default=2.0, help="slot arm length, mm")
    p.add_argument("--w", type=float, default=0.35, help="slot arm width, mm")
    p.add_argument("--t", type=float, default=0.25, help="shared wall thickness, mm")
    p.add_argument("--spacing", type=float, default=None, help="slot spacing per axis, mm (default lambda_g/4)")
    p.add_argument("--feed", type=float, default=1.0, help="feed length in guide wavelengths")
    p.add_argument("--out", default="design/design.json", help="where to write the design JSON")
    p.add_argument("--plot", default=None, help="optional PNG of the slot layout and ripple")
    args = p.parse_args(argv)

    inp = MorenoInputs(args.wg, args.f0, args.L, args.w, args.t, args.spacing, args.feed)
    r = compute(inp)
    print(summary_table(r))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(r, indent=2))
    print("\nWrote %s" % out)
    if args.plot:
        plot_layout(r, Path(args.plot))
        print("Wrote %s" % args.plot)
    return 0 if r["all_checks_ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
