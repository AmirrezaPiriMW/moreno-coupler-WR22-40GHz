"""First-order absolute coupling estimate (Bethe small-aperture theory).

For TE10 in two identical guides crossing at right angles, the magnetic-dipole
coupling of one slot at (x0, z0) into the coupled port is, for incident TE10
amplitude 1 (Pozar, Microwave Engineering, Sec. 4.8 and 7.4),

    |c| = pi * alpha_m * |sin(pi (x0 + z0) / a)| / (a^2 b)

The guide wave impedance cancels, so the magnetic coupling of crossed guides
is independent of frequency; only alpha_m, the wall attenuation and the
array factor of the slot pair vary.  This module adds three first-order
refinements:

* narrow-slot polarizability along each arm (Collin, Foundations for
  Microwave Engineering, Table 4.3):  alpha_m = pi L^3 / (24 (ln(4L/w) - 1));
  the cross is taken as in-plane isotropic, alpha_xx = alpha_zz = alpha_m;
* resonance enhancement 1 / (1 - (2 L / lambda0)^2) (Cohn's correction);
* attenuation through the wall as a below-cut-off guide of cut-off 2L.

Electric-dipole coupling is neglected (it limits directivity, not coupling).
The estimate is meant for orientation; full-wave tuning sets the final L.
See docs/report for its measured accuracy against openEMS.
"""
from __future__ import annotations

import math

from .design import free_space_wavelength_mm, pair_factor, slot_resonance_factor, slot_wall_loss_db


def slot_polarizability_mm3(L: float, w: float) -> float:
    return math.pi * L ** 3 / (24.0 * (math.log(4.0 * L / w) - 1.0))


resonance_factor = slot_resonance_factor


def coupling_estimate_db(f_ghz, a, b, d, spacing, L, w, t, resonance=True, wall=True) -> float:
    """Estimated coupling S31 of the slot pair, in positive dB."""
    alpha = slot_polarizability_mm3(L, w)
    mag = math.pi * alpha / (a * a * b) * pair_factor(f_ghz, a, d, spacing)
    if resonance:
        mag *= resonance_factor(f_ghz, L)
    db = -20.0 * math.log10(mag)
    if wall:
        db += slot_wall_loss_db(f_ghz, L, t)
    return db


def estimate_from_design(r: dict, f_ghz: float | None = None, **kw) -> float:
    return coupling_estimate_db(
        f_ghz or r["f0_ghz"], r["a_mm"], r["b_mm"], r["d_mm"], r["spacing_mm"],
        r["slot_length_mm"], r["slot_width_mm"], r["wall_mm"], **kw)


def directivity_estimate_db(r: dict, f_ghz: float | None = None, arms: int = 1) -> float:
    """First-order directivity limit set by the slots' electric-dipole coupling.

    On the diagonal the magnetic coupling to the isolated port vanishes; what remains
    is electric coupling through the normal E field, which is not directional and,
    with lambda_g/4 spacing, adds in phase at port 4. Narrow-slot electric
    polarizability alpha_e = pi w^2 L / 16 per arm (Collin); `arms` = 1 or 2 brackets
    how the two arms of a cross combine. Through the wall, the electric coupling
    decays like the slot's cut-off mode across its width (kc = pi/w), the magnetic
    one like the mode along its length (kc = pi/L).
    """
    f = f_ghz or r["f0_ghz"]
    a, d, L, w, t = r["a_mm"], r["d_mm"], r["slot_length_mm"], r["slot_width_mm"], r["wall_mm"]
    fc = r["fc_te10_ghz"]
    k0 = 2.0 * math.pi / free_space_wavelength_mm(f)
    beta = k0 * math.sqrt(1.0 - (fc / f) ** 2)
    p = beta * a / math.pi
    zr = 1.0 / (1.0 - (fc / f) ** 2)            # (Z_TE / eta0)^2
    X = math.pi * d / a
    alpha_e = arms * math.pi * w * w * L / 16.0
    ratio = zr * p * (alpha_e / slot_polarizability_mm3(L, w)) * math.sin(X) ** 2 / abs(math.sin(2 * X))

    def att(kc):
        return 8.686 * t * math.sqrt(max(kc * kc - k0 * k0, 0.0))

    return -20.0 * math.log10(ratio) + att(math.pi / w) - att(math.pi / L)
