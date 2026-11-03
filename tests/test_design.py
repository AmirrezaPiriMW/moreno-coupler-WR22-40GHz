import math

import pytest

from moreno.design import MorenoInputs, compute, guide_wavelength_mm, pair_factor, ripple


def test_wr22_40ghz_numbers():
    r = compute(MorenoInputs())
    assert r["fc_te10_ghz"] == pytest.approx(26.344, abs=1e-3)
    assert r["lambda_g_mm"] == pytest.approx(9.960, abs=2e-3)
    assert r["spacing_mm"] == pytest.approx(2.490, abs=1e-3)
    assert r["d_mm"] == pytest.approx(1.600, abs=1e-3)
    assert r["slot_centres_mm"][1][0] == pytest.approx(4.090, abs=1e-3)
    assert r["all_checks_ok"]


def test_wr90_10ghz_guide_wavelength():
    assert guide_wavelength_mm(10.0, 22.86) == pytest.approx(39.71, abs=0.02)


def test_quarter_wave_spacing_is_90_degrees():
    r = compute(MorenoInputs())
    assert r["electrical_spacing_deg"] == pytest.approx(90.0, abs=1e-9)


def test_below_cutoff_raises():
    with pytest.raises(ValueError):
        compute(MorenoInputs(waveguide="WR-22", f0_ghz=20.0))


def test_pair_peaks_with_quarter_wave_spacing():
    a = 5.69
    lamg = guide_wavelength_mm(40.0, a)
    best = pair_factor(40.0, a, (a - lamg / 4) / 2, lamg / 4)
    worse = pair_factor(40.0, a, (a - lamg / 8) / 2, lamg / 8)
    assert best > worse


def test_ripple_zero_at_f0():
    r = compute(MorenoInputs())
    data = ripple(r, n=201)
    f0_val = min(data, key=lambda p: abs(p[0] - 40.0))
    assert abs(f0_val[1]) < 0.05
    assert all(math.isfinite(v) for _, v in data)


def test_wall_clearance_is_arm_end_to_wall():
    r = compute(MorenoInputs())
    # arm end face at x = d - L/2 = 1.600 - 1.000
    assert r["clearance_mm"] == pytest.approx(0.600, abs=1e-3)
    # nearest features are arm ends offset diagonally: sqrt(2) * (Delta - (L + w) / 2)
    assert r["slot_gap_mm"] == pytest.approx(math.sqrt(2) * (r["spacing_mm"] - (2.0 + 0.35) / 2), abs=1e-12)


def test_overlapping_slots_fail_check():
    r = compute(MorenoInputs(slot_length_mm=4.8))
    names = {c["name"]: c["ok"] for c in r["checks"]}
    assert not names["slots_separate"] and not names["wall_clearance"]


def test_directivity_law_from_mode_fields():
    """Magnetic-dipole coupling equals p*sin(pi(x+z)/a) to port 3 and p*sin(pi(x-z)/a) to port 4."""
    import numpy as np

    a, f = 5.69, 40.0
    p = 2 * math.pi / guide_wavelength_mm(f, a) * a / math.pi
    rng = np.random.default_rng(0)
    for x, z in rng.uniform(0, a, (500, 2)):
        X, Z = math.pi * x / a, math.pi * z / a
        h1 = np.array([1j * p * math.sin(X), math.cos(X)])        # guide 1, +z
        h2_plus = np.array([math.cos(Z), 1j * p * math.sin(Z)])   # guide 2, +x
        h2_minus = np.array([-math.cos(Z), 1j * p * math.sin(Z)])  # guide 2, -x
        assert abs(h1 @ h2_plus) == pytest.approx(p * abs(math.sin(X + Z)), abs=1e-12)
        assert abs(h1 @ h2_minus) == pytest.approx(p * abs(math.sin(X - Z)), abs=1e-12)


def test_conductor_loss_matches_textbook_wr90():
    from moreno.design import conductor_loss_db_per_m

    # copper WR-90 at 10 GHz: alpha_c = 0.0125 Np/m = 0.108 dB/m (standard handbook value)
    assert conductor_loss_db_per_m(10.0, 22.86, 10.16) == pytest.approx(0.108, abs=0.002)


def test_slot_resonance_factor():
    from moreno.design import slot_resonance_factor

    assert slot_resonance_factor(40.0, 1e-6) == pytest.approx(1.0)
    assert slot_resonance_factor(40.0, 299.792458 / 40.0 / 2.0) == float("inf")


def test_crossed_guide_magnetic_coupling_is_frequency_independent():
    """|c| = pi alpha |sin(pi(x+z)/a)| / (a^2 b): no beta left after normalisation."""
    from moreno.estimate import coupling_estimate_db

    a, b, d, sp, L, w, t = 5.69, 2.845, 1.6, 2.49, 2.0, 0.35, 0.25
    c1 = coupling_estimate_db(36.0, a, b, d, sp, L, w, t, resonance=False, wall=False)
    c2 = coupling_estimate_db(44.0, a, b, d, sp, L, w, t, resonance=False, wall=False)
    # only the array factor sin(beta * Delta) differs between the two frequencies
    assert abs(c1 - c2) < 0.6
