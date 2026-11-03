"""Physical consistency checks on the committed openEMS results."""
import json
from pathlib import Path

import numpy as np
import pytest

from moreno.touchstone import coupler_metrics, read_touchstone, summarize

RES = Path(__file__).resolve().parents[1] / "results" / "openems"
REF = RES / "reference" / "moreno_openems.s4p"
FULL = RES / "full_matrix" / "moreno_openems.s4p"


@pytest.mark.skipif(not FULL.exists(), reason="no full-matrix result committed")
def test_full_matrix_is_lossless_and_reciprocal():
    _, S, _ = read_touchstone(FULL)
    # PEC model: every column carries unit power
    assert np.all(np.abs(np.sum(np.abs(S) ** 2, axis=1) - 1) < 0.01)
    assert np.max(np.abs(S - np.transpose(S, (0, 2, 1)))) < 0.02


@pytest.mark.skipif(not FULL.exists(), reason="no full-matrix result committed")
def test_full_matrix_symmetry_of_the_junction():
    """Rotation by 180 deg (ports 1<->2, 3<->4) and the x<->z mirror with the guides exchanged
    (ports 1<->3, 2<->4) require |S21| = |S43|, |S31| = |S42|, |S41| = |S32| and equal |Sii|."""
    _, S, _ = read_touchstone(FULL)
    m = np.abs(S)
    assert np.max(np.abs(m[:, 1, 0] - m[:, 3, 2])) < 0.01
    assert np.max(np.abs(m[:, 2, 0] - m[:, 3, 1])) < 0.01
    assert np.max(np.abs(m[:, 3, 0] - m[:, 2, 1])) < 0.01
    # rotation pairs of the reflections are exact; guide 1 versus guide 2 reflections (about -38 dB)
    # differ at the level of the numerical port reflection, see report Sec. 6.1 and 6.3
    assert np.max(np.abs(m[:, 0, 0] - m[:, 1, 1])) < 0.002
    assert np.max(np.abs(m[:, 2, 2] - m[:, 3, 3])) < 0.002
    assert np.max(np.ptp(np.stack([m[:, i, i] for i in range(4)]), axis=0)) < 0.02


@pytest.mark.skipif(not REF.exists(), reason="no reference result committed")
def test_reference_power_balance_and_port_assignment():
    f, S, _ = read_touchstone(REF)
    assert np.all(np.abs(np.sum(np.abs(S[:, :, 0]) ** 2, axis=1) - 1) < 0.01)
    s = summarize(coupler_metrics(f, S), 40.0)
    assert s["directivity_db_at_f0"] > 25          # port 3 is the coupled port
    assert 20 < s["coupling_db_at_f0"] < 35


@pytest.mark.skipif(not REF.exists(), reason="no reference result committed")
def test_reference_summary_matches_touchstone():
    f, S, _ = read_touchstone(REF)
    s = summarize(coupler_metrics(f, S), 40.0)
    stored = json.loads((RES / "reference" / "summary.json").read_text())
    for key in ("coupling_db_at_f0", "directivity_db_at_f0", "return_loss_db_at_f0"):
        assert s[key] == pytest.approx(stored[key], abs=1e-3)


STU = RES / "studies"
KEY = Path(__file__).resolve().parents[1] / "docs" / "report" / "key_numbers.json"


@pytest.mark.skipif(not (STU / "convergence.csv").exists() or not KEY.exists(), reason="no study results committed")
def test_headline_coupling_is_the_convergence_bracket():
    """The two edge treatments converge from opposite sides; the quoted value is the bracket centre."""
    import csv

    def series(name):
        rows = sorted(csv.DictReader(open(STU / name)), key=lambda x: -float(x["slot_cell_mm"]))
        return [float(x["coupling_db_at_f0"]) for x in rows]

    thirds, edges = series("convergence.csv"), series("staircase.csv")
    assert all(b > a for a, b in zip(thirds, thirds[1:]))      # rises on refinement
    assert all(b < a for a, b in zip(edges, edges[1:]))        # falls on refinement
    assert edges[-1] > thirds[-1]                              # still bracketing
    k = json.loads(KEY.read_text())
    assert k["coupling_extrapolated_db"] == pytest.approx((thirds[-1] + edges[-1]) / 2, abs=1e-9)
    assert k["coupling_slot_mesh_uncertainty_db"] == pytest.approx((edges[-1] - thirds[-1]) / 2, abs=1e-9)
    assert k["coupling_uncertainty_db"] >= k["coupling_slot_mesh_uncertainty_db"]
    readme = (Path(__file__).resolve().parents[1] / "README.md").read_text()
    assert "%.2f ± %.2f dB" % (k["coupling_extrapolated_db"], k["coupling_uncertainty_db"]) in readme
