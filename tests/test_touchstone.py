import numpy as np
import pytest

from moreno.touchstone import coupler_metrics, read_touchstone, summarize


def _write_s4p(path, freqs_ghz, S):
    lines = ["! synthetic", "# GHz S MA R 50"]
    for f, s in zip(freqs_ghz, S):
        vals = []
        for i in range(4):
            for j in range(4):
                vals += [abs(s[i, j]), np.degrees(np.angle(s[i, j]))]
        lines.append("%g " % f + " ".join("%.9g" % v for v in vals[:8]))
        for k in range(8, 32, 8):
            lines.append(" ".join("%.9g" % v for v in vals[k:k + 8]))
    path.write_text("\n".join(lines))


def test_four_port_roundtrip(tmp_path):
    f = np.array([39.0, 40.0, 41.0])
    S = np.zeros((3, 4, 4), complex)
    S[:, 0, 0] = 0.05
    S[:, 1, 0] = 0.999 * np.exp(-1j * 0.3)
    S[:, 2, 0] = 10 ** (-30 / 20)
    S[:, 3, 0] = 10 ** (-55 / 20)
    p = tmp_path / "c.s4p"
    _write_s4p(p, f, S)

    freq, S2, z0 = read_touchstone(p)
    assert z0 == 50
    assert np.allclose(freq, f * 1e9)
    assert np.allclose(S2, S, atol=1e-8)

    m = coupler_metrics(freq, S2)
    s = summarize(m, 40.0)
    assert s["coupling_db_at_f0"] == pytest.approx(30.0, abs=1e-6)
    assert s["directivity_db_at_f0"] == pytest.approx(25.0, abs=1e-6)


def test_two_port_column_order(tmp_path):
    p = tmp_path / "x.s2p"
    p.write_text("# GHz S RI R 50\n10 0.1 0 0.9 0 0.2 0 0.3 0\n")
    _, S, _ = read_touchstone(p)
    assert S[0, 1, 0] == pytest.approx(0.9)  # S21
    assert S[0, 0, 1] == pytest.approx(0.2)  # S12
