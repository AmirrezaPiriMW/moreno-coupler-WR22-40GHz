"""Minimal Touchstone v1 (.sNp) reader and coupler metrics. Needs only numpy."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

_FREQ_SCALE = {"HZ": 1.0, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9}


def read_touchstone(path):
    """Return (freq_hz, S, z0) with S shaped (n_freq, n_ports, n_ports)."""
    path = Path(path)
    m = re.search(r"\.s(\d+)p$", path.name, re.IGNORECASE)
    if not m:
        raise ValueError("Cannot infer port count from file name %s" % path.name)
    n = int(m.group(1))

    unit, fmt, z0 = "GHZ", "MA", 50.0
    numbers = []
    for raw in path.read_text().splitlines():
        line = raw.split("!", 1)[0].strip()
        if not line:
            continue
        if line.startswith("#"):
            tok = line[1:].upper().split()
            for i, tkn in enumerate(tok):
                if tkn in _FREQ_SCALE:
                    unit = tkn
                elif tkn in ("MA", "DB", "RI"):
                    fmt = tkn
                elif tkn in ("Y", "Z", "H", "G"):
                    raise ValueError("Only S-parameter files are supported")
                elif tkn == "R" and i + 1 < len(tok):
                    z0 = float(tok[i + 1])
            continue
        if line.startswith("["):  # Touchstone v2 keywords: not supported
            raise ValueError("Touchstone v2 keyword found; export as v1")
        numbers.extend(float(x) for x in line.split())

    per_freq = 1 + 2 * n * n
    if len(numbers) % per_freq:
        raise ValueError("Data length %d is not a multiple of %d" % (len(numbers), per_freq))
    arr = np.array(numbers).reshape(-1, per_freq)
    freq = arr[:, 0] * _FREQ_SCALE[unit]
    p1, p2 = arr[:, 1::2], arr[:, 2::2]

    if fmt == "MA":
        vals = p1 * np.exp(1j * np.deg2rad(p2))
    elif fmt == "DB":
        vals = 10 ** (p1 / 20.0) * np.exp(1j * np.deg2rad(p2))
    else:
        vals = p1 + 1j * p2

    S = vals.reshape(-1, n, n)
    if n == 2:  # v1 two-port order is S11 S21 S12 S22 (column-major)
        S = np.transpose(S, (0, 2, 1))
    return freq, S, z0


def db(x):
    return 20.0 * np.log10(np.maximum(np.abs(x), 1e-15))


def coupler_metrics(freq_hz, S, input_port=1, through=2, coupled=3, isolated=4):
    """Standard four-port coupler figures, all positive dB numbers."""
    i, t, c, iso = input_port - 1, through - 1, coupled - 1, isolated - 1
    coupling = -db(S[:, c, i])
    isolation = -db(S[:, iso, i])
    return {
        "freq_ghz": freq_hz / 1e9,
        "return_loss_db": -db(S[:, i, i]),
        "insertion_loss_db": -db(S[:, t, i]),
        "coupling_db": coupling,
        "isolation_db": isolation,
        "directivity_db": isolation - coupling,
    }


def summarize(metrics, f0_ghz):
    f = metrics["freq_ghz"]
    k = int(np.argmin(np.abs(f - f0_ghz)))
    out = {"f0_ghz": float(f[k])}
    for key in ("return_loss_db", "insertion_loss_db", "coupling_db", "isolation_db", "directivity_db"):
        out[key + "_at_f0"] = float(metrics[key][k])
    out["coupling_flatness_db"] = float(np.ptp(metrics["coupling_db"]))
    out["min_directivity_db"] = float(np.min(metrics["directivity_db"]))
    out["min_return_loss_db"] = float(np.min(metrics["return_loss_db"]))
    return out
