"""Sweep the slot arm length in openEMS to set the coupling level.

    python openems/sweep_slot_length.py --L 1.6 2.0 2.4
    python openems/sweep_slot_length.py --L 1.6 2.0 2.4 --target 20

Writes results/openems/sweep/slot_sweep.csv and slot_sweep.png. With --target it
interpolates the slot length that gives that coupling (dB) at f0.

The runs use the h = 0.035 mm slot mesh for speed. That mesh couples slightly
more strongly than the mesh-converged limit; the shift measured in the report's
convergence study (docs/report/key_numbers.json, "chart_mesh_shift_db") is added
to every value unless --no-shift is given.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import simulate_openems as sim  # noqa: E402
from moreno.design import MorenoInputs, compute  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--L", type=float, nargs="+", required=True, help="slot arm lengths, mm")
    p.add_argument("--target", type=float, help="target coupling at f0, dB (positive)")
    p.add_argument("--res", type=float, default=20)
    p.add_argument("--feed", type=float, default=1.5)
    p.add_argument("--fine", type=float, default=0.035, help="max cell size across slots and wall, mm")
    p.add_argument("--threads", type=int, default=0)
    p.add_argument("--no-shift", action="store_true", help="report raw h = 0.035 mm values")
    p.add_argument("--out", default=str(HERE.parent / "results" / "openems" / "sweep"))
    a = p.parse_args(argv)

    base = sim.load_design(None)["inputs"]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    run_args = sim.argparse.Namespace(res=a.res, feed=a.feed, fine=a.fine, thirds=True, points=201, max_steps=150000,
                                      end_criteria=1e-5, threads=a.threads, view=False, verbose=0,
                                      dump_fields=False, out=str(out))
    shift = 0.0
    key = HERE.parent / "docs" / "report" / "key_numbers.json"
    if not a.no_shift and key.exists() and abs(a.fine - 0.035) < 1e-9:
        import json

        shift = json.loads(key.read_text()).get("chart_mesh_shift_db", 0.0)
        print("Adding the mesh shift of %+.2f dB from the convergence study (h = 0.035 mm -> 0)" % shift)
    rows = []
    for L in a.L:
        inp = MorenoInputs(**{**base, "slot_length_mm": L})
        r = compute(inp)
        freq, col, _ = sim.run_once(r, 1, run_args, str(out / "sim" / ("L_%.3f" % L)))
        k = int(np.argmin(np.abs(freq - r["f0_ghz"] * 1e9)))
        c = -20 * np.log10(abs(col[2][k])) + shift
        d = 20 * np.log10(abs(col[2][k]) / abs(col[3][k]))
        rl = -20 * np.log10(abs(col[0][k]))
        rows.append((L, c, d, rl))
        print("L = %.3f mm: coupling %.2f dB, directivity %.1f dB, return loss %.1f dB" % (L, c, d, rl))

    with open(out / "slot_sweep.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["slot_length_mm", "coupling_db", "directivity_db", "return_loss_db"])
        w.writerows([["%.4f" % v for v in row] for row in rows])

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Ls, cs = [r_[0] for r_ in rows], [r_[1] for r_ in rows]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(Ls, cs, "o-", c="#2A5788")
    ax.set_xlabel("Slot arm length L (mm)")
    ax.set_ylabel("Coupling at f0 (dB)")
    ax.invert_yaxis()
    ax.grid(alpha=0.3)
    if a.target and len(rows) > 1:
        order = np.argsort(cs)
        L_t = float(np.interp(a.target, np.array(cs)[order], np.array(Ls)[order]))
        ax.axhline(a.target, ls="--", c="#A8461A")
        ax.axvline(L_t, ls="--", c="#A8461A")
        print("\nSlot length for %.1f dB coupling: about %.3f mm" % (a.target, L_t))
        print("Set it with: python -m moreno.design --L %.3f" % L_t)
    fig.tight_layout()
    fig.savefig(out / "slot_sweep.png", dpi=150)
    print("Wrote %s" % out)


if __name__ == "__main__":
    main()
