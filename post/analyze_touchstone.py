"""Analyse an HFSS-exported .s4p of the Moreno coupler.

Usage:
    python post/analyze_touchstone.py results/moreno.s4p --f0 40 --out results/
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from moreno.touchstone import coupler_metrics, read_touchstone, summarize  # noqa: E402


def plot(metrics, f0, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    f = metrics["freq_ghz"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    ax1.plot(f, -metrics["return_loss_db"], label="S11 (match)")
    ax1.plot(f, -metrics["insertion_loss_db"], label="S21 (through)")
    ax1.plot(f, -metrics["coupling_db"], label="S31 (coupled)", lw=2)
    ax1.plot(f, -metrics["isolation_db"], label="S41 (isolated)")
    ax1.set_ylabel("dB")
    ax1.legend(loc="lower right", fontsize=8)
    ax1.grid(alpha=0.3)
    ax1.axvline(f0, ls="--", c="grey", lw=1)

    ax2.plot(f, metrics["directivity_db"], c="#A8461A", lw=2)
    ax2.set_ylabel("Directivity (dB)")
    ax2.set_xlabel("Frequency (GHz)")
    ax2.grid(alpha=0.3)
    ax2.axvline(f0, ls="--", c="grey", lw=1)
    fig.suptitle("Moreno cross-guide coupler")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def analyze(touchstone, f0=40.0, out="results", ports=(1, 2, 3, 4), quiet=False):
    freq, S, _ = read_touchstone(touchstone)
    m = coupler_metrics(freq, S, *ports)
    s = summarize(m, f0)

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "coupler_metrics.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        keys = list(m)
        w.writerow(keys)
        for row in zip(*(m[k] for k in keys)):
            w.writerow(["%.6g" % v for v in row])
    (out / "summary.json").write_text(json.dumps(s, indent=2))
    plot(m, f0, out / "sparams.png")

    if not quiet:
        print("At %.3f GHz:" % s["f0_ghz"])
        print("  coupling     %6.2f dB" % s["coupling_db_at_f0"])
        print("  directivity  %6.2f dB" % s["directivity_db_at_f0"])
        print("  isolation    %6.2f dB" % s["isolation_db_at_f0"])
        print("  return loss  %6.2f dB" % s["return_loss_db_at_f0"])
        print("  through loss %6.3f dB" % s["insertion_loss_db_at_f0"])
        print("Across sweep: coupling flatness %.2f dB, min directivity %.2f dB"
              % (s["coupling_flatness_db"], s["min_directivity_db"]))
        if s["directivity_db_at_f0"] < 0:
            print("\nNOTE: port 4 is stronger than port 3. Swap the labels (--ports 1,2,4,3) "
                  "or mirror the slots onto the other diagonal.")
        print("Wrote %s" % out)
    return s


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("touchstone")
    p.add_argument("--f0", type=float, default=40.0)
    p.add_argument("--out", default="results")
    p.add_argument("--ports", default="1,2,3,4", help="input,through,coupled,isolated")
    args = p.parse_args(argv)
    analyze(args.touchstone, args.f0, args.out, tuple(int(x) for x in args.ports.split(",")))


if __name__ == "__main__":
    main()
