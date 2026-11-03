"""Gather an openEMS study campaign into results/openems/ (run once after the campaign).

    python docs/report/collect.py <campaign_dir> [<legacy_campaign_dir>]

Each campaign sub-folder is the --out folder of one simulate_openems.py run
(run_info.json, summary.json, coupler_metrics.csv, moreno_openems.s4p).
The optional legacy campaign holds runs of the superseded model (air carved out
of a metal block); only its mesh series is kept, for the comparison in the report.
"""
from __future__ import annotations

import csv
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
import os
OUT = Path(os.environ.get("MORENO_RESULTS", ROOT / "results" / "openems"))

GROUPS = {
    # slot/wall cell h refined, 1/3-2/3 edge rule, L = 2.00 mm
    "convergence": ["c050", "c035", "c025", "c0175", "c0125"],
    # same, mesh lines on the slot edges instead of the 1/3-2/3 rule
    "staircase": ["s035", "s025", "s0175"],
    # other numerical parameters, each changed alone from the h = 0.035 mm run
    "numerics": ["c035", "n_res30", "n_feed25", "n_e4"],
    # design chart at h = 0.035 mm
    "slot_length": ["L1.6", "L1.8", "c035", "L2.2", "L2.4", "L2.6"],
    # the same slot lengths at h = 0.025 mm, to check that the mesh offset does not depend on L
    "offset_check": ["o025_L1.6", "c025", "o025_L2.4"],
    # manufacturing tolerances at h = 0.035 mm
    "tolerance": ["c035", "t020", "t030", "w030", "w040"],
}
LEGACY = ["full_r20", "f060", "f045", "f035", "f025", "f0175"]
REFERENCE = "c0125"
FULL = "full"
FIELDS = "fields"

KEYS = ["coupling_db_at_f0", "directivity_db_at_f0", "isolation_db_at_f0", "return_loss_db_at_f0",
        "insertion_loss_db_at_f0", "coupling_flatness_db", "min_directivity_db", "min_return_loss_db"]


def load(camp: Path, name: str):
    p = camp / name / "run_info.json"
    return json.loads(p.read_text()) if p.exists() else None


def row(name, info):
    run = info["runs"][0]
    d = info["design"]
    out = {"run": name, "slot_cell_mm": round(run["fine_mm"], 5), "cells": run["cells"],
           "runtime_s": run["runtime_s"], "end_criteria": run["end_criteria"],
           "cells_per_lambda": run["res_cells_per_lambda"], "feed_lambda_g": run["feed_lambda_g"],
           "edge_rule": "thirds" if run.get("thirds_rule") else "lines on edges",
           "L_mm": d["slot_length_mm"], "w_mm": d["slot_width_mm"], "t_mm": d["wall_mm"],
           "power_balance_min": info["power_balance_port1"]["min"],
           "power_balance_max": info["power_balance_port1"]["max"]}
    out.update({k: info["summary"][k] for k in KEYS})
    return out


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def copy_run(camp, name, dest):
    dest.mkdir(parents=True, exist_ok=True)
    for f in ("moreno_openems.s4p", "coupler_metrics.csv", "summary.json", "run_info.json", "sparams.png",
              "fields_f0.npz"):
        src = camp / name / f
        if src.exists():
            shutil.copy(src, dest / f)


def main(camp, legacy=None):
    camp = Path(camp)
    studies = OUT / "studies"
    if studies.exists():
        shutil.rmtree(studies)
    studies.mkdir(parents=True)
    for group, names in GROUPS.items():
        rows = [row(n, i) for n in names if (i := load(camp, n))]
        write_csv(studies / ("%s.csv" % group), rows)
        print("%-13s %d/%d runs" % (group, len(rows), len(names)))
        for n in names:
            if (camp / n / "coupler_metrics.csv").exists():
                shutil.copy(camp / n / "coupler_metrics.csv", studies / ("%s__%s.csv" % (group, n)))
    # the raw result of every run (Touchstone file and settings), for tools/verify_results.py
    every = sorted({n for names in GROUPS.values() for n in names} | {REFERENCE, FULL, FIELDS})
    for n in every:
        if load(camp, n):
            dest = studies / "runs" / n
            dest.mkdir(parents=True, exist_ok=True)
            for f in ("moreno_openems.s4p", "run_info.json"):
                shutil.copy(camp / n / f, dest / f)
    print("runs          %d raw results kept in studies/runs/" % len(list((studies / "runs").iterdir())))
    if legacy:
        legacy = Path(legacy)
        rows = [row(n, i) for n in LEGACY if (i := load(legacy, n))]
        for r in rows:
            r["edge_rule"] = "lines on edges (air-priority model)"
        write_csv(studies / "legacy_air_priority_model.csv", rows)
        print("legacy        %d runs" % len(rows))

    for stale in ("moreno_openems.s4p", "coupler_metrics.csv", "summary.json", "sparams.png",
                  "slot_sweep.csv", "slot_sweep.png", "run_info.json"):
        (OUT / stale).unlink(missing_ok=True)
    for name, dest in ((REFERENCE, "reference"), (FULL, "full_matrix"), (FIELDS, "fields")):
        if (OUT / dest).exists():
            shutil.rmtree(OUT / dest)
        if load(camp, name):
            copy_run(camp, name, OUT / dest)
            (OUT / dest / "SOURCE.txt").write_text("campaign run '%s'\n" % name)
            print("%-13s <- %s" % (dest, name))


if __name__ == "__main__":
    main(*sys.argv[1:3])
