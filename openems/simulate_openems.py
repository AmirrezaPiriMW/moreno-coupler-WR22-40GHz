"""Full-wave FDTD simulation of the Moreno coupler with openEMS (free, no licence).

Builds the same geometry as the HFSS builders from design/design.json,
excites port 1 with the TE10 mode, and writes

    results/openems/moreno_openems.s4p   Touchstone (see note on columns)
    results/openems/sparams.png          plot + summary.json + coupler_metrics.csv

    python openems/simulate_openems.py                 # port 1 excitation (~4 min, 2 cores)
    python openems/simulate_openems.py --excite all    # full 4x4 S-matrix (4 runs)
    python openems/simulate_openems.py --fine 0.025    # finer slot mesh (convergence check)
    python openems/simulate_openems.py --staircase     # edges on mesh lines (for comparison)
    python openems/simulate_openems.py --view          # open geometry in AppCSXCAD

Install openEMS and its Python bindings first: see openems/README.md.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "post"))

from moreno.design import MorenoInputs, compute  # noqa: E402

C0 = 299792458.0

# Windows: let Python find the openEMS DLLs (set OPENEMS_INSTALL_PATH, e.g. C:\\opt\\openEMS)
if os.name == "nt" and os.environ.get("OPENEMS_INSTALL_PATH"):
    os.add_dll_directory(os.environ["OPENEMS_INSTALL_PATH"])


def load_design(path):
    p = Path(path) if path else ROOT / "design" / "design.json"
    if p.exists():
        return json.loads(p.read_text())
    return compute(MorenoInputs())


# ---------------------------------------------------------------------------
# Geometry and mesh (all lengths in mm)
# ---------------------------------------------------------------------------
class Geometry:
    def __init__(self, r, feed_lg):
        self.a, self.b, self.t = r["a_mm"], r["b_mm"], r["wall_mm"]
        self.L, self.w = r["slot_length_mm"], r["slot_width_mm"]
        self.d = r["d_mm"]
        self.lamg = r["lambda_g_mm"]
        self.feed = feed_lg * self.lamg        # guide length beyond the shared square
        self.y2 = self.b + self.t               # floor of guide 2
        self.top = 2 * self.b + self.t
        self.centres = [(self.d, self.d), (self.a - self.d, self.a - self.d)]

    def slot_arms(self):
        """The four slot arms as (x0, x1, z0, z1) rectangles in the common wall."""
        L, w = self.L, self.w
        arms = []
        for cx, cz in self.centres:
            arms.append((cx - L / 2, cx + L / 2, cz - w / 2, cz + w / 2))  # arm along x
            arms.append((cx - w / 2, cx + w / 2, cz - L / 2, cz + L / 2))  # arm along z
        return arms

    def metal_boxes(self, lo, hi):
        """All metal as boxes: the complement of the two guides and the four slot arms.

        The model is built from metal only, on a vacuum background, on purpose.
        openEMS makes an E-field edge PEC only when the highest-priority primitive
        at the edge's midpoint is a metal (Operator::CalcPEC_Range). If air boxes
        are carved out of a metal block with higher priority, the air also owns
        every edge lying ON a metal surface, so each conductor surface moves one
        cell into the metal: slots get wider, the wall thinner and the guides
        wider, by amounts that change with the mesh. With metal-only geometry the
        metal owns its surfaces and every wall sits exactly where it is drawn.
        """
        a, b, y2, top = self.a, self.b, self.y2, self.top
        boxes = [
            ([lo, 0, lo], [0, b, hi]), ([a, 0, lo], [hi, b, hi]),          # beside guide 1
            ([lo, y2, lo], [hi, top, 0]), ([lo, y2, a], [hi, top, hi]),    # beside guide 2
        ]
        # common wall (b <= y <= y2): grid it on every slot edge, keep the cells
        # that are not inside an arm, and merge them into runs along z
        arms = self.slot_arms()
        xs = sorted({lo, hi} | {v for r in arms for v in r[:2]})
        zs = sorted({lo, hi} | {v for r in arms for v in r[2:]})

        def is_air(xm, zm):
            return any(x0 < xm < x1 and z0 < zm < z1 for x0, x1, z0, z1 in arms)

        for x0, x1 in zip(xs[:-1], xs[1:]):
            run = None
            for z0, z1 in zip(zs[:-1], zs[1:]):
                if is_air(0.5 * (x0 + x1), 0.5 * (z0 + z1)):
                    if run:
                        boxes.append(([x0, b, run], [x1, y2, z0]))
                        run = None
                elif run is None:
                    run = z0
            if run is not None:
                boxes.append(([x0, b, run], [x1, y2, zs[-1]]))
        return boxes


def _clean(lines, tol):
    """Sort, round and drop lines closer than tol (avoids near-zero cells)."""
    out = []
    for v in sorted(np.round(np.asarray(lines, float), 6)):
        if not out or v - out[-1] > tol:
            out.append(float(v))
    return np.array(out)


def build_mesh(mesh, g, f_max_ghz, res, port_planes=(), fine_max=None, thirds=False):
    """Graded mesh: coarse in the feeds, fine around the slots and the wall.

    thirds=True applies the 1/3-2/3 rule to the slot edges in x and z: instead of a
    mesh line on each metal edge, lines are placed h/3 inside the metal and 2h/3
    inside the slot. This compensates the field singularity at thin conducting
    edges, which otherwise makes FDTD converge slowly (see docs/report).
    """
    lam_min = C0 / (f_max_ghz * 1e9) * 1e3
    coarse = lam_min / res
    fine = min(g.w / 3.0, g.t / 3.0, coarse / 3.0)
    if fine_max:
        fine = min(fine, fine_max)

    lo, hi = round(-g.feed, 4), round(g.a + g.feed, 4)
    slot_x = sorted({c for cx, _ in g.centres for c in (cx - g.L / 2, cx - g.w / 2, cx + g.w / 2, cx + g.L / 2)})

    x = [lo, 0.0, g.a, hi] + slot_x + list(port_planes)
    y = [0.0, g.b, g.top] + list(np.linspace(g.b, g.y2, max(4, int(np.ceil(g.t / fine)) + 1)))
    if thirds:
        anchors = []
        for cx, _ in g.centres:
            for e in (cx - g.L / 2, cx - g.w / 2):      # metal on the low side
                anchors += [e - fine / 3, e + 2 * fine / 3]
            for e in (cx + g.w / 2, cx + g.L / 2):      # metal on the high side
                anchors += [e - 2 * fine / 3, e + fine / 3]
        anchors.sort()
        x = [lo, 0.0, g.a, hi] + list(port_planes)
        for s0, s1 in zip(anchors[:-1], anchors[1:]):
            n = max(1, int(np.ceil((s1 - s0) / fine - 1e-9)))
            x += list(np.linspace(s0, s1, n + 1))
    else:
        # fine lines across each slot arm, with every slot edge on a mesh line
        for s0, s1 in zip(slot_x[:-1], slot_x[1:]):
            n = max(2, int(np.ceil((s1 - s0) / fine)))
            x += list(np.linspace(s0, s1, n + 1))
    z = list(x)  # the slot pattern is symmetric under x <-> z

    for ny, lines, res_ny in (("x", x, coarse), ("z", z, coarse), ("y", y, min(coarse, g.b / 8.0))):
        mesh.SetLines(ny, _clean(lines, fine / 4))
        mesh.SmoothMeshLines(ny, res_ny, 1.3)
        mesh.SetLines(ny, _clean(mesh.GetLines(ny), fine / 4))
    return coarse, fine


# ---------------------------------------------------------------------------
# One FDTD run with one port excited
# ---------------------------------------------------------------------------
def run_once(r, excite_port, args, sim_path):
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS

    g = Geometry(r, args.feed)
    f0 = r["f0_ghz"] * 1e9
    f_lo, f_hi = (v * 1e9 for v in r["sweep_ghz"])
    fc = 0.5 * (f_hi - f_lo) * 1.2  # a little wider than the sweep

    FDTD = openEMS(NrTS=args.max_steps, EndCriteria=args.end_criteria)
    FDTD.SetGaussExcite(f0, fc)
    # PEC on the y walls; PML where the four guides leave the box (x and z ends)
    FDTD.SetBoundaryCond(["PML_8", "PML_8", "PEC", "PEC", "PML_8", "PML_8"])

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    lo, hi = round(-g.feed, 4), round(g.a + g.feed, 4)
    coarse = C0 / f_hi * 1e3 / args.res
    # ports sit 10 cells in from the PML; start->stop points into the structure.
    # Their planes must lie exactly on mesh lines (same planes in x and z).
    off, span = round(10 * coarse, 4), round(4 * coarse, 4)
    planes = [lo + off, lo + off + span, hi - off, hi - off - span]
    coarse, fine = build_mesh(mesh, g, f_hi / 1e9, args.res, planes, getattr(args, "fine", None),
                              getattr(args, "thirds", False))

    # PEC walls on a vacuum background (see Geometry.metal_boxes for why the air
    # is not drawn). y = 0 and y = top are closed by the PEC boundary condition.
    metal = CSX.AddMetal("pec_walls")
    for start, stop in g.metal_boxes(lo, hi):
        metal.AddBox(start, stop, priority=10)

    a_m, b_m = g.a * 1e-3, g.b * 1e-3
    ports = [
        # guide 1: TE10 with 'a' along x (first transverse axis for z-propagation)
        FDTD.AddRectWaveGuidePort(0, [0, 0, lo + off], [g.a, g.b, lo + off + span], "z",
                                  a_m, b_m, "TE10", excite_port == 1, priority=100),
        FDTD.AddRectWaveGuidePort(1, [0, 0, hi - off], [g.a, g.b, hi - off - span], "z",
                                  a_m, b_m, "TE10", excite_port == 2, priority=100),
        # guide 2: transverse axes are (y, z); the broad wall runs along z -> TE01
        FDTD.AddRectWaveGuidePort(2, [lo + off, g.y2, 0], [lo + off + span, g.top, g.a], "x",
                                  b_m, a_m, "TE01", excite_port == 3, priority=100),
        FDTD.AddRectWaveGuidePort(3, [hi - off, g.y2, 0], [hi - off - span, g.top, g.a], "x",
                                  b_m, a_m, "TE01", excite_port == 4, priority=100),
    ]

    dump_planes = []
    if getattr(args, "dump_fields", False):
        # E-field phasor at f0 on the mid-height plane of each guide (HDF5)
        for name, yv in (("Ef_guide1", 0.5 * g.b), ("Ef_guide2", g.y2 + 0.5 * g.b)):
            dmp = CSX.AddDump(name, dump_type=10, file_type=1, frequency=[f0], dump_mode=1)
            dmp.AddBox([lo, yv, lo], [hi, yv, hi])
            dump_planes.append(name)

    if args.view:
        xml = os.path.join(sim_path, "geometry.xml")
        os.makedirs(sim_path, exist_ok=True)
        CSX.Write2XML(xml)
        viewer = shutil.which("AppCSXCAD") or shutil.which("AppCSXCAD.sh")
        if viewer:
            os.system('"%s" "%s"' % (viewer, xml))
        else:
            print("Geometry written to %s (open it with AppCSXCAD)" % xml)
        sys.exit(0)

    n = [len(mesh.GetLines(ax)) for ax in "xyz"]
    print("Port %d excited | mesh %d x %d x %d = %.2f M cells | coarse %.3f mm, fine %.3f mm"
          % (excite_port, n[0], n[1], n[2], np.prod(n) / 1e6, coarse, fine))

    os.makedirs(sim_path, exist_ok=True)
    t0 = time.time()
    FDTD.Run(sim_path, cleanup=True, verbose=args.verbose, numThreads=args.threads)
    print("  solved in %.0f s" % (time.time() - t0))

    freq = np.linspace(f_lo, f_hi, args.points)
    # De-embed: move every reference plane from its port to the edge of the
    # shared a x a square (z = 0, z = a, x = 0, x = a). ref_plane_shift is the
    # distance from each port's start plane, in drawing units.
    edges = [abs(0.0 - (lo + off)), abs((hi - off) - g.a), abs(0.0 - (lo + off)), abs((hi - off) - g.a)]
    for p, shift in zip(ports, edges):
        p.CalcPort(sim_path, freq, ref_plane_shift=shift)
    inc = ports[excite_port - 1].uf_inc
    column = np.array([p.uf_ref / inc for p in ports])  # S_k,exc
    if dump_planes:
        save_fields(sim_path, dump_planes, f0, g, args.out)
    info = {"excited_port": excite_port, "mesh_lines_xyz": n, "cells": int(np.prod(n)),
            "coarse_mm": coarse, "fine_mm": fine, "runtime_s": round(time.time() - t0, 1),
            "end_criteria": args.end_criteria, "res_cells_per_lambda": args.res,
            "feed_lambda_g": args.feed, "y_lines_in_wall": int(round(g.t / fine)),
            "thirds_rule": bool(getattr(args, "thirds", False))}
    return freq, column, info


def save_fields(sim_path, names, f0, g, out_dir):
    """Store the dumped E-field planes (complex, at f0) as one compact .npz."""
    from openEMS.utilities import HDF5Dump

    data = {"a_mm": g.a, "b_mm": g.b, "y2_mm": g.y2, "f0_hz": f0}
    for name in names:
        with HDF5Dump(os.path.join(sim_path, name + ".h5")) as h5:
            mesh = h5.GetMesh(region=True)          # lines in metres
            E = h5.GetFieldAtFrequency(f0)          # (3, Nx, Ny, Nz)
        data[name + "_x_mm"] = np.asarray(mesh["lines"][0]) * 1e3
        data[name + "_z_mm"] = np.asarray(mesh["lines"][2]) * 1e3
        data[name + "_Ey"] = np.squeeze(E[1]).astype(np.complex64)   # (Nx, Nz)
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    np.savez_compressed(Path(out_dir) / "fields_f0.npz", **data)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def write_s4p(path, freq, S, note):
    lines = ["! Moreno cross-guide coupler, openEMS FDTD", "! " + note,
             "! S-parameters are normalised to each port's TE10 wave impedance (generalised",
             "! S-matrix); 'R 50' below is only the placeholder the format requires.",
             "! Reference planes: edges of the common a x a square.",
             "# GHz S MA R 50"]
    for k, f in enumerate(freq):
        vals = []
        for i in range(4):
            for j in range(4):
                vals += [abs(S[k, i, j]), np.degrees(np.angle(S[k, i, j]))]
        rows = [vals[m:m + 8] for m in range(0, 32, 8)]
        lines.append("%.6f " % (f / 1e9) + " ".join("%.6e" % v for v in rows[0]))
        for row in rows[1:]:
            lines.append("  " + " ".join("%.6e" % v for v in row))
    Path(path).write_text("\n".join(lines) + "\n")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--design", help="design JSON (default design/design.json)")
    p.add_argument("--excite", default="1", help="'1' (default) or 'all' for the full S-matrix")
    p.add_argument("--res", type=float, default=20, help="cells per wavelength at f_max (default 20)")
    p.add_argument("--feed", type=float, default=1.5, help="feed length in guide wavelengths")
    p.add_argument("--fine", type=float, default=0.035,
                   help="max cell size across slots and wall, mm (default 0.035 = w/10; see docs/report "
                        "for the convergence study - coarser slot meshes make the coupling too strong)")
    p.add_argument("--thirds", dest="thirds", action="store_true", default=True,
                   help="1/3-2/3 mesh rule at the slot edges (default; much faster convergence)")
    p.add_argument("--staircase", dest="thirds", action="store_false",
                   help="put mesh lines on the slot edges instead (slow convergence; for comparison)")
    p.add_argument("--points", type=int, default=201)
    p.add_argument("--max-steps", type=int, default=200000)
    p.add_argument("--end-criteria", type=float, default=1e-5, help="energy decay to stop (1e-5 = -50 dB)")
    p.add_argument("--threads", type=int, default=0, help="0 = all cores")
    p.add_argument("--out", default=str(ROOT / "results" / "openems"))
    p.add_argument("--sim-dir", default=None, help="scratch folder for field data")
    p.add_argument("--view", action="store_true", help="only write geometry and open AppCSXCAD")
    p.add_argument("--dump-fields", action="store_true",
                   help="also save the E-field at f0 on the mid-planes of both guides (fields_f0.npz)")
    p.add_argument("--verbose", type=int, default=0)
    args = p.parse_args(argv)

    r = load_design(args.design)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    sim_root = Path(args.sim_dir) if args.sim_dir else out / "sim"

    excites = [1, 2, 3, 4] if args.excite == "all" else [int(args.excite)]
    S = None
    runs = []
    for e in excites:
        freq, col, info = run_once(r, e, args, str(sim_root / ("port%d" % e)))
        runs.append(info)
        if S is None:
            S = np.zeros((len(freq), 4, 4), complex)
        S[:, :, e - 1] = col.T

    if len(excites) == 4:
        note = "full 4x4 matrix, one run per port"
    else:
        e = excites[0] - 1
        S[:, e, :] = S[:, :, e]  # reciprocity fills the matching row
        note = ("only column %d (port %d excited) and, by reciprocity, row %d are simulated; "
                "other entries are zero" % (e + 1, e + 1, e + 1))
    s4p = out / "moreno_openems.s4p"
    write_s4p(s4p, freq, S, note)
    print("Wrote %s" % s4p)

    from analyze_touchstone import analyze

    summary = analyze(str(s4p), f0=r["f0_ghz"], out=str(out))
    # lossless (PEC) model: power leaving the excited port must sum to 1
    e = excites[0] - 1
    balance = np.sum(np.abs(S[:, :, e]) ** 2, axis=1)
    meta = {"design": {k: r[k] for k in ("a_mm", "b_mm", "wall_mm", "slot_length_mm", "slot_width_mm",
                                         "spacing_mm", "d_mm", "f0_ghz")},
            "runs": runs, "reference_planes": "edges of the shared a x a square",
            "power_balance_port%d" % (e + 1): {"min": float(balance.min()), "max": float(balance.max())},
            "summary": summary}
    if len(excites) == 4:
        meta["reciprocity_max_abs_Sij_minus_Sji"] = float(np.max(np.abs(S - np.transpose(S, (0, 2, 1)))))
        meta["power_balance_all_ports"] = {"min": float(np.min(np.sum(np.abs(S) ** 2, axis=1))),
                                           "max": float(np.max(np.sum(np.abs(S) ** 2, axis=1)))}
    (out / "run_info.json").write_text(json.dumps(meta, indent=2))
    print("Power balance (sum |S_k%d|^2): %.5f .. %.5f" % (e + 1, balance.min(), balance.max()))


if __name__ == "__main__":
    main()
