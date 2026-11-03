"""Build, solve and export the Moreno coupler in HFSS through PyAEDT.

Runs from a normal Python 3.10+ prompt (not inside AEDT). PyAEDT launches
or attaches to AEDT and picks the right solution-type name for the release,
which avoids the "unable to insert design" problems of raw scripts.

    pip install -r requirements-hfss.txt
    python hfss/build_pyaedt.py --version 2026.1 --student            # build only
    python hfss/build_pyaedt.py --version 2026.1 --student --analyze  # build + solve + post

Geometry, ports and variables match moreno/design.py exactly.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from moreno.design import MorenoInputs, compute  # noqa: E402


def load_design(path: str | None) -> dict:
    if path:
        return json.loads(Path(path).read_text())
    default = ROOT / "design" / "design.json"
    if default.exists():
        return json.loads(default.read_text())
    return compute(MorenoInputs())


def import_hfss():
    try:
        from ansys.aedt.core import Hfss  # PyAEDT >= 0.9
    except ImportError:
        try:
            from pyaedt import Hfss  # older PyAEDT
        except ImportError:
            sys.exit("PyAEDT is not installed: pip install -r requirements-hfss.txt")
    return Hfss


def build(hfss, r: dict) -> None:
    mm = lambda v: "%.4fmm" % v  # noqa: E731
    variables = {
        "a": mm(r["a_mm"]),
        "b": mm(r["b_mm"]),
        "t": mm(r["wall_mm"]),
        "Ls": mm(r["slot_length_mm"]),
        "ws": mm(r["slot_width_mm"]),
        "Delta": mm(r["spacing_mm"]),
        "d": "(a-Delta)/2",
        "Lext": mm(r["feed_mm"]),
    }
    for k, v in variables.items():
        hfss[k] = v

    hfss.modeler.model_units = "mm"
    m = hfss.modeler
    vac = "vacuum"
    m.create_box(["0mm", "0mm", "-Lext"], ["a", "b", "a+2*Lext"], name="Air", material=vac)
    m.create_box(["-Lext", "b+t", "0mm"], ["a+2*Lext", "b", "a"], name="G2", material=vac)
    # slot 1 at (d, d), slot 2 at (a-d, a-d); "+" shape = two crossed arms each
    m.create_box(["d-Ls/2", "b", "d-ws/2"], ["Ls", "t", "ws"], name="S1x", material=vac)
    m.create_box(["d-ws/2", "b", "d-Ls/2"], ["ws", "t", "Ls"], name="S1z", material=vac)
    m.create_box(["a-d-Ls/2", "b", "a-d-ws/2"], ["Ls", "t", "ws"], name="S2x", material=vac)
    m.create_box(["a-d-ws/2", "b", "a-d-Ls/2"], ["ws", "t", "Ls"], name="S2z", material=vac)
    # Each slot's two arms overlap, so they are united into one solid per slot. The
    # guides and slots only touch (no overlap), so they stay separate solids; this
    # lets a mesh operation target the slots alone.
    m.unite(["S1x", "S1z"])
    m.unite(["S2x", "S2z"])
    # Seed the slots at w/5; adaptive refinement does the rest. (A much finer forced
    # mesh adds elements without improving the converged answer.)
    hfss.mesh.assign_length_mesh(["S1x", "S2x"], inside_selection=True,
                                 maximum_length="%.4fmm" % (r["slot_width_mm"] / 5.0),
                                 maximum_elements=None, name="SlotMesh")

    g1, g2 = m["Air"], m["G2"]
    # every face touching the background is PEC by default -> copper walls
    port_faces = {
        "1": g1.bottom_face_z,  # input,    z = -Lext
        "2": g1.top_face_z,     # through,  z = a + Lext
        "3": g2.bottom_face_x,  # coupled,  x = -Lext
        "4": g2.top_face_x,     # isolated, x = a + Lext
    }
    for name, face in port_faces.items():
        # renormalize=False keeps the generalized S-matrix referenced to each port's own
        # TE10 wave impedance (Zpi). Renormalizing a waveguide port to 50 ohm would
        # distort every |S|.
        hfss.wave_port(face.id, modes=1, name=name, renormalize=False)

    setup = hfss.create_setup("Setup1")
    setup.props["Frequency"] = "%.4fGHz" % r["f0_ghz"]
    # |S31| is only about 0.03 (-30 dB), so the usual Delta-S of 0.02 would stop the
    # adaptive loop long before the coupling is resolved. 0.002 resolves |S31| to a
    # few percent (~0.5 dB). Two converged passes guard against a lucky pass.
    setup.props["MaximumPasses"] = 20
    setup.props["MinimumPasses"] = 3
    setup.props["MinimumConvergedPasses"] = 2
    setup.props["MaxDeltaS"] = 0.002
    setup.update()

    lo, hi = r["sweep_ghz"]
    hfss.create_linear_count_sweep(
        setup="Setup1", unit="GHz", start_frequency=lo, stop_frequency=hi,
        num_of_freq_points=101, name="Sweep", save_fields=False, sweep_type="Interpolating",
    )

    try:
        hfss.post.create_report(
            expressions=["dB(S(1,1))", "dB(S(2,1))", "dB(S(3,1))", "dB(S(4,1))"],
            setup_sweep_name="Setup1 : Sweep", primary_sweep_variable="Freq", plot_name="S-parameters",
        )
    except Exception as exc:  # report API differs between releases
        print("Skipping S-parameter report: %s" % exc)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--design", help="design JSON (default: design/design.json)")
    p.add_argument("--version", default="2026.1", help="AEDT version, e.g. 2025.2 or 2026.1")
    p.add_argument("--student", action="store_true", help="use AEDT Student")
    p.add_argument("--non-graphical", action="store_true")
    p.add_argument("--analyze", action="store_true", help="solve, export .s4p and post-process")
    p.add_argument("--cores", type=int, default=4)
    p.add_argument("--out", default=str(ROOT / "results"))
    p.add_argument("--keep-open", action="store_true", help="leave AEDT open when done")
    args = p.parse_args(argv)

    r = load_design(args.design)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    project = out / "Moreno_%s_%gGHz.aedt" % (r["inputs"]["waveguide"].replace("-", ""), r["f0_ghz"])

    Hfss = import_hfss()
    hfss = Hfss(
        project=str(project),
        design="Moreno_coupler",
        solution_type="Modal",
        version=args.version,
        student_version=args.student,
        non_graphical=args.non_graphical,
        new_desktop=True,
        close_on_exit=not args.keep_open,
    )
    try:
        build(hfss, r)
        hfss.save_project()
        print("Built %s" % project)

        if args.analyze:
            hfss.analyze_setup("Setup1", cores=args.cores)
            s4p = out / "moreno.s4p"
            hfss.export_touchstone(setup="Setup1", sweep="Sweep", output_file=str(s4p), renormalization=False)
            hfss.save_project()
            sys.path.insert(0, str(ROOT / "post"))
            from analyze_touchstone import analyze

            analyze(str(s4p), f0=r["f0_ghz"], out=str(out))
    finally:
        if not args.keep_open:
            hfss.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
