"""Build docs/report/Moreno_Coupler_WR22_40GHz_Report.pdf from the design and result files.

    python docs/report/collect.py <campaign>      # once, after a simulation campaign
    python docs/report/build_report.py

Every number in the report is computed here from design/design.json and
results/openems/**; nothing is typed in by hand. key_numbers.json (written
next to the PDF) holds the headline values for the README.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

import figures as FIG  # noqa: E402
from analysis import band_where, fit_power_law, read_csv, richardson  # noqa: E402
from layout import P, Report, bullets, figure, kpis, spacer, table  # noqa: E402
from reportlab.platypus import KeepTogether, NextPageTemplate, PageBreak  # noqa: E402
from scipy.interpolate import PchipInterpolator  # noqa: E402

from moreno.design import MorenoInputs, compute  # noqa: E402
from moreno.estimate import (coupling_estimate_db, directivity_estimate_db, resonance_factor,  # noqa: E402
                             slot_polarizability_mm3)
from moreno.touchstone import read_touchstone  # noqa: E402

import os  # noqa: E402
RES = Path(os.environ.get("MORENO_RESULTS", ROOT / "results" / "openems"))
STU = RES / "studies"
FIGDIR = Path(os.environ.get("MORENO_FIGDIR", HERE / "figures"))
OUTPDF = Path(os.environ.get("MORENO_PDF", HERE / "Moreno_Coupler_WR22_40GHz_Report.pdf"))
AUTHOR = "Amirreza Piri"
AFFIL = "Università degli Studi di Pavia"
DATE = "3 November 2026"
VERSION = "2.0"
FIGDIR.mkdir(exist_ok=True)


def fmt(x, n=2):
    return ("%%.%df" % n) % x


def sgn(x, n=2):
    v = round(float(x), n)
    v = 0.0 if v == 0 else v
    return (("%%+.%df" % n) % v).replace("-", "−")


def hmm(h):
    """Cell size without trailing zeros: 0.035, 0.0175, 0.0125."""
    return ("%.4f" % h).rstrip("0")


def rows(name):
    p = STU / ("%s.csv" % name)
    return read_csv(p) if p.exists() else []


def by_run(rs):
    return {r_["run"]: r_ for r_ in rs}


def curve(group, run):
    p = STU / ("%s__%s.csv" % (group, run))
    d = read_csv(p)
    return np.array([x["freq_ghz"] for x in d]), {k: np.array([x[k] for x in d]) for k in d[0] if k != "freq_ghz"}


# ---------------------------------------------------------------------------
# Load and analyse
# ---------------------------------------------------------------------------
r = json.loads((ROOT / "design" / "design.json").read_text())
assert r == json.loads(json.dumps(compute(MorenoInputs()))), "design/design.json is stale"
f0 = r["f0_ghz"]
a, b, t, L0, w0 = r["a_mm"], r["b_mm"], r["wall_mm"], r["slot_length_mm"], r["slot_width_mm"]
eta0 = 376.730313
Zte = eta0 / math.sqrt(1 - (r["fc_te10_ghz"] / f0) ** 2)
p_ratio = r["beta_rad_per_mm"] * a / math.pi

conv = sorted(rows("convergence"), key=lambda x: -x["slot_cell_mm"])        # coarse -> fine
stair = rows("staircase")
legacy = sorted(rows("legacy_air_priority_model"), key=lambda x: -x["slot_cell_mm"])
num = by_run(rows("numerics"))
sl = sorted(rows("slot_length"), key=lambda x: x["L_mm"])
off = sorted(rows("offset_check"), key=lambda x: x["L_mm"])
tol = by_run(rows("tolerance"))

fine3 = sorted(conv, key=lambda x: x["slot_cell_mm"])[:3]                    # h1 < h2 < h3
RICH = {k: richardson([x["slot_cell_mm"] for x in fine3], [x[k] for x in fine3])
        for k in ("coupling_db_at_f0", "directivity_db_at_f0", "return_loss_db_at_f0")}
RC = RICH["coupling_db_at_f0"]
finest = fine3[0]
h035 = next(x for x in conv if abs(x["slot_cell_mm"] - 0.035) < 1e-6)
# two independent cross-checks of the limit: the refine-by-2 triple and a power-law fit to all meshes
hs = {round(x["slot_cell_mm"], 5): x["coupling_db_at_f0"] for x in conv}
R2 = (richardson([0.0125, 0.025, 0.05], [hs[0.0125], hs[0.025], hs[0.05]])
      if all(k_ in hs for k_ in (0.0125, 0.025, 0.05)) else None)
FIT = fit_power_law(list(hs), list(hs.values()))
# The 1/3-2/3 series rises and the lines-on-edges series falls monotonically towards the same limit, so their
# finest values bracket it. With observed orders near 0.5, single-series extrapolation is not reliable; the
# bracket is. The quoted value is its midpoint, the uncertainty its half-width.
stair_all = sorted(rows("staircase"), key=lambda x: x["slot_cell_mm"])
thirds_seq = [x["coupling_db_at_f0"] for x in sorted(conv, key=lambda x: -x["slot_cell_mm"])]
stair_seq = [x["coupling_db_at_f0"] for x in sorted(stair_all, key=lambda x: -x["slot_cell_mm"])]
BRACKET = None
if stair_all and all(b_ > a_ for a_, b_ in zip(thirds_seq, thirds_seq[1:])) and \
        all(b_ < a_ for a_, b_ in zip(stair_seq, stair_seq[1:])) and stair_seq[-1] > thirds_seq[-1]:
    BRACKET = (thirds_seq[-1], stair_seq[-1])
if BRACKET:
    C_conv, U_conv = 0.5 * (BRACKET[0] + BRACKET[1]), 0.5 * (BRACKET[1] - BRACKET[0])
else:
    C_conv = RC["f_ext"]
    U_conv = max(RC["gci"], abs(C_conv - R2["f_ext"]) if R2 else 0.0, abs(C_conv - FIT["f0"]))
U_slot = U_conv
# the other numerical parameters (Sec. 6.3), each changed alone from the h = 0.035 mm run, add to the
# slot-mesh uncertainty; they are independent, so they are combined root-sum-square
NUM_D = {k_: num[k_]["coupling_db_at_f0"] - num["c035"]["coupling_db_at_f0"]
         for k_ in ("n_res30", "n_feed25", "n_e4") if k_ in num and "c035" in num}
U_conv = math.sqrt(U_slot ** 2 + sum(v * v for v in NUM_D.values()))
# directivity: the two edge treatments also approach it from opposite sides
D_th = [x["directivity_db_at_f0"] for x in sorted(conv, key=lambda x: -x["slot_cell_mm"])]
D_st = [x["directivity_db_at_f0"] for x in sorted(stair_all, key=lambda x: -x["slot_cell_mm"])]
D_BR = None
if D_st and all(b_ < a_ for a_, b_ in zip(D_th, D_th[1:])) and all(b_ > a_ for a_, b_ in zip(D_st, D_st[1:])) \
        and D_st[-1] < D_th[-1]:
    D_BR = (D_st[-1], D_th[-1])
D_conv = 0.5 * (D_BR[0] + D_BR[1]) if D_BR else D_th[-1]
U_D = 0.5 * (D_BR[1] - D_BR[0]) if D_BR else float("nan")
shift035 = C_conv - h035["coupling_db_at_f0"]           # correction applied to the h = 0.035 mm chart

# mesh shift between h = 0.035 and 0.025 mm at three slot lengths
sl_by_L = {round(x["L_mm"], 3): x for x in sl}
shift_check = []
for x in off:
    Lk = round(x["L_mm"], 3)
    if Lk in sl_by_L:
        shift_check.append((Lk, x["coupling_db_at_f0"] - sl_by_L[Lk]["coupling_db_at_f0"]))
shift_ref = dict(shift_check).get(round(L0, 3), None)
shift_spread = (max(s for _, s in shift_check) - min(s for _, s in shift_check)) if len(shift_check) > 1 else 0.0
# the h = 0.035 -> h = 0 correction scaled with the measured 0.035 -> 0.025 shift at each L
if shift_ref and len(shift_check) >= 2:
    ks = np.array([s / shift_ref for _, s in shift_check])
    kL = np.array([L_ for L_, _ in shift_check])
    corr_of_L = lambda L_: shift035 * np.interp(L_, kL, ks)  # noqa: E731
else:
    corr_of_L = lambda L_: shift035 + 0 * np.asarray(L_)  # noqa: E731

Ls = np.array([x["L_mm"] for x in sl])
C35 = np.array([x["coupling_db_at_f0"] for x in sl])
Ccorr = C35 + corr_of_L(Ls)
pchip = PchipInterpolator(Ls, Ccorr)
inv = PchipInterpolator(Ccorr[::-1], Ls[::-1])
U_chart = U_conv + abs(shift035) * (max(abs(k - 1) for k in ks) if shift_ref and len(shift_check) >= 2 else 0.0) \
    if shift_ref else U_conv
targets = []
for Ct in (20.0, 25.0, 30.0):
    targets.append((Ct, float(inv(Ct)) if Ccorr.min() <= Ct <= Ccorr.max() else None))
dCdL = float(pchip.derivative()(L0))                      # dB per mm (negative: longer = stronger)
est = np.array([coupling_estimate_db(f0, a, b, r["d_mm"], r["spacing_mm"], L_, w0, t) for L_ in Ls])
est_err = Ccorr - est
dir_est1 = directivity_estimate_db(r, arms=1)
dir_est2 = directivity_estimate_db(r, arms=2)

# tolerance sensitivities (central differences, dB per 0.01 mm)
dCdt = (tol["t030"]["coupling_db_at_f0"] - tol["t020"]["coupling_db_at_f0"]) / 0.10 * 0.01
dCdw = (tol["w040"]["coupling_db_at_f0"] - tol["w030"]["coupling_db_at_f0"]) / 0.10 * 0.01
dCdL01 = dCdL * 0.01

# reference (finest mesh) frequency response
ref = json.loads((RES / "reference" / "run_info.json").read_text())
refsum = ref["summary"]
fr, Sr, _ = read_touchstone(RES / "reference" / "moreno_openems.s4p")
fq = fr / 1e9
S1 = Sr[:, :, 0]
dB = lambda x: 20 * np.log10(np.maximum(np.abs(x), 1e-15))  # noqa: E731
Cf, If_, RLf, ILf = -dB(S1[:, 2]), -dB(S1[:, 3]), -dB(S1[:, 0]), -dB(S1[:, 1])
Df = If_ - Cf
at = lambda arr, fg: float(np.interp(fg, fq, arr))  # noqa: E731
d30 = band_where(fq, Df, 30.0)
d35 = band_where(fq, Df, 35.0)
bal_ref = np.sum(np.abs(S1) ** 2, axis=1)

# full matrix
ff, SF, _ = read_touchstone(RES / "full_matrix" / "moreno_openems.s4p")
full_info = json.loads((RES / "full_matrix" / "run_info.json").read_text())
recip = float(np.max(np.abs(SF - np.transpose(SF, (0, 2, 1)))))
colbal = np.sum(np.abs(SF) ** 2, axis=1)
mag = np.abs(SF)
sym = {
    "|S21| − |S43|": float(np.max(np.abs(mag[:, 1, 0] - mag[:, 3, 2]))),
    "|S31| − |S42|": float(np.max(np.abs(mag[:, 2, 0] - mag[:, 3, 1]))),
    "|S41| − |S32|": float(np.max(np.abs(mag[:, 3, 0] - mag[:, 2, 1]))),
    "|S11| − |S22|, |S33| − |S44|": float(max(np.max(np.abs(mag[:, 0, 0] - mag[:, 1, 1])),
                                              np.max(np.abs(mag[:, 2, 2] - mag[:, 3, 3])))),
}
diag_mirror = float(np.max(np.abs(mag[:, 0, 0] - mag[:, 2, 2])))
rl_g1 = -dB(SF[:, 0, 0])
rl_g2 = -dB(SF[:, 2, 2])
k40 = int(np.argmin(np.abs(ff - 40e9)))

KEY = {
    "coupling_extrapolated_db": C_conv, "coupling_uncertainty_db": U_conv, "finest": finest,
    "coupling_slot_mesh_uncertainty_db": U_slot, "directivity_converged_db": D_conv,
    "directivity_uncertainty_db": U_D,
    "finest_slot_cell_mm": finest["slot_cell_mm"], "reference_summary": refsum,
    "directivity_above_30_ghz": d30[:2], "chart": [{"L_mm": float(L_), "C_corrected": float(c)} for L_, c in zip(Ls, Ccorr)],
    "targets": targets, "chart_mesh_shift_db": float(shift035), "dC_dL_db_per_0.01mm": dCdL01, "dC_dw_db_per_0.01mm": dCdw, "dC_dt_db_per_0.01mm": dCdt,
    "copper_loss_db_per_m": r["copper_loss_db_per_m"], "reciprocity": recip,
    "power_balance_full": [float(colbal.min()), float(colbal.max())],
}
(OUTPDF.parent / "key_numbers.json").write_text(json.dumps(KEY, indent=2))

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
FIG.geometry(r, FIGDIR / "fig_geometry.png")
FIG.cross_section(r, FIGDIR / "fig_section.png")
FIG.coupling_law(r, FIGDIR / "fig_law.png")
FIG.convergence(conv, stair, legacy, {"f_ext": C_conv, "gci": U_slot}, FIGDIR / "fig_convergence.png")
FIG.sparams(RES / "reference" / "coupler_metrics.csv", f0, FIGDIR / "fig_sparams.png")
FIG.design_chart(list(Ls), list(C35), pchip, list(est), targets, FIGDIR / "fig_chart.png")
fam = []
for x in sl:
    fqq, cc = curve("slot_length", x["run"])
    fam.append(("L = %.1f mm" % x["L_mm"], fqq, cc["coupling_db"] + float(corr_of_L(x["L_mm"]))))
FIG.freq_family(fam, FIGDIR / "fig_family_c.png")
famd = [(lab, fqq, curve("slot_length", x["run"])[1]["directivity_db"]) for (lab, fqq, _), x in zip(fam, sl)]
FIG.freq_family(famd, FIGDIR / "fig_family_d.png", ylabel="directivity (dB)", invert=False)
has_fields = (RES / "fields" / "fields_f0.npz").exists()
if has_fields:
    FIG.field_maps(RES / "fields" / "fields_f0.npz", r, FIGDIR / "fig_fields.png")

# ---------------------------------------------------------------------------
# Document
# ---------------------------------------------------------------------------
st = []
A = st.append


def E(items):
    st.extend(items if isinstance(items, list) else [items])


def H(text, level=1):
    E(P(text, "h1" if level == 1 else "h2"))


# --- cover -------------------------------------------------------------------
E([spacer(70), P("TECHNICAL REPORT · VERSION %s" % VERSION, "bodyl")])
st[-1].style = st[-1].style.clone("kicker", fontName="Sans", fontSize=8.5, textColor="#A8461A")
E([spacer(6), P("Moreno Cross-Guide Directional Coupler", "title"), spacer(10),
   P("WR-22 rectangular waveguide, centre frequency 40 GHz: design, full-wave verification and "
     "manufacturing data", "subtitle"), spacer(24)])
A(table([["Author", "%s, %s" % (AUTHOR, AFFIL)],
         ["Date", DATE],
         ["Design", "Two crossed (+) slots in the common broad wall, on its diagonal, spaced λ<sub>g</sub>/4"],
         ["Verification", "openEMS FDTD, metal-only model, mesh-converged (two bracketing refinement series), full 4 × 4 S-matrix"],
         ["HFSS", "Parametric models supplied as scripts (PyAEDT and AEDT Run Script)"],
         ["Repository", "moreno-coupler: design code, HFSS and openEMS models, all results, this report"]],
        [0.2, 0.8], header=False, zebra=False, bold_first_col=True))
E([spacer(22), P(
    "<b>Abstract.</b> A cross-guide (Moreno-type) directional coupler is designed for WR-22 waveguide at 40 GHz. "
    "Two guides cross at right angles and share one broad wall; two crossed slots on the wall diagonal, spaced a "
    "quarter guide wavelength, couple magnetically into one end of the secondary guide only, where their "
    "contributions add. The slot positions follow in closed form from the TE<sub>10</sub> fields "
    "(d = %s mm from both walls, Δ = %s mm). The design was verified by full-wave FDTD simulation of a model "
    "whose own accuracy is demonstrated: power balance within %s %%, reciprocity error %s, the junction's "
    "symmetries reproduced, and two mesh-refinement series that bracket the result. With 2.00 × 0.35 mm slot arms in a "
    "0.25 mm wall the coupler gives <b>%s ± %s dB coupling</b> and <b>%s ± %s dB directivity</b> at 40 GHz, "
    "directivity above 30 dB over %s–%s GHz and input return loss of at least %s dB over 35–45 GHz. A design chart "
    "gives the slot length for other coupling levels (%s), and sensitivities to slot length, slot width and wall "
    "thickness set the manufacturing tolerances. During this review a modelling error in the earlier openEMS "
    "model was found and corrected; Appendix A lists it with every other correction."
    % (fmt(r["d_mm"], 3), fmt(r["spacing_mm"], 3), fmt(100 * max(abs(colbal.min() - 1), abs(colbal.max() - 1)), 2),
       "%.1e" % recip, fmt(C_conv, 2), fmt(U_conv, 2), fmt(D_conv, 1), fmt(U_D, 1),
       fmt(d30[0], 1), fmt(d30[1], 1), fmt(float(RLf.min()), 1),
       "; ".join("L ≈ %s mm for %g dB" % (fmt(L_, 2), C) for C, L_ in targets if L_ is not None)))])
A(NextPageTemplate("body"))
A(PageBreak())

# --- contents ------------------------------------------------------------------
H("Contents")
for line in ["1  Summary of results", "2  Scope and design inputs", "3  Theory of the cross-guide coupler",
             "4  Design calculation", "5  Full-wave model", "6  Verification of the numerical model",
             "7  Performance of the reference design", "8  Design chart and manufacturing tolerances",
             "9  HFSS model", "10  Limitations and recommendations", "11  Reproducing the results",
             "References", "Appendix A  Corrections made during this review", "Appendix B  Nomenclature"]:
    A(P(line, "bodyl"))
A(P("The PDF outline (bookmarks) links to every section and subsection.", "cap"))

# --- 1 summary -------------------------------------------------------------------
H("1  Summary of results")
A(kpis([("%s dB" % fmt(C_conv, 2), "coupling at 40 GHz<br/>converged, ± %s dB" % fmt(U_conv, 2)),
        ("%s dB" % fmt(D_conv, 1), "directivity at 40 GHz<br/>converged, ± %s dB" % fmt(U_D, 1)),
        ("%s dB" % fmt(refsum["return_loss_db_at_f0"], 1), "input return loss at<br/>40 GHz (lower bound)"),
        ("%s–%s" % (fmt(d30[0], 1), fmt(d30[1], 1)), "GHz with directivity<br/>above 30 dB"),
        ("%s dB" % fmt(abs(dCdL01), 2), "coupling change per<br/>0.01 mm of slot length")]))
A(spacer(8))
E(bullets([
    "<b>Geometry.</b> WR-22 (a = %s mm, b = %s mm), common wall t = %s mm, two crossed slots with arms "
    "L × w = %s × %s mm centred at (%s, %s) mm and (%s, %s) mm. Port 3, the end of the secondary guide nearest "
    "the corner (0, 0) where the slot diagonal starts, is the coupled port (Fig. 1)."
    % (fmt(a, 3), fmt(b, 3), fmt(t, 2), fmt(L0, 2), fmt(w0, 2), *[fmt(v, 3) for v in sum(r["slot_centres_mm"], [])]),
    "<b>Coupling.</b> Two mesh-refinement series with different edge treatments converge from opposite "
    "sides and bracket the converged coupling between %s and %s dB; with the smaller effects of the other "
    "numerical parameters the coupling is <b>%s ± %s dB</b>. As built, with L = 2.00 mm, the coupler is a "
    "%d dB coupler to within the numerical uncertainty."
    % (fmt(BRACKET[0], 2), fmt(BRACKET[1], 2), fmt(C_conv, 2), fmt(U_conv, 2), round(C_conv)),
    "<b>Directivity and match.</b> Directivity %s ± %s dB at 40 GHz (bracketed in the same way) and return "
    "loss %s dB; over 35–45 GHz the directivity on the finest mesh stays above %s dB and the return loss is at "
    "least %s dB (a lower bound, Sec. 6.3). The coupling rises by %s dB from 35 to "
    "45 GHz, mostly because the slot polarizability grows towards resonance (Sec. 3.5)."
    % (fmt(D_conv, 1), fmt(U_D, 1), fmt(refsum["return_loss_db_at_f0"], 1), fmt(float(Df.min()), 1),
       fmt(float(RLf.min()), 1), fmt(at(Cf, 35) - at(Cf, 45), 2)),
    "<b>Other coupling levels.</b> The slot length sets the level: %s (Fig. 7). Each 0.01 mm of extra slot "
    "length makes the coupling %s dB stronger; the same change in slot width or wall thickness moves it by "
    "%s dB and %s dB."
    % ("; ".join("%g dB needs L = %s mm" % (C, fmt(L_, 3)) for C, L_ in targets if L_ is not None),
       fmt(abs(dCdL01), 2), sgn(dCdw, 2), sgn(dCdt, 2)),
    "<b>Model validity.</b> Power balance %s–%s over all four excitations, reciprocity error %s, and the three "
    "transmission and coupling symmetries of the junction hold to within %s. Small-aperture theory predicts the coupling "
    "level within %s dB over the whole design chart and a directivity of %s–%s dB, close to the simulated value."
    % (fmt(colbal.min(), 4), fmt(colbal.max(), 4), "%.1e" % recip,
       "%.1e" % max(v for k_, v in sym.items() if "S11" not in k_),
       fmt(float(np.max(np.abs(est_err))), 1), fmt(min(dir_est1, dir_est2), 1), fmt(max(dir_est1, dir_est2), 1)),
    "<b>Correction.</b> The values given earlier in this project (24.3 dB coupling, 38.6 dB directivity) came "
    "from a model in which openEMS assigned every metal surface to the adjacent air and are withdrawn "
    "(Sec. 6.4, Appendix A).",
]))

# --- 2 scope -------------------------------------------------------------------
H("2  Scope and design inputs")
A(P("This report documents the design of a cross-guide directional coupler of the type introduced by Riblet and "
    "Saad [2] and tabulated in Moreno's design handbook [1]: two rectangular waveguides cross at right angles, one "
    "on top of the other, and couple through apertures in their common broad wall. The objectives were (i) to "
    "place the apertures so that the coupler is directional at 40 GHz, (ii) to verify the design with a full-wave "
    "simulation whose own accuracy is demonstrated, and (iii) to provide the data needed to choose the coupling "
    "level and to manufacture the part. Table 1 lists the inputs."))
A(KeepTogether([table([
    ["Quantity", "Symbol", "Value", "Source / note"],
    ["Waveguide", "–", "WR-22 (IEC R400)", "standard size for 33–50 GHz"],
    ["Broad / narrow inner dimension", "a, b", "%s / %s mm" % (fmt(a, 3), fmt(b, 3)), "0.224 × 0.112 in"],
    ["Recommended band", "–", "%g–%g GHz" % tuple(r["band_ghz"]), "single-mode TE<sub>10</sub>"],
    ["Centre frequency", "f<sub>0</sub>", "%s GHz" % fmt(f0, 1), "design requirement"],
    ["Common-wall thickness", "t", "%s mm" % fmt(t, 2), "thin for low slot loss, machinable"],
    ["Slot arm length × width", "L × w", "%s × %s mm" % (fmt(L0, 2), fmt(w0, 2)), "L sets the coupling (Sec. 8)"],
    ["Wall material (simulation)", "–", "perfect conductor", "copper loss estimated separately (Sec. 7.3)"],
    ["Analysis band", "–", "%g–%g GHz" % tuple(r["sweep_ghz"]), "f<sub>0</sub> ± 12.5 %"]],
    [0.31, 0.11, 0.22, 0.36]), P("Table 1. Design inputs.", "cap")]))

# --- 3 theory ------------------------------------------------------------------
H("3  Theory of the cross-guide coupler")
H("3.1  Fields at the common wall", 2)
A(P("Let the primary guide (guide 1) run along z and the secondary guide (guide 2) along x, both of width a; the "
    "common wall is the square 0 ≤ x, z ≤ a (Fig. 1). For a TE<sub>10</sub> wave of unit electric-field amplitude "
    "travelling in +z in guide 1, the magnetic field tangential to the wall is"))
A(P("<b>H</b><sub>1</sub> = (β/ωμ<sub>0</sub>) [ −sin(πx/a) <b>x̂</b> + (j/p) cos(πx/a) <b>ẑ</b> ],   "
    "p = βa/π         (1)", "eq"))
A(P("and the fields of the two waves that guide 2 can carry, towards −x (port 3) and towards +x (port 4), follow "
    "by exchanging x and z. A small aperture at (x<sub>0</sub>, z<sub>0</sub>) acts as a magnetic dipole "
    "<b>M</b> = α<sub>m</sub><b>H</b><sub>1</sub> (Bethe [3]); the amplitude of each wave it launches in guide 2 is "
    "jωμ<sub>0</sub> <b>H</b><sub>2</sub>·<b>M</b> / P<sub>10</sub> with the mode normalisation "
    "P<sub>10</sub> = ab/Z<sub>TE</sub> [4]. For an aperture whose magnetic polarizability is the same in every "
    "in-plane direction (true for a symmetric cross, α<sub>xx</sub> = α<sub>zz</sub> = α<sub>m</sub>) this gives"))
A(P("|c<sub>3</sub>| = (π α<sub>m</sub> / a<sup>2</sup>b) |sin π(x<sub>0</sub> + z<sub>0</sub>)/a|,   "
    "|c<sub>4</sub>| = (π α<sub>m</sub> / a<sup>2</sup>b) |sin π(x<sub>0</sub> − z<sub>0</sub>)/a|         (2)", "eq"))
A(P("Each scalar product has one transverse and one longitudinal term, so β cancels against the normalisation: "
    "<b>between crossed guides the magnetic coupling of an aperture does not depend on frequency</b> except through "
    "α<sub>m</sub>. Towards port 4 the two terms subtract and on the diagonal x<sub>0</sub> = z<sub>0</sub> they "
    "cancel exactly, so an aperture there couples magnetically only towards port 3. Fig. 2 maps both factors; "
    "equation (2) is checked against the explicit mode fields at 500 random positions in the test suite "
    "(tests/test_design.py)."))
A(figure(FIGDIR / "fig_law.png", "Figure 2. Relative magnetic coupling of a point aperture towards port 3 (left) "
         "and port 4 (right) versus its position on the common wall, from (2). Crosses: slot centres of this design.",
         0.82))
H("3.2  Two-slot array", 2)
A(P("Two slots are placed on the diagonal, at (d, d) and (a − d, a − d), symmetrically about the centre of the "
    "wall. From (2), their port-3 contributions have opposite sign, because sin[2π(a − d)/a] = −sin(2πd/a). The "
    "incident wave reaches slot 2 later by βΔ, with Δ = a − 2d, and the wave launched towards port 3 by slot 2 "
    "travels a further Δ to pass slot 1, so slot 2's contribution lags by 2βΔ in total. Choosing βΔ = π/2 makes "
    "this lag π, which cancels the sign difference, and the two contributions add at port 3:"))
A(P("Δ = λ<sub>g</sub>/4,   d = (a − λ<sub>g</sub>/4)/2         (3)", "eq"))
A(P("The pair then couples 2|sin(2πd/a)||sin βΔ| times as strongly as one slot at the best position d = a/4; here "
    "that is %s dB. Towards port 4 the path lengths z<sub>0</sub> + (a − x<sub>0</sub>) of the two slots are equal "
    "(= a), and the magnetic coupling of each is zero on the diagonal." % sgn(r["pair_vs_single_slot_db"], 2)))
A(P("The junction has two symmetries: a 180° rotation about the vertical axis through the wall centre, which maps "
    "port 1 ↔ 2 and 3 ↔ 4, and the mirror x ↔ z combined with an exchange of the two guides, which maps "
    "port 1 ↔ 3 and 2 ↔ 4. Together with reciprocity they require |S<sub>21</sub>| = |S<sub>43</sub>|, "
    "|S<sub>31</sub>| = |S<sub>42</sub>|, |S<sub>41</sub>| = |S<sub>32</sub>| and equal |S<sub>ii</sub>| at all four "
    "ports; Sec. 6.1 uses these as a check of the simulation."))
H("3.3  What limits the directivity", 2)
A(P("A slot also couples through the electric field normal to the wall, as an electric dipole of polarizability "
    "α<sub>e</sub> ≈ πw<sup>2</sup>L/16 per arm for a narrow slot [5]. This coupling is proportional to "
    "sin(πx<sub>0</sub>/a) sin(πz<sub>0</sub>/a), equal in both directions, and therefore not directional. With the "
    "λ<sub>g</sub>/4 spacing the electric contributions of the two slots cancel at port 3 (path difference "
    "λ<sub>g</sub>/2, equal amplitudes) and add at port 4 (equal paths): the electric leakage, not the magnetic "
    "coupling, sets the isolation. Relative to the magnetic coupling it is smaller by the factor "
    "(Z<sub>TE</sub>/η<sub>0</sub>)<sup>2</sup> p (α<sub>e</sub>/α<sub>m</sub>) sin<sup>2</sup>X / sin 2X, X = πd/a, "
    "and it decays faster through the wall, like the slot's cut-off mode across its width (k<sub>c</sub> = π/w) "
    "instead of along its length (k<sub>c</sub> = π/L). The resulting first-order directivity is "
    "<b>%s dB</b> counting one arm per cross and <b>%s dB</b> counting both, against %s ± %s dB simulated "
    "(Sec. 6.2)."
    % (fmt(dir_est1, 1), fmt(dir_est2, 1), fmt(D_conv, 1), fmt(U_D, 1))))
H("3.4  Thickness of the common wall", 2)
A(P("Each slot arm is a short rectangular guide of length t whose lowest mode has cut-off wavelength 2L. Below "
    "resonance (L &lt; λ<sub>0</sub>/2) the field decays through it, and the coupling falls by "
    "8.686 t [(π/L)<sup>2</sup> − k<sub>0</sub><sup>2</sup>]<sup>1/2</sup> dB, %s dB for this design. The wall "
    "must therefore be thin compared with L, but not so thin that it is hard to machine flat."
    % fmt(r["wall_loss_db"], 2)))
H("3.5  Closed-form estimate of the coupling level", 2)
alpha_m = slot_polarizability_mm3(L0, w0)
res_fac = resonance_factor(f0, L0)
e_thin = coupling_estimate_db(f0, a, b, r["d_mm"], r["spacing_mm"], L0, w0, t, resonance=False, wall=False)
e_res = coupling_estimate_db(f0, a, b, r["d_mm"], r["spacing_mm"], L0, w0, t, resonance=True, wall=False)
e_full = coupling_estimate_db(f0, a, b, r["d_mm"], r["spacing_mm"], L0, w0, t)
A(P("Equation (2) with the narrow-slot polarizability α<sub>m</sub> = πL<sup>3</sup> / [24 (ln(4L/w) − 1)] [5] "
    "(%s mm<sup>3</sup> here), Cohn's resonance factor 1/[1 − (2L/λ<sub>0</sub>)<sup>2</sup>] (%s) [6] and the "
    "wall attenuation of Sec. 3.4 gives for the slot pair at 40 GHz: %s dB (thin wall, static polarizability), "
    "%s dB (with resonance) and <b>%s dB</b> (with resonance and wall). Sec. 8 compares this estimate with the "
    "simulation over the slot length: it follows the slope of the design chart and differs from the "
    "mesh-converged coupling by %s to %s dB. That is closer than first-order small-aperture theory can promise "
    "for slots of half their resonant length, and the corrections partly compensate one another, so the "
    "estimate is used to choose a starting slot size and to read trends; the coupling level itself comes from "
    "the full-wave model."
    % (fmt(alpha_m, 3), fmt(res_fac, 3), fmt(e_thin, 1), fmt(e_res, 1), fmt(e_full, 1),
       sgn(float(np.min(-est_err)), 2), sgn(float(np.max(-est_err)), 2))))

# --- 4 design ----------------------------------------------------------------------
H("4  Design calculation")
A(P("Table 2 lists the calculation in order. All values come from moreno/design.py and are stored in "
    "design/design.json, which the HFSS and openEMS models both read."))
c1, c2 = r["slot_centres_mm"]
A(KeepTogether([table([
    ["Step", "Expression", "Value"],
    ["TE<sub>10</sub> / TE<sub>20</sub> cut-off", "f<sub>c</sub> = c/(2a), 2f<sub>c</sub>",
     "%s / %s GHz" % (fmt(r["fc_te10_ghz"], 4), fmt(r["fc_te20_ghz"], 4))],
    ["Free-space wavelength", "λ<sub>0</sub> = c/f<sub>0</sub>", "%s mm" % fmt(r["lambda0_mm"], 4)],
    ["Guide wavelength", "λ<sub>g</sub> = λ<sub>0</sub> / √(1 − (f<sub>c</sub>/f<sub>0</sub>)<sup>2</sup>)",
     "%s mm" % fmt(r["lambda_g_mm"], 4)],
    ["Phase constant", "β = 2π/λ<sub>g</sub>", "%s rad/mm" % fmt(r["beta_rad_per_mm"], 4)],
    ["TE<sub>10</sub> wave impedance", "Z<sub>TE</sub> = η<sub>0</sub>λ<sub>g</sub>/λ<sub>0</sub>", "%s Ω" % fmt(Zte, 1)],
    ["Slot spacing (per axis)", "Δ = λ<sub>g</sub>/4", "%s mm" % fmt(r["spacing_mm"], 4)],
    ["Spacing along the diagonal", "√2 Δ", "%s mm" % fmt(r["spacing_diagonal_mm"], 4)],
    ["Slot offset from both walls", "d = (a − Δ)/2", "%s mm" % fmt(r["d_mm"], 4)],
    ["Slot 1 / slot 2 centre (x, z)", "(d, d) / (a − d, a − d)",
     "(%s, %s) / (%s, %s) mm" % (fmt(c1[0], 4), fmt(c1[1], 4), fmt(c2[0], 4), fmt(c2[1], 4))],
    ["Electrical spacing", "βΔ", "%s°" % fmt(r["electrical_spacing_deg"], 2)],
    ["Pair versus one slot at a/4", "2|sin(2πd/a)||sin βΔ|", "%s dB" % sgn(r["pair_vs_single_slot_db"], 2)],
    ["Slot resonant length", "λ<sub>0</sub>/2", "%s mm (L = %s %%)" % (fmt(r["slot_resonant_length_mm"], 3),
                                                                        fmt(100 * L0 / r["slot_resonant_length_mm"], 0))],
    ["Arm end to side wall", "d − L/2", "%s mm" % fmt(r["clearance_mm"], 3)],
    ["Narrowest web between crosses", "√2 [Δ − (L + w)/2]", "%s mm" % fmt(r["slot_gap_mm"], 3)],
    ["Slot attenuation through wall", "8.686 t √((π/L)<sup>2</sup> − k<sub>0</sub><sup>2</sup>)",
     "%s dB" % fmt(r["wall_loss_db"], 3)],
    ["Copper loss of the feeds", "TE<sub>10</sub> α<sub>c</sub> [4, Sec. 3.3], σ = 5.8 × 10<sup>7</sup> S/m", "%s dB/m" % fmt(r["copper_loss_db_per_m"], 3)]],
    [0.31, 0.39, 0.30]),
    P("Table 2. Design calculation for f<sub>0</sub> = 40 GHz (c = 299.792458 mm·GHz, η<sub>0</sub> = 376.73 Ω).", "cap")]))
A(P("Every design check of moreno/design.py passes: f<sub>0</sub> lies between 1.25 f<sub>c</sub> and "
    "0.95 × 2f<sub>c</sub> (%s–%s GHz); the arm ends stay %s mm from the side walls and the two crosses %s mm "
    "apart; L is %s %% of the resonant length; βΔ = 90°; and the arm aspect ratio L/w = %s keeps each arm a "
    "narrow slot. Fig. 1 shows the layout to scale and Fig. 3 a section through slot 1."
    % (fmt(1.25 * r["fc_te10_ghz"], 2), fmt(1.9 * r["fc_te10_ghz"], 2), fmt(r["clearance_mm"], 3),
       fmt(r["slot_gap_mm"], 3), fmt(100 * L0 / r["slot_resonant_length_mm"], 0), fmt(L0 / w0, 1))))
A(figure(FIGDIR / "fig_geometry.png", "Figure 1. (a) Top view of the junction with the port assignment; guide 1 "
         "(input) lies below the common wall and guide 2 above it. (b) The common wall with the slot coordinates. "
         "Dimensions in mm, to scale."))
A(figure(FIGDIR / "fig_section.png", "Figure 3. Section in the x–y plane through slot 1 (z = d), to scale.", 0.62))

# --- 5 model -----------------------------------------------------------------------
H("5  Full-wave model")
A(P("The design was verified with openEMS [7], a finite-difference time-domain (FDTD) solver [8], because the HFSS "
    "Student installation available for this project refuses to open an HFSS design (a licensing fault outside the "
    "project, Sec. 9). The model, openems/simulate_openems.py, is built from design/design.json and has the same "
    "geometry, port assignment and S-parameter normalisation as the HFSS scripts; it differs only in what is "
    "specific to each solver (feed length, absorbing boundaries, mesh)."))
run0 = ref["runs"][0]
E(bullets([
    "<b>Geometry: metal only.</b> All metal is drawn as boxes, the complement of the two guides and the four slot "
    "arms (19 boxes), on a vacuum background; the top and bottom faces are closed by the PEC boundary condition. "
    "Air is deliberately not drawn: openEMS makes an E-field edge a perfect conductor only if the highest-priority "
    "primitive at the edge's midpoint is a metal, so air boxes carved out of a metal block with higher priority "
    "would also own every edge lying on a metal surface and shift each wall one cell into the metal (Sec. 6.4). "
    "The box decomposition is checked against the intended air region at 200 000 random points with no mismatch.",
    "<b>Walls.</b> Perfect electric conductor; the model is lossless (copper loss is estimated in Sec. 7.3).",
    "<b>Feeds and ports.</b> Each guide extends %s λ<sub>g</sub> (%s mm) beyond the common square into an 8-cell "
    "PML. TE<sub>10</sub> wave ports with analytic mode profiles sit 10 feed cells inside the PML-free region "
    "(guide 2's broad dimension lies along z, its second transverse axis, so its ports use TE<sub>01</sub> with a "
    "and b exchanged). S-parameters are normalised to each port's TE<sub>10</sub> wave impedance, and every "
    "reference plane is moved to the edge of the common square (z = 0, z = a, x = 0, x = a)."
    % (fmt(run0["feed_lambda_g"], 1), fmt(run0["feed_lambda_g"] * r["lambda_g_mm"], 2)),
    "<b>Excitation.</b> Gaussian pulse covering 35–45 GHz with margin; the run stops when the field energy has "
    "decayed by 50 dB. Frequency responses are evaluated at 201 points.",
    "<b>Mesh.</b> Graded rectilinear mesh, %d cells per wavelength at 45 GHz in the feeds; a fine cell h across the "
    "slot arms and through the wall; mesh lines exactly on every guide wall and wall face; the 1/3–2/3 rule at the "
    "slot edges (lines h/3 inside the metal and 2h/3 inside the slot), which places the effective conducting edge "
    "where the singular edge field needs it [7]. h is the critical numerical parameter (Sec. 6.2)."
    % run0["res_cells_per_lambda"],
    "<b>Reference run.</b> h = %s mm: %d × %d × %d lines, %s million cells, %s s on two CPU cores."
    % (fmt(run0["fine_mm"], 4), *run0["mesh_lines_xyz"], fmt(run0["cells"] / 1e6, 2), fmt(run0["runtime_s"], 0)),
]))

# --- 6 verification ----------------------------------------------------------------
H("6  Verification of the numerical model")
H("6.1  Power balance, reciprocity and symmetry", 2)
A(P("The full 4 × 4 scattering matrix was computed with four runs, one per excited port (h = %s mm, "
    "results/openems/full_matrix). For a lossless, reciprocal network every column must carry unit power and the "
    "matrix must be symmetric; the junction's two symmetries (Sec. 3.2) add further relations. Table 3 shows "
    "that over the whole band the power balance holds to %s %%, reciprocity to %s and every transmission "
    "and coupling symmetry to %s. The one relation that holds less well is the equality of the input "
    "reflections of the two guides: |S<sub>11</sub>| = |S<sub>22</sub>| and |S<sub>33</sub>| = |S<sub>44</sub>| "
    "hold exactly, but the return loss of guide 1 (%s–%s dB) and of guide 2 (%s–%s dB) differ by up to %s dB. "
    "These reflections are only %s–%s dB below the incident wave and partly numerical (Sec. 6.3), so the "
    "difference reflects the port discretisation rather than the junction."
    % (fmt(full_info["runs"][0]["fine_mm"], 3),
       fmt(100 * max(abs(colbal.min() - 1), abs(colbal.max() - 1)), 2), "%.1e" % recip,
       "%.1e" % max(v for k_, v in sym.items() if "S11" not in k_), fmt(rl_g1.min(), 1), fmt(rl_g1.max(), 1), fmt(rl_g2.min(), 1),
       fmt(rl_g2.max(), 1), fmt(float(np.max(np.abs(rl_g1 - rl_g2))), 1),
       "%d" % math.floor(min(rl_g1.min(), rl_g2.min())), "%d" % math.ceil(max(rl_g1.max(), rl_g2.max())))))
A(KeepTogether([table([
    ["Check (35–45 GHz, 201 points)", "Exact value", "Simulated"],
    ["Σ<sub>k</sub>|S<sub>kj</sub>|<sup>2</sup>, all four columns", "1",
     "%s … %s" % (fmt(colbal.min(), 4), fmt(colbal.max(), 4))],
    ["max |S<sub>ij</sub> − S<sub>ji</sub>| (reciprocity)", "0", "%.1e" % recip]] +
    [["max %s (symmetry)" % k_.replace("S21", "S<sub>21</sub>").replace("S43", "S<sub>43</sub>")
      .replace("S31", "S<sub>31</sub>").replace("S42", "S<sub>42</sub>").replace("S41", "S<sub>41</sub>")
      .replace("S32", "S<sub>32</sub>").replace("S11", "S<sub>11</sub>").replace("S22", "S<sub>22</sub>")
      .replace("S33", "S<sub>33</sub>").replace("S44", "S<sub>44</sub>"), "0",
      "%.1e" % v] for k_, v in sym.items()] +
    [["max ||S<sub>11</sub>| − |S<sub>33</sub>|| (mirror, reflections)", "0", "%.1e (see text)" % diag_mirror]],
    [0.52, 0.18, 0.30]),
    P("Table 3. Consistency checks of the full scattering matrix (magnitudes are linear).", "cap")]))
all_runs = list({x["run"]: x for x in conv + stair + sl + list(tol.values()) + list(num.values()) + off}.values())
dev = lambda x: max(abs(x["power_balance_min"] - 1), abs(x["power_balance_max"] - 1))  # noqa: E731
regular = [x for x in all_runs if x["run"] != "n_e4"]
A(P("The single-excitation study runs (Tables 4, 5, 7 and 8) satisfy the same power balance: across the %d runs "
    "that stop at −50 dB energy decay, Σ<sub>k</sub>|S<sub>k1</sub>|<sup>2</sup> deviates from 1 by at most %s %%. The one run deliberately "
    "stopped early, at −40 dB (Sec. 6.3), deviates by up to %s %%, as expected for a truncated time signal. Ports, "
    "mode normalisation and de-embedding are therefore consistent."
    % (len(regular), fmt(100 * max(dev(x) for x in regular), 2), fmt(100 * dev(num["n_e4"]), 1) if "n_e4" in num else "–")))

H("6.2  Mesh convergence at the slots", 2)
A(P("Table 4 and Fig. 4 show the results at 40 GHz as the slot and wall cell h is refined from %s mm to %s mm "
    "(from %d to %d cells across the slot width), everything else fixed."
    % (fmt(conv[0]["slot_cell_mm"], 3), fmt(conv[-1]["slot_cell_mm"], 4),
       round(w0 / conv[0]["slot_cell_mm"]), round(w0 / conv[-1]["slot_cell_mm"]))))
trow = [["Slot cell h (mm)", "Cells (M)", "Run time (s)", "Coupling (dB)", "Directivity (dB)", "Return loss (dB)"]]
for x in conv:
    trow.append([hmm(x["slot_cell_mm"]), fmt(x["cells"] / 1e6, 2), fmt(x["runtime_s"], 0),
                 fmt(x["coupling_db_at_f0"], 2), fmt(x["directivity_db_at_f0"], 2), fmt(x["return_loss_db_at_f0"], 2)])
trow.append(["converged (bracket)", "–", "–", "<b>%s ± %s</b>" % (fmt(C_conv, 2), fmt(U_slot, 2)),
             "<b>%s ± %s</b>" % (fmt(D_conv, 2), fmt(U_D, 2)), "–"])
A(KeepTogether([table(trow, [0.2, 0.12, 0.14, 0.18, 0.18, 0.18], align_right_cols=(1, 2, 3, 4, 5)),
                P("Table 4. Slot-mesh refinement at 40 GHz (L = 2.00 mm, w = 0.35 mm, t = 0.25 mm, 1/3–2/3 edge "
                  "rule). Last row: converged coupling and directivity, bracketed by this series and the "
                  "lines-on-edges series (slot-mesh uncertainty only; see text).", "cap")]))
A(P("With the 1/3–2/3 rule the coupling converges monotonically towards weaker coupling: successive refinements "
    "change it by %s dB. The same model with mesh lines exactly on the slot edges (green in Fig. 4a) converges "
    "monotonically from the other side: %s. The difference between the two edge treatments shrinks steadily, "
    "from %s dB at h = %s mm to %s dB at h = %s mm, so both series approach one limit and their finest values "
    "bound it: <b>%s dB ≤ C ≤ %s dB</b>, i.e. %s ± %s dB as far as the slot mesh is concerned. Adding the "
    "other numerical parameters of Sec. 6.3 root-sum-square (%s) gives the coupling quoted throughout, "
    "<b>%s ± %s dB</b>."
    % (", ".join(fmt(thirds_seq[i + 1] - thirds_seq[i], 2) for i in range(len(thirds_seq) - 1)),
       "; ".join("%s dB at h = %s mm" % (fmt(x["coupling_db_at_f0"], 2), hmm(x["slot_cell_mm"])) for x in stair_all[::-1]),
       fmt(stair_all[-1]["coupling_db_at_f0"] - hs[round(stair_all[-1]["slot_cell_mm"], 5)], 2),
       hmm(stair_all[-1]["slot_cell_mm"]),
       fmt(stair_all[0]["coupling_db_at_f0"] - hs[round(stair_all[0]["slot_cell_mm"], 5)], 2),
       hmm(stair_all[0]["slot_cell_mm"]),
       fmt(BRACKET[0], 2), fmt(BRACKET[1], 2), fmt(C_conv, 2), fmt(U_slot, 2),
       ", ".join("%s dB" % sgn(v, 2) for v in NUM_D.values()), fmt(C_conv, 2), fmt(U_conv, 2))))
RS = (richardson([x["slot_cell_mm"] for x in stair_all[:3]], [x["coupling_db_at_f0"] for x in stair_all[:3]])
      if len(stair_all) >= 3 else None)
ests = [RC["f_ext"], FIT["f0"]] + ([R2["f_ext"]] if R2 else []) + ([RS["f_ext"]] if RS else [])
A(P("Both series converge slowly (observed orders p ≈ %s and %s), as expected where the field is singular at "
    "the re-entrant corners of the crosses. At such low orders Richardson extrapolation of a single series [9] is "
    "unreliable: it gives %s dB from the three finest 1/3–2/3 meshes (grid-convergence index %s dB), %s dB from "
    "the refine-by-two sequence h = 0.05, 0.025, 0.0125 mm, %s dB from a power-law fit to all five meshes and %s dB "
    "from the lines-on-edges series. These estimates scatter over %s dB; the bracket does not depend on any "
    "assumed order of convergence and is therefore used instead."
    % (fmt(RC["p"], 2), fmt(RS["p"], 2) if RS else "–", fmt(RC["f_ext"], 2), fmt(RC["gci"], 2),
       fmt(R2["f_ext"], 2) if R2 else "–", fmt(FIT["f0"], 2), fmt(RS["f_ext"], 2) if RS else "–",
       fmt(max(ests) - min(ests), 2))))
A(P("The directivity behaves in the same way: it decreases on refinement with the 1/3–2/3 rule (%s dB) and "
    "increases with lines on the slot edges (%s dB), so it is bracketed as %s dB ≤ D ≤ %s dB, i.e. <b>%s ± %s dB</b>. "
    "It depends on a residual coupling to port 4 that is %s–%s dB below the incident wave, and it stays well above "
    "30 dB on every mesh. The return loss changes by less than %s dB on the slot-mesh refinement; it is governed by "
    "the feed mesh instead (Sec. 6.3)."
    % (" → ".join(fmt(v, 2) for v in D_th), " → ".join(fmt(v, 2) for v in D_st), fmt(D_BR[0], 2), fmt(D_BR[1], 2),
       fmt(D_conv, 2), fmt(U_D, 2), "%d" % math.floor(If_.min()), "%d" % math.ceil(If_.max()),
       fmt(np.ptp([x["return_loss_db_at_f0"] for x in conv]), 1))))
A(figure(FIGDIR / "fig_convergence.png", "Figure 4. (a) Coupling at 40 GHz versus slot cell size h for the "
         "metal-only model with the 1/3–2/3 edge rule (blue) and with mesh lines on the slot edges (green); star and "
         "band: the bracketed converged value and its uncertainty. (b) The same blue series against the superseded air-priority model "
         "(grey, Sec. 6.4). (c) Directivity and return loss of the metal-only model."))

H("6.3  Other numerical parameters", 2)
base = num.get("c035")
nrows = [["Change from the h = 0.035 mm run", "Coupling (dB)", "Δ (dB)", "Directivity (dB)", "Return loss (dB)"]]
labels = {"c035": "none (baseline: 20 cells/λ, feed 1.5 λ<sub>g</sub>, −50 dB stop)",
          "n_res30": "feed mesh 30 cells/λ", "n_feed25": "feeds 2.5 λ<sub>g</sub> long",
          "n_e4": "stop at −40 dB energy decay"}
for k_ in ("c035", "n_res30", "n_feed25", "n_e4"):
    if k_ in num:
        x = num[k_]
        nrows.append([labels[k_], fmt(x["coupling_db_at_f0"], 2), sgn(x["coupling_db_at_f0"] - base["coupling_db_at_f0"], 2),
                      fmt(x["directivity_db_at_f0"], 2), fmt(x["return_loss_db_at_f0"], 2)])
A(KeepTogether([table(nrows, [0.43, 0.14, 0.11, 0.16, 0.16], align_right_cols=(1, 2, 3, 4)),
                P("Table 5. Sensitivity to the remaining numerical parameters, each changed alone.", "cap")]))
dn = [abs(num[k_]["coupling_db_at_f0"] - base["coupling_db_at_f0"]) for k_ in ("n_res30", "n_feed25", "n_e4") if k_ in num]
rl_d = num["n_res30"]["return_loss_db_at_f0"] - base["return_loss_db_at_f0"] if "n_res30" in num else 0.0
A(P("None of these parameters moves the coupling by more than %s dB, below the slot-mesh uncertainty; they are "
    "nevertheless included in the uncertainty quoted for the coupling (Sec. 6.2). The return loss is different: a finer "
    "feed mesh raises it by %s dB, because part of the computed reflection is numerical (the analytic port mode "
    "and the discretised guide differ slightly in dispersion). The simulated return loss is therefore a lower "
    "bound on that of the real junction." % (fmt(max(dn), 2), fmt(rl_d, 1))))

H("6.4  A modelling error found and corrected", 2)
lg_spread = (legacy[-1]["coupling_db_at_f0"] - legacy[0]["coupling_db_at_f0"]) if legacy else float("nan")
A(P("The first openEMS model drew a solid metal block and carved the guides and slots out of it as air boxes of "
    "higher priority. openEMS sets an E-field edge to PEC only if the highest-priority primitive containing the "
    "edge's midpoint is a metal (Operator::CalcPEC_Range in the openEMS source). Edges lying on a metal surface are "
    "contained in both the metal block and the air box, so the air won and the tangential electric field on every "
    "conducting surface was left free: each wall effectively moved one cell into the metal. The slots appeared "
    "longer and wider, the common wall thinner and the guides wider, all by amounts that depend on the local mesh. "
    "The symptom was a coupling that kept changing on refinement: %s dB over the six meshes of the superseded "
    "model (grey in Fig. 4a), including the 24.3 dB reported earlier in this project. In the metal-only model of "
    "Sec. 5 the metal owns its surfaces, and the coupling changes by only %s dB over its whole refinement "
    "(h = %s to %s mm). All "
    "results in this report come from the corrected model; the superseded series is kept in "
    "results/openems/studies/legacy_air_priority_model.csv for reference only."
    % (fmt(abs(lg_spread), 1), fmt(abs(conv[0]["coupling_db_at_f0"] - conv[-1]["coupling_db_at_f0"]), 2),
       fmt(conv[0]["slot_cell_mm"], 3), fmt(conv[-1]["slot_cell_mm"], 4))))

# --- 7 performance -----------------------------------------------------------------
H("7  Performance of the reference design")
A(P("Fig. 5 shows the response of the reference design on the finest mesh (h = %s mm) and Table 6 summarises it. "
    "Reference planes are at the edges of the common square." % fmt(finest["slot_cell_mm"], 4)))
A(figure(FIGDIR / "fig_sparams.png", "Figure 5. Scattering parameters for excitation at port 1 (top) and "
         "directivity 20 log(|S<sub>31</sub>|/|S<sub>41</sub>|) (bottom), finest mesh (h = 0.0125 mm). Dashed line: "
         "40 GHz."))
A(KeepTogether([table([
    ["Quantity", "35 GHz", "40 GHz", "45 GHz", "Over 35–45 GHz"],
    ["Coupling −20 log|S<sub>31</sub>| (dB)", fmt(at(Cf, 35), 2), fmt(at(Cf, 40), 2), fmt(at(Cf, 45), 2),
     "%s–%s (range %s)" % (fmt(Cf.min(), 2), fmt(Cf.max(), 2), fmt(np.ptp(Cf), 2))],
    ["Coupling, mesh-converged (dB)", "–", "%s ± %s" % (fmt(C_conv, 2), fmt(U_conv, 2)), "–", "–"],
    ["Directivity (dB)", fmt(at(Df, 35), 1), fmt(at(Df, 40), 1), fmt(at(Df, 45), 1), "min %s" % fmt(Df.min(), 1)],
    ["Directivity, mesh-converged (dB)", "–", "%s ± %s" % (fmt(D_conv, 1), fmt(U_D, 1)), "–", "–"],
    ["Isolation −20 log|S<sub>41</sub>| (dB)", fmt(at(If_, 35), 1), fmt(at(If_, 40), 1), fmt(at(If_, 45), 1),
     "min %s" % fmt(If_.min(), 1)],
    ["Return loss −20 log|S<sub>11</sub>| (dB)", fmt(at(RLf, 35), 1), fmt(at(RLf, 40), 1), fmt(at(RLf, 45), 1),
     "min %s" % fmt(RLf.min(), 1)],
    ["Through loss −20 log|S<sub>21</sub>| (dB)", fmt(at(ILf, 35), 3), fmt(at(ILf, 40), 3), fmt(at(ILf, 45), 3),
     "max %s" % fmt(ILf.max(), 3)]],
    [0.32, 0.12, 0.15, 0.12, 0.29], align_right_cols=(1, 2, 3)),
    P("Table 6. Performance of the reference design (lossless walls, finest mesh). Directivity above 30 dB from "
      "%s to %s GHz%s; above 35 dB from %s to %s GHz%s."
      % (fmt(d30[0], 1), fmt(d30[1], 1), " (whole band)" if d30[2] else "", fmt(d35[0], 1), fmt(d35[1], 1),
         " (whole band)" if d35[2] else ""), "cap")]))
H("7.1  Through loss and coupled power", 2)
A(P("In the lossless model the through loss is set entirely by the power taken out by the coupler and by "
    "reflection: −10 log(1 − |S<sub>31</sub>|<sup>2</sup> − |S<sub>41</sub>|<sup>2</sup> − |S<sub>11</sub>|<sup>2</sup>) = "
    "%s dB at 40 GHz, against %s dB read directly from S<sub>21</sub>; the small difference is within the "
    "power-balance accuracy of Sec. 6.1."
    % (fmt(-10 * np.log10(1 - 10 ** (-at(Cf, 40) / 10) - 10 ** (-at(If_, 40) / 10) - 10 ** (-at(RLf, 40) / 10)), 4),
       fmt(at(ILf, 40), 4))))
if has_fields:
    H("7.2  Field distribution", 2)
    FD = np.load(RES / "fields" / "fields_f0.npz")
    xf, zf, Ef = FD["Ef_guide2_x_mm"], FD["Ef_guide2_z_mm"], np.abs(FD["Ef_guide2_Ey"])
    iz = int(np.argmin(np.abs(zf - a / 2)))
    lvl3 = 20 * np.log10(Ef[(xf > -12) & (xf < -3), iz])
    lvl4 = 20 * np.log10(Ef[(xf > 9) & (xf < 17), iz])
    fsum = json.loads((RES / "fields" / "summary.json").read_text())
    A(P("Fig. 6 shows the magnitude of the vertical electric field at 40 GHz in the mid-plane of each guide (run with "
        "h = 0.035 mm). In guide 2 the coupled wave leaves the junction towards port 3 (−x); along the centre line "
        "of the guide the field towards port 4 is %s dB lower, the same as the directivity of this run "
        "(%s dB). Along the feed to port 3 the field magnitude varies by only %s dB, the signature of a single "
        "travelling TE<sub>10</sub> wave with no standing-wave component."
        % (fmt(lvl3.mean() - lvl4.mean(), 1), fmt(fsum["directivity_db_at_f0"], 1), fmt(np.ptp(lvl3), 2))))
    A(figure(FIGDIR / "fig_fields.png", "Figure 6. |E<sub>y</sub>| at 40 GHz for excitation at port 1, in dB relative "
             "to the incident peak, on the mid-planes of (a) guide 1 and (b) guide 2, each drawn along its own axis of "
             "propagation; light blue: slot outlines. PML regions omitted (h = 0.035 mm run)."))
H("7.3  Conductor loss", 2)
A(P("Copper walls (σ = 5.8 × 10<sup>7</sup> S/m, smooth) add the TE<sub>10</sub> conductor loss of %s dB/m at "
    "40 GHz [4, Sec. 3.3] along the through path, e.g. %s dB for a 30 mm long part. Surface roughness typically "
    "raises this by tens of percent at Q band; plating or polishing the guide interiors keeps it small. Over the "
    "few millimetres of the junction itself the conductor loss is of order 0.005 dB, so it does not change the "
    "coupling level measurably."
    % (fmt(r["copper_loss_db_per_m"], 3), fmt(r["copper_loss_db_per_m"] * 0.03, 3))))

# --- 8 chart and tolerances --------------------------------------------------------
H("8  Design chart and manufacturing tolerances")
H("8.1  Choosing the coupling level", 2)
A(P("The slot length sets the coupling level. The chart (Fig. 7, Table 7) was computed at h = 0.035 mm and "
    "corrected to the mesh-converged level by the shift found at L = 2.00 mm (%s dB). The shift was checked at two "
    "other slot lengths on the h = 0.025 mm mesh: the change between h = 0.035 and 0.025 mm is %s, so the "
    "correction is scaled with L accordingly, and the chart uncertainty is ± %s dB."
    % (sgn(shift035, 2), ", ".join("%s dB at L = %s mm" % (sgn(s, 2), fmt(L_, 1)) for L_, s in shift_check),
       fmt(U_chart, 2))))
A(figure(FIGDIR / "fig_chart.png", "Figure 7. Coupling at 40 GHz versus slot arm length. Grey: simulated at "
         "h = 0.035 mm; blue: mesh-corrected (monotone cubic interpolation); dashed: closed-form estimate of "
         "Sec. 3.5. Diamonds: slot length for standard coupling levels."))
crow = [["L (mm)", "L / (λ<sub>0</sub>/2)", "Coupling (dB)", "Directivity (dB)", "Return loss (dB)", "Estimate (dB)"]]
for x, cc, ee in zip(sl, Ccorr, est):
    crow.append([fmt(x["L_mm"], 2), fmt(x["L_mm"] / r["slot_resonant_length_mm"], 3), fmt(cc, 2),
                 fmt(x["directivity_db_at_f0"], 1), fmt(x["return_loss_db_at_f0"], 1), fmt(ee, 2)])
A(KeepTogether([table(crow, [0.13, 0.17, 0.17, 0.18, 0.18, 0.17], align_right_cols=(0, 1, 2, 3, 4, 5)),
                P("Table 7. Design chart at 40 GHz (coupling mesh-corrected; directivity and return loss at "
                  "h = 0.035 mm, where the directivity at L = 2.00 mm is %s dB above its converged value). Slot "
                  "width 0.35 mm, wall 0.25 mm, slot positions as in Table 2."
                  % fmt(h035["directivity_db_at_f0"] - D_conv, 1), "cap")]))
dir_L = {round(x["L_mm"], 2): x["directivity_db_at_f0"] for x in sl}
A(P("Interpolating the chart: %s. Fig. 8 shows how the frequency response changes with L: longer slots couple more "
    "strongly and vary more across the band because they approach resonance (L = λ<sub>0</sub>/2 = %s mm). The "
    "directivity at 40 GHz falls with increasing L above 2.0 mm (%s dB at L = 2.0 mm, %s dB at %s mm) as the "
    "electric leakage of Sec. 3.3 grows; for shorter slots it is set by residual couplings to port 4 that partly "
    "cancel and varies non-monotonically (%s). After choosing L, set it with "
    "<font face='Mono'>python -m moreno.design --L &lt;value&gt;</font> and re-run the simulation to confirm."
    % ("; ".join("%g dB → L = %s mm" % (C, fmt(L_, 3)) for C, L_ in targets if L_ is not None),
       fmt(r["slot_resonant_length_mm"], 3), fmt(dir_L[2.0], 1), fmt(dir_L[round(Ls.max(), 2)], 1),
       fmt(Ls.max(), 1), ", ".join("%s dB at %s mm" % (fmt(dir_L[k_], 1), fmt(k_, 1)) for k_ in sorted(dir_L) if k_ < 2.0))))
for C, L_ in targets:
    if L_ is None:
        A(P("A %g dB coupler lies just outside the chart (%s dB at L = %s mm): it needs L slightly above %s mm, where "
            "the directivity of this two-slot layout is down to about %s dB. For tighter coupling with good "
            "directivity, use more apertures along the diagonal or the multi-aperture layouts of [1, 2] rather "
            "than longer slots."
            % (C, fmt(Ccorr[-1], 2), fmt(Ls[-1], 1), fmt(Ls[-1], 1), fmt(dir_L[round(Ls.max(), 2)], 0))))
from reportlab.platypus import Image as RLImage, Table as RLTable  # noqa: E402
from layout import TEXT_W  # noqa: E402
from PIL import Image as PILImage  # noqa: E402


def img(path, width):
    w_, h_ = PILImage.open(path).size
    return RLImage(str(path), width=width, height=width * h_ / w_)


A(KeepTogether([RLTable([[img(FIGDIR / "fig_family_c.png", TEXT_W / 2 - 4), img(FIGDIR / "fig_family_d.png", TEXT_W / 2 - 4)]],
                        colWidths=[TEXT_W / 2] * 2),
                P("Figure 8. Coupling (mesh-corrected, left) and directivity (h = 0.035 mm, right; about %s dB "
                  "optimistic at L = 2.00 mm, Table 7) versus frequency for the slot lengths of Table 7."
                  % fmt(h035["directivity_db_at_f0"] - D_conv, 1), "cap")]))

H("8.2  Tolerances", 2)
trows = [["Case (one parameter changed)", "Coupling (dB)", "Δ coupling (dB)", "Directivity (dB)", "Return loss (dB)"]]
names = {"c035": "nominal (L 2.00, w 0.35, t 0.25 mm)", "t020": "wall t = 0.20 mm", "t030": "wall t = 0.30 mm",
         "w030": "slot width w = 0.30 mm", "w040": "slot width w = 0.40 mm"}
for k_ in ("c035", "t020", "t030", "w030", "w040"):
    x = tol[k_]
    trows.append([names[k_], fmt(x["coupling_db_at_f0"], 2), sgn(x["coupling_db_at_f0"] - tol["c035"]["coupling_db_at_f0"], 2),
                  fmt(x["directivity_db_at_f0"], 1), fmt(x["return_loss_db_at_f0"], 1)])
A(KeepTogether([table(trows, [0.38, 0.15, 0.16, 0.16, 0.15], align_right_cols=(1, 2, 3, 4)),
                P("Table 8. Tolerance runs at 40 GHz (h = 0.035 mm; differences are not affected by the mesh shift).",
                  "cap")]))
A(P("Per 0.01 mm, the coupling changes by %s dB with slot length, %s dB with slot width and %s dB with wall "
    "thickness (negative: stronger coupling). With a typical machining tolerance of ± 0.01 mm on each, applied "
    "independently and combined in quadrature, the coupling is held to ± %s dB; ± 0.005 mm halves this. "
    "In the ± 0.05 mm cases of Table 8 the directivity stays above %s dB on the h = 0.035 mm mesh; allowing for "
    "the %s dB by which that mesh overstates the directivity (Table 7), the worst case is about %s dB. The slot "
    "positions matter less: the directivity "
    "condition x<sub>0</sub> = z<sub>0</sub> is first-order insensitive to a common shift along the diagonal, and a "
    "0.01 mm error in Δ changes βΔ by only %s°."
    % (sgn(dCdL01, 2), sgn(dCdw, 2), sgn(dCdt, 2), fmt(math.sqrt(dCdL01 ** 2 + dCdw ** 2 + dCdt ** 2), 2),
       fmt(min(x["directivity_db_at_f0"] for x in tol.values()), 1), fmt(h035["directivity_db_at_f0"] - D_conv, 1),
       fmt(min(x["directivity_db_at_f0"] for x in tol.values()) - (h035["directivity_db_at_f0"] - D_conv), 1),
       fmt(math.degrees(r["beta_rad_per_mm"] * 0.01), 2))))
H("8.3  Manufacturing notes", 2)
E(bullets([
    "Make the common wall a separate plate (t = %s mm) between the two guide blocks; wire EDM or laser cutting "
    "holds the slot outline to a few micrometres, and the plate can be lapped flat." % fmt(t, 2),
    "Milled slots have rounded ends. With full-radius ends (r = w/2) each end loses 2(1 − π/4)r<sup>2</sup> of "
    "area, which by area equivalence shortens the arm by about 0.11 w per end, %s mm in total here; from the chart "
    "slope that weakens the coupling by roughly %s dB. Either cut square ends (wire EDM) or model the radius and "
    "lengthen L to compensate." % (fmt(2 * (1 - math.pi / 4) * w0 / 2, 3), fmt(abs(dCdL) * 2 * (1 - math.pi / 4) * w0 / 2, 2)),
    "Align the two guide blocks to the plate with dowel pins; a lateral offset of the upper guide shifts both "
    "slots relative to its walls and couples into the directivity condition.",
    "Use standard UG-383/U (WR-22) flanges on all four ports; terminate port 4 in a matched load.",
]))

# --- 9 HFSS ---------------------------------------------------------------------------
H("9  HFSS model")
A(P("Two equivalent HFSS builders are provided; both read design/design.json and create the same parametric model "
    "(design variables a, b, t, Ls, ws, Delta, d = (a − Delta)/2, Lext)."))
E(bullets([
    "<b>hfss/build_pyaedt.py</b> runs from a normal Python prompt through PyAEDT, builds, solves, exports "
    "results/moreno.s4p and runs the same post-processing as the openEMS model. All PyAEDT calls were checked "
    "against the PyAEDT 1.9.0 signatures.",
    "<b>hfss/build_in_aedt.py</b> runs inside AEDT (Automation › Run Script) with no installation and tries every "
    "solution-type name AEDT has used for Driven Modal designs.",
    "<b>Model.</b> Vacuum solids for the two guides and the two slots; every face touching the background is a "
    "perfect conductor by default. Wave ports with one mode, <b>not renormalised</b> (the generalised S-matrix "
    "refers to each port's TE<sub>10</sub> impedance, like the openEMS model; renormalising a waveguide port to "
    "50 Ω would distort every |S|).",
    "<b>Solver settings.</b> Adaptive meshing at 40 GHz with ΔS = 0.002 and two converged passes, because "
    "|S<sub>31</sub>| ≈ 0.03 cannot be resolved with the common ΔS = 0.02; a length-based mesh seed of w/5 in the "
    "slots; an interpolating sweep 35–45 GHz. Isolation near −65 dB is below what a ΔS criterion controls, so "
    "the HFSS directivity should be confirmed with a discrete point at 40 GHz on the final mesh.",
    "<b>Status.</b> The scripts could not be run in this project because the AEDT Student 2025 R2 installation "
    "rejects any HFSS design with the licence error 'No such feature exists: hfss_gui' (FlexNet −5,357). This "
    "is a known fault of that Student package; a current Student release or a full licence is needed.",
]))
A(P("The HFSS result should agree with Table 6 to within the mesh uncertainties of the two solvers; a difference "
    "of more than about 0.5 dB in coupling would point to a modelling difference worth investigating."))

# --- 10 limitations ------------------------------------------------------------------
H("10  Limitations and recommendations")
E(bullets([
    "The walls are perfect conductors; copper loss is estimated analytically (Sec. 7.3), not simulated.",
    "Slot ends are square. Model the cutter radius of the actual process before freezing L (Sec. 8.3).",
    "The coupling level is known to ± %s dB from the numerical study (Sec. 6.2–6.3); a second solver (HFSS) "
    "would give an independent check." % fmt(U_conv, 2),
    "Flanges, the matched load on port 4 and any transitions are not modelled; their mismatch adds to the "
    "directivity error of a measured coupler and usually dominates at this directivity level.",
    "For a coupler with flatter coupling over the band, use more than two slots with a binomial or Chebyshev "
    "distribution along the diagonal, or the multi-hole cross-guide designs of the literature [1, 2].",
]))

# --- 11 reproducing --------------------------------------------------------------------
H("11  Reproducing the results")
A(P("All files are in the repository; every number in this report is recomputed from them by "
    "docs/report/build_report.py.", "bodyl"))
A(P("pip install -r requirements.txt<br/>python tools/verify_results.py                     # re-check every number here<br/>"
    "bash tools/run_campaign.sh                         # re-run all 24 simulations (~2 h)<br/>"
    "python docs/report/collect.py campaign             # gather them into results/<br/>"
    "python -m moreno.design                            # design/design.json<br/>"
    "python openems/simulate_openems.py --fine 0.0125   # reference run (results/openems/)<br/>"
    "python openems/simulate_openems.py --excite all    # full scattering matrix<br/>"
    "python openems/sweep_slot_length.py --L 1.6 2.0 2.4 --target 25<br/>"
    "python hfss/build_pyaedt.py --version 2026.1 --student --analyze<br/>"
    "python docs/report/build_report.py                 # this report<br/>pytest -q", "code"))
A(P("The test suite checks the design equations, the coupling law (2) against explicit mode fields, the Touchstone "
    "reader, and the committed results (power balance, reciprocity, symmetry, port assignment, the convergence "
    "bracket). tools/verify_results.py recomputes every number of this report from the raw result of each run "
    "(results/openems/studies/runs) without using the report code.", "bodyl"))

# --- references -----------------------------------------------------------------------
H("References")
for i, ref_ in enumerate([
    "T. Moreno, <i>Microwave Transmission Design Data</i>. New York: McGraw-Hill, 1948 (Dover reprint, 1958).",
    "H. J. Riblet and T. S. Saad, “A new type of waveguide directional coupler,” <i>Proc. IRE</i>, vol. 36, "
    "no. 1, pp. 61–64, Jan. 1948.",
    "H. A. Bethe, “Theory of diffraction by small holes,” <i>Phys. Rev.</i>, vol. 66, no. 7–8, pp. 163–182, 1944.",
    "D. M. Pozar, <i>Microwave Engineering</i>, 4th ed. Hoboken, NJ: Wiley, 2012, Secs. 3.3, 4.8 and 7.4.",
    "R. E. Collin, <i>Foundations for Microwave Engineering</i>, 2nd ed. New York: McGraw-Hill, 1992 "
    "(IEEE Press reprint, 2001), small-aperture coupling.",
    "S. B. Cohn, “Microwave coupling by large apertures,” <i>Proc. IRE</i>, vol. 40, no. 6, pp. 696–699, 1952.",
    "T. Liebig et al., openEMS – a free and open electromagnetic field solver using the EC-FDTD method, "
    "https://openems.de (source used: github.com/thliebig/openEMS, built Nov. 2026).",
    "A. Taflove and S. C. Hagness, <i>Computational Electrodynamics: The Finite-Difference Time-Domain Method</i>, "
    "3rd ed. Norwood, MA: Artech House, 2005.",
    "P. J. Roache, “Perspective: a method for uniform reporting of grid refinement studies,” <i>J. Fluids "
    "Eng.</i>, vol. 116, no. 3, pp. 405–413, 1994.",
], 1):
    A(P("[%d]  %s" % (i, ref_), "ref"))

# --- appendix A -----------------------------------------------------------------------
H("Appendix A  Corrections made during this review")
A(P("Every file of the project was re-checked for this report. The following errors were found and corrected; "
    "results given earlier in the project that depend on them are withdrawn.", "bodyl"))
A(table([
    ["Item", "Was", "Now", "Effect"],
    ["openEMS geometry", "air carved out of metal with higher priority", "metal-only boxes on vacuum",
     "every wall was shifted one cell into the metal; coupling 24.3 dB → %s dB, directivity 38.6 → %s dB"
     % (fmt(C_conv, 2), fmt(D_conv, 1))],
    ["openEMS mesh at the slots", "h = 0.083 mm, lines on edges", "h ≤ 0.035 mm, 1/3–2/3 rule, convergence study",
     "coupling now known to ± %s dB" % fmt(U_conv, 2)],
    ["Wall clearance check", "d − L/2 − w/2 = 0.425 mm", "d − L/2 = %s mm" % fmt(r["clearance_mm"], 3),
     "the arm end, not its side, is nearest the wall"],
    ["Separation of the crosses", "not checked", "exact shortest distance, %s mm" % fmt(r["slot_gap_mm"], 3),
     "new design check"],
    ["Closed-form frequency trend", "array factor and wall loss only", "plus slot resonance (Cohn)",
     "trend now follows the simulation in sign and order of magnitude"],
    ["PyAEDT wave ports", "renormalised to 50 Ω (PyAEDT default)", "not renormalised",
     "|S| from HFSS would have been wrong"],
    ["HFSS convergence", "ΔS = 0.02, 1 converged pass", "ΔS = 0.002, 2 passes", "coupling of −30 dB now resolved"],
    ["Touchstone header", "'R 50' without comment", "normalisation stated in the file", "no misreading as 50 Ω data"],
    ["Symmetry statement", "|S<sub>21</sub>| = |S<sub>43</sub>| from rotation", "from the x ↔ z mirror",
     "correct reasoning; relations verified (Table 3)"],
], [0.2, 0.25, 0.25, 0.30]))

# --- appendix B -----------------------------------------------------------------------
H("Appendix B  Nomenclature")
A(table([
    ["Symbol", "Meaning"],
    ["a, b", "broad and narrow inner dimensions of the waveguide"],
    ["t", "thickness of the common wall"],
    ["L, w", "length and width of each slot arm"],
    ["d", "offset of the slot centres from the side walls, slot 1 at (d, d), slot 2 at (a − d, a − d)"],
    ["Δ", "slot spacing per axis, a − 2d"],
    ["f<sub>c</sub>, λ<sub>0</sub>, λ<sub>g</sub>, β", "TE<sub>10</sub> cut-off, free-space and guide wavelength, phase constant"],
    ["Z<sub>TE</sub>, η<sub>0</sub>", "TE<sub>10</sub> wave impedance, free-space impedance"],
    ["α<sub>m</sub>, α<sub>e</sub>", "magnetic and electric polarizability of a slot arm"],
    ["h", "FDTD cell size across the slots and through the wall"],
    ["C, D, I", "coupling −20 log|S<sub>31</sub>|, directivity I − C, isolation −20 log|S<sub>41</sub>|"],
], [0.22, 0.78]))

doc = Report(str(OUTPDF), "Moreno Cross-Guide Directional Coupler – WR-22, 40 GHz",
             "Moreno cross-guide coupler · WR-22 · 40 GHz · Technical report v%s" % VERSION, AUTHOR, DATE)
doc.build(st)
print("Wrote %s" % OUTPDF)
