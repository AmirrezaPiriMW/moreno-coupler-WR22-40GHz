"""Independent re-derivation of every headline number of the report from the raw results.

    python tools/verify_results.py                   # uses results/openems/studies/runs
    python tools/verify_results.py <campaign_dir>    # or a fresh campaign (tools/run_campaign.sh)

Uses only numpy and the raw .s4p / run_info.json of each run; it does not import the
report code, so it checks docs/report/key_numbers.json and the study tables independently.
Exits non-zero if anything disagrees.
"""
import csv
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

R = Path(__file__).resolve().parents[1]
C = Path(sys.argv[1]) if len(sys.argv) > 1 else R / "results" / "openems" / "studies" / "runs"
fails = []


def check(name, got, want, tol):
    ok = abs(got - want) <= tol
    print("%-58s %12.5f %12.5f %s" % (name, got, want, "ok" if ok else "MISMATCH"))
    if not ok:
        fails.append(name)


def s4p(path):
    nums, unit = [], 1e9
    for line in Path(path).read_text().splitlines():
        line = line.split("!")[0].strip()
        if not line or line.startswith("#"):
            continue
        nums += [float(v) for v in line.split()]
    a = np.array(nums).reshape(-1, 33)
    f = a[:, 0] * unit
    mag, ang = a[:, 1::2], a[:, 2::2]
    S = (mag * np.exp(1j * np.radians(ang))).reshape(-1, 4, 4)
    return f, S


def metrics(name, f0=40e9):
    f, S = s4p(C / name / "moreno_openems.s4p")
    k = np.argmin(abs(f - f0))
    assert abs(f[k] - f0) < 1, "40 GHz is not a sample point"
    d = lambda x: -20 * np.log10(abs(x))  # noqa: E731
    c, i, rl, il = d(S[k, 2, 0]), d(S[k, 3, 0]), d(S[k, 0, 0]), d(S[k, 1, 0])
    return dict(C=c, D=i - c, I=i, RL=rl, IL=il, f=f, S=S)


info = lambda n: json.loads((C / n / "run_info.json").read_text())  # noqa: E731

print("== run parameters really are what the report says")
expect = {"c050": (0.05, 2.0, .35, .25, True), "c035": (0.035, 2.0, .35, .25, True), "c025": (0.025, 2.0, .35, .25, True),
          "c0175": (0.0175, 2.0, .35, .25, True), "c0125": (0.0125, 2.0, .35, .25, True),
          "s035": (0.035, 2.0, .35, .25, False), "s025": (0.025, 2.0, .35, .25, False), "s0175": (0.0175, 2.0, .35, .25, False),
          "L1.6": (0.035, 1.6, .35, .25, True), "L1.8": (0.035, 1.8, .35, .25, True), "L2.2": (0.035, 2.2, .35, .25, True),
          "L2.4": (0.035, 2.4, .35, .25, True), "L2.6": (0.035, 2.6, .35, .25, True),
          "t020": (0.035, 2.0, .35, .20, True), "t030": (0.035, 2.0, .35, .30, True),
          "w030": (0.035, 2.0, .30, .25, True), "w040": (0.035, 2.0, .40, .25, True),
          "o025_L1.6": (0.025, 1.6, .35, .25, True), "o025_L2.4": (0.025, 2.4, .35, .25, True),
          "n_res30": (0.035, 2.0, .35, .25, True), "n_feed25": (0.035, 2.0, .35, .25, True), "n_e4": (0.035, 2.0, .35, .25, True),
          "full": (0.035, 2.0, .35, .25, True), "fields": (0.035, 2.0, .35, .25, True)}
for n, (h, L, w, t, thirds) in expect.items():
    I = info(n)
    run, d = I["runs"][0], I["design"]
    good = (abs(run["fine_mm"] - h) < 1e-9 and abs(d["slot_length_mm"] - L) < 1e-9 and abs(d["slot_width_mm"] - w) < 1e-9
            and abs(d["wall_mm"] - t) < 1e-9 and run["thirds_rule"] == thirds and abs(d["f0_ghz"] - 40) < 1e-9)
    extra = {"n_res30": run["res_cells_per_lambda"] == 30, "n_feed25": run["feed_lambda_g"] == 2.5,
             "n_e4": run["end_criteria"] == 1e-4}.get(n, run["res_cells_per_lambda"] == 20 and run["feed_lambda_g"] == 1.5
                                                       and run["end_criteria"] == 1e-5)
    # d and spacing must follow lambda_g/4 for the design's f0
    sp_ok = abs(d["spacing_mm"] - 2.48997) < 1e-4 and abs(d["d_mm"] - (5.69 - d["spacing_mm"]) / 2) < 1e-9
    print("%-10s h=%-7g L=%-4g w=%-4g t=%-4g thirds=%-5s %s" % (n, h, L, w, t, thirds,
          "ok" if good and extra and sp_ok else "MISMATCH"))
    if not (good and extra and sp_ok):
        fails.append("params " + n)

print("\n== every committed study row equals the raw .s4p value")
for g in ["convergence", "staircase", "numerics", "slot_length", "offset_check", "tolerance"]:
    for row in csv.DictReader(open(R / ("results/openems/studies/%s.csv" % g))):
        m = metrics(row["run"])
        for key, mk in (("coupling_db_at_f0", "C"), ("directivity_db_at_f0", "D"), ("return_loss_db_at_f0", "RL")):
            if abs(float(row[key]) - m[mk]) > 1e-6:
                fails.append("%s %s %s" % (g, row["run"], key))
print("study tables vs raw files: %s" % ("ok" if not any(x.startswith(("conv", "stair", "num", "slot", "off", "tol")) for x in fails) else "MISMATCH"))

print("\n== convergence bracket")
th = [metrics(n)["C"] for n in ["c050", "c035", "c025", "c0175", "c0125"]]
st = [metrics(n)["C"] for n in ["s035", "s025", "s0175"]]
print("thirds", np.round(th, 4), "monotone up", all(np.diff(th) > 0))
print("lines-on-edges", np.round(st, 4), "monotone down", all(np.diff(st) < 0))
lo, hi = th[-1], st[-1]
Cc, U = (lo + hi) / 2, (hi - lo) / 2
key = json.loads((R / "docs/report/key_numbers.json").read_text())
check("converged coupling (bracket centre)", key["coupling_extrapolated_db"], Cc, 1e-9)
check("slot-mesh uncertainty (half width)", key["coupling_slot_mesh_uncertainty_db"], U, 1e-9)
dn = [metrics(n)["C"] - metrics("c035")["C"] for n in ("n_res30", "n_feed25", "n_e4")]
check("total uncertainty (RSS with numerics)", key["coupling_uncertainty_db"], math.sqrt(U ** 2 + sum(v * v for v in dn)), 1e-9)
Dth = [metrics(n)["D"] for n in ["c050", "c035", "c025", "c0175", "c0125"]]
Dst = [metrics(n)["D"] for n in ["s035", "s025", "s0175"]]
print("directivity thirds", np.round(Dth, 3), "edges", np.round(Dst, 3))
check("directivity bracket centre", key["directivity_converged_db"], (Dth[-1] + Dst[-1]) / 2, 1e-9)
check("directivity half width", key["directivity_uncertainty_db"], (Dth[-1] - Dst[-1]) / 2, 1e-9)

print("\n== chart and targets (independent: linear correction + numpy cubic through points)")
Ls = np.array([1.6, 1.8, 2.0, 2.2, 2.4, 2.6])
C35 = np.array([metrics(n)["C"] for n in ["L1.6", "L1.8", "c035", "L2.2", "L2.4", "L2.6"]])
shift = Cc - metrics("c035")["C"]
sh_check = {1.6: metrics("o025_L1.6")["C"] - C35[0], 2.0: metrics("c025")["C"] - C35[2], 2.4: metrics("o025_L2.4")["C"] - C35[4]}
print("0.035->0.025 shift by L:", {k: round(v, 4) for k, v in sh_check.items()})
ratio = np.interp(Ls, [1.6, 2.0, 2.4], [sh_check[1.6] / sh_check[2.0], 1.0, sh_check[2.4] / sh_check[2.0]])
Ccorr = C35 + shift * ratio
for c, row in zip(Ccorr, key["chart"]):
    check("chart L=%.1f" % row["L_mm"], row["C_corrected"], c, 1e-9)
poly = np.polyfit(Ls, Ccorr, 3)
for Ct, Lk in key["targets"]:
    if Lk is None:
        print("target %g dB outside chart; chart range %.2f..%.2f" % (Ct, Ccorr.min(), Ccorr.max()))
        continue
    roots = [r_.real for r_ in np.roots(poly - np.array([0, 0, 0, Ct])) if abs(r_.imag) < 1e-9 and 1.6 <= r_.real <= 2.6]
    check("L for %g dB (cubic fit vs PCHIP)" % Ct, Lk, roots[0], 0.003)
dCdL = np.polyval(np.polyder(poly), 2.0) * 0.01
check("dC/dL per 0.01 mm", key["dC_dL_db_per_0.01mm"], dCdL, 0.003)
check("dC/dw per 0.01 mm", key["dC_dw_db_per_0.01mm"], (metrics("w040")["C"] - metrics("w030")["C"]) / 10, 1e-9)
check("dC/dt per 0.01 mm", key["dC_dt_db_per_0.01mm"], (metrics("t030")["C"] - metrics("t020")["C"]) / 10, 1e-9)

print("\n== full matrix")
f, S = s4p(C / "full" / "moreno_openems.s4p")
bal = np.sum(abs(S) ** 2, axis=1)
check("power balance min", key["power_balance_full"][0], bal.min(), 1e-9)
check("power balance max", key["power_balance_full"][1], bal.max(), 1e-9)
check("reciprocity", key["reciprocity"], np.max(abs(S - S.transpose(0, 2, 1))), 1e-9)
m = abs(S)
print("|S21|-|S43| %.2e  |S31|-|S42| %.2e  |S41|-|S32| %.2e  S11-S22 %.2e S33-S44 %.2e S11-S33 %.2e" % (
    np.max(abs(m[:, 1, 0] - m[:, 3, 2])), np.max(abs(m[:, 2, 0] - m[:, 3, 1])), np.max(abs(m[:, 3, 0] - m[:, 2, 1])),
    np.max(abs(m[:, 0, 0] - m[:, 1, 1])), np.max(abs(m[:, 2, 2] - m[:, 3, 3])), np.max(abs(m[:, 0, 0] - m[:, 2, 2]))))

print("\n== reference (finest) response")
ref = metrics("c0125")
f, S1 = ref["f"], ref["S"][:, :, 0]
d = lambda x: -20 * np.log10(abs(x))  # noqa: E731
Cf, If, RLf, ILf = d(S1[:, 2]), d(S1[:, 3]), d(S1[:, 0]), d(S1[:, 1])
Df = If - Cf
fg = f / 1e9
for fx in (35, 40, 45):
    k = np.argmin(abs(fg - fx))
    print("%g GHz: C %.2f D %.1f I %.1f RL %.1f IL %.3f" % (fx, Cf[k], Df[k], If[k], RLf[k], ILf[k]))
print("C range %.2f..%.2f  Dmin %.2f  Imin %.2f  RLmin %.2f  ILmax %.3f" % (Cf.min(), Cf.max(), Df.min(), If.min(), RLf.min(), ILf.max()))
ok30 = fg[Df >= 30]
ok35 = fg[Df >= 35]
print("D>=30 from %.2f to %.2f (contiguous=%s); D>=35 from %.2f to %.2f" % (ok30.min(), ok30.max(), len(ok30) == len(fg), ok35.min(), ok35.max()))
k = np.argmin(abs(fg - 40))
pw = -10 * np.log10(1 - abs(S1[k, 2]) ** 2 - abs(S1[k, 3]) ** 2 - abs(S1[k, 0]) ** 2)
print("through loss from power: %.4f dB, from S21: %.4f dB" % (pw, ILf[k]))

print("\n== design numbers")
c0 = 299.792458
a, b = 5.69, 2.845
fc = c0 / (2 * a)
lg = (c0 / 40) / math.sqrt(1 - (fc / 40) ** 2)
print("fc %.4f  lam0 %.4f  lg %.4f  beta %.4f  ZTE %.1f  Delta %.4f  d %.4f" % (
    fc, c0 / 40, lg, 2 * math.pi / lg, 376.730313 / math.sqrt(1 - (fc / 40) ** 2), lg / 4, (a - lg / 4) / 2))
print("\nFAILURES:", fails if fails else "none")
sys.exit(1 if fails else 0)
