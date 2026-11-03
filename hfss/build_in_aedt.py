# -*- coding: utf-8 -*-
"""Build the Moreno coupler from INSIDE AEDT (no pip installs needed).

AEDT 2025 R2 or 2026 R1 (incl. Student):  Automation > Run Script... > this file.

Reads ../design/design.json when it can find it; otherwise uses the
built-in WR-22 / 40 GHz values below (identical to moreno/design.py defaults).
Saves Moreno_coupler.aedt in your AEDT project directory.
"""
import json
import math
import os

# ---------------- desktop handle ----------------
try:
    oDesktop  # noqa: F821  (injected by Run Script)
except NameError:
    import ScriptEnv
    ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")


def msg(text, level=0):
    oDesktop.AddMessage("", "", level, "[Moreno] " + text)


# ---------------- design values ----------------
def default_design():
    c0, f0, a = 299.792458, 40.0, 5.690
    lamg = (c0 / f0) / math.sqrt(1.0 - (c0 / (2 * a) / f0) ** 2)
    sp = lamg / 4.0
    return {"a_mm": a, "b_mm": 2.845, "wall_mm": 0.25, "slot_length_mm": 2.0,
            "slot_width_mm": 0.35, "spacing_mm": sp, "d_mm": (a - sp) / 2,
            "feed_mm": round(lamg, 3), "f0_ghz": f0, "sweep_ghz": [35.0, 45.0],
            "lambda_g_mm": lamg}


r = None
try:
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "..", "design", "design.json")
    if os.path.exists(path):
        with open(path) as fh:
            r = json.load(fh)
        msg("Loaded %s" % os.path.normpath(path))
except Exception:
    r = None
if r is None:
    r = default_design()
    msg("Using built-in WR-22 / 40 GHz values")

a, b, t = r["a_mm"], r["b_mm"], r["wall_mm"]
Lext = r["feed_mm"]

# ---------------- project / design ----------------
oDesktop.RestoreWindow()
oProject = oDesktop.NewProject()
design_name = "Moreno_coupler"


def get_design(name):
    try:
        return oProject.GetDesign(name)
    except Exception:
        return None


oDesign = None
for sol in ["HFSS Modal Network", "Modal", "DrivenModal", ""]:
    try:
        oProject.InsertDesign("HFSS", design_name, sol, "")
    except Exception as e:
        msg("InsertDesign('%s') raised: %s" % (sol, e), 1)
    oDesign = get_design(design_name)
    if oDesign is not None:
        msg("Design inserted with solution type '%s'" % sol)
        break
if oDesign is None:
    msg("Could not insert an HFSS design. If Project > Insert HFSS Design also fails by hand, "
        "the HFSS licence in this install is the problem, not the script.", 2)
    raise SystemExit("HFSS design insertion failed")

oDesign = oProject.SetActiveDesign(design_name)
for sol in ["HFSS Modal Network", "Modal", "DrivenModal"]:
    try:
        oDesign.SetSolutionType(sol, ["NAME:Options", "EnableAutoOpen:=", False])
        break
    except Exception:
        pass

oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.SetModelUnits(["NAME:Units Parameter", "Units:=", "mm", "Rescale:=", False])


def addvar(name, value):
    oDesign.ChangeProperty(
        ["NAME:AllTabs",
         ["NAME:LocalVariableTab",
          ["NAME:PropServers", "LocalVariables"],
          ["NAME:NewProps",
           ["NAME:" + name, "PropType:=", "VariableProp", "UserDef:=", True, "Value:=", value]]]])


for name, val in [("a", a), ("b", b), ("t", t), ("Ls", r["slot_length_mm"]),
                  ("ws", r["slot_width_mm"]), ("Delta", r["spacing_mm"]), ("Lext", Lext)]:
    addvar(name, "%.4fmm" % val)
addvar("d", "(a-Delta)/2")


def box(name, pos, size):
    oEditor.CreateBox(
        ["NAME:BoxParameters",
         "XPosition:=", pos[0], "YPosition:=", pos[1], "ZPosition:=", pos[2],
         "XSize:=", size[0], "YSize:=", size[1], "ZSize:=", size[2]],
        ["NAME:Attributes", "Name:=", name, "Flags:=", "", "Color:=", "(132 170 210)",
         "Transparency:=", 0.7, "PartCoordinateSystem:=", "Global", "UDMId:=", "",
         "MaterialValue:=", "\"vacuum\"", "SurfaceMaterialValue:=", "\"\"",
         "SolveInside:=", True, "ShellElement:=", False, "ShellElementThickness:=", "0mm",
         "IsMaterialEditable:=", True, "UseMaterialAppearance:=", False, "IsLightweight:=", False])


box("Air", ["0mm", "0mm", "-Lext"], ["a", "b", "a+2*Lext"])
box("G2", ["-Lext", "b+t", "0mm"], ["a+2*Lext", "b", "a"])
box("S1x", ["d-Ls/2", "b", "d-ws/2"], ["Ls", "t", "ws"])
box("S1z", ["d-ws/2", "b", "d-Ls/2"], ["ws", "t", "Ls"])
box("S2x", ["a-d-Ls/2", "b", "a-d-ws/2"], ["Ls", "t", "ws"])
box("S2z", ["a-d-ws/2", "b", "a-d-Ls/2"], ["ws", "t", "Ls"])
# Unite each slot's two arms; guides and slots only touch, so they stay separate
# solids and the slots can carry their own mesh operation.
for arms in ("S1x,S1z", "S2x,S2z"):
    oEditor.Unite(["NAME:Selections", "Selections:=", arms],
                  ["NAME:UniteParameters", "KeepOriginals:=", False, "TurnOnNBodyBoolean:=", True])
oDesign.GetModule("MeshSetup").AssignLengthOp(
    ["NAME:SlotMesh", "RefineInside:=", True, "Enabled:=", True, "Objects:=", ["S1x", "S2x"],
     "RestrictElem:=", False, "NumMaxElem:=", "1000", "RestrictLength:=", True,
     "MaxLength:=", "%.4fmm" % (r["slot_width_mm"] / 5.0)])


def face_at(body, x, y, z):
    return int(oEditor.GetFaceByPosition(
        ["NAME:FaceParameters", "BodyName:=", body,
         "XPosition:=", "%.5fmm" % x, "YPosition:=", "%.5fmm" % y, "ZPosition:=", "%.5fmm" % z]))


y2 = b + t + b / 2.0
ports = [("1", face_at("Air", a / 2, b / 2, -Lext)), ("2", face_at("Air", a / 2, b / 2, a + Lext)),
         ("3", face_at("G2", -Lext, y2, a / 2)), ("4", face_at("G2", a + Lext, y2, a / 2))]
oBnd = oDesign.GetModule("BoundarySetup")
for name, fid in ports:
    try:
        oBnd.AssignWavePort(
            ["NAME:" + name, "Faces:=", [fid], "NumModes:=", 1,
             "UseLineModeAlignment:=", False, "DoDeembed:=", False, "RenormalizeAllTerminals:=", True,
             ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", False, "CharImp:=", "Zpi"]],
             "ShowReporterFilter:=", False, "ReporterFilter:=", [True], "UseAnalyticAlignment:=", False])
    except Exception as e:
        msg("Wave port %s failed on face %s: %s" % (name, fid, e), 2)
        raise

oAna = oDesign.GetModule("AnalysisSetup")
oAna.InsertSetup("HfssDriven",
                 ["NAME:Setup1", "SolveType:=", "Single", "Frequency:=", "%.4fGHz" % r["f0_ghz"],
                  "MaxDeltaS:=", 0.002, "MaximumPasses:=", 20, "MinimumPasses:=", 3,
                  "MinimumConvergedPasses:=", 2, "PercentRefinement:=", 30,
                  "IsEnabled:=", True, "BasisOrder:=", 1])
lo, hi = r["sweep_ghz"]
oAna.InsertFrequencySweep("Setup1",
                          ["NAME:Sweep", "IsEnabled:=", True, "RangeType:=", "LinearCount",
                           "RangeStart:=", "%.4fGHz" % lo, "RangeEnd:=", "%.4fGHz" % hi, "RangeCount:=", 101,
                           "Type:=", "Interpolating", "SaveFields:=", False, "SaveRadFields:=", False,
                           "InterpTolerance:=", 0.5, "InterpMaxSolns:=", 50, "InterpMinSolns:=", 0,
                           "InterpMinSubranges:=", 1, "ExtrapToDC:=", False])
try:
    oDesign.GetModule("ReportSetup").CreateReport(
        "S Parameters", "Modal Solution Data", "Rectangular Plot", "Setup1 : Sweep",
        ["Domain:=", "Sweep"], ["Freq:=", ["All"]],
        ["X Component:=", "Freq",
         "Y Component:=", ["dB(S(1,1))", "dB(S(2,1))", "dB(S(3,1))", "dB(S(4,1))"]])
except Exception:
    msg("S-parameter plot not created; add it via Results > Create Modal Solution Data Report.", 1)

proj = oDesktop.GetProjectDirectory().rstrip("/\\") + "/Moreno_coupler.aedt"
oProject.SaveAs(proj, True)
msg("Built: lambda_g = %.3f mm, Delta = %.3f mm, d = %.3f mm. Saved to %s"
    % (r.get("lambda_g_mm", 0.0), r["spacing_mm"], r["d_mm"], proj))
