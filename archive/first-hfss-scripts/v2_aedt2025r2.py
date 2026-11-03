# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Moreno cross-guide coupler, WR-22, f0 = 40 GHz
# Target: Ansys Electronics Desktop 2025 R2 (incl. Student), CPython scripting
#
# Run:  Automation tab (or Tools) > Run Script...  > pick this .py file
#
# Matches the design canvas:
#   WR-22  a = 5.690 mm, b = 2.845 mm, shared wall t = 0.25 mm
#   "+" crossed slots, arm 2.00 x 0.35 mm
#   lambda_g(40 GHz) = 9.960 mm, Delta = lambda_g/4 = 2.490 mm, d = 1.600 mm
#   Slot centres (x, z) = (1.600, 1.600) and (4.090, 4.090) mm
#   Port 1 input (-z), 2 through (+z), 3 coupled (-x), 4 isolated (+x)
#
# Air-only model: every face that touches the background is PEC by default.
# Feeds are 1 lambda_g long to keep the mesh small for the Student licence.
# ----------------------------------------------------------------------------
import math

# ---------------- get the desktop handle (works for CPython in 2025 R2) -----
try:
    oDesktop  # injected by Run Script
except NameError:
    import ScriptEnv
    ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")

# ---------------- design values (mm, GHz) ----------------
C0 = 299.792458
f0 = 40.0
a, b, t = 5.690, 2.845, 0.25
Ls, ws = 2.0, 0.35

fc = C0 / (2.0 * a)
lam0 = C0 / f0
lamg = lam0 / math.sqrt(1.0 - (fc / f0) ** 2)
Delta = lamg / 4.0
d = (a - Delta) / 2.0
Lext = round(lamg, 3)

# ---------------- project / design ----------------
oDesktop.RestoreWindow()
oProject = oDesktop.NewProject()
design_name = "Moreno_WR22_40GHz"

def msg(text, level=0):
    oDesktop.AddMessage("", "", level, "[Moreno] " + text)


def get_design(name):
    try:
        return oProject.GetDesign(name)
    except Exception:
        return None


# Insert the design: try each solution-type spelling, stop at the first that
# actually produces a design (a failed insert does not always raise).
oDesign = None
for sol in ["HFSS Modal Network", "Modal", "DrivenModal", ""]:
    try:
        oProject.InsertDesign("HFSS", design_name, sol, "")
    except Exception as e:
        msg("InsertDesign with solution type '%s' raised: %s" % (sol, e), 1)
    oDesign = get_design(design_name)
    if oDesign is not None:
        msg("Design inserted (solution type tried: '%s')." % sol)
        break

if oDesign is None:
    msg("Could not insert an HFSS design. Try Project > Insert HFSS Design by hand: "
        "if that also fails, the HFSS licence in this Student install is the problem, not the script.", 2)
    raise SystemExit("HFSS design insertion failed")

oDesign = oProject.SetActiveDesign(design_name)

# 2024 R1+ may create a hybrid design by default: force plain modal network
for sol in ["HFSS Modal Network", "Modal", "DrivenModal"]:
    try:
        oDesign.SetSolutionType(sol, ["NAME:Options", "EnableAutoOpen:=", False])
        msg("Solution type set to '%s'." % sol)
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


addvar("a", "%.3fmm" % a)
addvar("b", "%.3fmm" % b)
addvar("t", "%.3fmm" % t)
addvar("Ls", "%.3fmm" % Ls)
addvar("ws", "%.3fmm" % ws)
addvar("Delta", "%.3fmm" % Delta)
addvar("d", "(a-Delta)/2")
addvar("Lext", "%.3fmm" % Lext)


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


# Guide 1 (lower, along z) and guide 2 (upper, along x)
box("Air", ["0mm", "0mm", "-Lext"], ["a", "b", "a+2*Lext"])
box("G2", ["-Lext", "b+t", "0mm"], ["a+2*Lext", "b", "a"])
# Slot 1 at (d, d)
box("S1x", ["d-Ls/2", "b", "d-ws/2"], ["Ls", "t", "ws"])
box("S1z", ["d-ws/2", "b", "d-Ls/2"], ["ws", "t", "Ls"])
# Slot 2 at (a-d, a-d)
box("S2x", ["a-d-Ls/2", "b", "a-d-ws/2"], ["Ls", "t", "ws"])
box("S2z", ["a-d-ws/2", "b", "a-d-Ls/2"], ["ws", "t", "Ls"])

oEditor.Unite(["NAME:Selections", "Selections:=", "Air,G2,S1x,S1z,S2x,S2z"],
              ["NAME:UniteParameters", "KeepOriginals:=", False, "TurnOnNBodyBoolean:=", True])


# ---------------- wave ports ----------------
def face_at(x, y, z):
    return int(oEditor.GetFaceByPosition(
        ["NAME:FaceParameters", "BodyName:=", "Air",
         "XPosition:=", "%.5fmm" % x, "YPosition:=", "%.5fmm" % y, "ZPosition:=", "%.5fmm" % z]))


y2 = b + t + b / 2.0
ports = [
    ("1", face_at(a / 2.0, b / 2.0, -Lext)),       # input
    ("2", face_at(a / 2.0, b / 2.0, a + Lext)),    # through
    ("3", face_at(-Lext, y2, a / 2.0)),            # coupled
    ("4", face_at(a + Lext, y2, a / 2.0)),         # isolated
]

oBnd = oDesign.GetModule("BoundarySetup")
for name, fid in ports:
    oBnd.AssignWavePort(
        ["NAME:" + name, "Faces:=", [fid], "NumModes:=", 1,
         "UseLineModeAlignment:=", False, "DoDeembed:=", False,
         "RenormalizeAllTerminals:=", True,
         ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", False, "CharImp:=", "Zpi"]],
         "ShowReporterFilter:=", False, "ReporterFilter:=", [True], "UseAnalyticAlignment:=", False])

# ---------------- solution setup (kept light for the Student licence) -------
oAna = oDesign.GetModule("AnalysisSetup")
oAna.InsertSetup("HfssDriven",
                 ["NAME:Setup1", "SolveType:=", "Single", "Frequency:=", "40GHz",
                  "MaxDeltaS:=", 0.02, "MaximumPasses:=", 12, "MinimumPasses:=", 2,
                  "MinimumConvergedPasses:=", 1, "PercentRefinement:=", 30,
                  "IsEnabled:=", True, "BasisOrder:=", 1])
oAna.InsertFrequencySweep("Setup1",
                          ["NAME:Sweep", "IsEnabled:=", True, "RangeType:=", "LinearCount",
                           "RangeStart:=", "35GHz", "RangeEnd:=", "45GHz", "RangeCount:=", 101,
                           "Type:=", "Interpolating", "SaveFields:=", False, "SaveRadFields:=", False,
                           "InterpTolerance:=", 0.5, "InterpMaxSolns:=", 50, "InterpMinSolns:=", 0,
                           "InterpMinSubranges:=", 1, "ExtrapToDC:=", False])

# ---------------- S-parameter plot ----------------
try:
    oRep = oDesign.GetModule("ReportSetup")
    oRep.CreateReport("S Parameters", "Modal Solution Data", "Rectangular Plot", "Setup1 : Sweep",
                      ["Domain:=", "Sweep"], ["Freq:=", ["All"]],
                      ["X Component:=", "Freq",
                       "Y Component:=", ["dB(S(1,1))", "dB(S(2,1))", "dB(S(3,1))", "dB(S(4,1))"]])
except Exception:
    oDesktop.AddMessage("", "", 1, "S-parameter plot not created; add it via Results > Create Modal Solution Data Report.")

# ---------------- save ----------------
proj_dir = oDesktop.GetProjectDirectory().rstrip("/\\")
oProject.SaveAs(proj_dir + "/Moreno_WR22_40GHz.aedt", True)
oDesktop.AddMessage("", "", 0,
                    "Moreno coupler built: lambda_g = %.3f mm, Delta = %.3f mm, d = %.3f mm. Saved to %s"
                    % (lamg, Delta, d, proj_dir + "/Moreno_WR22_40GHz.aedt"))
