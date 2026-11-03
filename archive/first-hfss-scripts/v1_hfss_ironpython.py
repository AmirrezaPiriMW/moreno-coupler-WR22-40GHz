# ----------------------------------------------------------------------------
# Moreno cross-guide coupler, WR-22, f0 = 40 GHz
# Builds a parametric HFSS (Driven Modal) model and saves it as an .aedt project.
#
# Run in AEDT:  Tools > Run Script...  (IronPython or CPython, AEDT 2019 R1+)
#
# Geometry (air volume only; every face on the background is PEC by default):
#   Guide 1 (lower) runs along +z, x in [0,a], y in [0,b]
#   Shared wall      y in [b, b+t]
#   Guide 2 (upper) runs along x, z in [0,a], y in [b+t, 2b+t]
#   Two "+" crossed slots through the wall, centred at (d,d) and (a-d,a-d)
#   with d = (a - Delta)/2 and Delta = lambda_g/4 at f0.
#
# Ports: 1 input (z = -Lext), 2 through (z = a+Lext),
#        3 coupled (x = -Lext), 4 isolated (x = a+Lext)
# ----------------------------------------------------------------------------
import math

try:
    import ScriptEnv
    ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
except Exception:
    pass  # newer AEDT already provides oDesktop

# ---------------- design values (mm, GHz) ----------------
C0 = 299.792458          # mm*GHz
f0 = 40.0
a = 5.690                # WR-22 broad wall
b = 2.845                # WR-22 narrow wall
t = 0.25                 # shared wall thickness
Ls = 2.0                 # slot arm length
ws = 0.35                # slot arm width

fc = C0 / (2.0 * a)
lam0 = C0 / f0
lamg = lam0 / math.sqrt(1.0 - (fc / f0) ** 2)
Delta = lamg / 4.0
d = (a - Delta) / 2.0
Lext = round(2.0 * lamg, 3)  # feed length beyond the shared square

# ---------------- project / design ----------------
oDesktop.RestoreWindow()
oProject = oDesktop.NewProject()
design_name = "Moreno_WR22_40GHz"
oProject.InsertDesign("HFSS", design_name, "DrivenModal", "")
oDesign = oProject.SetActiveDesign(design_name)
oEditor = oDesign.SetActiveEditor("3D Modeler")
oEditor.SetModelUnits(["NAME:Units Parameter", "Units:=", "mm", "Rescale:=", False])


def addvar(name, value):
    oDesign.ChangeProperty(
        ["NAME:AllTabs",
         ["NAME:LocalVariableTab",
          ["NAME:PropServers", "LocalVariables"],
          ["NAME:NewProps",
           ["NAME:" + name, "PropType:=", "VariableProp", "UserDef:=", True, "Value:=", value]]]])


addvar("a", "%.4fmm" % a)
addvar("b", "%.4fmm" % b)
addvar("t", "%.4fmm" % t)
addvar("Ls", "%.4fmm" % Ls)
addvar("ws", "%.4fmm" % ws)
addvar("Delta", "%.4fmm" % Delta)      # lambda_g/4 at 40 GHz
addvar("d", "(a-Delta)/2")
addvar("Lext", "%.3fmm" % Lext)


def box(name, pos, size):
    oEditor.CreateBox(
        ["NAME:BoxParameters",
         "XPosition:=", pos[0], "YPosition:=", pos[1], "ZPosition:=", pos[2],
         "XSize:=", size[0], "YSize:=", size[1], "ZSize:=", size[2]],
        ["NAME:Attributes", "Name:=", name, "Flags:=", "", "Color:=", "(132 170 210)",
         "Transparency:=", 0.7, "PartCoordinateSystem:=", "Global", "UDMId:=", "",
         "MaterialValue:=", "\"vacuum\"", "SolveInside:=", True])


# Guides
box("Air", ["0mm", "0mm", "-Lext"], ["a", "b", "a+2*Lext"])
box("G2", ["-Lext", "b+t", "0mm"], ["a+2*Lext", "b", "a"])

# Slot 1 at (d, d)
box("S1x", ["d-Ls/2", "b", "d-ws/2"], ["Ls", "t", "ws"])
box("S1z", ["d-ws/2", "b", "d-Ls/2"], ["ws", "t", "Ls"])
# Slot 2 at (a-d, a-d)
box("S2x", ["a-d-Ls/2", "b", "a-d-ws/2"], ["Ls", "t", "ws"])
box("S2z", ["a-d-ws/2", "b", "a-d-Ls/2"], ["ws", "t", "Ls"])

oEditor.Unite(["NAME:Selections", "Selections:=", "Air,G2,S1x,S1z,S2x,S2z"],
              ["NAME:UniteParameters", "KeepOriginals:=", False])

# ---------------- wave ports ----------------
def face_at(x, y, z):
    return oEditor.GetFaceByPosition(
        ["NAME:FaceParameters", "BodyName:=", "Air",
         "XPosition:=", "%.5fmm" % x, "YPosition:=", "%.5fmm" % y, "ZPosition:=", "%.5fmm" % z])


y2 = b + t + b / 2.0
port_faces = [
    ("1", face_at(a / 2.0, b / 2.0, -Lext)),          # input
    ("2", face_at(a / 2.0, b / 2.0, a + Lext)),       # through
    ("3", face_at(-Lext, y2, a / 2.0)),               # coupled
    ("4", face_at(a + Lext, y2, a / 2.0)),            # isolated
]

oBnd = oDesign.GetModule("BoundarySetup")
for name, fid in port_faces:
    oBnd.AssignWavePort(
        ["NAME:" + name, "Faces:=", [fid], "NumModes:=", 1,
         "RenormalizeAllTerminals:=", True, "UseLineModeAlignment:=", False,
         "DoDeembed:=", False,
         ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", False, "CharImp:=", "Zpi"]],
         "ShowReporterFilter:=", False, "ReporterFilter:=", [True], "UseAnalyticAlignment:=", False])

# ---------------- solution setup ----------------
oAna = oDesign.GetModule("AnalysisSetup")
oAna.InsertSetup("HfssDriven",
                 ["NAME:Setup1", "Frequency:=", "40GHz", "MaxDeltaS:=", 0.01,
                  "MaximumPasses:=", 20, "MinimumPasses:=", 3, "MinimumConvergedPasses:=", 2,
                  "PercentRefinement:=", 30, "IsEnabled:=", True, "BasisOrder:=", 1])
oAna.InsertFrequencySweep("Setup1",
                          ["NAME:Sweep", "IsEnabled:=", True, "RangeType:=", "LinearCount",
                           "RangeStart:=", "35GHz", "RangeEnd:=", "45GHz", "RangeCount:=", 201,
                           "Type:=", "Interpolating", "SaveFields:=", False, "SaveRadFields:=", False,
                           "InterpTolerance:=", 0.5, "InterpMaxSolns:=", 250, "InterpMinSolns:=", 0,
                           "InterpMinSubranges:=", 1, "ExtrapToDC:=", False])

# ---------------- S-parameter report ----------------
try:
    oRep = oDesign.GetModule("ReportSetup")
    oRep.CreateReport("S Parameters", "Modal Solution Data", "Rectangular Plot", "Setup1 : Sweep",
                      ["Domain:=", "Sweep"], ["Freq:=", ["All"]],
                      ["X Component:=", "Freq",
                       "Y Component:=", ["dB(S(1,1))", "dB(S(2,1))", "dB(S(3,1))", "dB(S(4,1))"]], [])
except Exception:
    pass  # report API differs slightly between releases; add the plot by hand if needed

# ---------------- save ----------------
path = oDesktop.GetProjectDirectory().rstrip("/\\") + "/Moreno_WR22_40GHz.aedt"
oProject.SaveAs(path, True)
oDesktop.AddMessage("", "", 0, "Moreno coupler built: d = %.3f mm, Delta = %.3f mm, lambda_g = %.3f mm. Saved to %s"
                    % (d, Delta, lamg, path))
