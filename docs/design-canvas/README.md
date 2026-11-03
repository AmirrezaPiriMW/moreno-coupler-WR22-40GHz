# Design canvas source

Source of the interactive design sheet made alongside this project (a Claude
Design canvas):

| File | Content |
|---|---|
| `Main.dc.html` | Design calculator: waveguide presets (WR-137 … WR-22), inputs f0, L, w, t, Δ; slot layout to scale; computed geometry; the same design checks as `moreno/design.py`; closed-form band trend |
| `Assembly.dc.html` | Port map, top view (WR-22, 40 GHz) |
| `canvas.json` | Canvas layout |

The files need the canvas runtime (`support.js`) to render, so they do not
open as plain web pages; they are kept so the calculator logic can be read and
reused. All of it is in the `<script type="text/x-dc">` block of
`Main.dc.html`, in plain JavaScript that mirrors `moreno/design.py` (same
formulas, same checks, same thresholds). `moreno/design.py` remains the
reference implementation.
