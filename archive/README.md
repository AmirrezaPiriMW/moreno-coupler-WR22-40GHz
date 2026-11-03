# Archive: superseded files

Kept for the history of the project only. **Do not use them**; use `hfss/`
instead. Both scripts were written before the review documented in the
technical report (Appendix A).

| File | What it was | Why it was replaced |
|---|---|---|
| `first-hfss-scripts/v1_hfss_ironpython.py` | First AEDT Run Script (IronPython era), 2 λg feeds, guides and slots united into one body | Inserts the design with the old solution-type name only, which AEDT 2025 R2 rejects ("unable to insert design"); adaptive ΔS = 0.01 cannot resolve a −30 dB coupling; no slot mesh seed |
| `first-hfss-scripts/v2_aedt2025r2.py` | Rewrite for AEDT 2025 R2 (CPython), 1 λg feeds | Fixed values instead of `design/design.json`; ΔS = 0.02 (too loose, see report Sec. 9); no slot mesh seed; superseded by `hfss/build_in_aedt.py` |

The first PyAEDT builder (not kept) left PyAEDT's default of renormalising
the wave ports to 50 Ω, which distorts every waveguide S-parameter; the current
`hfss/build_pyaedt.py` sets `renormalize=False`.
