# Free full-wave simulation with openEMS

[openEMS](https://www.openems.de) is a free, open-source FDTD solver with a
Python interface. These scripts build the same coupler as the HFSS builders,
from `design/design.json`, and write a Touchstone file that
`post/analyze_touchstone.py` reads unchanged. No licence needed.

## Run

```bash
python openems/simulate_openems.py                    # port 1 excited, ~4 min on 2 cores
python openems/simulate_openems.py --excite all       # full 4x4 S-matrix (4 runs)
python openems/simulate_openems.py --fine 0.0125      # reference-quality slot mesh (~25 min)
python openems/simulate_openems.py --view             # write geometry, open in AppCSXCAD
python openems/sweep_slot_length.py --L 2.0 2.2 2.4 --target 25
```

Outputs go to `results/openems/` (or `--out`): `moreno_openems.s4p`,
`sparams.png`, `summary.json`, `coupler_metrics.csv` and `run_info.json`
(mesh, run time, power balance). Scratch field data goes to
`results/openems/sim/`, which git ignores.

With the default single excitation, the `.s4p` contains column 1 (port 1
excited) and, by reciprocity, row 1. That is everything the coupler metrics
need; the file header says which entries are simulated. Use `--excite all`
for a complete matrix.

## How the model works

* **Metal only.** All metal is drawn as boxes, the complement of the two
  guides and the four slot arms, on a vacuum background; top and bottom are
  closed by the PEC boundary condition. Air is deliberately *not* drawn:
  openEMS makes a cell edge a perfect conductor only when the highest-priority
  primitive at the edge's midpoint is a metal (`Operator::CalcPEC_Range`), so
  air boxes carved out of a metal block with a higher priority would own every
  edge on a metal surface and move each wall one cell into the metal. That
  error made the slots look bigger and the wall thinner by mesh-dependent
  amounts; the first version of this model had it (see the report, Sec. 6.4).
* PML on the four guide ends.
* TE10 wave ports 1.5 λg from the common square. Guide 2's broad wall runs
  along z, its second transverse axis in openEMS, so its ports use `TE01` with
  a and b exchanged. S-parameters are normalised to each port's TE10 wave
  impedance; reference planes are moved to the edges of the square.
* Graded mesh: λ/20 in the feeds and a separate fine cell h (`--fine`, default
  0.035 mm) across the slot arms and through the wall. Guide walls, wall faces
  and port planes lie on mesh lines; the slot edges follow the **1/3–2/3 rule**
  (a line h/3 inside the metal and 2h/3 inside the slot; `--staircase` puts
  lines on the edges instead, for comparison). Runs stop at −50 dB
  field-energy decay.
* `--dump-fields` also stores the E-field at f0 on the mid-plane of each guide
  (`fields_f0.npz`).

**h is the only numerical parameter that matters.** The convergence study
(two edge treatments that bracket the converged coupling) and the
mesh-corrected design chart are in the
[technical report](../docs/report/Moreno_Coupler_WR22_40GHz_Report.pdf);
`sweep_slot_length.py` applies the measured h = 0.035 mm → 0 shift
automatically.

## Install openEMS

### Linux, or Windows through WSL (Ubuntu 24.04) — tested

```bash
sudo apt install build-essential cmake git libhdf5-dev libvtk9-dev libcgal-dev \
    libtinyxml-dev libboost-thread-dev libboost-system-dev libboost-date-time-dev \
    libboost-serialization-dev libboost-chrono-dev libboost-program-options-dev cython3
P=$HOME/opt/openEMS
for r in fparser CSXCAD openEMS; do
  git clone --depth 1 https://github.com/thliebig/$r.git
  (cd $r && mkdir -p build && cd build &&
   cmake .. -DCMAKE_INSTALL_PREFIX=$P -DFPARSER_ROOT_DIR=$P -DCSXCAD_ROOT_DIR=$P -DCMAKE_BUILD_TYPE=Release &&
   make -j"$(nproc)" && make install)
done
(cd CSXCAD/python  && CSXCAD_INSTALL_PATH=$P pip install .)
(cd openEMS/python && OPENEMS_INSTALL_PATH=$P CSXCAD_INSTALL_PATH=$P pip install .)
pip install --upgrade --ignore-installed h5py   # the h5py build must match your numpy
export LD_LIBRARY_PATH=$P/lib:$LD_LIBRARY_PATH   # add to ~/.bashrc
```

Check: `python -c "from openEMS import openEMS; print('ok')"`.

### Native Windows

Download the Windows build from the
[openEMS releases page](https://github.com/thliebig/openEMS-Project/releases),
unzip it (for example to `C:\opt\openEMS`), and install its Python modules
following the openEMS documentation for your release. Then

```bat
set OPENEMS_INSTALL_PATH=C:\opt\openEMS
python openems\simulate_openems.py
```

The script adds that folder to the DLL search path. The Windows build also
contains `AppCSXCAD.exe` for `--view`.

## Reference results

`results/openems/reference/` (finest mesh), `full_matrix/`, `fields/` and
`studies/` (every run of the verification campaign), summarised in the
technical report. Values quoted in earlier versions of this project
(24.3 dB coupling, 38.6 dB directivity, "L ≈ 2.23 mm for 20 dB") came from the
superseded air-priority model and are withdrawn.
