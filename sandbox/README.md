# BeamBeam Sandbox

This sandbox contains OPALX BeamBeam inputs, diagnostics, visualization tools,
and the manufactured witness-field note.  The current layout keeps generated
data separate from source scripts and notes:

- `python/`: reusable Python analysis and visualization tools.
- `note/`: LaTeX note and generated table inputs.
- `data/`: OPALX-generated diagnostics, figures, movies, and temporary outputs.
- `track-e-p/`: gamma-gamma pair tracking inputs, FROMFILE distributions, and
  generated OPALX element/particle visualization scripts.
- `oldpy/`: older one-off helper scripts kept for reference.

## Python Environment

Assume Python 3.11.  The sandbox environment is the repository-local
`.venv-h6` venv:

```bash
cd /path/to/opalx-beambeam
python3.11 -m venv --system-site-packages .venv-h6
source .venv-h6/bin/activate
python -m pip install -r sandbox/requirements-h6.txt
```

The scripts expect NumPy, pandas, matplotlib, Pillow, h5py, SciPy, PyVista/VTK,
and imageio.  Tkinter is needed for the native GUI but is not installed by pip.
For a MacPorts Python 3.11 environment, install Tk support with:

```bash
sudo port install \
  python311 \
  py311-tkinter
```

Current `.venv-h6` snapshot on this checkout:

```text
Python executable: .venv-h6/bin/python
Python version:    3.11.15
Tk/Tcl:            8.6 / 8.6

numpy              2.4.3
pandas             3.0.2
matplotlib         3.10.8
Pillow             9.5.0
h5py               3.16.0
scipy              1.17.1
pyvista            0.48.2
vtk                9.6.1
imageio            2.37.3
imageio-ffmpeg     0.6.0
```

Installed distributions reported by `python -m pip list --format=freeze`:

```text
attrs==26.1.0
beautifulsoup4==4.14.3
Brotli==1.2.0
certifi==2026.4.22
charset-normalizer==3.4.7
contourpy==1.3.3
cycler==0.12.1
cyclopts==4.11.2
docstring_parser==0.18.0
docutils==0.22.4
fonttools==4.62.1
h5py==3.16.0
idna==3.14
ImageIO==2.37.3
imageio-ffmpeg==0.6.0
kiwisolver==1.5.0
lxml==6.0.2
markdown-it-py==4.2.0
matplotlib==3.10.8
mdurl==0.1.2
numpy==2.4.3
oldest-supported-numpy==0.1
olefile==0.47
packaging==26.0
pandas==3.0.2
Pillow==9.5.0
pip==24.0
platformdirs==4.9.6
pooch==1.9.0
pycairo==1.29.0
Pygments==2.20.0
pyparsing==3.3.2
python-dateutil==2.9.0.post0
pytz==2026.1.post1
pyvista==0.48.2
requests==2.33.1
rich==15.0.0
rich-rst==1.3.2
scipy==1.17.1
scooby==0.11.2
setuptools==79.0.1
six==1.17.0
soupsieve==2.8.3
typing_extensions==4.15.0
unicodedata2==17.0.1
urllib3==2.7.0
vtk==9.6.1
zopfli==0.4.1
```

Quick check:

```bash
export MPLCONFIGDIR="${TMPDIR:-/tmp}/opalx-mpl"
python - <<'PY'
import tkinter as tk
import numpy
import pandas
import matplotlib
import PIL
import h5py
import scipy
import pyvista
import vtk
import imageio
import imageio_ffmpeg

print("tk", tk.TkVersion, tk.TclVersion)
print("numpy", numpy.__version__)
print("pandas", pandas.__version__)
print("matplotlib", matplotlib.__version__)
print("PIL", PIL.__version__)
print("h5py", h5py.__version__)
print("scipy", scipy.__version__)
print("pyvista", pyvista.__version__)
print("vtk", vtk.vtkVersion.GetVTKVersion())
print("imageio", imageio.__version__)
print("imageio-ffmpeg", imageio_ffmpeg.__version__)
PY
```

## Main Tools

- `python/beambeam_analysis.py`: combined GUI and CLI front-end for collision
  window, manufactured solution, and OPALX comparison views.
- `python/checkCollWin.py`: collision-window scalar dump viewer; writes PNG or
  GIF from OPALX `RHO_scalar-collwin_vis` dumps.
- `python/read_beambeam_h5.py`: BeamBeam diagnostics HDF5 overview and simple
  line-density/gallery plots.
- `python/beam-beam-manufactured-solution.py`: manufactured BeamBeam field
  comparison helper.
- `python/boosted_gaussian_witness.py`: Lorentz-transformed boosted-Gaussian
  manufactured witness calculation used by the LaTeX note.
- `note/boosted_gaussian_witness.tex`: current physics note.

## Combined GUI

```bash
cd /path/to/opalx-beambeam
source .venv-h6/bin/activate
python sandbox/python/beambeam_analysis.py gui
```

If Tk is unavailable, the tool falls back to a local browser UI.  For the native
Tk GUI, run from a normal terminal with the `.venv-h6` environment active.

The GUI stores its last inputs in:

```text
sandbox/.beambeam_analysis_state.json
```

## Collision-Window Movies

`checkCollWin.py` renders OPALX collision-window scalar dumps.  Give a seed dump
and a step range; a multi-step range writes a GIF.

For the current `BeamBeam-static-1V.in` example, OPALX normally writes only the
two active-window dumps at steps 5 and 6.  If the input is run from the repository
root as `sandbox/BeamBeam-static-1V.in`, those files are under `data/sandbox/`:

```bash
cd /path/to/opalx-beambeam
source .venv-h6/bin/activate

python sandbox/python/checkCollWin.py \
  data/sandbox/BeamBeam-static-1V-RHO_scalar-collwin_vis-000005.dat \
  --start 5 --end 6 \
  --physical \
  --output data/sandbox/BeamBeam-static-1V-collwin.gif \
  --input sandbox/BeamBeam-static-1V.in \
  --duration 120 \
  --save-frames
```

Useful options:

- `--physical`: plot physical coordinates in mm.  This is the preferred view
  for checking the BeamBeam element and collision-window positions.
- `--grid`: plot raw grid indices `i,j,k`; useful for debugging the mesh dump
  itself, but less useful for checking the longitudinal window location.
- `--no-geometry`: hide BeamBeam/window overlays.
- `--save-frames`: also write one PNG per movie frame.
- `--input path/to/input.in`: explicitly provide the OPALX input file used to
  infer BeamBeam geometry.

The collision-window figure shows:

- Three central slices of the dumped charge density: `x-y`, `x-s`, and `y-s`.
- Green dashed lines and the gold marker: the interaction point.
- Orange dashed rectangle: the full BeamBeam element span.
- Crimson shaded band and solid horizontal lines: the active collision-window
  interval.  The lower and upper lines are labelled `window start` and
  `window end`.  If the window is the full BeamBeam element, the label states
  `collision window = BeamBeam`.
- Cyan `b1` marker: the simulated bunch center.  Magenta `b2` marker: the
  mirrored bunch used for the symmetric collision-window visualization.

The title records the OPALX step, time, dump state, active-window flag, bunch
center, IP position, BeamBeam `s` range, collision-window `s` range, mesh
spacing, and integrated dumped charge.  If the dump range contains before/after
snapshots with the same OPALX `global_step`, both frames are kept in the GIF so
the transition into the collision window remains visible.

The combined front-end can run the same path:

```bash
python sandbox/python/beambeam_analysis.py collwin \
  data/sandbox/BeamBeam-static-1V-RHO_scalar-collwin_vis-000005.dat \
  --start 5 --end 6 \
  --physical \
  --output data/sandbox/BeamBeam-static-1V-collwin.gif \
  --input sandbox/BeamBeam-static-1V.in
```

If OPALX is run from inside `sandbox/` instead, replace `data/sandbox/` in the
commands above with `sandbox/data/` or with the directory where the dumps were
written.

## ASCII Diagnostics

The current BeamBeam diagnostics workflow is based on ASCII scalar/vector dumps,
not on a required HDF5 diagnostics file.  Use
`beam-beam-manufactured-solution.py` with one matched dump triplet from the same
OPALX step:

- `RHO_scalar-beambeam_rho_pre-*.dat`: binned source charge density.
- `PHI_scalar-beambeam_phi-*.dat`: solved electrostatic potential.
- `EF_vector-beambeam_e-*.dat`: solved electric field.

Example:

```bash
cd /path/to/opalx-beambeam
source .venv-h6/bin/activate

python sandbox/python/beam-beam-manufactured-solution.py \
  --compare-rho-dump data/sandbox/BeamBeam-2-RHO_scalar-beambeam_rho_pre-000005.dat \
  --compare-phi-dump data/sandbox/BeamBeam-2-PHI_scalar-beambeam_phi-000005.dat \
  --compare-e-dump data/sandbox/BeamBeam-2-EF_vector-beambeam_e-000005.dat \
  --output data/sandbox/BeamBeam-2-ascii-diagnostics.png
```

The script prints scalar error metrics and writes the requested output image plus
five derived line-profile figures:

- `BeamBeam-2-ascii-diagnostics.png`: central `x-z` slice comparison for
  `rho` and `phi`.  Each row shows analytic, OPALX, and OPALX-minus-analytic
  panels with max, L2, and relative L2 difference annotations.
- `BeamBeam-2-ascii-diagnostics-rho-z-axis.png`: `rho(z)` on the central beam
  axis and the corresponding difference curve.
- `BeamBeam-2-ascii-diagnostics-phi-z-axis.png`: `phi(z)` on the central beam
  axis and the corresponding difference curve.
- `BeamBeam-2-ascii-diagnostics-ez-z-axis.png`: `E_z(z)` on the central beam
  axis and the corresponding difference curve.  This is the main longitudinal
  field diagnostic near the interaction point.
- `BeamBeam-2-ascii-diagnostics-ex-x-axis.png`: `E_x(x)` through the grid plane
  nearest the interaction point and the corresponding difference curve.
- `BeamBeam-2-ascii-diagnostics-ey-y-axis.png`: `E_y(y)` through the grid plane
  nearest the interaction point and the corresponding difference curve.

For active BeamBeam-window snapshots, the script reconstructs the mirrored
Gaussian source from the dump metadata and uses the local interaction-point
position as the symmetry plane.

## Optional Diagnostics HDF5

Recent simple BeamBeam examples may not write a `*-beambeam_diagnostics.h5`
file.  The `read_beambeam_h5.py` helper is still useful for older runs or inputs
where that diagnostics file is explicitly produced, but it is not required for
the collision-window dump workflow above.

When such an HDF5 file exists, use:

```bash
python sandbox/python/read_beambeam_h5.py \
  --overview --table \
  path/to/example-beambeam_diagnostics.h5
```

Line-density GIF:

```bash
python sandbox/python/read_beambeam_h5.py \
  --line-density-z \
  --start 6 --end 57 \
  --bins 64 \
  --gif path/to/example-line-density-z.gif \
  path/to/example-beambeam_diagnostics.h5
```

Static projected gallery:

```bash
python sandbox/python/read_beambeam_h5.py \
  --gallery-xz \
  --start 6 --end 20 \
  --gallery-cols 5 \
  path/to/example-beambeam_diagnostics.h5
```

## Manufactured Witness Note

The LaTeX note is kept in this repository under `sandbox/note/` and is also
available on Overleaf:

```text
https://www.overleaf.com/project/6a00de23024e8d1800affb90
```

Regenerate the tables and figures used by the LaTeX note:

```bash
cd /path/to/opalx-beambeam
source .venv-h6/bin/activate

python sandbox/python/boosted_gaussian_witness.py

python sandbox/python/boosted_gaussian_witness.py \
  --pair-demo \
  --output-prefix sandbox/data/boosted_gaussian_witness_physical
```

Build the note:

```bash
cd sandbox/note
latexmk -pdf boosted_gaussian_witness.tex
```

The script writes generated LaTeX tables into `sandbox/note/` and pair-demo CSV
and figures into `sandbox/data/`.

## Gamma-Gamma Pair Tracking

The current gamma-gamma pair sandbox lives in `track-e-p/`.  It contains the
OPALX inputs, low-energy `e-/e+` FROMFILE distributions, and generated
visualization output.

Main files:

- `track-e-p/gamma_gamma_pairs.in`: baseline run without witness space-charge
  kicks.
- `track-e-p/gamma_gamma_pairs-2.in`: comparison run with witness-container
  space-charge sampling enabled.
- `track-e-p/fort98.txt`: original pair list from the gamma-gamma process.
- `track-e-p/gamma_gamma_electrons.fromfile`: electron FROMFILE distribution.
- `track-e-p/gamma_gamma_positrons.fromfile`: positron FROMFILE distribution.

Regenerate the FROMFILE distributions:

```bash
cd sandbox/track-e-p
source ../../.venv-h6/bin/activate

python convert_fort98_to_fromfile.py fort98.txt
```

Run OPALX from `track-e-p/` so generated output lands under `track-e-p/data/`:

```bash
cd sandbox/track-e-p

/path/to/opalx gamma_gamma_pairs.in
```

For MPI:

```bash
cd sandbox/track-e-p

mpirun -np 2 /path/to/opalx gamma_gamma_pairs.in
```

OPALX writes generated element visualization scripts next to the output, for
example:

```text
track-e-p/data/gamma_gamma_pairs_ElementPositions.py
track-e-p/data/gamma_gamma_pairs-2_ElementPositions.py
```

These generated scripts can visualize the lattice, BeamBeam aperture, H5
particles, and write particle movies.  Do not edit them for persistent changes:
update `src/Structure/MeshGenerator.cpp`, rebuild OPALX, and rerun the input.

Interactive particle view with a step slider:

```bash
cd sandbox/track-e-p
source ../../.venv-h6/bin/activate

python data/gamma_gamma_pairs-2_ElementPositions.py \
  --show \
  --part-vis \
  --part-slider \
  --part-max-primary 20000 \
  --part-max-pairs 10000
```

Static latest-step view:

```bash
python data/gamma_gamma_pairs-2_ElementPositions.py \
  --show \
  --part-vis \
  --part-step latest
```

Specific global step:

```bash
python data/gamma_gamma_pairs-2_ElementPositions.py \
  --show \
  --part-vis \
  --part-step 100
```

Write a particle GIF:

```bash
python data/gamma_gamma_pairs-2_ElementPositions.py \
  --part-movie pairs.gif \
  --movie-fps 12 \
  --movie-step-stride 2 \
  --part-max-primary 20000 \
  --part-max-pairs 10000
```

Notes:

- `--show` is required for interactive/static visualization.  Without it, the
  generated script prints help and exits.
- `--part-slider` is interactive only and is disabled while writing a movie.
- `--part-step latest` is the default for static particle views.
- `--movie-fps` sets playback rate.
- `--movie-step-stride` keeps every Nth stored H5 step to reduce movie size.
- `--part-debug` prints H5 particle-step discovery details.
- GIF output uses PyVista plus imageio; MP4-style filenames use PyVista's movie
  writer and may also need `imageio-ffmpeg`.

Momentum diagnostics:

```bash
cd sandbox/track-e-p
source ../../.venv-h6/bin/activate

python analyze_pair_momenta.py --plot
```

By default the momentum script reads `gamma_gamma_electrons.fromfile` and
`gamma_gamma_positrons.fromfile`, then writes figures and CSV summaries into
`track-e-p/data/`.
