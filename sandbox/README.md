# Sandbox Environment

This directory contains local BeamBeam analysis tools and test inputs.

## Recommended Environment

Use the Tk-capable MacPorts Python 3.11 environment:

```bash
cd /Users/adelmann/git/opalx-beambeam

/opt/local/bin/python3.11 -m venv --system-site-packages .venv-h6
source .venv-h6/bin/activate
```

The following MacPorts packages are expected to be installed:

```bash
sudo port install \
  python311 \
  py311-tkinter \
  py311-numpy \
  py311-matplotlib \
  py311-Pillow \
  py311-h5py
```

Quick verification:

```bash
python - <<'PY'
import tkinter as tk
import numpy
import matplotlib
import PIL
import h5py

print("tk", tk.TkVersion, tk.TclVersion)
print("numpy", numpy.__version__)
print("matplotlib", matplotlib.__version__)
print("PIL", PIL.__version__)
print("h5py", h5py.__version__)
PY
```

Expected:
- `tk 8.6 8.6`

## Launch

Combined GUI front-end:

```bash
cd /Users/adelmann/git/opalx-beambeam
source .venv-h6/bin/activate
python sandbox/beambeam_analysis.py gui
```

On systems where Tk is unavailable, or when the process is launched from a
sandbox that cannot create macOS GUI windows, the tool falls back to a local
browser UI. In this checkout `.venv-h6` is based on MacPorts Python 3.11 and
imports Tk 8.6 correctly; use a normal terminal session for the native Tk GUI.

## Notes

- The GUI stores the last inputs in:
  - [`.beambeam_analysis_state.json`](/Users/adelmann/git/opalx/sandbox/.beambeam_analysis_state.json)
- For movie views in the Tk GUI:
  - if no output path is given, the tool renders to a temporary file
  - the viewer window provides `Start` / `Stop` and a frame slider
- If `matplotlib` complains about its cache path, use:

```bash
export MPLCONFIGDIR=/tmp/mpl
```

## Main Tools

- [beambeam_analysis.py](/Users/adelmann/git/opalx/sandbox/beambeam_analysis.py)
  - combined front-end
- [checkCollWin.py](/Users/adelmann/git/opalx/sandbox/checkCollWin.py)
  - collision-window viewer
- [beam-beam-manufactured-solution.py](/Users/adelmann/git/opalx/sandbox/beam-beam-manufactured-solution.py)
  - manufactured-solution comparison
- [read_beambeam_h5.py](/Users/adelmann/git/opalx/sandbox/read_beambeam_h5.py)
  - HDF5 diagnostics reader

## Still Used: read_beambeam_h5.py

`read_beambeam_h5.py` is still part of the normal workflow. The combined GUI does not replace it for all HDF5 diagnostics tasks.

Typical usage:

```bash
cd /Users/adelmann/git/opalx-beambeam
source .venv-h6/bin/activate

python sandbox/read_beambeam_h5.py \
  --overview --table \
  data/sandbox/beambeam-1-beambeam_diagnostics.h5
```

```bash
python sandbox/read_beambeam_h5.py \
  --gallery-xz \
  --start 20 --end 35 \
  data/sandbox/beambeam-1-beambeam_diagnostics.h5
```

```bash
python sandbox/read_beambeam_h5.py \
  --line-density-z --bins 32 \
  --start 20 --end 35 \
  data/sandbox/beambeam-1-beambeam_diagnostics.h5
```
