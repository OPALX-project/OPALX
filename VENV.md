# Python Environment

This repository uses a local Tk-capable Python virtual environment for the BeamBeam analysis tools.

## Prerequisites

Install the required MacPorts packages:

```bash
sudo port install \
  python311 \
  py311-tkinter \
  py311-numpy \
  py311-matplotlib \
  py311-Pillow \
  py311-h5py
```

## Create The Virtual Environment

From the repository root:

```bash
cd /Users/adelmann/git/opalx
/opt/local/bin/python3.11 -m venv --system-site-packages .venv-h6
```

## Activate

```bash
cd /Users/adelmann/git/opalx
source .venv-h6/bin/activate
export MPLCONFIGDIR=/tmp/mpl
```

## Verify

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

## Main Commands

Manufactured solution script:

```bash
python sandbox/beam-beam-manufactured-solution.py --help
```

Combined GUI:

```bash
python sandbox/beambeam_analysis.py gui
```

HDF diagnostics reader:

```bash
python sandbox/read_beambeam_h5.py --help
```
