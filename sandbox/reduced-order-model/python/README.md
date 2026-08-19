# Python reference model

`rigid_two_gaussian_fields.py` evaluates the lab-frame electric and magnetic
fields of the two rigid counter-propagating primary Gaussian beams. It produces
snapshots when their centroids are 3, 2, 1, and 0 longitudinal rms sizes apart.

Run from the repository root:

```bash
~/.venv-h6/bin/python \
  sandbox/reduced-order-model/python/rigid_two_gaussian_fields.py
```

The script writes one CSV field table per separation, a numerical summary, a
JSON record of assumptions and derived parameters, and a combined x-z field
figure under `sandbox/reduced-order-model/outputs/fields/`. Every CSV contains
all three components of E and B in SI units; the plot shows the nonzero
components Ex, Ez, and By in the symmetry plane y=0.

Run the focused checks with:

```bash
~/.venv-h6/bin/python -m unittest \
  sandbox/reduced-order-model/python/test_rigid_two_gaussian_fields.py
```

The historical manufactured-solution material remains in
`sandbox/analytic-model/`; it is intentionally not moved while the reduced model
interface and comparison quantities are being fixed.
