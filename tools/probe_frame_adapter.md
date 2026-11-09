# Probe coordinate-frame adapter

`probe_frame_adapter.py` expresses one station's OPALX laboratory PROBE records
in an explicitly defined Cartesian frame. It works with cyclotron and FFA output
and has no lattice, ring-radius or example-specific assumptions. Native output
remains canonical. This is a standalone analysis utility, separate from native
C++ unit tests. Use Python with NumPy, pandas and h5py.

## Frame definition

Supply the target origin **o** in laboratory metres and a matrix **B** whose rows
are the target frame's unit axes expressed in laboratory coordinates:

```text
r_target = B (r_lab - o)
p_target = B p_lab
```

The default origin is (0,0,0) and the default basis is the identity. Arbitrary
translations and orthonormal bases are supported, including axis permutations,
sign changes and either handedness. Both position and momentum are polar vectors;
only position is translated. The basis must be finite and orthonormal within an
absolute 1e-12 input-validation tolerance (zero relative tolerance), allowing
rounded decimal axes while rejecting scaling or shear. The utility does not
transform fields, lattice placements or probe detector geometry.

For example, convert a frame whose origin is laboratory (-4.42,0,0) m and whose
X,Y,Z axes are laboratory X,Z,Y:

```sh
~/.venv-h6/bin/python tools/probe_frame_adapter.py \
    native/PROBE1.h5 analysis/PROBE1.loss \
    --origin-m -4.42 0 0 --basis 1 0 0 0 0 1 0 1 0
```

This gives (x+4.42,z,y) and (px,pz,py). The nine basis entries are supplied
row-major. A 90-degree rotated frame about laboratory Y instead uses
`--basis 0 0 -1 0 1 0 1 0 0`; the origin is independent of orientation.
No origin or orientation is inferred from particle trajectories.

## Records and precision

Native ASCII `.loss` and HDF5 `.h5` input are supported. Output uses the standard
ten ASCII columns `x,y,z,px,py,pz,id,turn,bunchNumber,time`. Positions remain in
metres, momenta in p/(mc), and time in seconds. HDF5 charge and mass datasets are
outside this ten-column format and are not exported.

Particle IDs, bunch numbers, reference-turn labels and times are preserved;
there is no relabelling option. A coordinate transform cannot reconcile different
tracking-counter definitions. Match independent runs by station, particle ID,
time and chronological passage, checking their turn definitions separately.
IDs need not be contiguous. Rows are sorted stably by `bunchNumber,id,time`,
without assuming one crossing per turn. HDF5 `Step#` groups are flush batches,
not turns; all batches are concatenated.

Input values must be finite and identities exact nonnegative signed-64-bit
integers. Native ASCII normally has six significant digits; conversion cannot
recover precision already lost. Output uses 17 significant digits, while HDF5
input retains its source precision subject to floating-point transformation.

The output directory must exist. Files are created exclusively to protect
existing data. Run once per station. Missing files are errors, including no-hit
runs for which the native writer created no file. Header-only ASCII input
produces header-only output; no hit records are invented.

## Python API

- `read_probe(path)` returns a pandas DataFrame.
- `transform_frame(frame, origin_m=(0,0,0), basis=None)` returns a transformed copy.
- `write_probe_ascii(frame, path)` creates a standard ten-column ASCII file.

The DataFrame uses the ten columns listed above. To invert a transform, use
basis `B.T` and origin `-B @ o` expressed in the transformed frame.
