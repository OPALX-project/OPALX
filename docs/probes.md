# Cyclotron and FFA probes

`PROBE` is a passive finite detector segment for cyclotron-sector and native `RING` particle tracking. It records crossings in both directions without changing fields, particle/reference states, accepted timesteps or stopping conditions. Each declared station is global and is included once; magnetic ring symmetry does not replicate detectors.

## Geometry and units

`XSTART`/`XEND` specify laboratory X, and the legacy names `YSTART`/`YEND` specify laboratory Z. The horizontal plane is X-Z; laboratory Y is vertical. Endpoints and `STEP` are in millimetres. A segment has finite horizontal extent and no vertical clipping. `WIDTH` is retained but unused. Detector coordinates do not inherit magnetic element placement.

```opal
P1: PROBE, XSTART=-4420, XEND=-4420,
           YSTART=3420, YEND=5420, OUTFN="PROBE1";
// Include P1 once in the RING together with its magnetic elements.
```

This example is a radial station for a ring centred at laboratory (X,Z)=(-4.42,0) m. Translated or rotated lattices require explicitly transformed detector endpoints.

## Native output

The standard loss writer produces `OUTFN.loss` with `ASCIIDUMP=TRUE`, otherwise `OUTFN.h5`. Output is in the fixed OPALX laboratory frame, with positions in metres, momenta in p/(mc), and time in seconds. The ten ASCII columns are:

```text
x y z px py pz id turn bunchNumber time
```

HDF5 stores the corresponding named datasets in `Step#` groups. Groups contain flushed batches, not individual turns. ASCII uses the existing loss writer's six significant digits. Particle IDs remain stable through sorting, compaction and MPI migration; file row order is not a particle identity or chronological ordering guarantee. Join records using station, particle ID, bunch and time.

`turn` is the zero-based shared reference-return count at the start of the accepted step, independently maintained for each particle container. It is not a per-particle crossing number. A fixed launch plane normal to the initial reference momentum defines the reference return. A step crossing that plane retains its step-start label. Other trackers may use different turn-counter definitions. A coordinate transform does not reconcile them; compare station, particle ID, time and chronological passage while retaining the native reference-return labels.

Cyclotron output additionally contains lab-origin radial `.hist` counts and unweighted `.peaks` turn-mean radii. FFA produces phase-space records only: the laboratory origin need not coincide with its ring centre. `STEP` remains finite and positive, default 1 mm; FFA does not use it to control tracking or sample frequency. No-hit FFA stations produce no phase-space file, matching the standard loss writer. Use fresh output basenames; existing probe products are rejected to prevent mixing runs.

## Crossing and preparation semantics

Only accepted tracking endpoints are observed, after physical coordinate restoration and loss processing. For signed plane distances d0 and d1, the bracketed crossing fraction is alpha=-d0/(d1-d0). Position, momentum and time use the same linear interpolation fraction, followed by finite-segment clipping. Initial points on a plane are suppressed, and directed hysteresis prevents duplicate endpoint/jitter hits.

For smooth, resolved transverse motion, linear interpolation has second-order event-location error in the step duration. Hidden intra-step recrossings, tangencies and RF discontinuities are not reconstructed. Old OPAL used extrapolation from the current state, so exact file equality is not expected; compare matched phase space and timestep convergence. Reducing diagnostic sampling frequency can miss hits and is not a performance substitute.

COF accepts and skips probes. Closed-orbit trials, finite-difference rays and OrbitThreader preparation do not allocate probe recorders or write crossings. Spectral-tune trial tracking with PROBE, checkpoint restart and tracking continuation remain unsupported.

Histories are registered particle attributes containing fixed laboratory endpoints. They migrate with stable IDs, including when a space-charge solver temporarily changes coordinates or ownership. They must not be transformed into a solver's temporary coordinate frame. Recording and output remain independent of the selected supported transport/space-charge backend.

Pending hits are flushed collectively every 100 global tracking steps and at normal completion; implemented RF failure paths also flush accepted hits. General out-of-field-map shutdown and diagnostic preservation remain tracked in issue #672.

## Coordinate-frame conversion

Keep native laboratory output canonical. The optional [probe-frame adapter](../tools/probe_frame_adapter.md) converts records into any explicitly supplied Cartesian frame, independently of the lattice or example. Given a target origin o in laboratory metres and an orthonormal matrix B whose rows are target axes expressed in laboratory coordinates, it applies `r_target = B (r_lab - o)` and `p_target = B p_lab`. It supports arbitrary translations, rotations and axis conventions, preserving particle IDs, bunch numbers, reference-turn labels and times. This analysis utility does not change magnetic fields, lattice placement or detector geometry, and is not a native unit-test dependency.

Validation uses self-contained native C++ crossing/migration checks and separate old-OPAL integration examples with matched complete six-dimensional initial states. Probe-enabled/disabled runs should have identical physical output. Performance optimization remains deferred.
