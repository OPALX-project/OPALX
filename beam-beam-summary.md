BeamBeam Summary

Physics Model

- A single physical bunch (`b1`) is tracked through the lattice.
- Inside the BeamBeam window, the opposing bunch (`b2`) is approximated by a mirrored copy of `b1`.
- The mirror is taken about the interaction point:
  - `z_mirror = 2 z_IP - z_actual`
- Outside BeamBeam tracking there is no mirrored bunch contribution.
- The BeamBeam window is the BeamBeam element itself:
  - `beginS = element begin`
  - `endS = element end`
  - `interactionPointS = 0.5 * (beginS + endS)`
- The intended runtime transition is:
  - before entering the BeamBeam element: normal single-bunch tracking
  - inside the BeamBeam element: real bunch plus mirrored copy for the field solve
  - after the first particle exits the BeamBeam element: back to normal single-bunch tracking
- The field-mesh model is:
  - outside BeamBeam: bunch-following in `z`
  - inside BeamBeam: fixed in the lab / BeamBeam frame over the full BeamBeam element
  - `x/y` remain bunch-following

Implementation

- BeamBeam uses the element stack:
  - `BeamBeam`
  - `BeamBeamRep`
  - `OpalBeamBeam`
- Shared BeamBeam type definitions live in [BeamBeamDefinitions.h](/Users/adelmann/git/opalx/src/AbsBeamline/BeamBeamDefinitions.h):
  - `BEAMBEAM::Config`
  - `BEAMBEAM::ActualGeometry`
  - `BEAMBEAM::WindowState`
  - `BEAMBEAM::Runtime`
  - `BEAMBEAM::Diagnostics`
- The runtime state machine is:
  - `Inactive`
  - `Active`
  - `Completed`
- `ParallelTracker` owns:
  - BeamBeam entry / exit detection
  - actual lab-frame BeamBeam geometry
  - the high-level BeamBeam state machine
- `PartBunch` owns:
  - the fixed BeamBeam mesh handoff
  - mirrored deposition
  - restoration of the normal co-moving field domain

Runtime Logic

- Entry detection occurs before the self-field solve of the first active BeamBeam step.
- BeamBeam activation now starts when the full bunch is inside the BeamBeam element.
- While `Active`:
  - the longitudinal mesh is fixed to the BeamBeam element in the lab frame
  - the bunch moves through that fixed mesh
  - the particle container is not periodically re-wrapped each active step
- Mirrored deposition is implemented at scatter time:
  - zero `rho`
  - scatter the real bunch
  - if `COPY=TRUE`, mirror particle `z` about the BeamBeam IP and scatter again
  - then solve on that resulting density
- Exit occurs when the first particle exits the BeamBeam element.
- After exit:
  - copying stops
  - the field-domain state from before BeamBeam is restored
  - tracking returns to the ordinary co-moving mesh

Diagnostics

- ASCII diagnostics:
  - `data/collision_ascii_frames.txt`
  - `BeamBeamWindowAnimation`
- HDF5 diagnostics:
  - [H5BeamBeamDiagnosticsWriter.h](/Users/adelmann/git/opalx/src/Structure/H5BeamBeamDiagnosticsWriter.h)
  - [H5BeamBeamDiagnosticsWriter.cpp](/Users/adelmann/git/opalx/src/Structure/H5BeamBeamDiagnosticsWriter.cpp)
- Current user-facing HDF states are:
  - `normal tracking`
  - `beambeam tracking`
- The HDF path writes one file per run and currently stores:
  - `rho`
  - `phi`
  - `Ex`, `Ey`, `Ez`
  - mesh origin / spacing
  - BeamBeam metadata
  - particle mean position and charge metadata

Python Tools

- [checkCollWin.py](/Users/adelmann/git/opalx/sandbox/checkCollWin.py)
  - collision-window / BeamBeam scalar-dump viewer
- [read_beambeam_h5.py](/Users/adelmann/git/opalx/sandbox/read_beambeam_h5.py)
  - HDF5 overview, line density, and `x,z` gallery tool
- [beam-beam-manufactured-solution.py](/Users/adelmann/git/opalx/sandbox/beam-beam-manufactured-solution.py)
  - manufactured-solution generator and OPALX comparison tool
- [beambeam_analysis.py](/Users/adelmann/git/opalx/sandbox/beambeam_analysis.py)
  - combined front-end with Tk GUI / browser fallback

Manufactured Solution

- The current manufactured solution is an exact symmetric 3D isotropic Gaussian pair.
- It compares OPALX against analytic:
  - `rho`
  - `phi`
- Supported manufactured-compare views include:
  - `rho+phi`
  - `rho`
  - `phi`
  - `rho-z-axis`
  - `phi-z-axis`
- Movie mode uses fixed axes across the selected range.

Current Design Choices

- The BeamBeam window is the BeamBeam element.
- The IP remains the midpoint of that element.
- The raw HDF schema still stores internal snapshot kinds, but user-facing tools map them to:
  - `normal tracking`
  - `beambeam tracking`
- `read_beambeam_h5.py` remains part of the workflow and is not replaced by the combined GUI.
