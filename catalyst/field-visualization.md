# External Field Visualization — Design Notes

## Context

The `catalyst-viz-v2-jb` branch adds Catalyst in-situ visualization to OPALX.
Particle channels (`MUON`, `BEAM1`, ...) and a self-field mesh channel
(`bunch_mesh`) are already working. This document captures what we learned
about the OPALX field architecture on `origin/master` (PR #555,
`45e126f3a8e7`) and how to extend the Catalyst layer to visualize **external**
element fields (field maps and analytic models).

---

## 1. Field architecture after the PR #555 refactor

### 1.1 Self-fields (already exported)

Self-fields live on a PIC mesh owned by `CartesianPIC3DFieldStorage`
(in `src/SpaceCharge/CartesianPIC3D/`). The ownership chain is:

```
TrackRun
  ├── bunch_m (PartBunch — owns CartesianDomain: mesh + field layout)
  ├── spaceChargeSolver_m (SpaceChargeSolver)
  │     └── algorithm_m (CartesianPIC3DAlgorithm)
  │           └── fieldStorage_m (CartesianPIC3DFieldStorage)
  │                 ├── electricField()   → ippl::Field<Vector<T,Dim>, ...>
  │                 ├── chargeDensity()   → ippl::Field<T, ...>
  │                 ├── potential()       → ippl::Field<T, ...> (CG only)
  │                 ├── accumulatedElectricField()
  │                 ├── accumulatedMagneticField()
  │                 └── flippedZSlabField()
  └── itsTracker_m (ParallelTracker — borrows spaceChargeSolver_m)
```

Flow: particles → scatter → `rho` (mesh) → Poisson → `E`/`phi` (mesh) →
gather → `pc.E`, `pc.B`.

**Export status:** Works. `ParallelTracker.cpp` registers `bunch_mesh` with
`addMeshChannel<3>()` and adds `E_self`, `B_self`, `rho`, `phi` as arrays.
**Needs adaptation** after rebasing onto `origin/master` — `FieldContainer` is
deleted; access must go through `spaceChargeSolver_m` → algorithm → field
storage (see section 4).

### 1.2 External fields (not yet exported)

External fields are computed **at each particle's position** by the element's
field model. There is no intermediate PIC mesh for external fields.

**Two categories:**

#### a) Field-map elements (grid-based)

Elements that load a field map from disk at `goOnline()` time. The grid data
is stored in `Kokkos::DualView<double*>` members on the `Fieldmap` subclass.

| Fieldmap class          | MapType              | Grid dims | Time-varying?         |
|-------------------------|----------------------|-----------|-----------------------|
| `FM2DMagnetoStatic`     | T2DMagnetoStatic     | (r, z)    | No — static B         |
| `FM2DElectroStatic`     | T2DElectroStatic     | (r, z)    | No — static E         |
| `FM2DDynamic`           | T2DDynamic           | (r, z)    | RF scalar × cos(ωt+φ)|
| `Astra1DDynamic`        | TAstraDynamic        | 1D (z)    | RF scalar × cos(ωt+φ)|
| `Astra1DMagnetoStatic`  | TAstraMagnetoStatic  | 1D (z)    | No — static B         |
| `G4BL2DMagnetoStatic`   | TG4BL2DMagnetoStatic | (x, z)    | No — static B         |

**Key point:** The grid is always static. "Dynamic" fieldmaps have a
time-varying scalar amplitude (`electricScale = V_peak × cos(ωt + φ)`,
`magneticScale` similarly), but the grid data, coordinates, and interpolation
kernel are fixed for the entire run. The `applyField()` / `applyRFField()`
methods run a Kokkos kernel over all particles, doing bilinear interpolation
from the grid at `R(i)` and accumulating into `pc.E` / `pc.B`.

Grid metadata available on each `Fieldmap` subclass:
- `zbegin_m`, `zend_m` — longitudinal bounds
- `rbegin_m`, `rend_m` — radial bounds (or x-bounds for G4BL)
- `hr_m`, `hz_m` — grid spacing
- `num_gridpr_m`, `num_gridpz_m` — grid point counts
- `FieldstrengthBz_m`, `FieldstrengthBr_m`, etc. — `DualView<double*>` with
  the actual field values (host-accessible via `.view_host()`)

#### b) Analytic elements (no grid)

Elements like `Multipole` compute fields from a closed-form formula
(multipole expansion, hard-edge model, etc.) at each particle position. No
intermediate grid exists.

```cpp
// Multipole::apply(pc) runs:
Kokkos::parallel_for("Multipole::apply()", nLocal, KOKKOS_LAMBDA(const size_t i) {
    computeField(Rview(i), Ef, Bf, normalComponents, skewComponents, ...);
    Eview(i) += Ef;
    Bview(i) += Bf;
});
```

### 1.3 How fields reach the pusher

```
Self-fields:   particles → scatter → mesh (rho) → Poisson → E/B on mesh
               → gather → pc.E, pc.B

External:      element.apply(pc)  →  Kokkos kernel at R(i)
               →  accumulate into pc.E, pc.B

Pusher:        reads pc.E, pc.B (already merged) → Boris kick
```

Both contributions are merged into `pc.E` / `pc.B` before the pusher runs.
The old `FieldContainer::getExternalE()` / `getExternalB()` (separate
`ippl::Field`s on the PIC mesh) no longer exist.

---

## 2. Proposed design: per-element Catalyst registration

### 2.1 Overview

Each element registers its field with the Catalyst adaptor during
`goOnline()` (or a dedicated `registerVizChannel()` method called from
`ParallelTracker` after the beamline is prepared). The `UpdatePolicy`
mechanism already in the IPPL Catalyst API handles the static vs. dynamic
distinction with zero per-step cost for static meshes.

### 2.2 Static fieldmaps (magnetostatic, electrostatic)

- Register once at start with `UpdatePolicy::Static`.
- Export: grid coordinates (from `zbegin`, `zend`, `hr`, `hz`,
  `num_gridpr`, `num_gridpz`), field arrays (`Bz`, `Br`, etc. from
  `DualView::view_host()`).
- The element's local-to-lab coordinate transform (`getCSTrafoLab2Local()`)
  is exported as a `TransformData` on the channel, identical to the particle
  bunch transform.
- Zero per-step cost — the `Static` policy skips the deep-copy on every
  `Execute()`.

Channel naming: `element_<name>_<field_type>`, e.g.
`element_Q1_magnetic`, `element_CAV1_EZ`.

### 2.3 RF fieldmaps (dynamic)

Two options:

**Option A — Export grid once, scalar per step:**
- Register the grid with `UpdatePolicy::Static` (same as above).
- Also register a small `UpdatePolicy::Always` channel carrying the scalar
  amplitude (`V_peak × cos(ωt + φ)`). The Python script applies this scale
  when rendering.

**Option B — Export scaled grid per step:**
- Register the grid with `UpdatePolicy::Always`.
- Each step: multiply the static grid by the current scalar and deep-copy.
- Higher per-step cost but simpler script-side logic.

Option A is preferred — it reuses the `Static` policy for the grid (zero
cost) and the scalar is tiny.

### 2.4 Analytic elements (no grid)

These require a **user-defined sampling mesh**:

- **Particle-bounded mesh:** Create a grid covering the current particle
  bounding box. Sample `getFieldstrength(R, E, B)` at each grid point.
  Requires `UpdatePolicy::Always` (mesh bounds change each step) or
  `UpdatePolicy::Manual` with periodic refresh.

- **Fixed region mesh:** User specifies a region (e.g. around a specific
  element or a slice through the apparatus). Grid is static; field values
  are sampled once (static elements) or per step (time-varying elements).

- **Hybrid:** A fixed coarse mesh covering the full apparatus, refined
  locally around the bunch. This is the most complex but gives the best
  visual coverage.

For a first implementation, the particle-bounded mesh with
`UpdatePolicy::Always` is the simplest. The sampling loop calls
`element->getFieldstrength(gridPoint, E, B)` on the host (this method is
already virtual on `Fieldmap` and available on analytic elements via
`ElementBase::apply(R, P, t, E, B)`).

### 2.5 Registration flow

```
ParallelTracker::execute()
  ├── (existing) addParticleChannel("MUON", ...)
  ├── (existing) addMeshChannel("bunch_mesh", ...) [self-fields]
  │
  └── (new) for each element in itsOpalBeamline_m:
        if element has fieldmap:
          element->registerVizChannel(cat_vis);  // static grid, Static policy
        elif element has analytic field and user requests sampling:
          element->registerVizChannel(cat_vis, samplingMeshConfig);
```

Elements that carry no field (drifts, markers, monitors) are skipped.

### 2.6 What `ElementBase` needs

A new virtual method (default no-op):

```cpp
virtual void registerVizChannel(
    ippl::CatalystAdaptor& adaptor,
    const VizMeshConfig& config = {}) {}
```

Subclasses with field maps override it to:
1. Read grid dims, spacing, bounds from their `Fieldmap*` member.
2. Get the host view of field arrays.
3. Call `adaptor.addMeshChannel<Dim>(name, mesh, layout, nghost)`.
4. Call `meshCh.addArray(fieldName, fieldOrView, UpdatePolicy::Static)`.
5. Set the element-to-lab transform.

For analytic elements, the override samples `getFieldstrength()` onto the
provided `config.mesh` and registers with the appropriate policy.

### 2.7 Conduit channel structure

Each registered element becomes a Conduit channel of type `"mesh"` (for 2D
field maps) or `"multimesh"` (if we combine multiple field components).
The Python state script (`catalyst-muon-combined.py` or a new
`catalyst-elements.py`) binds to these channels by `registrationName` and
applies the element-to-lab transform for correct placement in the scene.

For 2D field maps (r, z), the channel carries:
- `coordsets`: uniform grid in (r, z)
- `topologies`: `unstructured` with the field map's cell layout
- `fields`: `Bz`, `Br` (or `Ez`, `Er`) as point/cell data
- `transform`: element-local → lab (for ParaView placement)

---

## 3. What needs to change in the existing code

### 3.1 After rebasing `catalyst-viz-v2-jb` onto `origin/master`

| File | Change |
|------|--------|
| `ParallelTracker.h` | `SelfFieldSystem*` → `SpaceChargeSolver*` (already done on master); add `#include` for solver header |
| `ParallelTracker.cpp` (self-field export) | Replace `itsBunch_m->getFieldContainer()` with access through `spaceChargeSolver_m` → algorithm → `CartesianPIC3DFieldStorage`. Need public accessors on `SpaceChargeSolver` and `CartesianPIC3DAlgorithm`. |
| `ParallelTracker.cpp` (external field export) | New: iterate `itsOpalBeamline_m.getElements()`, call `element->registerVizChannel(cat_vis)` |
| `SpaceCharge/SpaceChargeSolver.h` | Add `CartesianPIC3DFieldStorage& fieldStorage()` public accessor (or a dedicated `vizFields()` method) |
| `SpaceCharge/CartesianPIC3D/CartesianPIC3DAlgorithm.h` | Add `CartesianPIC3DFieldStorage& fieldStorage()` public accessor |
| `AbsBeamline/ElementBase.h` | Add `virtual void registerVizChannel(...)` (default no-op) |
| `AbsBeamline/Multipole.h`, `RFCavity.h`, etc. | Override `registerVizChannel()` for field-map and analytic elements |
| `catalyst/catalyst-muon-combined.py` | Add element field channel bindings; drop `E_external`/`B_external` from `bunch_mesh` |

### 3.2 Self-field adaptation (minimal)

The current code:
```cpp
auto fcontainer = itsBunch_m->getFieldContainer();
auto meshCh = cat_vis.addMeshChannel<3>("bunch_mesh",
    fcontainer->getMesh(), fcontainer->getFL(), fcontainer->getE().getNghost());
meshCh.addArray("E_self", fcontainer->getE());
...
```

After rebase, becomes something like:
```cpp
// Mesh + layout from PartBunch's CartesianDomain
auto& domain = itsBunch_m->cartesianDomain();
auto meshCh = cat_vis.addMeshChannel<3>("bunch_mesh",
    domain.mesh(), domain.layout(), /* nghost from field storage */);

// Fields from the solver's field storage
auto& fs = spaceChargeSolver_m->fieldStorage();
meshCh.addArray("E_self", fs.electricField());
meshCh.addArray("rho", fs.chargeDensity());
if (solverType == CG) meshCh.addArray("phi", fs.potential());
// B_self from accumulatedMagneticField() if needed
```

### 3.3 External field export (new)

For static fieldmaps:
```cpp
void RFCavity::registerVizChannel(ippl::CatalystAdaptor& adaptor) {
    if (!fieldmap_m) return;
    // Get grid metadata + host views from fieldmap_m
    // Register as Static mesh channel with element-to-lab transform
}
```

For analytic elements (if requested by user):
```cpp
void Multipole::registerVizChannel(
    ippl::CatalystAdaptor& adaptor, const VizMeshConfig& config) {
    // Create sampling grid from config
    // For each grid point: getFieldstrength(R, E, B)
    // Register as Always or Manual channel
}
```

---

## 4. Open questions

1. **Access to `CartesianPIC3DFieldStorage`:** The solver → algorithm →
   field storage chain is private. We need public accessors, or
   alternatively a `vizFields()` callback interface that the Catalyst
   adaptor calls. The latter avoids exposing internal solver types.

2. **2D field maps in ParaView:** The field map grids are 2D (r, z). ParaView
   can render 2D datasets, but they need to be placed in 3D using the
   element's coordinate transform. The existing `TransformData` mechanism
   (used for particle bunches) should work, but needs verification with
   non-particle mesh channels.

3. **RF cavity time-varying amplitude:** If we export the grid as `Static`
   and the amplitude as a separate small channel, the Python script needs
   to combine them. Alternatively, the Catalyst script can read the
   frequency/phase from the input file and compute the scale itself.

4. **Performance:** Static field maps are zero-cost per step. Analytic
   element sampling on a user-defined grid has a per-step cost proportional
   to the grid size × number of elements. This should be optional and
   coarsely configurable.

5. **Multi-rank:** Field maps are replicated on all ranks (each element owns
   its own copy). The Catalyst adaptor only exports from rank 0, so this is
   not an issue. The self-field PIC mesh is decomposed and already handled
   by the existing `MergeBlocks` logic in the Python script.

---

## 5. References

- PR #555 (`45e126f3a8e7`): "Big FieldSolver refactor" — deleted
  `FieldContainer.hpp`, `BinnedFieldSolver`, `FieldSolver`, moved PIC
  fields into `SpaceCharge/CartesianPIC3D/CartesianPIC3DFieldStorage`.
- `catalyst-editing.md`: Workflow for building and iterating on Catalyst
  state scripts.
- `catalyst-viz-plan.md`: Overall plan for the mesh-centric channel API.
- IPPL `catalyst-viz-dev` branch: `Channel.h/.hpp` (mesh/particle channels,
  `UpdatePolicy`), `CatalystAdaptor.h/.hpp` (registration + execute).
