# External Field Visualization — Design Notes

## Context

The `catalyst-viz-v2-jb` branch adds Catalyst in-situ visualization to OPALX.
Each particle container is published as a self-contained `particle_bunch_N`
multimesh channel with `particles`, `boundingbox`, `transform`, and `mesh`
blocks. This document captures the OPALX field architecture on `origin/master`
(PR #555, `45e126f3a8e7`) and how the Catalyst layer visualizes **external**
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

**Export status:** Works. `ParallelTracker.cpp` registers the self-field
mesh via `addMeshBlock("particle_bunch_0", "mesh", ...)` and adds
`E_self`, `rho` as arrays. Access goes through `spaceChargeSolver_m` →
`algorithm()` → `dynamic_cast<CartesianPIC3DAlgorithm*>` → `fieldStorage()`.

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

All 2D fieldmaps (`FM2DMagnetoStatic`, `FM2DDynamic`, `G4BL2DMagnetoStatic`)
are **cylindrically symmetric** — confirmed from their `isInside()`
implementations (`sqrt(r(0)² + r(1)²) < rend_m`) and `computeField()`
kernels (`RR = sqrt(R(0)² + R(1)²)`, radial component distributed as
`Bx = Br * R(0)/RR`, `By = Br * R(1)/RR`).

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

## 2. User-facing options

Two independent, composable options control how external element fields
are visualized. Both can be enabled simultaneously; each produces its own
set of Catalyst channels.

| Option | What it does | When to use |
|--------|-------------|-------------|
| **ExternalFieldsOnBunchMesh** | Sample all external fields onto the existing bunch PIC mesh. One mesh, resampled every step. | First implementation. Simple, one mesh, works for all element types. |
| **ExternalFieldsOnElementMesh** | Generate a per-element mesh for each element, place fields on it. Fieldmap elements export native (r, z) grid + transform; analytic elements get a local sampling mesh. | Higher fidelity per element; native fieldmap resolution; 2D cylindrical rendering. |

### 2.1 ExternalFieldsOnBunchMesh (implement first)

**One mesh, all fields.** The existing `bunch_mesh` channel (3D uniform
Cartesian, decomposed across ranks) is reused. External fields are
sampled onto it at each step, so they are always correct regardless of
bunch motion or mesh repartition.

**How it works:**

1. **Registration (one-time, before tracking loop):**
   - `ElementVizRegistry::createFields()` iterates
     `itsOpalBeamline_m.getElements()`.
   - For each element with a field (fieldmap or analytic), the registry
     creates two `VectorField` objects (`E` and `B`) on the bunch mesh +
     layout. These are owned by the registry and live for the duration
     of tracking.
   - `forEachField()` is called from `ParallelTracker.cpp` to register
     each field as an array (`E_ext_<name>`, `B_ext_<name>`) on the
     `mesh` block via `addMeshBlock()`.
   - `UpdatePolicy::Always` for all external field arrays, since they
     must be resampled every step to stay correct as the bunch moves.

2. **Per-step sampling (inside tracking loop, before `cat_vis.Execute()`):**
   - `ElementVizRegistry::sampleDynamic()` iterates all registered
     elements, zeroing both E and B before accumulating.
   - For each element, for each local grid point `R`:
     ```
     E = 0, B = 0
     element->apply(localR, P, t, localE, localB)  // bounds-checked
     refE = localToRefCSTrafo.rotateTo(localE)
     refB = localToRefCSTrafo.rotateTo(localB)
     field(i,j,k) = refE  // or refB for B fields
     ```
   - Elements whose longitudinal extent doesn't overlap the bunch
     bounding box are skipped (coarse cull).
   - `element->apply(R, P, t, E, B)` has internal bounds checking —
     points outside the element return zero.

3. **What the Python state script sees:**
   - `particle_bunch_0` multimesh with `particles`, `boundingbox`,
     `transform`, and `mesh` blocks.
   - `mesh` block cell data: `E_self`, `rho` (self-fields, if solver
     active) + `E_ext_<name>`, `B_ext_<name>` (one pair per element
     with a field).

**Pros:**
- One mesh, already decomposed, already in the `MergeBlocks` pipeline.
- No per-element grid generation, no per-element mesh ownership.
- Works for all element types (fieldmap, analytic, static, dynamic).
- External fields are always correct (resampled every step).

**Cons:**
- Resolution is limited to the bunch mesh resolution (may be coarser
  than a fieldmap's native grid).
- Per-step sampling cost for all elements near the bunch (proportional
  to `N_active_elements × N_local_grid_points`).
- Fieldmap elements are sampled even though their native grid is already
  available (wasteful but correct).

### 2.2 ExternalFieldsOnElementMesh (planned, not yet implemented)

**Per-element meshes, native fidelity.** Each element gets its own
Catalyst channel with a mesh appropriate to its field type.

#### a) Fieldmap elements (2D cylindrical, static grid)

- Export the native `(r, z)` grid directly via `MeshChannelT<2>`.
- Grid metadata: `zbegin`, `zend`, `hr`, `hz`, `num_gridpr`,
  `num_gridpz` from the `Fieldmap` subclass.
- Field arrays: `Bz`, `Br` (or `Ez`, `Er`, `Bt`) from
  `DualView::view_host()` — zero per-step cost for magnetostatic.
- Element-to-lab transform: `getCSTrafoLab2Local()` exported as
  `TransformData` on the channel, same mechanism as particle bunches.
- `UpdatePolicy::Static` for magnetostatic / electrostatic fieldmaps.
- For RF fieldmaps: grid `Static` + a small `Always` channel carrying
  the scalar amplitude (`V_peak × cos(ωt + φ)`). The Python script
  applies the scale when rendering. Alternatively, export the scaled
  grid per step (`Always`) — simpler script, higher cost.

**Python-side rendering with `vtkCylindricalTransform`:**
- ParaView's `vtkCylindricalTransform` converts `(r, φ, z)` → `(x, y, z)`.
- The state script applies `RotationalExtrusion` (angle = 360°, axis = z)
  to the 2D `(r, z)` dataset, producing a 3D volume.
- The element's lab transform (`TransformData`) places the extruded
  volume in the correct position and orientation in the scene.
- All 2D fieldmaps are cylindrically symmetric (see section 1.2a).

**Channel naming:** `element_<name>_<field_type>`, e.g.
`element_Q1_magnetic`, `element_CAV1_EZ`.

#### b) Analytic elements (no native grid)

- The registry generates a mesh in the element's local coordinate space,
  bounded by the element's field extent.
- For static analytic elements (e.g. DC multipole): sample once at
  registration, `UpdatePolicy::Static`. The registry owns a frozen
  `UniformCartesian<double, 3>` mesh (origin/spacing fixed at
  registration time, immune to bunch mesh motion).
- For time-varying analytic elements: sample every step,
  `UpdatePolicy::Always`. The registry resamples using
  `element->getFieldstrength(R, t, E, B)` at each local grid point.
- `GeometryPolicy::Static` (new, see section 3) prevents per-step
  geometry updates for the frozen mesh.

#### c) Open questions for this option

- **2D fieldmap extrusion:** Export 2D + transform and extrude in Python
  (`RotationalExtrusion`) vs. extrude in C++ and export 3D. Python
  extrusion is preferred for the first implementation.
- **Analytic element mesh resolution:** Match the bunch mesh spacing?
  User-configurable? Element-specific?
- **Performance:** Per-element sampling for time-varying analytic elements
  is `N_grid × N_elements` per step. Should be optional and coarsely
  configurable.

---

## 3. IPPL Catalyst API additions

### 3.1 GeometryPolicy

A new viz-only policy on `MeshChannelT`, separate from the field-data
`UpdatePolicy`:

```cpp
enum class GeometryPolicy {
    Dynamic,  // Re-read mesh origin/spacing every Execute (default,
              //   current behavior — bunch may have moved)
    Static,   // Mesh origin/spacing set once at init, never re-read.
              //   Zero per-step geometry cost.
};
```

**Why viz-only:** The OPALX solver has no analogous concept — it always
recomputes everything. `GeometryPolicy` lives entirely in the IPPL
Catalyst layer (`Channel.h` / `Channel.hpp`) and controls only whether
`MeshChannelT::execute()` re-reads `mesh_m.getOrigin()` and
`mesh_m.getMeshSpacing()` on each step (lines 357-393 of `Channel.hpp`).

**Combination matrix:**

| GeometryPolicy | UpdatePolicy | Use case |
|----------------|-------------|----------|
| Static  | Static  | Static fieldmap native grid (zero cost) |
| Static  | Always  | Frozen mesh, resampled field (analytic static element) |
| Dynamic | Static  | Bunch mesh, static field (not typical) |
| Dynamic | Always  | Bunch mesh, resampled field (ExternalFieldsOnBunchMesh) |

### 3.2 TransformData on mesh channels

`TransformData` (already defined in `Channel.h:45-50`) is currently used
only by `ParticleChannel`. For `ExternalFieldsOnElementMesh`, mesh
channels also need a transform (element-local → lab).

**Option A:** Add `setTransform()` to `MeshChannelT` (mirrors
`ParticleChannel::setTransform()`). The Conduit node gets a
`block_transform` block, same as particle channels.

**Option B:** The Python state script reads the element transform from
a separate metadata channel and applies it in the pipeline.

Option A is cleaner and consistent with the particle channel approach.

---

## 4. ElementVizRegistry

### 4.1 Purpose

A dedicated OPALX-side class that owns visualization fields for external
element visualization. Lives in `src/Algorithms/ElementVizRegistry.h/.cpp`
(not in IPPL, not in `ElementBase` subclasses).

```cpp
class ElementVizRegistry {
public:
    using Mesh_t      = ippl::UniformCartesian<double, 3>;
    using Layout_t    = ippl::FieldLayout<3>;
    using Centering_t = Mesh_t::DefaultCentering;
    using VectorField = ippl::Field<ippl::Vector<double, 3>, 3, Mesh_t, Centering_t>;

    /// One-time creation of E and B fields on the bunch mesh.
    /// Caller registers each field via forEachField().
    void createFields(OpalBeamline& beamline, PartBunch_t& bunch);

    /// Per-step sampling of external fields.
    void sampleDynamic(
        OpalBeamline& beamline, PartBunch_t& bunch,
        const Vector_t<double, 3>& bunchBBoxMin,
        const Vector_t<double, 3>& bunchBBoxMax);

    /// Iterate over all registered element fields.
    /// Callback receives (name, field) for each E_ext and B_ext.
    template <typename Fn>
    void forEachField(Fn&& fn);

    size_t size() const;

private:
    struct ElementField {
        std::string name;               // element name (e.g. "Q01", "SOL")
        std::shared_ptr<ElementBase> element;
        VectorField E;                  // electric field on bunch mesh
        VectorField B;                  // magnetic field on bunch mesh
    };
    std::vector<ElementField> elementFields_m;
};
```

### 4.2 Why the registry, not ElementBase

- **ElementBase** is a physics class. Adding `registerVizChannel()`
  pollutes it with visualization concerns and creates a dependency from
  `AbsBeamline` to `ippl::CatalystAdaptor`.
- **PartBunch** already owns the self-field domain, but stretching it
  to own per-element sampling meshes is a leaky abstraction.
- **ElementVizRegistry** keeps viz code in one place, owns the
  `VectorField` objects (required by `MeshFieldArrayImpl` which holds
  `const Field_type&`), and can query the orbit tracker's active set
  to skip elements far from the bunch.
- **CatalystAdaptor.h** cannot be included in `ElementVizRegistry.cpp`
  (non-template method definitions cause multiple definition errors), so
  registration is done from `ParallelTracker.cpp` via `forEachField()`.

### 4.3 Integration with ParallelTracker

```
ParallelTracker::execute()
  ├── prepareSections()
  ├── OrbitThreader::execute()
  ├── globalBoundingBox computed
  │
  ├── #ifdef IPPL_ENABLE_CATALYST
  │     cat_vis declared
  │     ├── addParticleChannel("particle_bunch_0", pc).setTransform(td)
  │     ├── addMeshBlock("particle_bunch_0", "mesh", domain.mesh(), ...)
  │     │     .addArray("E_self", fs.electricField())
  │     │     .addArray("rho", fs.chargeDensity())
  │     ├── elementVizRegistry.createFields(itsOpalBeamline_m, *itsBunch_m)
  │     ├── elementVizRegistry.forEachField([&](name, field) {
  │     │     meshCh.addArray(name, field);
  │     │   })
  │     └── cat_vis.Initialize(steerReg)
  │
  └── main tracking loop
        └── (before cat_vis.Execute()):
              elementVizRegistry.sampleDynamic(
                itsOpalBeamline_m, *itsBunch_m, rmin, rmax)
              cat_vis.Execute(step, step)
```

### 4.4 Orbit tracker integration

`OrbitThreader` (`src/Algorithms/OrbitThreader.cpp:126`) computes which
elements are active at a given reference position via
`IndexMap::value_t activeSet`. The registry queries this to skip
elements far from the bunch — their external field arrays stay zero,
which is correct and free.

For the first implementation (ExternalFieldsOnBunchMesh), the simplest
approach is to sample all elements whose bounds intersect the bunch
bounding box, using the existing `globalBoundingBox` computed at
line 395-400. Tighter orbit-tracker-based culling can be added later.

---

## 5. Channel restructuring

### 5.1 Background: why mesh and particles were separate

The `e2230386699a` refactor (Sep 24) created `MeshChannelT<Dim>` and
`ParticleChannel` as independent classes with different Conduit layouts.
The mesh was kept in a separate channel as a **code structure artifact**,
not by design intent. Three weeks earlier, commit `1a89ca349ca6` (Sep 4)
had deliberately **embedded the transform into the particle multimesh**
with the stated goal: "This keeps each bunch self-contained: particles,
helper box, and transform travel together as one dataset collection."

The mesh was never merged because `MeshChannelT` was a new general-purpose
class and nobody implemented the merge path — not because it was
considered and rejected. The IPPL demos (PenningTrap by SKlapproth
`6776302453d3`, FEL by SKlapproth `c559dcae4f20`) both used the old flat
registry API with no separate mesh channel concept at all.

### 5.2 New structure: one multimesh per bunch (implemented)

Each particle container is published as a self-contained
`particle_bunch_N` multimesh channel with `particles`, `boundingbox`,
`transform`, and `mesh` blocks.

**Channel naming:** `particle_bunch_0`, `particle_bunch_1`, ...
(indexed by particle container, not simulation-specific names like
`MUON`).

**Block naming** (dropping the redundant `block_` prefix):

| Old | New | Content |
|-----|-----|---------|
| `block_main` | `particles` | Particle positions, IDs, attributes |
| `block_help` | `boundingbox` | 2×2×2 uniform grid covering particle domain |
| `block_transform` | `transform` | Single-point mesh with origin/rotation vectors |
| *(separate channel)* | `mesh` | Self-field PIC grid (E_self, B_self, rho, phi, ...) |

**Conduit structure after restructuring:**
```
catalyst/channels/particle_bunch_0
  type = "multimesh"
  assembly/particles   = "particles"
  assembly/boundingbox = "boundingbox"
  assembly/transform   = "transform"
  assembly/mesh        = "mesh"

  data/particles/
    coordsets/p_explicit_coords/   (explicit points: x, y, z)
    topologies/p_unstructured_topo (point cells)
    fields/position, RankID, ParticleIDs, <user attrs>

  data/boundingbox/
    coordsets/bound_helper_coords/ (uniform, 2x2x2)
    topologies/bound_helper_topo

  data/transform/
    coordsets/tf_coords/           (single point)
    topologies/tf_topo
    fields/BunchOrigin, BunchRotRow0, ..., BunchInvRotRow2

  data/mesh/                       ← was separate "bunch_mesh" channel
    coordsets/cart_uniform_coords/ (uniform grid, per-rank partition)
    topologies/fmesh_topo
    fields/E_self, B_self, rho, phi, <per-element ext fields>, RankID
```

### 5.3 IPPL API: `MeshChannelT` basePath (implemented)

`MeshChannelT` has a `basePath` parameter. When non-empty, it writes
into `catalyst/channels/<basePath>/data/<name>` instead of a top-level
channel. `MeshChannelT::init()` sets `assembly/<name> = <name>` on the
parent channel so the block appears in the DataAssembly.

```cpp
// ParallelTracker.cpp (implemented)
auto handle = cat_vis.addParticleChannel("particle_bunch_0", pc).setTransform(td);

auto meshCh = cat_vis.addMeshBlock<3>(
    "particle_bunch_0", "mesh",
    domain.mesh(), domain.layout(), /*nghost=*/1, /*useGhostMasks=*/false);

meshCh.addArray("E_self", fs.electricField());
meshCh.addArray("rho", fs.chargeDensity());

elementVizRegistry_m->forEachField(
    [&meshCh](const std::string& name, auto& field) {
        meshCh.addArray(name, field);
    });
```

### 5.4 `pipeline_default.py` dispatch (updated)

Dispatches on `ctype`:
- `"multimesh"` → `ExtractBlock` for `//particles`, `//boundingbox`,
  `//transform`; also `//mesh` if present (extract mesh block, apply
  `MergeBlocks` + `CellDatatoPointData`).

### 5.5 VTK extract format (implemented)

- `particle_bunch_0` (type `multimesh`) → one `.vtpc` for everything
- The mesh block is inside the `.vtpc`, not a separate `.vtpd`
- Offline workflow: one `.vtpc` per bunch per step

### 5.6 Parallel data flow (verified)

Each rank contributes its particles + its mesh partition as separate
blocks in the multimesh. ParaView's Catalyst connector creates:

```
particle_bunch_0 (vtkPartitionedDataSetCollection):
  block "particles"   → vtkPDS, N partitions (vtkUnstructuredGrid — points)
  block "boundingbox" → vtkPDS, N partitions (vtkImageData — 2x2x2)
  block "transform"   → vtkPDS, N partitions (vtkUnstructuredGrid — point)
  block "mesh"        → vtkPDS, N partitions (vtkImageData — local PIC grid)
```

`MergeBlocks(MergePartitionsOnly=1)` flattens per-rank partitions into
a single dataset for visualization. Multi-rank (2 ranks) verified:
50 steps, 0 Catalyst errors, all four assembly entries present in
every snapshot.

### 5.7 What does NOT break

- IPPL core mechanism — persistent mirrors, `UpdatePolicy`,
  `GeometryPolicy` (to be added) work regardless of channel grouping
- `ParticleAttribBase::prepareConduitNode/updateConduitData` — particle
  attribute serialization is independent of channel structure
- `OrbitThreader` active set — independent of visualization
- Self-field access chain (after rebase) — `spaceChargeSolver_m →
  algorithm → fieldStorage` is independent of channel grouping
- `ElementVizRegistry` design — can work with either separate or merged
  channels
- Conduit blueprint compliance — multimesh with mixed topology types
  (unstructured points + uniform grid) is valid; each block is an
  independent mesh

---

## 6. Implementation status

### Phase 1: Rebase + restructure + ExternalFieldsOnBunchMesh (DONE)

**IPPL branch (`catalyst-viz-dev`):**
1. ✅ Block renaming: `block_main` → `particles`, `block_help` →
   `boundingbox`, `block_transform` → `transform`.
2. ✅ `basePath` parameter on `MeshChannelT`; writes into
   `catalyst/channels/<basePath>/data/<name>`.
3. ✅ `addMeshBlock(parentName, blockName, mesh, layout, nghost,
   useGhostMasks, geometryPolicy)` on `CatalystAdaptor`.
4. ✅ `GeometryPolicy` enum (`Dynamic`/`Static`); skips geometry update
   in `execute()` when `Static`.
5. ✅ Assembly entries (`particles`, `boundingbox`, `transform`, `mesh`)
   set on both `init()` and re-init paths (after particle channel reset).
6. ✅ `MeshChannelT::execute()` detects wiped node (missing `type`)
   and rebuilds from scratch, including the assembly entry.

**OPALX branch (`catalyst-viz-v2-jb`):**
1. ✅ Rebased onto `origin/master` (`45e126f3a8e7`).
2. ✅ Self-field export via `spaceChargeSolver_m` → `algorithm()` →
   `dynamic_cast<CartesianPIC3DAlgorithm*>` → `fieldStorage()`.
3. ✅ Public accessors on `SpaceChargeSolver` and
   `CartesianPIC3DAlgorithm`.
4. ✅ Removed `E_external`/`B_external` (no longer exist on PIC mesh).
5. ✅ Channel names: `particle_bunch_0` (indexed by container).
6. ✅ Self-field mesh merged into multimesh via `addMeshBlock()`.
7. ✅ `ElementVizRegistry` (`src/Algorithms/ElementVizRegistry.h/.cpp`):
   - `createFields()`: creates `E` and `B` VectorFields per element.
   - `forEachField()`: registers `E_ext_<name>` and `B_ext_<name>`.
   - `sampleDynamic()`: per-step sampling via `element->apply()`.
8. ✅ Python state scripts updated (selectors, registrationName).
9. ✅ Single-rank smoke test: 50 steps, 0 Catalyst errors.
10. ✅ Multi-rank (2 ranks): 50 steps, 0 Catalyst errors.

### Phase 2 (future): ExternalFieldsOnElementMesh

1. **Fieldmap element export** — native (r, z) grid + transform.
   - `ElementVizRegistry::registerFieldmapElement()`: create
     `MeshChannelT<2>` with grid from `Fieldmap` metadata, add arrays
     from `DualView::view_host()`, set transform from
     `getCSTrafoLab2Local()`.
   - `UpdatePolicy::Static` for magnetostatic / electrostatic.
   - RF fieldmaps: grid `Static` + scalar `Always` channel (preferred),
     or grid `Always` with scaled values.

2. **Python-side cylindrical rendering**
   - State script applies `vtkCylindricalTransform` +
     `RotationalExtrusion` (angle = 360°) to 2D datasets.
   - Element lab transform places extruded volume in scene.

3. **Analytic element per-element mesh**
   - Registry generates local mesh bounded by element field extent.
   - `GeometryPolicy::Static` for time-independent analytic elements
     (sample once).
   - `UpdatePolicy::Always` for time-varying analytic elements.

4. **User configuration** — input file or environment variables to
   select which elements to visualize, mesh resolution, etc.

---

## 7. What changed in the existing code (Phase 1)

| File | Change |
|------|--------|
| IPPL `Channel.h` | Block renaming, `GeometryPolicy` enum, `MeshChannelT` basePath parameter. |
| IPPL `Channel.hpp` | Block renaming, skip geometry update when `GeometryPolicy::Static`; re-init mesh block + assembly entry when parent channel resets. |
| IPPL `CatalystAdaptor.h` | `addMeshBlock(parentName, blockName, ...)` declaration. |
| IPPL `CatalystAdaptor.hpp` | `addMeshBlock()` implementation; block renaming in particle init/execute fns; `assembly/transform` added to new API path. |
| IPPL `pipeline_default.py` | Block selectors updated: `//particles`, `//boundingbox`, `//transform`, `//mesh`. |
| `ParallelTracker.h` | Forward declaration of `opalx::ElementVizRegistry`; `std::unique_ptr<opalx::ElementVizRegistry> elementVizRegistry_m;` |
| `ParallelTracker.cpp` (channel names) | `MUON` → `particle_bunch_0`; merge mesh via `addMeshBlock()`. |
| `ParallelTracker.cpp` (self-field) | Replace `getFieldContainer()` with `spaceChargeSolver_m` → `algorithm()` → `dynamic_cast<CartesianPIC3DAlgorithm*>` → `fieldStorage()`. |
| `ParallelTracker.cpp` (external fields) | `ElementVizRegistry::createFields()` + `forEachField()` before tracking loop; `sampleDynamic()` inside loop before `cat_vis.Execute()`. |
| `SpaceCharge/SpaceChargeSolver.h` | Added `SpaceChargeAlgorithm& algorithm()` public accessor. |
| `SpaceCharge/CartesianPIC3D/CartesianPIC3DAlgorithm.h` | Added `FieldStorage& fieldStorage()` public accessor. |
| **NEW** `src/Algorithms/ElementVizRegistry.h` | Registry class declaration. |
| **NEW** `src/Algorithms/ElementVizRegistry.cpp` | Registry class implementation — field ownership, sampling. |
| `catalyst/catalyst-muon-combined.py` | `registrationName='particle_bunch_0'`, selectors `//particles`/`//boundingbox`/`//transform`/`//mesh`. |
| `catalyst/catalyst-muon-test-1.py` | Same renaming. |
| `catalyst/catalyst-muon-fields.py` | Updated to extract mesh block from multimesh. |

---

## 8. Open questions

1. **Self-field access:** Resolved. `SpaceChargeSolver::algorithm()` and
   `CartesianPIC3DAlgorithm::fieldStorage()` are public accessors.
   `ParallelTracker.cpp` uses `dynamic_cast` to access the storage.

2. **External fields on mesh block:** Resolved. External element fields
   (`E_ext_<name>`, `B_ext_<name>`) are added as arrays on the same
   `mesh` block as the self-fields. One mesh block, all field arrays.

3. **RF fieldmap amplitude (Phase 2):** If we export the grid as `Static`
   and the amplitude as a separate small channel, the Python script
   needs to combine them. Alternatively, the Catalyst script can read
   the frequency/phase from the input file and compute the scale itself.

4. **Multi-rank:** Resolved. Verified with 2 ranks, 50 steps, 0 Catalyst
   errors. All four assembly entries present in every snapshot.

5. **Sampling cost:** For ExternalFieldsOnBunchMesh, the per-step cost
   is `N_active_elements × N_local_grid_points` calls to
   `element->apply()`. For a typical bunch mesh of ~100³ local points
   and ~10 active elements, this is ~10⁷ evaluations per step —
   comparable to the existing particle push. Should be acceptable for
   visualization but worth measuring with a real solver.

---

## 9. References

- PR #555 (`45e126f3a8e7`): "Big FieldSolver refactor" — deleted
  `FieldContainer.hpp`, `BinnedFieldSolver`, `FieldSolver`, moved PIC
  fields into `SpaceCharge/CartesianPIC3D/CartesianPIC3DFieldStorage`.
- `catalyst-editing.md`: Workflow for building and iterating on Catalyst
  state scripts.
- `catalyst-viz-plan.md`: Overall plan for the mesh-centric channel API.
- IPPL `catalyst-viz-dev` branch: `Channel.h/.hpp` (mesh/particle channels,
  `UpdatePolicy`), `CatalystAdaptor.h/.hpp` (registration + execute).
- `FM2DMagnetoStatic.h:113-114`, `FM2DDynamic.h:88-89`,
  `G4BL2DMagnetoStatic.h:113-114`: All 2D fieldmaps are cylindrically
  symmetric (`isInside()` checks `sqrt(r(0)² + r(1)²) < rend_m`).
- `OrbitThreader.cpp:126`, `IndexMap.cpp:115`: Active element set
  computation for culling elements far from the bunch.
- IPPL commit `1a89ca349ca6` (Sep 4): "Embed bunch transform as a block
  in the particle multimesh" — design goal: "This keeps each bunch
  self-contained: particles, helper box, and transform travel together
  as one dataset collection."
- IPPL commit `e2230386699a` (Sep 24): "mesh-centric channel API" —
  created `MeshChannelT` and `ParticleChannel` as separate classes.
  Mesh kept separate as code structure artifact, not by design intent.
- IPIL demos: PenningTrap `6776302453d3` (SKlapproth), FEL
  `c559dcae4f20` (SKlapproth) — both used old flat registry API,
  migrated to channel API in `e2230386699a`.
