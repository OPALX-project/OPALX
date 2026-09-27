# External Field Visualization — Design Notes

## Context

The `catalyst-viz-v2-jb` branch adds Catalyst in-situ visualization to OPALX.
Particle channels (`MUON`, `BEAM1`, ...) and a self-field mesh channel
(`bunch_mesh`) are already working. This document captures what we learned
about the OPALX field architecture on `origin/master` (PR #555,
`45e126f3a8e7`) and how to extend the Catalyst layer to visualize **external**
element fields (field maps and analytic models).

**Planned restructuring** (section 5): merge the self-field mesh into the
particle multimesh, rename channels to `particle_bunch_0/1/...`, and rename
blocks to `particles`, `boundingbox`, `transform`, `mesh`. This makes each
bunch fully self-contained — consistent with the design goal stated in
IPPL commit `1a89ca349ca6`.

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
storage (see section 6).

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
   - `ElementVizRegistry` iterates `itsOpalBeamline_m.getElements()`.
   - For each element with a field (fieldmap or analytic), the registry
     creates an `ippl::Field<Vector_t, 3>` on the bunch mesh + layout.
     These `Field` objects are owned by the registry and live for the
     duration of tracking.
   - Each field is registered as a new array on the existing `bunch_mesh`
     channel (or a separate `external_fields` channel sharing the same
     mesh — TBD during implementation).
   - `UpdatePolicy::Always` for all external field arrays, since they
     must be resampled every step to stay correct as the bunch moves.

2. **Per-step sampling (inside tracking loop, before `cat_vis.Execute()`):**
   - The registry iterates the orbit tracker's active set
     (`IndexMap::value_t`) to skip elements far from the bunch.
   - For each active element, for each local grid point `R`:
     ```
     E = 0, B = 0
     element->getFieldstrength(R, E, B)  // or apply(R, P, t, E, B)
     // If R is outside element bounds, E and B stay zero
     externalField(i) = E + B  (as appropriate)
     ```
   - Inactive elements contribute zero (their arrays stay zero from the
     last clear).
   - A single `Kokkos::parallel_for` over all local grid points × active
     elements is possible but per-element calls are simpler initially.

3. **What the Python state script sees:**
   - `bunch_mesh` channel with `E_self`, `B_self`, `rho`, `phi` (existing)
   - `bunch_mesh` channel with `E_ext_Q1`, `B_ext_Q1`, `E_ext_CAV1`, ...
     (new — one pair per element)
   - Or a separate `external_fields` channel on the same mesh (cleaner
     for ParaView pipeline separation — decide during implementation).

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

A dedicated OPALX-side class that owns visualization meshes and fields
for external element visualization. Lives in the OPALX Catalyst layer
(not in IPPL, not in `ElementBase` subclasses).

```cpp
class ElementVizRegistry {
public:
    /// One-time registration of all external fields.
    /// Called from ParallelTracker::execute() before the tracking loop.
    void registerElements(
        OpalBeamline& beamline,
        PartBunch_t& bunch,
        ippl::CatalystAdaptor& cat_vis);

    /// Per-step sampling of dynamic external fields.
    /// Called from ParallelTracker inside the tracking loop, before
    /// cat_vis.Execute().
    void sampleDynamic(
        OpalBeamline& beamline,
        const BoundingBox& bunchBBox,
        double time);

private:
    /// Owned ippl::Field objects for external fields on the bunch mesh.
    /// Key: element name, Value: sampled E/B field on bunch mesh.
    std::unordered_map<std::string, Field3D_t> externalFields_m;

    /// For ExternalFieldsOnElementMesh (future):
    /// Per-element owned meshes and fields.
    // std::unordered_map<std::string, ElementVizData> elementMeshes_m;

    /// Active element tracking (from orbit tracker / IndexMap)
    // IndexMap::value_t lastActiveSet_m;
};
```

### 4.2 Why the registry, not ElementBase

- **ElementBase** is a physics class. Adding `registerVizChannel()`
  pollutes it with visualization concerns and creates a dependency from
  `AbsBeamline` to `ippl::CatalystAdaptor`.
- **PartBunch** already owns the self-field domain, but stretching it
  to own per-element sampling meshes is a leaky abstraction.
- **ElementVizRegistry** keeps viz code in one place, owns the `ippl::Field`
  objects (required by `MeshFieldArrayImpl` which holds
  `const Field_type&`), and can query the orbit tracker's active set to
  skip elements far from the bunch.

### 4.3 Integration with ParallelTracker

```
ParallelTracker::execute()
  ├── prepareSections()                           // line 236
  ├── OrbitThreader::execute()                    // line 385
  ├── globalBoundingBox computed                  // line 395-400
  │
  ├── #ifdef IPPL_ENABLE_CATALYST
  │     cat_vis declared
  │     ├── addMeshChannel<3>("bunch_mesh", ...)  // self-fields (existing)
  │     ├── addParticleChannel("MUON", ...)       // particles (existing)
  │     └── elementVizRegistry.registerElements(  // NEW
  │           itsOpalBeamline_m, itsBunch_m, cat_vis)
  │
  └── main tracking loop                          // line 527+
        └── (before cat_vis.Execute()):
              elementVizRegistry.sampleDynamic(
                itsOpalBeamline_m, globalBoundingBox, time)
              cat_vis.Execute()
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

### 5.2 New structure: one multimesh per bunch

Merging the self-field mesh into the particle multimesh so each bunch is
fully self-contained — particles, bounding box, transform, and mesh all
travel together as one Conduit channel.

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

### 5.3 IPPL API change: `MeshChannelT` basePath

**Approach B** (simplest): add a `basePath` parameter to `MeshChannelT`
so it writes into `catalyst/channels/particle_bunch_0/data/mesh/` instead
of `catalyst/channels/mesh/`. `MeshChannelT` and `ParticleChannel` stay
independent classes. The particle channel's init fn adds
`assembly/mesh = "mesh"` to the Conduit node. Init order: particle
channel first (creates the multimesh node), then the mesh block writes
into `data/mesh/`.

```cpp
// ParallelTracker.cpp (after restructuring)
auto pcHandle = cat_vis.addParticleChannel("particle_bunch_0", pc);
pcHandle.setTransform(td);

// Mesh block writes into the same multimesh
auto meshCh = cat_vis.addMeshBlock("particle_bunch_0", "mesh",
                                    domain.mesh(), domain.layout(), nghost);
meshCh.addArray("E_self", fs.electricField());
meshCh.addArray("rho", fs.chargeDensity());
// ...
```

### 5.4 `pipeline_default.py` dispatch change

Currently dispatches on `ctype`:
- `"multimesh"` → `ExtractBlock` for particles, helper, transform
- `"mesh"` → `MergeBlocks` + `CellDatatoPointData` + `ResampleToImage`

After restructuring, all channels are `"multimesh"`. The script needs a
new path: "multimesh that contains a mesh block → extract the mesh block,
apply mesh pipeline (`MergeBlocks` + `Cell2Point`); extract particle
block, apply particle pipeline."

### 5.5 VTK extract format change

Currently:
- `bunch_mesh` (type `mesh`) → `.vtpd` (Partitioned Dataset)
- `MUON` (type `multimesh`) → `.vtpc` (Partitioned Dataset Collection)

After:
- `particle_bunch_0` (type `multimesh`) → one `.vtpc` for everything
- The mesh block is inside the `.vtpc`, not a separate `.vtpd`
- Offline workflow: instead of separate `.vti` + `.vtpc` snapshots,
  there's one `.vtpc` per bunch per step

### 5.6 Parallel data flow

Each rank contributes its particles + its mesh partition as separate
blocks in the multimesh. ParaView's Catalyst connector collects these
and creates:

```
particle_bunch_0 (vtkPartitionedDataSetCollection):
  block "particles"   → vtkPDS, N partitions (vtkUnstructuredGrid — points)
  block "boundingbox" → vtkPDS, N partitions (vtkImageData — 2x2x2)
  block "transform"   → vtkPDS, N partitions (vtkUnstructuredGrid — point)
  block "mesh"        → vtkPDS, N partitions (vtkImageData — local PIC grid)
```

`MergeBlocks(MergePartitionsOnly=1)` already flattens per-rank
partitions into a single dataset for visualization. Multi-rank testing
is required to verify this works correctly with the merged structure.

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

## 6. Implementation phases

### Phase 1 (now): Rebase + restructure + ExternalFieldsOnBunchMesh

**Step 1: IPPL branch (`catalyst-viz-dev`) — independent of OPALX**

1. Block renaming: `block_main` → `particles`, `block_help` →
   `boundingbox`, `block_transform` → `transform`.
2. Add `basePath` parameter to `MeshChannelT` so it can write into a
   multimesh `data/mesh/` path instead of a top-level channel.
3. Add `addMeshBlock(parentName, blockName, mesh, layout, nghost)` to
   `CatalystAdaptor` — registers a mesh block within an existing
   particle channel's multimesh.
4. Add `GeometryPolicy` enum (`Dynamic`/`Static`) to `Channel.h`,
   `MeshChannelT` constructor parameter. Skip geometry update in
   `execute()` when `Static`. Default `Dynamic` preserves current
   behavior.

**Step 2: OPALX rebase onto `origin/master` (`45e126f3a8e7`)**
   - Fix self-field export: replace `getFieldContainer()` with access
     through `spaceChargeSolver_m` → algorithm → field storage.
   - Add public accessors on `SpaceChargeSolver` and
     `CartesianPIC3DAlgorithm` (or a `vizFields()` callback).
   - Drop `E_external` / `B_external` from the mesh (they no longer
     exist as separate `ippl::Field`s on the PIC mesh).

**Step 3: OPALX channel restructuring + ElementVizRegistry**
   - Rename channel `MUON` → `particle_bunch_0` (indexed by container).
   - Merge self-field mesh into the multimesh as `data/mesh/` block
     via `addMeshBlock()`.
   - Add `ElementVizRegistry`:
     - Owns `ippl::Field<Vector_t, 3>` objects on the bunch mesh, one
       per element with a field.
     - `registerElements()`: iterate beamline, create fields, register
       as arrays on the mesh block.
     - `sampleDynamic()`: iterate active elements, sample
       `getFieldstrength(R, E, B)` at each local grid point, zero
       outside bounds. Called every step.

**Step 4: Python state scripts**
   - `catalyst-muon-combined.py`: update `registrationName` to
     `particle_bunch_0`, selectors to `//particles`, `//boundingbox`,
     `//transform`, `//mesh`.
   - Drop `E_external`/`B_external` bindings.
   - Add per-element external field bindings.
   - Update `pipeline_default.py` dispatch for merged multimesh
     (extract mesh block, apply mesh pipeline).

**Step 5: Multi-rank testing**
   - Verify merged structure works correctly across N ranks.
   - Verify `MergeBlocks(MergePartitionsOnly=1)` flattens both particle
     and mesh partitions.
   - Verify bounding box and transform are consistent across ranks.

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

## 7. What needs to change in the existing code (Phase 1)

| File | Change |
|------|--------|
| IPPL `Channel.h` | Block renaming, `GeometryPolicy` enum, `MeshChannelT` basePath parameter. |
| IPPL `Channel.hpp` | Block renaming, skip geometry update when `GeometryPolicy::Static`. |
| IPPL `CatalystAdaptor.h` | `addMeshBlock(parentName, blockName, ...)` declaration. |
| IPPL `CatalystAdaptor.hpp` | `addMeshBlock()` implementation; block renaming in particle init/execute fns. |
| IPPL `pipeline_default.py` | Dispatch for merged multimesh (extract mesh block + mesh pipeline). |
| `ParallelTracker.h` | Add `ElementVizRegistry elementVizRegistry_m;` Add `#include` for registry header. |
| `ParallelTracker.cpp` (channel names) | `MUON` → `particle_bunch_0`; merge mesh via `addMeshBlock()`. |
| `ParallelTracker.cpp` (self-field) | Replace `getFieldContainer()` with `spaceChargeSolver_m` → algorithm → `CartesianPIC3DFieldStorage`. Need public accessors on `SpaceChargeSolver` and `CartesianPIC3DAlgorithm`. |
| `ParallelTracker.cpp` (external fields) | Remove `E_external`/`B_external` from mesh block. Call `elementVizRegistry_m.registerElements(...)` before tracking loop. Call `elementVizRegistry_m.sampleDynamic(...)` inside loop before `cat_vis.Execute()`. |
| `SpaceCharge/SpaceChargeSolver.h` | Add `CartesianPIC3DFieldStorage& fieldStorage()` public accessor (or `vizFields()` callback). |
| `SpaceCharge/CartesianPIC3D/CartesianPIC3DAlgorithm.h` | Add `CartesianPIC3DFieldStorage& fieldStorage()` public accessor. |
| **NEW** `catalyst/ElementVizRegistry.h` | Registry class declaration. |
| **NEW** `catalyst/ElementVizRegistry.cpp` | Registry class implementation — field ownership, registration, sampling. |
| `catalyst/catalyst-muon-combined.py` | `registrationName='particle_bunch_0'`, selectors `//particles`/`//boundingbox`/`//transform`/`//mesh`. Drop `E_external`/`B_external`. Add per-element external field bindings. |
| `catalyst/catalyst-muon-test-1.py` | Same renaming. |
| `catalyst/catalyst-muon-fields.py` | Update or merge into combined script. |
| `catalyst/catalyst-editing.md` | Update channel names, block names, snapshot instructions. |

---

## 8. Open questions

1. **Access to `CartesianPIC3DFieldStorage`:** The solver → algorithm →
   field storage chain is private. We need public accessors, or
   alternatively a `vizFields()` callback interface that the Catalyst
   adaptor calls. The latter avoids exposing internal solver types.

2. **External fields on mesh block vs. separate arrays:** External
   element fields are sampled onto the bunch mesh and added as arrays
   on the `data/mesh/` block. Alternatively, a separate `data/external/`
   block could be used for cleaner ParaView pipeline separation. Decide
   during implementation.

3. **RF fieldmap amplitude (Phase 2):** If we export the grid as `Static`
   and the amplitude as a separate small channel, the Python script
   needs to combine them. Alternatively, the Catalyst script can read
   the frequency/phase from the input file and compute the scale itself.

4. **Multi-rank:** Field maps are replicated on all ranks (each element
   owns its own copy). The Catalyst adaptor only exports from rank 0, so
   this is not an issue. The self-field PIC mesh is decomposed and
   already handled by the existing `MergeBlocks` logic in the Python
   script. External fields on the bunch mesh are also decomposed and
   follow the same `MergeBlocks` path. **Multi-rank testing is required**
   to verify the merged multimesh structure (particles + mesh in one
   channel) works correctly across N ranks.

5. **Sampling cost:** For ExternalFieldsOnBunchMesh, the per-step cost
   is `N_active_elements × N_local_grid_points` calls to
   `getFieldstrength()`. For a typical bunch mesh of ~100³ local points
   and ~10 active elements, this is ~10⁷ evaluations per step —
   comparable to the existing particle push. Should be acceptable for
   visualization but worth measuring.

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
