# Catalyst Visualization Integration Plan

This plan tracks the work to evolve OPALX's Catalyst in-situ visualization from a
hardcoded, per-step-allocating prototype into a reusable, memory-lean subsystem
built on a mesh-centric export API in IPPL.

Companion document: `~/src/ippl/catalyst/catalyst-viz-performance-todo.md`
(untracked) contains the original per-item performance assessment of the
`catalyst-viz-v2.0.0` code. Its P0/P1 items are folded into Phase 2 below.

---

## Phase 1 — Port and commit current work

Status: **completed**

- [x] Fix Catalyst discovery: remove `CATALYST_HINT_PATH` from
  `cmake/OPALXOptions.cmake`, use standard `Catalyst_DIR` directly.
- [x] Update `README.md` to document `Catalyst_DIR` instead of
  `CATALYST_HINT_PATH`.
- [x] Build and smoke-test OPALX with the `Catalyst_DIR` change.
- [x] Commit IPPL `catalyst-viz-dev` fixes:
  - `VIZ: style: apply clang-format to Catalyst files`
  - `VIZ: fix(conduit): map bool particle attributes to conduit uint8`
  - `VIZ: fix(catalyst): make InitVizChannel/ExecVizChannel public`
  - `VIZ: fix(catalyst): restore requires(!is_vector_v<T>) on scalar field InitVizChannel`
  - `VIZ: fix(catalyst): copy field device view to host before extracting non-ghost interior`
  - `VIZ: feat(catalyst): add TransformData bunch-transform channel`
- [x] Commit OPALX work:
  - Field export via `OPALX_CATALYST_FIELDS`
  - External field sampling onto the bunch mesh
  - Numeric snapshot dumping via `OPALX_CATALYST_DUMP_SNAPSHOT`
  - Bunch-transform channel support
  - CMake Catalyst forwarding
  - New scripts: `catalyst-muon-fields.py`, `catalyst-editing.md`

---

## Findings motivating the Phase 2 redesign

Verified in the current code (file:line as of the `catalyst-viz-dev` /
`catalyst-viz-v2-jb` heads):

### Memory / per-step allocation (the dominant cost)

1. **Cache destroyed every step.**
   `ippl/src/Stream/InSitu/CatalystAdaptor.hpp:1254-1256` clears
   `viewRegistry_m`, `ghostMaskCache_m`, and resets `node_m` at the end of every
   `Execute()`. All host mirrors, ghost masks, RankID/iota arrays, and the whole
   Conduit node tree are freed and rebuilt each step. Root cause of items 2-8.
2. **Inverted refresh semantics that work by accident.** The early return
   `refreshOnly && !forceHostCopy_m[label]` (`CatalystAdaptor.hpp:181-183`,
   `:603-605`) would serve stale data if the caches were *not* cleared;
   the per-step clear is what makes data flow. `rememberNow()`
   (`:1106-1134`) is a flag save/restore hack to refresh fields without
   double-copying within one step.
3. **Double copy + double buffer per field.** Non-ghost path
   (`CatalystAdaptor.hpp:271-294`): full device view (ghosts included) mirrored
   to host, then interior deep-copied into a second fresh buffer. ~2x data
   movement and 2x peak allocation per field per step.
4. **RankID duplicated per field.** `rank_id_view_cells`
   (`CatalystAdaptor.hpp:330-356`) is a full `int` array per field per step;
   identical content for all fields sharing a mesh. Should be one array per
   mesh, allocated once.
5. **Ghost masks re-derived every step.** With `IPPL_CATALYST_GHOST_MASKS=ON`,
   the per-step cache clear forces a fresh device `Field<unsigned char>`
   allocation, kernel launch, and D2H copy per field (`:461-512`). Masks depend
   only on (mesh, layout, nghost) — static until repartition.
6. **Particle attribute mirrors allocated per step.**
   `ippl/src/Particle/ParticleAttrib.hpp:602-608`:
   `signConduitBlueprintNode` does `create_mirror_view_and_copy` per attribute
   per call; fresh `iota` (int64) and `rank_id` per step
   (`CatalystAdaptor.hpp:649-650`).
7. **OPALX `sampleExternalFieldsToMesh` compounds it.**
   `src/Algorithms/ParallelTracker.cpp:1123-1128`: two full-mesh host mirrors
   allocated per step, a serial host triple loop over cells x elements, then
   copy back to device — after which the Catalyst layer copies device->host
   again. Same bytes cross the bus twice per step.
8. **Conduit node tree rebuilt per step.** Channels use `set_external`
   (zero-copy), so rebuilding the tree each step is pure allocation overhead.

### Abstraction (wrong seam)

- **Channel = registered object.** One Conduit channel per field duplicates the
  mesh description (coordset/topology), RankID, and ghost arrays per field, and
  shows N copies of "the mesh" in ParaView's pipeline browser.
- **`VisRegistryRuntime` is a flat labeled list.** It cannot express "fields of
  mesh M", "static array", or refresh policy; semantics are smuggled through
  string-keyed side tables (`forceHostCopy_m`, `transformChannels_m`).
- **`ViewRegistry` conflates three jobs** — staging-buffer ownership, refresh
  cache, and mirror lookup (with unnamed entries that exist only to keep memory
  alive, `ParticleAttrib.hpp:673`).
- **The bunch transform is embedded as a block inside the particle multimesh**
  (`CatalystAdaptor.hpp:873-908`) — a workaround for the channel-per-object
  model. Bunch-frame metadata is not particle data.
- **OPALX side**: ~170 lines of `#ifdef IPPL_ENABLE_CATALYST` inside
  `ParallelTracker::execute()` with env vars parsed in static lambdas;
  registration policy, refresh policy, and per-step plumbing interleaved with
  tracking logic.

### Multi-bunch requirement

`FodoCell-multibeam-fromfile.in` (regression tests) tracks two beams
simultaneously (`BEAM1` charge -1, `BEAM2` charge +1). `PicManager`
(`ippl/src/Manager/PicManager.h:72-120`) holds them as
`pcontainers_m[0]` and `[1]`, each with its own particle set and its own
`getToLabTrafo()`. The current Catalyst code in `ParallelTracker.cpp:466`
only registers container 0.

The mesh-centric API must support multiple independent particle channels,
each with its own transform, sharing one mesh channel for fields:

```cpp
auto meshCh = cat_vis.addMeshChannel("bunch", fc->getMesh(), fc->getFL());
meshCh->addArray("rho",    fc->getRho());
meshCh->addArray("E_self", fc->getE());

for (size_t i = 0; i < bunch.getNumParticleContainers(); ++i) {
    auto pc   = bunch.getParticleContainer(i);
    auto name = bunch.getParticleName(i);  // "BEAM1", "BEAM2"
    auto ch   = cat_vis.addParticleChannel(name, pc);
    ch->setTransform(rowsFromMatrix(pc->getToLabTrafo()...));
}
```

### Why Phase 2 and old Phase 3 are merged

Grouping channels by mesh and redesigning the registration API touch the same
lines (channel creation, visitor path, script channel naming). Old Phase 2 as
written proposed *auto-detecting* mesh identity inside `ExecVizChannel` via
pointer comparison — the wrong seam: the mesh relationship is known by the
caller at registration time (OPALX's `FieldContainer` hands over six fields on
one mesh). Inferring it at exec time is magic that fails edge cases (same mesh,
different layouts; same mesh deliberately split) and forces the adaptor to
accumulate grouping state mid-`Execute`. The merged phase makes grouping
explicit at registration (option B) and uses the new API in OPALX immediately.

---

## Phase 2 — Mesh-centric export API, persistent staging, OPALX adoption

Status: **pending** (merges old Phase 2 "group field channels by mesh" and old
Phase 3 "registry cleanup")

### Design

Core model: mesh-centric channels with explicit ownership and update policy.

```cpp
// IPPL: construction-time, declarative, grouped by mesh.
// One mesh channel per unique mesh; fields are arrays on the shared topology.
auto meshCh = cat_vis.addMeshChannel("bunch", fc->getMesh(), fc->getLayout());
meshCh->addArray("rho",    fc->getRho());                  // dynamic by default
meshCh->addArray("E_self", fc->getE());
meshCh->addArray("B_self", *fc->getTempBField());
meshCh->addArray("E_ext",  fc->getExternalE()).policy(UpdatePolicy::Static);

// One particle channel per particle container (supports multi-bunch).
for (size_t i = 0; i < bunch.getNumParticleContainers(); ++i) {
    auto pc   = bunch.getParticleContainer(i);
    auto name = bunch.getParticleName(i);                  // "BEAM1", "BEAM2", ...
    auto ch   = cat_vis.addParticleChannel(name, pc);
    ch->setTransform(rowsFromMatrix(pc->getToLabTrafo()...));
}
```

Design rules:

1. **One Conduit channel per mesh; arrays are fields on the shared topology.**
   Coordset, topology, RankID, ghost arrays: allocated once per mesh, shared by
   all arrays. Mesh identity is declared by the caller, never inferred.
2. **Persistent staging buffers.** Each exported array owns one host mirror,
   allocated at registration, reallocated only when the local domain changes.
   Per-step update = `deep_copy` into the existing buffer. Remove
   `viewRegistry_m.clear()` / `ghostMaskCache_m.clear()` / `node_m.reset()` from
   `Execute()`. Build the Conduit channel tree once; per step only update
   `state/cycle` and call `catalyst_execute`.
3. **`UpdatePolicy`: `Always | Static | Manual`.** `Always` (default) copies
   every `Execute` — fixes the inverted refresh semantics. `Static` copies once
   at registration (external magnet fields, non-moving meshes, ghost masks,
   RankID). `Manual` copies only on `refresh(label)` — replaces `rememberNow`'s
   flag save/restore.
4. **Single-buffer field path.** Prefer publishing *with* ghosts +
   `vtkGhostType` (machinery already exists,
   `CatalystAdaptor.hpp:434-546`): one full mirror, one D2H copy, no interior
   extraction. If ghost publishing stays off by default, at least cache the
   full mirror and the interior buffer separately and persistently.
5. **Zero-copy into Conduit stays** (`set_external` with strides) — already
   correct.
6. **Channel metadata to scripts.** Generalize `--channel_names` to
   per-channel metadata (channel -> mesh, arrays, types) so
   `catalyst-muon-*.py` and the IPPL PNG extractors select arrays generically.
   This absorbs old Phase 3's "dynamic discovery" item, since grouping renames
   every channel anyway.

### IPPL checklist

- [ ] `addMeshChannel` / `addParticleChannel` / `addTransformChannel` API with
      per-array `UpdatePolicy`
- [ ] Persistent host staging; remove per-`Execute()` clears; build Conduit
      tree once (P0 items 1, 1b of the performance todo)
- [ ] Per-mesh shared RankID + ghost masks, allocated once (P0 items 2, 5)
- [ ] Single-buffer field publishing (ghost-mode preferred)
- [ ] Per-attribute persistent particle mirrors; drop per-step `iota`/`rank_id`
      allocation (P0 item 3)
- [ ] Channel metadata passed to scripts (generic array discovery)
- [ ] Compat shims: keep `rememberNow()` working for the PenningTrap demo;
      mark deprecated
- [ ] Update PNG extractors (`png_ext_sfield.py` / `png_ext_vfield.py`) and
      `pipeline_default.py` for the new channel layout
- [ ] Build and unit-test IPPL standalone (Landau damping / PenningTrap demos)

### OPALX checklist

- [ ] New `InSituManager` (or `CatalystSession`) class owning registration and
      the per-step `execute()`; `ParallelTracker` shrinks to constructing it
      and calling `insitu.execute(step, t)`
- [ ] Move env vars (`OPALX_CATALYST_FIELDS`, `OPALX_CATALYST_LIVE_FREQ`,
      `OPALX_CATALYST_DUMP_SNAPSHOT`) into the manager as validated options
- [ ] Mark verified-static arrays `Static`; keep time-dependent element fields
      (`RFCavity`, `TravelingWave`, `VariableRFCavity`) dynamic
- [ ] `sampleExternalFieldsToMesh`: reuse host buffers across steps (P1 item);
      longer term move the sampling loop to a Kokkos kernel or sample directly
      into the staging buffer to eliminate the host->device->host round trip
- [ ] Update `catalyst-muon-fields.py` / `catalyst-muon-test-1.py` to bind
      mesh channels and select arrays by name
- [ ] Update `catalyst/catalyst-editing.md` (grouped channels, new names,
      migration notes)
- [ ] Tests: particle pipeline, field pipeline, snapshot dumping, live
      steering; memory check (RSS stable across steps)
- [ ] Commit sequence: IPPL API first, then OPALX adoption

---

## Phase 3 — Cleanup

Status: **pending** (slimmed; only after Phase 2 lands)

- [ ] Remove compat shims (`rememberNow`) and dead refresh code
- [ ] Optional: runtime channel add/remove (register-before-Initialize is
      sufficient for OPALX and PenningTrap)
- [ ] Optional: device-side external-field sampling in OPALX
- [ ] Migrate demos (PenningTrap / FEL) if affected

---

## Risks

- **Channel renames break every existing script/extractor.** Requires a
  migration section in `catalyst-editing.md` and updated PNG extractors in
  IPPL. Snapshot dumps written by old scripts remain readable offline.
- **Stable-node assumption.** Confirm Catalyst v2 / ParaView does not require
  channels to be re-set on each `catalyst_execute` (standard flow re-reads
  external arrays; verify explicitly with live + steering enabled).
- **`Static` correctness.** Known static element types in OPALX: `Multipole`
  (incl. `Quadrupole`) and `Solenoid` (constant coefficients, ignore time).
  Known dynamic: `RFCavity`, `TravelingWave`, `VariableRFCavity`. Default
  consumers to `Always`; mark `Static` only where verified.
- **`results_m` growth.** Verified that the default IPPL pipeline overwrites
  fixed paths (no unbounded growth), but document the contract or reset
  defensively for custom pipelines (P1 item 7 of the performance todo).

---

## Branches

- OPALX: `catalyst-viz-v2-jb` (worktree `~/src/opalx/insitu`)
- IPPL: `catalyst-viz-dev` (`~/src/ippl/catalyst`)

## Build directory

`/home/biddisco/build/opalx-catalyst` (executable: `bin/opalx`)
