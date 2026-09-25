# Editing the Catalyst visualization scripts

Tips for iterating on the in-situ pipelines (e.g. `catalyst/catalyst-muon-test-1.py`)
without rebuilding the pipeline by hand in ParaView every time.

## 1. Record a trace of your tweaks, then transplant the `Set()` blocks

The fastest loop for changing viz parameters (scalar bars, LUTs, camera, display
properties):

1. In ParaView GUI, load the same data your script shows (for geometry: the
   `.vtm` from `tools/element_positions/element_positions.py`; for particles: a
   `.vtp`/`.vti` dumped by `IPPL_CATALYST_VTK=ON`).
2. Build the same pipeline once (readers + programmable filter + Show). Save it
   as a normal ParaView state (`.pvsm`) so you can reload this scratch session
   anytime - this is your stand-in for the catalyst pipeline.
3. `Tools > Start Trace`, make your tweaks (move/resize scalar bars, change
   colors, camera), `Tools > Stop Trace`.
4. Copy the generated `GetScalarBar(...).Set(...)`,
   `GetColorTransferFunction(...).Set(...)`, `renderView1.Set(CameraPosition=...)`
   blocks into the catalyst script verbatim.

The trace's `FindSource(...)` / `SetActiveSource(...)` lines only identify which
proxy is being modified - drop them and keep the `Set()` calls, making sure the
variable names match the ones in the catalyst script. The addendum at the end of
a trace (layout size, camera placement) is also copyable - the camera block is
often all you want.

This is exactly how the `etype` categorical LUT and lower-left horizontal scalar
bar in `catalyst-muon-test-1.py` were produced: tweaked in the GUI, recorded,
pasted verbatim.

## 2. Bridging file-loaded and in-situ pipelines

The core annoyance: in the GUI your channels come from file readers, in a
catalyst run the same data arrives in-situ as `PVTrivialProducer` channels, so
naively the script needs hand-editing every time you move between the two. Three
techniques remove that editing entirely.

### 2a. Export from the GUI: File > Save Catalyst State

Build the state so the export works as-is (zero hand-edits - verified with
this ParaView build):

1. Load placeholder file data for every catalyst channel, naming each
   stand-in reader exactly like the channel (e.g. `MUON` for the particle
   bunch, `bunch_mesh` for the field grid) - the runtime substitution is
   automatic (2b). The `.vtm` geometry reader keeps its normal name.
2. Build the pipeline and looks you want. Two don'ts:
   - Never touch the stand-in reader's block-selection / array checkboxes:
     defaults load everything, and only *modified* properties are traced.
     Touched ones would be applied to the substituted channel producer and
     only produce `'...' is not available on Conduit` warnings.
   - Don't set the Point Gaussian Intensity/Radius/Type-scalars properties
     (GUI-only; they warn in batch). Use a fixed Gaussian Radius instead.

   `ExtractBlock` selectors can stay exactly as the GUI records them
   (`/Root/block_main` with Assembly 'Hierarchy') - verified to work in-situ
   too; `//block_name` also works everywhere if you prefer.
3. Add a PNG extractor: make the render view active, then the **Extractors**
   menu (top-level, next to Filters) -> PNG. In its properties set the
   filename template (absolute path with a `{timestep:06d}` token; there is
   also `{camera}` for camera orbits), ImageResolution, and the Trigger
   (TimeStep, frequency 1). The exported block looks like:

       pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
       pNG1.Writer.FileName = '/path/to/catalyst_viz_{timestep:06d}.png'
       pNG1.Writer.ImageResolution = [1024, 768]
       pNG1.Writer.Format = 'PNG'
       pNG1.Trigger.Frequency = 1

4. `File > Save Catalyst State`. The dialog only shows Catalyst options
   (extracts directory, frequency...); this ParaView version has no
   per-source channel marking. The exported script renders PNGs on every
   catalyst timestep with no `catalyst_execute` at all - extractors are the
   export's render mechanism (a state without extractors exports to a script
   that does nothing per timestep).

What you give up without a `catalyst_execute` (add one only if needed):

- Per-step `bLUT.RescaleTransferFunctionToDataRange()` - without it the LUT
  keeps the range baked at export time (from the snapshot). Fine while the
  field range barely moves during the run.
- The custom time-stamp text (the timestep is in the PNG filename anyway;
  try the Annotate Time filter, or keep a minimal execute for the label).
- The `IPPL_CATALYST_PNG` switch - extractor output is governed by its Trigger
  instead. The snapshot dump (2c) is controlled by `OPALX_CATALYST_DUMP_SNAPSHOT`.

PNG cadence is the catalyst execute cadence either way (every 5 steps with
the current ippl settings).

### 2b. Dual-mode input: one script for both worlds (automatic)

No conditional adapter is needed. Under Catalyst, `paraview.simple` replaces
**any** source created with a channel's `registrationName` by the live
channel producer (a `Conduit` source, because ippl publishes conduit
blueprint channels). So the same line

    muon_particles = XMLMultiBlockDataReader(
        registrationName='MUON',
        FileName=['/home/biddisco/temp/muon/MUON_snapshot/MUON_snapshot.vtm'])

reads the snapshot file when run in the GUI/pvpython (the offline stand-in)
and silently binds the live channel during a catalyst run. The substituted
producer ignores unknown properties (`IgnoreUnknownSetRequests`), which is
why reader-only `Set()` calls end up as warnings - avoid them by not
touching the reader widgets (2a don'ts). Everything
downstream (`ExtractBlock` selectors like `//block_main` and
`//block_transform`, the programmable filter, displays, LUTs, bars) is shared
verbatim, because the selectors only match block *names*.

The same mechanism works for field channels: a reader with
`registrationName='bunch_mesh'` is replaced by the live mesh channel.
The combined script (`catalyst-muon-combined.py`) uses this to bind both
`MUON` (particles) and `bunch_mesh` (fields) in the same state.

So: build and iterate in the GUI on the snapshot + `.vtm`, then run the
identical script in-situ. Note the geometry side needs no stand-in at all:
the `.vtm` is a file in both contexts, which is why
`XMLMultiBlockDataReader` already works in-situ.

### 2c. Capturing a faithful in-situ snapshot as your offline stand-in

The adapter above is only as good as the stand-in file. Do not try to
reconstruct the channel data by hand - capture it from a real catalyst run.
The state scripts register data extractors controlled by
`OPALX_CATALYST_DUMP_SNAPSHOT`:

- Set `OPALX_CATALYST_DUMP_SNAPSHOT=N` to write the stand-in on iteration `N`
  (and every `N` iterations after that). Set it to `0` or unset it to disable.
- The particle extractor (in `catalyst-muon-test-1.py` and
  `catalyst-muon-combined.py`) writes
  `MUON_snapshot/{timestep:06d}.vtpc` (a `vtkPartitionedDataSetCollection`
  with blocks `block_main`, `block_help`, `block_transform`).
- The field extractor (in `catalyst-muon-combined.py`, requires
  `OPALX_CATALYST_FIELDS=ON`) writes
  `bunch_mesh_snapshot/{timestep:06d}.vtpd` (a `vtkPartitionedDataSet`
  containing the uniform grid with all field arrays as cell data).

The scripts execute the pipeline every timestep so that `N` is interpreted as
the simulation iteration rather than the pipeline invocation count. Move or
rename the produced files as needed for your offline stand-in paths.

**Generating snapshots for offline editing of the combined script:**

Run a short simulation with both channels and snapshot enabled. Use
`pipeline_default.py` (the generic IPPL script) so the state scripts don't
conflict with the snapshot extractors:

    cd /home/biddisco/temp/muon

    OPALX_CATALYST_FIELDS=ON \
    OPALX_CATALYST_DUMP_SNAPSHOT=64 \
    CATALYST_PIPELINE_PATH=/home/biddisco/build/opalx-catalyst/_deps/ippl-build/share/ippl/catalyst_scripts/pipeline_default.py \
    IPPL_CATALYST_VIS=ON IPPL_CATALYST_LIVE=OFF IPPL_CATALYST_STEER=OFF \
    IPPL_CATALYST_PNG=OFF IPPL_CATALYST_VTK=OFF IPPL_CATALYST_VERBOSITY=5 \
    IPPL_CATALYST_PROXY_OPTION=OFF \
    CATALYST_IMPLEMENTATION_PATHS=/home/biddisco/build/paraview/lib/catalyst \
    CATALYST_IMPLEMENTATION_NAME=paraview \
    mpirun -n 1 /home/biddisco/build/opalx-catalyst/bin/opalx ./mue4_analytical.in --info 2

This writes both `MUON_snapshot/MUON_snapshot_000064.vtpc` and
`bunch_mesh_snapshot/bunch_mesh_snapshot_000064.vtpd` (plus their
subdirectories). For a particle-only snapshot, omit `OPALX_CATALYST_FIELDS=ON`.

**Creating a VTM stand-in from a VTPC snapshot:**

`catalyst-muon-test-1.py` and `catalyst-muon-combined.py` use
`XMLMultiBlockDataReader` (`.vtm`) as the offline stand-in because the PVSM
was originally authored with that reader type. To create a `.vtm` wrapper
around a `.vtpc` snapshot, copy the VTPC file and its subdirectory to a
stable location, then write a small `.vtm` that references the block files
by name. The block names (`block_main`, `block_help`, `block_transform`)
match exactly between the old `vtkMultiBlockDataSet` and the new
`vtkPartitionedDataSetCollection` formats, so all `ExtractBlock` selectors
(`//block_transform`, `/Root/block_main`, `//block_main`) work unchanged.

**Note on VTPD field snapshots:**

The `bunch_mesh` channel arrives as `vtkImageData` (a uniform grid). The
VTPD snapshot is a `vtkPartitionedDataSet` containing one image data
partition per MPI rank. The combined script uses `MergeBlocks` to flatten
these into a single dataset before slicing.

## 3. Catalyst Live: tweak while the simulation runs

`catalyst-muon-test-1.py` already sets `options.EnableCatalystLive = 1` with
`CatalystLiveTrigger = 'TimeStep'`. Run opalx with `IPPL_CATALYST_LIVE=ON`, then
in ParaView GUI: `Catalyst > Connect` (default port 22222). You can then change
display properties, scalar bars and camera live while the run streams timesteps.
Record a trace of the changes you like and transplant them (workflow 1).

## 4. Headless testing with pvpython (no GUI, no simulation)

Use the ParaView build's own python so imports and VTK features match what
catalyst will use:

    /home/biddisco/build/paraview/bin/pvpython my_test.py

(or `pvbatch` - it is what catalyst itself launches headlessly.)

A test script can rebuild the interesting part of the pipeline in ~20 lines:
create a `RenderView`, load the `.vtm` (and a small `.vtp` carrying
`BunchOrigin`/`BunchInvOrigin`... arrays as a fake transform block), create the
`ProgrammableFilter` with the same `Script` string the state file uses, `Show`
it, apply the LUT/scalar-bar settings, then `Render()` + `SaveScreenshot()` to a
PNG. Fetch data with
`proxy.GetClientSideObject().GetOutputDataObject(0)` to assert things
programmatically (block counts, transforms applied, arrays present).

Two gotchas learned the hard way:

- `GetFieldData().GetArray(name)` only sees numeric arrays; string arrays need
  `GetAbstractArray(name)`.
- Readers expose `FieldArrayStatus` separately from `CellArrayStatus`; if a
  field-data array does not show up, enable it explicitly on the reader proxy.

## 5. Anatomy of the state scripts - where each knob lives

Three state scripts are provided in `catalyst/`:

| Script | Channels | Description |
|--------|----------|-------------|
| `catalyst-muon-test-1.py` | `MUON` (particles only) | Geometry + particles, the original workflow |
| `catalyst-muon-fields.py` | `bunch_mesh` (fields only) | E/B field slices through the mesh |
| `catalyst-muon-combined.py` | `MUON` + `bunch_mesh` | Both particles and fields in one scene |

All three follow the same structure:

- Camera / view: `renderView1.Set(CameraPosition=..., CameraFocalPoint=...,
  CameraViewUp=..., CameraParallelScale=...)` near the top. `CameraParallelScale`
  is the zoom (smaller = closer). The `_DisableFirstRenderCameraReset()` call
  keeps ParaView from reframing on the first Show.
- Particle source: `XMLMultiBlockDataReader(registrationName='MUON')` is the
  catalyst channel input; `ExtractBlock` selectors pull `//block_main`
  (bunch) and `//block_transform` (pose) out of it.
- Field source (combined/fields only): `XMLImageDataReader(registrationName='bunch_mesh')`
  reads the uniform grid; `MergeBlocks` flattens per-rank partitions; `Slice`
  creates 2D cuts through the mesh for visualization.
- Geometry: the `XMLMultiBlockDataReader` on the `.vtm` (not a Catalyst
  channel; reads from file in both contexts).
- Filter wiring: `ProgrammableFilter(Input=[geometry, transform_extract])`.
- Display properties: the `Show(...)` + `Display.Set(...)` blocks
  (Representation, ColorArrayName, Opacity, GaussianRadius...).
- Color: `GetColorTransferFunction(name)` (LUT) and
  `GetOpacityTransferFunction(name)` (opacity curve). The LUT proxy is keyed by
  array name, so GUI and script settings converge on the same object.
- Scalar bars: `GetScalarBar(LUT, renderView1)`.
- Per-timestep behavior: `catalyst_execute(info)` at the bottom - rescaling LUTs
  (`bLUT.RescaleTransferFunctionToDataRange()`), updating the time annotation,
  and `SaveScreenshot(...)` per `info.cycle`.

## 6. Keep filter bodies in separate files (the loader pattern)

Programmable filter scripts embedded as `Script="""..."""` strings are painful to
edit (escaped quotes, no syntax highlighting, no reuse). Instead keep the body in
a real file and load it:

    INVERSE = True   # parameter passed into the filter body
    import os
    _filter_script = os.environ.get(
        'BUNCH_TRANSFORM_FILTER',
        '/home/biddisco/src/opalx/insitu/catalyst/bunch_transform_filter.py')
    _filter_ns = {'self': self, 'INVERSE': INVERSE}
    with open(_filter_script) as _fh:
        exec(compile(_fh.read(), _filter_script, 'exec'), _filter_ns)

Now you edit `catalyst/bunch_transform_filter.py` with a normal editor and every
state script that loads it gets the change - no pipeline rebuilding, no
copy-paste drift. The body falls back to the `BUNCH_TRANSFORM_INVERSE`
environment variable (then a default) when `INVERSE` is not passed in, so it
still works when pasted directly into a GUI Programmable Filter.

## 7. Scalar bar cookbook

    bar = GetScalarBar(lut, renderView1)
    bar.Set(
        AutoOrient=0,                     # keep manual Orientation
        Orientation='Horizontal',         # or 'Vertical'
        WindowLocation='Lower Left Corner',
        # or free placement (normalized viewport coords):
        # WindowLocation='Any Location', Position=[0.50, 0.055], Position2=[0.44, 0.07],
        Title='element_type',
        ComponentTitle='',
        AllowOverlappingLabels=1,
        FontSize=12,
    )
    bar.Visibility = 1

- Corner locations do not stack; to place a bar above another annotation use
  `Any Location` with explicit `Position` (bottom-left of the bar) and
  `Position2` (its width/height), both in fractions of the viewport.
- `ScalarBarLength` shrinks the bar along its axis (try/except it - older builds
  only have `Position2`).

Categorical (discrete) coloring, as used for `etype`:

    lut = GetColorTransferFunction('element_type')
    lut.Set(
        InterpretValuesAsCategories=1,
        AnnotationsInitialized=1,
        Annotations=['DRIFT', 'DRIFT', 'MONITOR', 'MONITOR', ...],  # value, label pairs
        IndexedColors=[r, g, b, r, g, b, ...],                       # one color per category
        IndexedOpacities=[1.0, ...],
        EnableOpacityMapping=1,
    )

Note: string arrays live in per-block *field data*; color by them with
`ColorArrayName=['FIELD', 'name']` (`['POINTS', ...]` / `['CELLS', ...]` for the
other associations). String cell data does not survive the `.vtm` XML round trip
in every VTK version - field data does.

## 8. Environment-variable switches

Parameters the scripts read at load time, so you can change behavior without
editing:

- `CATALYST_PIPELINE_PATH` - the state script catalyst loads. Choose
  `catalyst-muon-test-1.py` (particles only), `catalyst-muon-fields.py`
  (fields only), or `catalyst-muon-combined.py` (both).
- `OPALX_CATALYST_FIELDS` - set to `ON`/`1`/`TRUE`/`YES` to have OPALX
  register the `bunch_mesh` field channel in addition to particle channels.
  Required for `catalyst-muon-fields.py` and `catalyst-muon-combined.py`.
- `OPALX_CATALYST_DUMP_SNAPSHOT` - set to `N` to write channel snapshots
  on iteration `N` (and every `N` after). See section 2c for details.
- `IPPL_CATALYST_PNG` - ippl's C++ side only accepts exactly `ON`; the state
  script's own guard accepts `1/ON/TRUE/YES`. Watch the spelling.
- `IPPL_CATALYST_VTK=ON` - dump channel data to files; great for building
  placeholder pipelines (workflow 1/2) and post-mortem debugging.
- `BUNCH_TRANSFORM_FILTER` / `BUNCH_TRANSFORM_INVERSE` - override the filter
  body location / transform direction.

Run template (particles only):

    export CATALYST_PIPELINE_PATH=/home/biddisco/src/opalx/insitu/catalyst/catalyst-muon-test-1.py
    IPPL_CATALYST_VIS=ON IPPL_CATALYST_LIVE=OFF IPPL_CATALYST_STEER=OFF \
    IPPL_CATALYST_PNG=1 IPPL_CATALYST_VTK=OFF IPPL_CATALYST_VERBOSITY=0 \
    IPPL_CATALYST_PROXY_OPTION=OFF \
    CATALYST_IMPLEMENTATION_PATHS=/home/biddisco/build/paraview/lib/catalyst \
    CATALYST_IMPLEMENTATION_NAME=paraview \
    /home/biddisco/build/opalx-catalyst/bin/opalx ./mue4_analytical.in --info 2

Run template (fields + particles combined):

    export CATALYST_PIPELINE_PATH=/home/biddisco/src/opalx/insitu/catalyst/catalyst-muon-combined.py
    OPALX_CATALYST_FIELDS=ON \
    IPPL_CATALYST_VIS=ON IPPL_CATALYST_LIVE=OFF IPPL_CATALYST_STEER=OFF \
    IPPL_CATALYST_PNG=OFF IPPL_CATALYST_VTK=OFF IPPL_CATALYST_VERBOSITY=5 \
    IPPL_CATALYST_PROXY_OPTION=OFF \
    CATALYST_IMPLEMENTATION_PATHS=/home/biddisco/build/paraview/lib/catalyst \
    CATALYST_IMPLEMENTATION_NAME=paraview \
    /home/biddisco/build/opalx-catalyst/bin/opalx ./mue4_analytical.in --info 2

## 9. Trace gotchas

- Traces open with `# script-version: 2.0` and `paraview.compatibility`
  pinning; keep those lines in state scripts but drop them when copying
  fragments.
- `GenerateRGBPoints(...)` is available in state scripts and traces; it builds a
  default ramp without spelling out RGBPoints tuples.
- Scripts embedded inside `Script="""..."""` properties have their quotes
  escaped (`\"`); when moving code between an embedded string and a real file,
  re-escape/un-escape accordingly (another reason to prefer the loader pattern).
- A trace never shows properties you did not touch - if a `Set()` does not
  appear, the value was already current when recording started.

## 10. Field channels

When `OPALX_CATALYST_FIELDS=ON` is set, OPALX exports a single mesh channel
named `bunch_mesh` (type `"mesh"`) in addition to the particle channel(s).
All field arrays live as cell data on this one uniform grid:

| Array        | Type   | Contents                                            |
|--------------|--------|-----------------------------------------------------|
| `E_self`     | vector | Space-charge electric field (`FieldContainer::E_m`) |
| `B_self`     | vector | Space-charge magnetic field (`Btmp_m`)              |
| `rho`        | scalar | Charge density                                      |
| `phi`        | scalar | Electrostatic potential (CG solver only)            |
| `E_external` | vector | Beamline element electric field sampled on the mesh |
| `B_external` | vector | Beamline element magnetic field sampled on the mesh |

Two state scripts consume this channel:

- `catalyst/catalyst-muon-fields.py` - fields only (E/B slices).
- `catalyst/catalyst-muon-combined.py` - fields + particles in one scene.

Use either with `IPPL_CATALYST_PNG=OFF`; otherwise IPPL's built-in per-channel
PNG extractors grab the same channels and the state script cannot bind them.

Run template for the combined pipeline:

    export CATALYST_PIPELINE_PATH=/home/biddisco/src/opalx/insitu/catalyst/catalyst-muon-combined.py
    OPALX_CATALYST_FIELDS=ON \
    IPPL_CATALYST_VIS=ON IPPL_CATALYST_LIVE=OFF IPPL_CATALYST_STEER=OFF \
    IPPL_CATALYST_PNG=OFF IPPL_CATALYST_VTK=OFF IPPL_CATALYST_VERBOSITY=5 \
    IPPL_CATALYST_PROXY_OPTION=OFF \
    CATALYST_IMPLEMENTATION_PATHS=/home/biddisco/build/paraview/lib/catalyst \
    CATALYST_IMPLEMENTATION_NAME=paraview \
    /home/biddisco/build/opalx-catalyst/bin/opalx ./mue4_analytical.in --info 2

The field channel arrives as `vtkImageData`. For offline development, dump a
stand-in snapshot from any run that registers the channels:

    OPALX_CATALYST_DUMP_SNAPSHOT=64 ...

This writes `bunch_mesh_snapshot/{timestep:06d}.vtpd` once the run reaches
iteration 64. The `.vtpd` is a `vtkPartitionedDataSet` containing one image
data partition per MPI rank; `MergeBlocks` in the combined script flattens
these into a single dataset for slicing.
