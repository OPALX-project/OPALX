# OPAL-X element geometry generator

This tool converts the element-description JSON written by OPAL-X into a
ParaView-readable `.vtm` multiblock file. Each accelerator element becomes a
separate mesh (drift tube, bend chamber, solenoid envelope, monitor disc, etc.)
that can be overlaid with particle traces or field data.

## Requirements

- Python 3
- `pyvista`
- `numpy`

Install with:

```bash
pip install pyvista numpy
```

## Quick start

OPAL-X writes `<base>_ElementPositions.json` when run with the appropriate
options. Convert it to a `.vtm` file with:

```bash
python element_positions.py data/mue4_analytical_ElementPositions.json
```

This writes `data/mue4_analytical_ElementPositions.vtm`. Open the `.vtm` in
ParaView; each element category appears as a separate block.

To preview interactively:

```bash
python element_positions.py data/mue4_analytical_ElementPositions.json --show
```

To write a screenshot without opening a window:

```bash
python element_positions.py data/mue4_analytical_ElementPositions.json \
    --screenshot geometry.png
```

## Command-line options

| Option | Description |
|--------|-------------|
| `json` | Input `<base>_ElementPositions.json` from OPAL-X. |
| `--out <path>` | Output `.vtm` path (default: `<json stem>.vtm`). |
| `--clamp-radius <m>` | Radius used for elements with unset/infinite apertures. Defaults to `reference_transverse` from the JSON. |
| `--clamp-threshold <m>` | Semi-axes larger than this are treated as unset (default: 10). |
| `--no-clamp` | Draw unset apertures at their raw size. |
| `--simple-geometry` | Use schematic geometry: circular bend tubes and solid solenoids. |
| `--capped` | Close the ends of drift tubes, bends, and generic elements. Default is open tubes so the beam vacuum stays visible. |
| `--hide <category>` | Skip a category (repeatable). Categories: `beamline`, `bending`, `focusing`, `rf`, `diagnostic`. |
| `--show` | Open an interactive PyVista window. |
| `--screenshot <path>` | Render offscreen and save an image. |

## Architecture

The converter is built around small `OpalX...Source` classes:

```text
OpalXElementSource (abstract)
├── OpalXBendSource          # exact aperture swept along design_path
├── OpalXSolenoidSource      # hollow tube (default) or solid magnet
├── OpalXMultipoleSource     # four poles for order-2, else cylinder
├── OpalXDriftSource         # cylinder or rectangular box
├── OpalXCavitySource        # stacked cell profile
├── OpalXMonitorSource       # thin annular disc
└── OpalXGenericSource       # fallback cylinder
```

Each source receives:

- `element`: the JSON description of one beamline element.
- `args`: the parsed command-line namespace.
- `reference`: the `(semi_x, semi_y)` reference transverse extent from the JSON.

It returns a `pyvista.PolyData` mesh in **element-local coordinates**. The caller
then transforms it to global coordinates with `place()` and attaches category
cell data.

The registry at the bottom of the file maps OPAL-X element type strings to
source classes:

```python
SOURCE_REGISTRY = {
    'SOLENOID': OpalXSolenoidSource,
    'SBEND': OpalXBendSource,
    ...
}
```

## Adding a new element source

1. Subclass `OpalXElementSource`.
2. Implement `execute()` to return a `pv.PolyData` in local coordinates, or
   `None` if the element cannot be represented.
3. Add the mapping to `SOURCE_REGISTRY`.

Example:

```python
class OpalXMyElementSource(OpalXElementSource):
    def execute(self):
        (semi_x, semi_y), clamped = self._effective_transverse()
        length = self.element.get('length', 0.0)
        if length <= 0.0:
            return None
        mesh = cylinder(length, semi_x, semi_y, capped=self.args.capped)
        mesh = place(mesh, self.element)
        return self._add_metadata(mesh, clamped)
```

## Geometry primitives

Wherever possible the sources reuse VTK/PyVista primitives:

- `pv.Cylinder` for straight drift/solenoid/multipole/generic tubes.
- `pv.Disc` extruded for monitors.
- `pv.Box` for rectangular apertures and quadrupole poles.
- `pv.lines_from_points(...).tube(...)` for simple bend tubes.
- `sweep_profile()` only for bends (curved path) and the hollow solenoid bore
  (needs both inner and outer surfaces).

## Output format

The output is a VTK multiblock dataset (`.vtm`). Each block corresponds to an
element category, and each block contains one named `PolyData` per element.
The following field data are attached to every element mesh:

- `element`: element name.
- `etype`: OPAL-X element type string.
- `clamped`: whether the aperture was clamped because it was unset/infinite.

Cell data `etype` contains a category ID so elements can be colored by category
in ParaView.
