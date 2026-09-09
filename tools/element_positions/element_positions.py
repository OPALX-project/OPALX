#!/usr/bin/env python3
"""Build a multiblock VTK dataset from an opalx element-position description file.

opalx writes ``<base>_ElementPositions.json`` (see MeshGenerator::writeDescriptions)
next to its other outputs. That file contains one description per beamline element
(name, type, length, global pose, transverse extent, and for bends the local
design path) -- no geometry. This script turns those descriptions into a
``vtkMultiBlockDataSet`` written as ``.vtm`` (plus a directory of ``.vtp`` parts
that must travel with it):

    Beamline
    ├── Solenoids
    │   ├── SOL
    │   └── ...
    ├── Drifts
    │   ├── STOP_0250
    │   └── ...
    ├── Monitors
    ├── Quadrupoles
    └── Bends

Each element is an individually named block, so ParaView's pipeline browser can
show/hide (or color, clip, extract) any element or whole category independently.

Elements whose transverse extent is unset or absurd (the 1e6 m "infinite
aperture" default, e.g. for MONITORs) are clamped to a representative radius:
the user-supplied --clamp-radius if given, otherwise the reference_transverse
radius that opalx derived from the beamline (half the transverse support of the
first element that has one), otherwise 0.1 m. Clamped elements carry a
``clamped`` field-data flag.

Requires numpy and pyvista.

Examples
--------
    python element_positions.py data/mue4_analytical_ElementPositions.json
    python element_positions.py data/mue4_analytical_ElementPositions.json --show
    python element_positions.py ... --clamp-radius 0.15 --hide Monitors --hide Bends
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pyvista as pv

# Element-type visualization categories: (category, color)
CATEGORY_COLORS = {
    'Solenoids': (1.0, 0.541, 0.0),
    'Drifts': (0.0, 0.0, 1.0),
    'Monitors': (0.55, 0.55, 0.55),
    'Quadrupoles': (1.0, 0.0, 0.0),
    'Multipoles': (0.537, 0.745, 0.525),
    'Bends': (1.0, 0.847, 0.0),
    'Collimators': (0.5, 0.5, 0.0),
    'Cavities': (1.0, 1.0, 0.0),
    'Other': (0.5, 0.5, 0.5),
}


def category_of(element):
    etype = element['type']
    if etype == 'SOLENOID':
        return 'Solenoids'
    if etype == 'DRIFT':
        return 'Drifts'
    if etype == 'MONITOR':
        return 'Monitors'
    if etype == 'MULTIPOLE':
        if element.get('multipole_order') == 2:
            return 'Quadrupoles'
        return 'Multipoles'
    if etype in ('SBEND', 'RBEND', 'RBEND3D'):
        return 'Bends'
    if etype == 'COLLIMATOR':
        return 'Collimators'
    if etype in ('RFCAVITY', 'TRAVELINGWAVE', 'CONSTANTEFIELDCAVITY'):
        return 'Cavities'
    return 'Other'


def pose_matrix(element):
    """4x4 homogeneous matrix mapping element-local to global coordinates."""
    T = np.eye(4)
    T[:3, :3] = np.asarray(element['rotation'], dtype=float)
    T[:3, 3] = np.asarray(element['position'], dtype=float)
    return T


def place(mesh, element):
    mesh.transform(pose_matrix(element), inplace=True)
    return mesh


def cylinder(length, semi_x, semi_y, capped=False, z0=0.0, scale=1.0):
    """Cylinder along local z from z0 to z0+length, elliptical via semi_x/semi_y."""
    height = max(length, 1e-9)
    mesh = pv.Cylinder(
        center=(0.0, 0.0, z0 + 0.5 * height),
        direction=(0.0, 0.0, 1.0),
        height=height,
        radius=max(scale * semi_x, 1e-9),
        resolution=36,
        capping=capped,
    )
    if semi_x > 0.0 and semi_y != semi_x:
        mesh.scale([1.0, semi_y / semi_x, 1.0], inplace=True)
    return mesh


def box(length, semi_x, semi_y):
    return pv.Box(bounds=[-semi_x, semi_x, -semi_y, semi_y, 0.0, length])


def build_solenoid(length, semi_x, semi_y):
    body = cylinder(length, semi_x, semi_y)
    collar_length = min(0.18 * length, 0.8 * min(semi_x, semi_y))
    if collar_length <= 1e-12 or 2.0 * collar_length >= length:
        return body
    mesh = body
    for z0 in (0.0, length - collar_length):
        mesh = mesh + cylinder(collar_length, semi_x, semi_y, z0=z0, scale=1.18)
    return mesh


def build_quadrupole(length, semi_x, semi_y):
    poles = []
    hw, hh = 0.45 * semi_x, 0.75 * semi_y
    for sign in (+1, -1):
        xc = sign * (semi_x + hw)
        poles.append(pv.Box(bounds=[xc - hw, xc + hw, -hh, hh, 0.0, length]))
    hw, hh = 0.75 * semi_x, 0.45 * semi_y
    for sign in (+1, -1):
        yc = sign * (semi_y + hh)
        poles.append(pv.Box(bounds=[-hw, hw, yc - hh, yc + hh, 0.0, length]))
    mesh = poles[0]
    for pole in poles[1:]:
        mesh = mesh + pole
    return mesh


def build_cells(length, semi_x, semi_y, profile):
    mesh = None
    segment = length / len(profile)
    for i, scale in enumerate(profile):
        part = cylinder(segment, semi_x, semi_y, z0=i * segment, scale=scale)
        mesh = part if mesh is None else mesh + part
    return mesh


def build_bend(element, radius):
    points = np.asarray(element['design_path'], dtype=float)
    if len(points) < 2:
        return None
    points = (pose_matrix(element) @ np.hstack([points, np.ones((len(points), 1))]).T).T[:, :3]
    return pv.lines_from_points(points).tube(radius=radius, n_sides=24)


def effective_transverse(element, args, reference):
    """Return ((semi_x, semi_y), clamped) applying the clamp policy."""
    transverse = element.get('transverse')
    absurd = (
        transverse is None
        or transverse['semi_x'] > args.clamp_threshold
        or transverse['semi_y'] > args.clamp_threshold
    )
    if not absurd:
        return (transverse['semi_x'], transverse['semi_y']), False
    if args.no_clamp and transverse is not None:
        return (transverse['semi_x'], transverse['semi_y']), False
    if args.clamp_radius is not None:
        radius = args.clamp_radius
    elif reference is not None:
        radius = 0.5 * (reference[0] + reference[1])
    else:
        radius = 0.1
    return (radius, radius), True


def build_element(element, args, reference):
    """Build the local geometry of one element, or None if it has none."""
    (semi_x, semi_y), clamped = effective_transverse(element, args, reference)
    etype = element['type']
    length = element.get('length', 0.0)

    if etype in ('SBEND', 'RBEND') and element.get('design_path'):
        mesh = build_bend(element, max(semi_x, semi_y))
    elif etype == 'RBEND3D' and element.get('design_path'):
        mesh = build_bend(element, max(semi_x, semi_y))
    elif length <= 0.0:
        return None
    elif etype == 'SOLENOID':
        mesh = build_solenoid(length, semi_x, semi_y)
    elif etype == 'MULTIPOLE' and element.get('multipole_order') == 2:
        mesh = build_quadrupole(length, semi_x, semi_y)
    elif etype == 'RFCAVITY':
        mesh = build_cells(length, semi_x, semi_y, [0.78, 1.15, 0.82, 1.25, 0.82, 1.15, 0.78])
    elif etype == 'TRAVELINGWAVE':
        mesh = build_cells(length, semi_x, semi_y, [0.92, 1.05] * 4)
    elif element.get('transverse', {}).get('kind') == 'rectangular':
        mesh = box(length, semi_x, semi_y)
    elif etype in ('MONITOR', 'MARKER', 'PROBE'):
        mesh = cylinder(length, semi_x, semi_y, capped=True)
    else:
        mesh = cylinder(length, semi_x, semi_y)
    if mesh is None:
        return None

    if 'design_path' not in element or etype not in ('SBEND', 'RBEND', 'RBEND3D'):
        place(mesh, element)
    mesh.field_data['element'] = [element['name']]
    mesh.field_data['etype'] = [etype]
    mesh.field_data['clamped'] = [np.uint8(clamped)]
    return mesh


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('json', help='path to <base>_ElementPositions.json written by opalx')
    parser.add_argument('--out', help='output .vtm path (default: <json stem>.vtm)')
    parser.add_argument(
        '--clamp-radius',
        type=float,
        help='radius [m] for elements without a physical transverse extent '
        '(default: reference_transverse from the json, else 0.1)',
    )
    parser.add_argument(
        '--clamp-threshold',
        type=float,
        default=10.0,
        help='semi-axes larger than this [m] count as unset/infinite (default: 10)',
    )
    parser.add_argument(
        '--no-clamp', action='store_true', help='draw unset apertures at their raw size'
    )
    parser.add_argument('--hide', action='append', default=[], help='category to skip (repeatable)')
    parser.add_argument('--show', action='store_true', help='open an interactive pyvista window')
    parser.add_argument('--screenshot', help='write an offscreen render to this path')
    args = parser.parse_args(argv)

    with open(args.json) as fh:
        data = json.load(fh)
    reference = data.get('reference_transverse')
    if reference is not None:
        reference = (reference['semi_x'], reference['semi_y'])

    blocks = {key: pv.MultiBlock() for key in CATEGORY_COLORS}
    for element in data['elements']:
        category = category_of(element)
        if category in args.hide:
            continue
        mesh = build_element(element, args, reference)
        if mesh is None:
            continue
        blocks[category][element['name']] = mesh

    beamline = pv.MultiBlock()
    used = 0
    for category, group in blocks.items():
        if len(group):
            beamline[category] = group
            used += 1
    if used == 0:
        sys.exit('no element geometry built; nothing to write')

    out = Path(args.out) if args.out else Path(args.json).with_suffix('.vtm')
    beamline.save(out)
    n_elements = sum(len(beamline[key]) for key in beamline.keys())
    print(f'wrote {out} ({n_elements} elements in {len(beamline.keys())} categories)')

    if args.show or args.screenshot:
        plotter = pv.Plotter(off_screen=bool(args.screenshot) and not args.show)
        shown = []
        for category in beamline.keys():
            color = CATEGORY_COLORS.get(category, (0.5, 0.5, 0.5))
            plotter.add_mesh(beamline[category], color=color, smooth_shading=True)
            shown.append((category, color))
        plotter.add_legend(shown, bcolor=(1.0, 1.0, 1.0), face='circle', border=True)
        plotter.add_axes()
        if args.screenshot:
            plotter.screenshot(args.screenshot)
            print(f'wrote {args.screenshot}')
        if args.show:
            plotter.show()


if __name__ == '__main__':
    main()
