#!/usr/bin/env python3
"""Validate map-3 fields at logged reference positions, not at nominal s.

OrbitThreader logs its midpoint position and field to eight significant digits.
The comparison bound propagates that text-rounding uncertainty; it is not an
integration/map accuracy tolerance. Samples whose rounded position is ambiguous
at a hard quadrupole face are reported separately, not treated as field errors.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re

os.environ.setdefault("MPLCONFIGDIR", "/tmp/opalx-map3-mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/opalx-map3-cache")

from analytic_fields import bend_field, layout, load_parameters, np, pd
from run_best_cases import METHODS, REPO, ROOT, digest, parse_map, validate_physics


def rotation_y(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def placed_elements(p):
    """Independent nominal entrance poses, matching ELEMEDGE placement."""
    position, angle, previous_end = np.zeros(3), 0.0, 0.0
    for row in layout(p).itertuples():
        rotation = rotation_y(-angle)  # Local-to-global: positive bend turns to -x.
        origin = position + rotation @ np.array([0, 0, row.body_start_m - previous_end])
        yield row, origin, rotation
        length = row.body_end_m - row.body_start_m
        if row.kind == "bend":
            theta = np.deg2rad(p.body_angle_deg)
            position = origin + rotation @ np.array([p.radius_m * (np.cos(theta) - 1),
                                                     0, p.radius_m * np.sin(theta)])
            angle += theta
        else:
            position = origin + rotation @ np.array([0, 0, length])
        previous_end = row.body_end_m


def arc_coordinates(local, p):
    """Sector body with straight entrance/exit tangent extensions; source chart."""
    result = np.array(local, copy=True)
    theta = np.deg2rad(p.body_angle_deg)
    offset = local - [p.radius_m * (np.cos(theta) - 1), 0, p.radius_m * np.sin(theta)]
    exit_local = offset @ rotation_y(-theta)
    past_exit = (local[:, 2] > 0) & (exit_local[:, 2] >= 0)
    body = (local[:, 2] > 0) & ~past_exit
    result[past_exit, 0] = exit_local[past_exit, 0]
    result[past_exit, 2] = p.body_length + exit_local[past_exit, 2]
    result[body, 0] = np.hypot(local[body, 0] + p.radius_m, local[body, 2]) - p.radius_m
    result[body, 2] = p.radius_m * np.arctan2(local[body, 2], local[body, 0] + p.radius_m)
    return result


def global_field(positions, p):
    """Native near-axis law in global coordinates; only longitudinal masks here.

    Logged references are on y=0. This validator is not a general aperture or
    material model; the run's IndexMap is authoritative for actual membership.
    """
    positions = np.asarray(positions, float).reshape((-1, 3))
    field = np.zeros_like(positions)
    for row, origin, rotation in placed_elements(p):
        local = (positions - origin) @ rotation
        if row.kind == "bend":
            arc = arc_coordinates(local, p)
            bf = bend_field(arc[:, 2], arc[:, 0], arc[:, 1], p)
            phi = np.clip(arc[:, 2], 0, p.body_length) / p.radius_m
            # Same tangent-to-entrance rotation as GeometryHelper.
            entry = np.column_stack((np.cos(phi) * bf[:, 0] - np.sin(phi) * bf[:, 2],
                                     bf[:, 1], np.sin(phi) * bf[:, 0] + np.cos(phi) * bf[:, 2]))
            field += entry @ rotation.T
        elif row.kind == "quadrupole":
            active = (local[:, 2] >= 0) & (local[:, 2] < p.quadrupole_length_m)
            bf = np.column_stack((p.quad_gradient * local[:, 1], p.quad_gradient * local[:, 0],
                                  np.zeros(len(local))))
            field += np.where(active[:, None], bf @ rotation.T, 0)
    return field


def read_path(path):
    rows = []
    names = ["s", "x", "y", "z", "ux", "uy", "uz", "Ex", "Ey", "Ez", "Bx", "By", "Bz", "kinetic_MeV", "t_ns"]
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split(maxsplit=15)
        rows.append([*map(float, parts[:15]), parts[15].strip() if len(parts) > 15 else ""])
    data = pd.DataFrame(rows, columns=[*names, "elements"])
    if len(data) < 10 or not np.isfinite(data[names]).all().all():
        raise ValueError("Missing/nonfinite design-path samples")
    if not data.s.is_monotonic_increasing:
        raise ValueError("Nonmonotone reference path")
    return data


def rounding_bound(values):
    values = np.asarray(values, float)
    magnitudes = np.where(values == 0, 1, abs(values))
    return np.where(values == 0, 0, 0.51 * 10.0**(np.floor(np.log10(magnitudes)) - 7))


def validate_case(case, p):
    data = read_path(case / "data/map-3-dba_DesignPath.dat")
    positions = data[["x", "y", "z"]].to_numpy()
    measured = data[["Bx", "By", "Bz"]].to_numpy()
    expected = global_field(positions, p)
    position_bound = rounding_bound(positions)
    ambiguous = np.zeros(len(data), dtype=bool)
    for row, origin, rotation in placed_elements(p):
        if row.kind != "quadrupole":
            continue
        local = (positions - origin) @ rotation
        z_bound = position_bound @ abs(rotation[:, 2]) + 1e-12
        ambiguous |= (abs(local[:, 2]) <= z_bound) | (abs(local[:, 2] - p.quadrupole_length_m) <= z_bound)
    bound = rounding_bound(measured) + 1e-12
    for dim in range(3):
        delta = np.zeros_like(positions)
        delta[:, dim] = 1e-6
        gradient = (global_field(positions + delta, p) - global_field(positions - delta, p)) / 2e-6
        bound += 1.2 * abs(gradient) * position_bound[:, dim, None]
    errors = abs(measured - expected)
    valid = ~ambiguous
    ratios = errors[valid] / bound[valid]
    drifts = data.elements.str.fullmatch(r"D[12],\s*")
    drift_zero = bool((measured[drifts.to_numpy()] == 0).all())
    comparison = data.copy()
    for i, component in enumerate(("Bx", "By", "Bz")):
        comparison[f"{component}_analytic_T"] = expected[:, i]
        comparison[f"{component}_error_T"] = errors[:, i]
        comparison[f"{component}_rounding_bound_T"] = bound[:, i]
    comparison["ambiguous_quad_boundary"] = ambiguous
    comparison.to_csv(case / "field-validation.csv", index=False, float_format="%.12e")
    timing = (case / "timing.dat").read_text()
    orbit = re.search(r"^OrbThreader\.+\s+1\s+(\S+)", timing, flags=re.M)
    main = re.search(r"^mainTimer\.+\s+1\s+(\S+)", timing, flags=re.M)
    result = dict(samples=len(data), compared_samples=int(valid.sum()),
                  ambiguous_quad_boundary_samples=int(ambiguous.sum()),
                  max_field_difference_T=float(errors[valid].max()),
                  max_rounding_bound_ratio=float(ratios.max()),
                  field_matches_within_text_rounding=bool((ratios <= 1).all()),
                  drift_samples=int(drifts.sum()), drift_fields_exactly_zero=drift_zero,
                  max_abs_y_m=float(abs(data.y).max()), max_abs_uy=float(abs(data.uy).max()),
                  final_logged_s_m=float(data.s.iloc[-1]),
                  orbit_threader_wall_s=float(orbit.group(1)) if orbit else None,
                  total_wall_s=float(main.group(1)) if main else None,
                  source_hashes={Path(__file__).name: digest(__file__),
                                 "analytic_fields.py": digest(ROOT / "analytic_fields.py"),
                                 "DesignPath.dat": digest(case / "data/map-3-dba_DesignPath.dat")})
    (case / "field-validation.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "opalx-best-settings")
    args = parser.parse_args()
    p = load_parameters()
    provenance = json.loads((args.output / "provenance.json").read_text())
    for name, sha in provenance["source_hashes"].items():
        assert digest(REPO / name) == sha, name
    assert digest(provenance["executable"]) == provenance["executable_sha256"]
    for method in METHODS:
        case = args.output / method.lower()
        result = json.loads((case / "result.json").read_text())
        assert result["status"] == "OK" and result["ranks"] == 1
        for name, sha in result["artifacts_sha256"].items():
            assert digest(case / name) == sha, (method, name)
        source = (case / "map-3-dba.in").read_text()
        validate_physics(source, p)
        stdout = (case / "map-3-dba.out").read_text()
        line = next(line for line in stdout.splitlines() if "Linear-map starting perturbations" in line)
        steps = np.array([float(value) for value in line.split(":", 1)[1].split()])
        np.testing.assert_allclose(steps, [result["epsilon"]] * 6, rtol=1e-12, atol=0)
        np.testing.assert_array_equal(parse_map(stdout), np.loadtxt(case / "matrix.txt"))
    records = [dict(integration_method=method, **validate_case(args.output / method.lower(), p)) for method in METHODS]
    print(pd.DataFrame(records).drop(columns="source_hashes").to_string(index=False))
    pd.DataFrame(records).drop(columns="source_hashes").to_csv(args.output / "field-validation-summary.csv", index=False)
    return 0 if all(row["field_matches_within_text_rounding"] and row["drift_fields_exactly_zero"]
                    and row["drift_samples"] > 0 for row in records) else 1


if __name__ == "__main__":
    raise SystemExit(main())
