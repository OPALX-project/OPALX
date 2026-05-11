#!/usr/bin/env python3

from __future__ import annotations

import argparse
import io
import re
from pathlib import Path

import numpy as np
from PIL import Image


HEADER_RE = re.compile(r"([A-Za-z0-9_.\[\]]+)\s*=\s*(\([^)]*\)|\[[^\]]*\]|[^#\s]+)")
STEP_FILE_RE = re.compile(r"^(.*-)(\d{6})(\.[^.]+)$")
IP_DEF_RE = re.compile(
    r"^\s*(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*:\s*(?:BEAMBEAM|IP)\s*,(?P<body>.*);"
)
ATTR_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([^,;]+)")


def parse_vector(text: str) -> tuple[float, ...]:
    parts = re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", text)
    return tuple(float(part) for part in parts)


def parse_dump(path: Path) -> dict:
    metadata: dict[str, str] = {}
    rows: list[tuple[int, int, int, float, float, float, float]] = []

    with path.open() as stream:
        for line in stream:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                for key, value in HEADER_RE.findall(line[1:].strip()):
                    metadata[key] = value.strip()
                continue

            parts = line.split()
            if len(parts) < 7:
                continue

            rows.append(
                (
                    int(parts[0]),
                    int(parts[1]),
                    int(parts[2]),
                    float(parts[3]),
                    float(parts[4]),
                    float(parts[5]),
                    float(parts[6]),
                )
            )

    if not rows:
        raise ValueError(f"no rho samples found in {path}")

    data = np.asarray(rows, dtype=float)
    i_vals = data[:, 0].astype(int)
    j_vals = data[:, 1].astype(int)
    k_vals = data[:, 2].astype(int)

    imin, jmin, kmin = i_vals.min(), j_vals.min(), k_vals.min()
    nx = i_vals.max() - imin + 1
    ny = j_vals.max() - jmin + 1
    nz = k_vals.max() - kmin + 1

    xs = np.full(nx, np.nan)
    ys = np.full(ny, np.nan)
    zs = np.full(nz, np.nan)
    rho = np.zeros((nx, ny, nz), dtype=float)

    for i_raw, j_raw, k_raw, x, y, z, rho_val in rows:
        i = i_raw - imin
        j = j_raw - jmin
        k = k_raw - kmin
        xs[i] = x
        ys[j] = y
        zs[k] = z
        rho[i, j, k] = rho_val

    if np.isnan(xs).any() or np.isnan(ys).any() or np.isnan(zs).any():
        raise ValueError(f"incomplete rho grid in {path}")

    return {
        "path": path,
        "metadata": metadata,
        "is": np.arange(imin, imin + nx),
        "js": np.arange(jmin, jmin + ny),
        "ks": np.arange(kmin, kmin + nz),
        "xs": xs,
        "ys": ys,
        "zs": zs,
        "rho": rho,
    }


def step_from_filename(path: Path) -> int | None:
    match = STEP_FILE_RE.match(path.name)
    if match is None:
        return None
    return int(match.group(2))


def resolve_step_paths(seed: Path, start: int, end: int) -> list[Path]:
    match = STEP_FILE_RE.match(seed.name)
    if match is None:
        raise ValueError("filename must end in -XXXXXX.<ext> so step ranges can be resolved")

    prefix, _step, suffix = match.groups()
    lo = min(start, end)
    hi = max(start, end)
    paths = [seed.with_name(f"{prefix}{step:06d}{suffix}") for step in range(lo, hi + 1)]
    missing = [path for path in paths if not path.exists()]
    if missing:
        preview = ", ".join(str(path) for path in missing[:5])
        raise FileNotFoundError(f"missing step files: {preview}")
    return paths


def collect_step_paths(seed: Path) -> list[Path]:
    match = STEP_FILE_RE.match(seed.name)
    if match is None:
        raise ValueError("filename must end in -XXXXXX.<ext> so step ranges can be resolved")

    prefix, _step, suffix = match.groups()
    candidates = []
    for path in seed.parent.glob(f"{prefix}*{suffix}"):
        match = STEP_FILE_RE.match(path.name)
        if match is None or match.group(1) != prefix or match.group(3) != suffix:
            continue
        candidates.append(path)

    return sorted(candidates)


def default_output_path(seed: Path, start: int, end: int) -> Path:
    match = STEP_FILE_RE.match(seed.name)
    if match is None:
        suffix = ".png" if start == end else ".gif"
        return seed.with_suffix(suffix)

    prefix, _step, _suffix = match.groups()
    if start == end:
        return seed.with_suffix(".png")

    lo = min(start, end)
    hi = max(start, end)
    return seed.with_name(f"{prefix}{lo:06d}-{hi:06d}.gif")


def normalize_output_path(output: Path, start: int, end: int) -> Path:
    if start == end:
        return output

    if output.suffix.lower() != ".gif":
        return output.with_suffix(".gif")

    return output


def filter_paths_by_unique_global_step(paths: list[Path]) -> list[Path]:
    filtered: list[Path] = []
    seen_steps: set[int] = set()
    for path in paths:
        dump = parse_dump(path)
        step_text = dump["metadata"].get("global_step")
        if step_text is None:
            filtered.append(path)
            continue
        step = int(parse_numeric(step_text))
        if step in seen_steps:
            continue
        seen_steps.add(step)
        filtered.append(path)
    return filtered


def axis_interval_from_centers(values: np.ndarray) -> tuple[float, float]:
    if len(values) == 1:
        center = float(values[0])
        return center - 0.5, center + 0.5

    lower_spacing = 0.5 * float(values[1] - values[0])
    upper_spacing = 0.5 * float(values[-1] - values[-2])
    return float(values[0] - lower_spacing), float(values[-1] + upper_spacing)


def axis_extent(a: np.ndarray, b: np.ndarray) -> list[float]:
    ax0, ax1 = axis_interval_from_centers(a)
    bx0, bx1 = axis_interval_from_centers(b)
    return [ax0, ax1, bx0, bx1]


def longitudinal_display_delta(dump: dict, use_grid: bool) -> float:
    if use_grid:
        return 0.0

    metadata = dump["metadata"]
    path_length_s = metadata_float(metadata, "path_length_s")
    if path_length_s is None:
        return 0.0

    return path_length_s


def color_limits(rho: np.ndarray) -> tuple[str, float, float]:
    rho_min = float(np.min(rho))
    rho_max = float(np.max(rho))

    if rho_max <= 0.0:
        return "Blues_r", rho_min, 0.0
    if rho_min >= 0.0:
        return "Reds", 0.0, rho_max

    vmax = max(abs(rho_min), abs(rho_max))
    return "seismic", -vmax, vmax


def format_float(value: float, digits: int = 6) -> str:
    return f"{float(value):.{digits}f}"


def format_sci(value: float, digits: int = 3) -> str:
    return f"{float(value):.{digits}e}"


def format_range(values: tuple[float, float] | None, digits: int = 6) -> str | None:
    if values is None:
        return None
    return f"[{format_float(values[0], digits)}, {format_float(values[1], digits)}]"


def parse_numeric(text: str) -> float:
    match = re.search(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", text)
    if match is None:
        raise ValueError(f"cannot parse numeric value from {text!r}")
    return float(match.group(0))


def find_default_input(seed: Path) -> Path | None:
    candidates: list[Path] = []
    data_dir = seed.parent
    if data_dir.name == "data":
        case_dir = data_dir.parent
        candidates.extend(
            [
                case_dir / f"{case_dir.name}.in",
                case_dir / "ip-1.in",
            ]
        )
    candidates.extend(
        [
            seed.parent / "ip-1.in",
            seed.parent / f"{seed.parent.name}.in",
        ]
    )

    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def parse_ip_geometry_from_input(path: Path) -> dict[str, float]:
    with path.open() as stream:
        for line in stream:
            match = IP_DEF_RE.match(line)
            if match is None:
                continue

            attributes = {
                key.upper(): value.strip()
                for key, value in ATTR_RE.findall(match.group("body"))
            }
            if "ELEMEDGE" not in attributes or "L" not in attributes:
                raise ValueError(f"IP definition in {path} is missing ELEMEDGE or L")

            element_begin_s = parse_numeric(attributes["ELEMEDGE"])
            element_length = parse_numeric(attributes["L"])
            element_end_s = element_begin_s + element_length
            interaction_point_s = element_begin_s + 0.5 * element_length

            geometry = {
                "interaction_point_s": interaction_point_s,
                "ip_begin_s": element_begin_s,
                "ip_end_s": element_end_s,
                "window_begin_s": element_begin_s,
                "window_end_s": element_end_s,
            }

            return geometry

    raise ValueError(f"no BEAMBEAM element definition found in {path}")


def build_input_geometry(args, seed: Path) -> dict[str, float]:
    geometry: dict[str, float] = {}

    if args.input is not None:
        geometry.update(parse_ip_geometry_from_input(Path(args.input)))
    else:
        default_input = find_default_input(seed)
        if default_input is not None:
            geometry.update(parse_ip_geometry_from_input(default_input))

    overrides = {
        "interaction_point_s": args.bb_s,
        "ip_begin_s": args.bb_begin_s,
        "ip_end_s": args.bb_end_s,
        "window_begin_s": args.window_begin_s,
        "window_end_s": args.window_end_s,
    }
    geometry.update({key: value for key, value in overrides.items() if value is not None})
    return geometry


def infer_inside_ip_paths(seed: Path, geometry: dict[str, float]) -> list[Path]:
    ip_begin_s = geometry.get("ip_begin_s")
    ip_end_s = geometry.get("ip_end_s")
    if ip_begin_s is None or ip_end_s is None:
        return []

    inside_paths: list[Path] = []
    for path in collect_step_paths(seed):
        dump = parse_dump(path)
        bounds_rmin = metadata_vector(dump["metadata"], "bunch_bounds_rmin")
        bounds_rmax = metadata_vector(dump["metadata"], "bunch_bounds_rmax")
        if bounds_rmin is not None and bounds_rmax is not None:
            bunch_tail_s = float(bounds_rmin[2])
            bunch_head_s = float(bounds_rmax[2])
            if bunch_head_s >= ip_begin_s and bunch_tail_s <= ip_end_s:
                inside_paths.append(path)
            continue

        path_length_s = metadata_float(dump["metadata"], "path_length_s")
        if path_length_s is not None and ip_begin_s <= path_length_s <= ip_end_s:
            inside_paths.append(path)

    return inside_paths


def metadata_float(metadata: dict[str, str], key: str) -> float | None:
    value = metadata.get(key)
    if value is None:
        return None
    return parse_numeric(value)


def metadata_vector(metadata: dict[str, str], key: str) -> tuple[float, ...] | None:
    value = metadata.get(key)
    if value is None:
        return None
    return tuple(parse_vector(value))


def primary_center_from_bounds(metadata: dict[str, str]) -> np.ndarray | None:
    rmin = metadata_vector(metadata, "bunch_bounds_rmin")
    rmax = metadata_vector(metadata, "bunch_bounds_rmax")
    if rmin is None or rmax is None or len(rmin) < 3 or len(rmax) < 3:
        return None
    return np.asarray(
        tuple(0.5 * (lo + hi) for lo, hi in zip(rmin[:3], rmax[:3])),
        dtype=float,
    )


def primary_center_from_mean(metadata: dict[str, str]) -> np.ndarray | None:
    mean_r = metadata_vector(metadata, "particle_mean_r")
    if mean_r is None or len(mean_r) < 3:
        return None
    return np.asarray(mean_r[:3], dtype=float)


def primary_center_from_rho(dump: dict) -> np.ndarray | None:
    rho = np.asarray(dump["rho"], dtype=float)
    weights = np.maximum(-rho, 0.0)
    weight_sum = float(weights.sum())
    if weight_sum <= 0.0:
        return None

    xs = np.asarray(dump["xs"], dtype=float)
    ys = np.asarray(dump["ys"], dtype=float)
    zs = np.asarray(dump["zs"], dtype=float)

    x_proj = weights.sum(axis=(1, 2))
    y_proj = weights.sum(axis=(0, 2))
    z_proj = weights.sum(axis=(0, 1))

    return np.asarray(
        (
            float((x_proj * xs).sum() / weight_sum),
            float((y_proj * ys).sum() / weight_sum),
            float((z_proj * zs).sum() / weight_sum),
        ),
        dtype=float,
    )


def primary_center_from_active_rho(
    dump: dict,
    interaction_point_z: float,
    primary_on_positive_z: bool,
) -> np.ndarray | None:
    rho = np.asarray(dump["rho"], dtype=float)
    zs = np.asarray(dump["zs"], dtype=float)

    if primary_on_positive_z:
        z_mask = zs >= interaction_point_z
    else:
        z_mask = zs <= interaction_point_z

    if not np.any(z_mask):
        return None

    rho_half = rho[:, :, z_mask]
    weights = np.maximum(-rho_half, 0.0)
    weight_sum = float(weights.sum())
    if weight_sum <= 0.0:
        return None

    xs = np.asarray(dump["xs"], dtype=float)
    ys = np.asarray(dump["ys"], dtype=float)
    zs_half = zs[z_mask]

    x_proj = weights.sum(axis=(1, 2))
    y_proj = weights.sum(axis=(0, 2))
    z_proj = weights.sum(axis=(0, 1))

    return np.asarray(
        (
            float((x_proj * xs).sum() / weight_sum),
            float((y_proj * ys).sum() / weight_sum),
            float((z_proj * zs_half).sum() / weight_sum),
        ),
        dtype=float,
    )


def convert_primary_center_to_local_z(
    primary_center: np.ndarray | None,
    metadata: dict[str, str],
    geometry: dict[str, float],
) -> np.ndarray | None:
    if primary_center is None:
        return None

    primary_local = np.asarray(primary_center, dtype=float).copy()
    interaction_window_active = int(parse_numeric(metadata.get("interaction_window_active", "0")))
    if interaction_window_active == 0:
        return primary_local

    ip_s = geometry.get("interaction_point_s")
    if ip_s is None:
        return primary_local

    interaction_point_local_z = (
        metadata_float(metadata, "interaction_point_beam_z")
        if "interaction_point_beam_z" in metadata
        else metadata_float(metadata, "interaction_point_local_z[m]")
    )
    path_length_s = metadata_float(metadata, "path_length_s")
    if path_length_s is None and interaction_point_local_z is not None:
        path_length_s = ip_s - interaction_point_local_z
    if path_length_s is None:
        return primary_local

    primary_local[2] = primary_local[2] - path_length_s
    return primary_local


def derive_overlay_geometry(dump: dict, geometry: dict[str, float]) -> dict[str, object]:
    metadata = dump["metadata"]
    overlay: dict[str, object] = {
        "interaction_point": None,
        "interaction_point_s": None,
        "ip_element_z_range": None,
        "collision_window_z_range": None,
        "ip_element_s_range": None,
        "collision_window_s_range": None,
        "primary_center": None,
        "mirrored_center": None,
    }

    if "interaction_point" in metadata:
        overlay["interaction_point"] = np.asarray(parse_vector(metadata["interaction_point"]), dtype=float)
        if "interaction_point_s" in metadata:
            overlay["interaction_point_s"] = parse_numeric(metadata["interaction_point_s"])
        if "ip_element_z_range" in metadata:
            overlay["ip_element_z_range"] = tuple(parse_vector(metadata["ip_element_z_range"]))
        if "collision_window_z_range" in metadata:
            overlay["collision_window_z_range"] = tuple(parse_vector(metadata["collision_window_z_range"]))
        elif overlay["ip_element_z_range"] is not None:
            overlay["collision_window_z_range"] = overlay["ip_element_z_range"]
        if "ip_element_s_range" in metadata:
            overlay["ip_element_s_range"] = tuple(parse_vector(metadata["ip_element_s_range"]))
        if "collision_window_s_range" in metadata:
            overlay["collision_window_s_range"] = tuple(parse_vector(metadata["collision_window_s_range"]))
        elif overlay["ip_element_s_range"] is not None:
            overlay["collision_window_s_range"] = overlay["ip_element_s_range"]
        interaction_window_active = int(parse_numeric(metadata.get("interaction_window_active", "0")))
        mean_center = primary_center_from_mean(metadata)
        bounds_center = primary_center_from_bounds(metadata)
        if interaction_window_active == 1:
            overlay["primary_center"] = mean_center
            if overlay["primary_center"] is None:
                overlay["primary_center"] = convert_primary_center_to_local_z(
                    bounds_center,
                    metadata,
                    geometry,
                )
            if overlay["primary_center"] is None:
                ip_s = overlay["interaction_point_s"]
                primary_on_positive_z = True
                if ip_s is not None and bounds_center is not None:
                    primary_on_positive_z = bool(bounds_center[2] >= ip_s)
                overlay["primary_center"] = primary_center_from_active_rho(
                    dump,
                    float(overlay["interaction_point"][2]),
                    primary_on_positive_z,
                )
        else:
            overlay["primary_center"] = mean_center
            if overlay["primary_center"] is None:
                overlay["primary_center"] = primary_center_from_rho(dump)
        active_flag = int(parse_numeric(metadata.get("interaction_window_active", "1")))
        if overlay["primary_center"] is not None and active_flag == 1:
            ip_point = np.asarray(overlay["interaction_point"], dtype=float)
            primary = np.asarray(overlay["primary_center"], dtype=float)
            distance = abs(primary[2] - ip_point[2])
            if primary[2] <= ip_point[2]:
                mirror_z = ip_point[2] + distance
            else:
                mirror_z = ip_point[2] - distance
            overlay["mirrored_center"] = np.asarray((primary[0], primary[1], mirror_z), dtype=float)
        return overlay

    ip_s = geometry.get("interaction_point_s")
    ip_begin_s = geometry.get("ip_begin_s")
    ip_end_s = geometry.get("ip_end_s")
    window_begin_s = geometry.get("window_begin_s")
    window_end_s = geometry.get("window_end_s")

    if ip_s is None:
        return overlay

    interaction_point_local_z = (
        metadata_float(metadata, "interaction_point_beam_z")
        if "interaction_point_beam_z" in metadata
        else metadata_float(metadata, "interaction_point_local_z[m]")
    )
    path_length_s = metadata_float(metadata, "path_length_s")
    interaction_window_active = int(parse_numeric(metadata.get("interaction_window_active", "0")))

    if interaction_point_local_z is not None:
        interaction_point_z = interaction_point_local_z
    elif path_length_s is not None:
        interaction_point_z = ip_s - path_length_s
    elif interaction_window_active == 0:
        interaction_point_z = ip_s
    else:
        return overlay

    overlay["interaction_point"] = np.asarray((0.0, 0.0, interaction_point_z), dtype=float)
    overlay["interaction_point_s"] = ip_s

    if ip_begin_s is not None and ip_end_s is not None:
        overlay["ip_element_s_range"] = (ip_begin_s, ip_end_s)
        if interaction_point_local_z is not None or path_length_s is not None:
            overlay["ip_element_z_range"] = (
                ip_begin_s - ip_s + interaction_point_z,
                ip_end_s - ip_s + interaction_point_z,
            )
        else:
            overlay["ip_element_z_range"] = (ip_begin_s, ip_end_s)

    if window_begin_s is not None and window_end_s is not None:
        overlay["collision_window_s_range"] = (window_begin_s, window_end_s)
        if interaction_point_local_z is not None or path_length_s is not None:
            overlay["collision_window_z_range"] = (
                window_begin_s - ip_s + interaction_point_z,
                window_end_s - ip_s + interaction_point_z,
            )
        else:
            overlay["collision_window_z_range"] = (window_begin_s, window_end_s)

    mean_center = primary_center_from_mean(metadata)
    bounds_center = primary_center_from_bounds(metadata)
    if interaction_window_active == 1:
        overlay["primary_center"] = mean_center
        if overlay["primary_center"] is None:
            overlay["primary_center"] = convert_primary_center_to_local_z(
                bounds_center,
                metadata,
                geometry,
            )
        if overlay["primary_center"] is None:
            primary_on_positive_z = True
            if ip_s is not None and bounds_center is not None:
                primary_on_positive_z = bool(bounds_center[2] >= ip_s)
            overlay["primary_center"] = primary_center_from_active_rho(
                dump,
                interaction_point_z,
                primary_on_positive_z,
            )
    else:
        overlay["primary_center"] = mean_center
        if overlay["primary_center"] is None:
            overlay["primary_center"] = primary_center_from_rho(dump)
    if (
        overlay["primary_center"] is not None
        and overlay["interaction_point"] is not None
        and interaction_window_active == 1
    ):
        ip_point = np.asarray(overlay["interaction_point"], dtype=float)
        primary = np.asarray(overlay["primary_center"], dtype=float)
        distance = abs(primary[2] - ip_point[2])
        if primary[2] <= ip_point[2]:
            mirror_z = ip_point[2] + distance
        else:
            mirror_z = ip_point[2] - distance
        overlay["mirrored_center"] = np.asarray((primary[0], primary[1], mirror_z), dtype=float)

    return overlay


def empty_overlay() -> dict[str, object]:
    return {
        "interaction_point": None,
        "interaction_point_s": None,
        "ip_element_z_range": None,
        "collision_window_z_range": None,
        "ip_element_s_range": None,
        "collision_window_s_range": None,
        "primary_center": None,
        "mirrored_center": None,
    }


def freeze_active_window_geometry(
    dumps: list[dict],
    overlays: list[dict[str, object]],
) -> list[dict[str, object]]:
    return overlays


def shift_overlay_z(overlay: dict[str, object], delta_z: float) -> dict[str, object]:
    shifted = dict(overlay)

    if shifted["interaction_point"] is not None:
        point = np.asarray(shifted["interaction_point"], dtype=float).copy()
        point[2] += delta_z
        shifted["interaction_point"] = point

    if shifted["primary_center"] is not None:
        point = np.asarray(shifted["primary_center"], dtype=float).copy()
        point[2] += delta_z
        shifted["primary_center"] = point

    if shifted["mirrored_center"] is not None:
        point = np.asarray(shifted["mirrored_center"], dtype=float).copy()
        point[2] += delta_z
        shifted["mirrored_center"] = point

    if shifted["ip_element_z_range"] is not None:
        z0, z1 = shifted["ip_element_z_range"]
        shifted["ip_element_z_range"] = (z0 + delta_z, z1 + delta_z)

    if shifted["collision_window_z_range"] is not None:
        z0, z1 = shifted["collision_window_z_range"]
        shifted["collision_window_z_range"] = (z0 + delta_z, z1 + delta_z)

    return shifted


def align_longitudinal_display_frame(
    dumps: list[dict],
    overlays: list[dict[str, object]],
    use_grid: bool,
) -> tuple[list[dict], list[dict[str, object]]]:
    if use_grid:
        return dumps, overlays

    reference_ip_z: float | None = None
    for dump, overlay in zip(dumps, overlays):
        active_flag = int(parse_numeric(dump["metadata"].get("interaction_window_active", "0")))
        if active_flag != 1 or overlay["interaction_point"] is None:
            continue
        reference_ip_z = float(np.asarray(overlay["interaction_point"], dtype=float)[2])
        break

    if reference_ip_z is None:
        for overlay in overlays:
            if overlay["interaction_point"] is not None:
                reference_ip_z = float(np.asarray(overlay["interaction_point"], dtype=float)[2])
                break

    if reference_ip_z is None:
        return dumps, overlays

    aligned_dumps: list[dict] = []
    aligned_overlays: list[dict[str, object]] = []

    for dump, overlay in zip(dumps, overlays):
        active_flag = int(parse_numeric(dump["metadata"].get("interaction_window_active", "0")))
        if overlay["interaction_point"] is None:
            aligned_dumps.append(dump)
            aligned_overlays.append(overlay)
            continue

        if active_flag == 1:
            aligned_dumps.append(dump)
            aligned_overlays.append(overlay)
            continue

        current_ip_z = float(np.asarray(overlay["interaction_point"], dtype=float)[2])
        delta_z = reference_ip_z - current_ip_z
        if abs(delta_z) < 1.0e-15:
            aligned_dumps.append(dump)
            aligned_overlays.append(overlay)
            continue

        shifted_dump = dict(dump)
        shifted_dump["zs"] = np.asarray(dump["zs"], dtype=float) + delta_z
        aligned_dumps.append(shifted_dump)
        aligned_overlays.append(shift_overlay_z(overlay, delta_z))

    return aligned_dumps, aligned_overlays


def convert_point_to_grid(
    point: np.ndarray,
    physical_coords: tuple[np.ndarray, np.ndarray, np.ndarray],
    index_coords: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> np.ndarray:
    converted = np.zeros(3, dtype=float)
    for axis, values in enumerate(physical_coords):
        if len(values) == 1:
            converted[axis] = float(index_coords[axis][0])
            continue
        spacing = values[1] - values[0]
        converted[axis] = index_coords[axis][0] + (point[axis] - values[0]) / spacing
    return converted


def convert_z_range_to_grid(
    z_range: tuple[float, float] | None,
    dump: dict,
) -> tuple[float, float] | None:
    if z_range is None:
        return None
    z0 = convert_point_to_grid(
        np.asarray((0.0, 0.0, z_range[0])),
        (dump["xs"], dump["ys"], dump["zs"]),
        (dump["is"], dump["js"], dump["ks"]),
    )[2]
    z1 = convert_point_to_grid(
        np.asarray((0.0, 0.0, z_range[1])),
        (dump["xs"], dump["ys"], dump["zs"]),
        (dump["is"], dump["js"], dump["ks"]),
    )[2]
    return (float(z0), float(z1))


def add_ip_overlay(ax, x_value: float, y_value: float, x_label: str, y_label: str) -> None:
    ax.axvline(x_value, color="limegreen", linestyle="--", linewidth=1.2, alpha=0.9)
    ax.axhline(y_value, color="limegreen", linestyle="--", linewidth=1.2, alpha=0.9)
    ax.plot([x_value], [y_value], marker="o", markersize=6, color="black", markerfacecolor="gold")
    ax.text(
        x_value,
        y_value,
        " IP",
        color="black",
        fontsize=9,
        va="bottom",
        ha="left",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.7, "pad": 0.2},
    )
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)


def add_bunch_center(
    ax,
    x_value: float,
    y_value: float,
    *,
    marker: str,
    edgecolor: str,
    facecolor: str,
    label: str,
    x_offset: float,
) -> None:
    ax.plot(
        [x_value],
        [y_value],
        marker=marker,
        markersize=7,
        markeredgewidth=1.3,
        markeredgecolor=edgecolor,
        markerfacecolor=facecolor,
        color=edgecolor,
        linestyle="None",
    )
    ax.text(
        x_value + x_offset,
        y_value,
        f" {label}",
        color=edgecolor,
        fontsize=8,
        va="top",
        ha="left",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.6, "pad": 0.2},
    )


def expand_limits_to_include_point(
    ax,
    extent: list[float],
    x_value: float,
    y_value: float,
) -> None:
    xmin, xmax, ymin, ymax = extent
    xpad = 0.05 * max(xmax - xmin, abs(x_value), 1.0)
    ypad = 0.05 * max(ymax - ymin, abs(y_value), 1.0)

    ax.set_xlim(min(xmin, x_value) - xpad, max(xmax, x_value) + xpad)
    ax.set_ylim(min(ymin, y_value) - ypad, max(ymax, y_value) + ypad)


def draw_domain_box(ax, extent: list[float]) -> None:
    import matplotlib.patches as patches

    xmin, xmax, ymin, ymax = extent
    rect = patches.Rectangle(
        (xmin, ymin),
        xmax - xmin,
        ymax - ymin,
        fill=False,
        edgecolor="black",
        linewidth=1.0,
        linestyle=":",
        alpha=0.9,
    )
    ax.add_patch(rect)


def draw_longitudinal_band(
    ax,
    x_extent: list[float],
    z_range: tuple[float, float],
    *,
    edgecolor: str,
    linestyle: str,
    label: str,
) -> None:
    import matplotlib.patches as patches

    xmin, xmax, _, _ = x_extent
    zmin, zmax = z_range
    rect = patches.Rectangle(
        (xmin, zmin),
        xmax - xmin,
        zmax - zmin,
        fill=False,
        edgecolor=edgecolor,
        linewidth=1.4,
        linestyle=linestyle,
        alpha=0.95,
    )
    ax.add_patch(rect)
    ax.text(
        xmax,
        zmax,
        f" {label}",
        color=edgecolor,
        fontsize=9,
        va="bottom",
        ha="left",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.6, "pad": 0.2},
    )


def same_z_range(
    left: tuple[float, float] | None,
    right: tuple[float, float] | None,
    tol: float = 1.0e-12,
) -> bool:
    if left is None or right is None:
        return False
    return abs(left[0] - right[0]) <= tol and abs(left[1] - right[1]) <= tol


def transformed_overlay(
    dump: dict,
    overlay: dict[str, object],
    use_grid: bool,
) -> tuple[np.ndarray | None, tuple[float, float] | None, tuple[float, float] | None]:
    interaction_point = overlay["interaction_point"]
    ip_element_z_range = overlay["ip_element_z_range"]
    collision_window_z_range = overlay["collision_window_z_range"]
    primary_center = overlay["primary_center"]
    mirrored_center = overlay["mirrored_center"]

    if interaction_point is None:
        return None, ip_element_z_range, collision_window_z_range, primary_center, mirrored_center

    interaction_point = np.asarray(interaction_point, dtype=float)
    if use_grid:
        interaction_point = convert_point_to_grid(
            interaction_point,
            (dump["xs"], dump["ys"], dump["zs"]),
            (dump["is"], dump["js"], dump["ks"]),
        )
        ip_element_z_range = convert_z_range_to_grid(ip_element_z_range, dump)
        collision_window_z_range = convert_z_range_to_grid(collision_window_z_range, dump)
        if primary_center is not None:
            primary_center = convert_point_to_grid(
                np.asarray(primary_center, dtype=float),
                (dump["xs"], dump["ys"], dump["zs"]),
                (dump["is"], dump["js"], dump["ks"]),
            )
        if mirrored_center is not None:
            mirrored_center = convert_point_to_grid(
                np.asarray(mirrored_center, dtype=float),
                (dump["xs"], dump["ys"], dump["zs"]),
                (dump["is"], dump["js"], dump["ks"]),
            )
    else:
        delta_z = longitudinal_display_delta(dump, use_grid)
        if abs(delta_z) > 0.0:
            interaction_point = interaction_point.copy()
            interaction_point[2] += delta_z
            if ip_element_z_range is not None:
                ip_element_z_range = (
                    ip_element_z_range[0] + delta_z,
                    ip_element_z_range[1] + delta_z,
                )
            if collision_window_z_range is not None:
                collision_window_z_range = (
                    collision_window_z_range[0] + delta_z,
                    collision_window_z_range[1] + delta_z,
                )
            if primary_center is not None:
                primary_center = np.asarray(primary_center, dtype=float).copy()
                primary_center[2] += delta_z
            if mirrored_center is not None:
                mirrored_center = np.asarray(mirrored_center, dtype=float).copy()
                mirrored_center[2] += delta_z

    return interaction_point, ip_element_z_range, collision_window_z_range, primary_center, mirrored_center


def compute_global_limits(
    dumps: list[dict],
    overlays: list[dict[str, object]],
    use_grid: bool,
) -> dict[str, object]:
    rho_min = min(float(np.min(dump["rho"])) for dump in dumps)
    rho_max = max(float(np.max(dump["rho"])) for dump in dumps)
    cmap, vmin, vmax = color_limits(np.asarray([rho_min, rho_max]))
    if vmin == vmax:
        vmax = vmin + 1.0

    xmins: list[float] = []
    xmaxs: list[float] = []
    ymins: list[float] = []
    ymaxs: list[float] = []
    zmins: list[float] = []
    zmaxs: list[float] = []

    for dump, overlay in zip(dumps, overlays):
        xs = dump["is"] if use_grid else dump["xs"]
        ys = dump["js"] if use_grid else dump["ys"]
        zs = dump["ks"] if use_grid else (np.asarray(dump["zs"], dtype=float) +
                                          longitudinal_display_delta(dump, use_grid))
        xmins.append(float(xs[0]))
        xmaxs.append(float(xs[-1]))
        ymins.append(float(ys[0]))
        ymaxs.append(float(ys[-1]))
        zmins.append(float(zs[0]))
        zmaxs.append(float(zs[-1]))

        (
            interaction_point,
            ip_element_z_range,
            collision_window_z_range,
            primary_center,
            mirrored_center,
        ) = transformed_overlay(dump, overlay, use_grid)
        if interaction_point is not None:
            xmins.append(float(interaction_point[0]))
            xmaxs.append(float(interaction_point[0]))
            ymins.append(float(interaction_point[1]))
            ymaxs.append(float(interaction_point[1]))
            zmins.append(float(interaction_point[2]))
            zmaxs.append(float(interaction_point[2]))
        for point in (primary_center, mirrored_center):
            if point is None:
                continue
            xmins.append(float(point[0]))
            xmaxs.append(float(point[0]))
            ymins.append(float(point[1]))
            ymaxs.append(float(point[1]))
            zmins.append(float(point[2]))
            zmaxs.append(float(point[2]))
        if ip_element_z_range is not None:
            zmins.append(float(ip_element_z_range[0]))
            zmaxs.append(float(ip_element_z_range[1]))
        if collision_window_z_range is not None:
            zmins.append(float(collision_window_z_range[0]))
            zmaxs.append(float(collision_window_z_range[1]))

    def padded(lo: float, hi: float) -> tuple[float, float]:
        span = hi - lo
        pad = 0.05 * max(span, abs(lo), abs(hi), 1.0)
        return lo - pad, hi + pad

    xlim = padded(min(xmins), max(xmaxs))
    ylim = padded(min(ymins), max(ymaxs))
    zlim = padded(min(zmins), max(zmaxs))
    if use_grid:
        global_i = np.asarray([idx for dump in dumps for idx in dump["is"]], dtype=float)
        global_j = np.asarray([idx for dump in dumps for idx in dump["js"]], dtype=float)
        global_k = np.asarray([idx for dump in dumps for idx in dump["ks"]], dtype=float)
        xlim = axis_interval_from_centers(global_i)
        ylim = axis_interval_from_centers(global_j)
        zlim = axis_interval_from_centers(global_k)
    else:
        xlim = (-0.01, 0.01)
        ylim = (-0.01, 0.01)
    return {"cmap": cmap, "vmin": vmin, "vmax": vmax, "xlim": xlim, "ylim": ylim, "zlim": zlim}


def build_title(dump: dict, overlay: dict[str, object]) -> str:
    metadata = dump["metadata"]
    details: list[str] = []

    if "global_step" in metadata:
        details.append(f"frame={int(parse_numeric(metadata['global_step']))}")

    field_origin, mesh_spacing, total_charge = compute_mesh_diagnostics(dump, use_grid=False)
    if field_origin is not None:
        details.append(
            "origin=(%s, %s, %s)"
            % (
                format_sci(field_origin[0]),
                format_sci(field_origin[1]),
                format_sci(field_origin[2]),
            )
        )
    if mesh_spacing is not None:
        details.append(
            "h=(%s, %s, %s)"
            % (
                format_sci(mesh_spacing[0]),
                format_sci(mesh_spacing[1]),
                format_sci(mesh_spacing[2]),
            )
        )
    details.append("Q=%s C" % format_sci(total_charge))

    return " | ".join(details)


def bunch_position_lab(dump: dict, center: np.ndarray | None) -> float | None:
    if center is None:
        return None

    metadata = dump["metadata"]
    delta_z = longitudinal_display_delta(dump, use_grid=False)
    path_length_s = metadata_float(metadata, "path_length_s")
    point = np.asarray(center, dtype=float)
    if path_length_s is not None:
        return float(point[2] + path_length_s)
    return float(point[2] + delta_z)


def format_bunch_position(dump: dict, center: np.ndarray | None) -> str | None:
    position = bunch_position_lab(dump, center)
    if position is None:
        return None
    return format_float(position, 4)


def compute_mesh_diagnostics(
    dump: dict,
    use_grid: bool,
) -> tuple[tuple[float, float, float] | None, tuple[float, float, float] | None, float]:
    metadata = dump["metadata"]

    field_origin = metadata_vector(metadata, "field_domain_rmin")
    mesh_spacing = metadata_vector(metadata, "mesh_spacing")
    if field_origin is None:
        field_origin = metadata_vector(metadata, "mesh_origin")

    if field_origin is not None and not use_grid:
        delta_z = longitudinal_display_delta(dump, use_grid=False)
        field_origin = (field_origin[0], field_origin[1], field_origin[2] + delta_z)

    if mesh_spacing is None:
        xs = np.asarray(dump["xs"], dtype=float)
        ys = np.asarray(dump["ys"], dtype=float)
        zs = np.asarray(dump["zs"], dtype=float)
        mesh_spacing = (
            float(xs[1] - xs[0]) if len(xs) > 1 else 0.0,
            float(ys[1] - ys[0]) if len(ys) > 1 else 0.0,
            float(zs[1] - zs[0]) if len(zs) > 1 else 0.0,
        )

    rho = np.asarray(dump["rho"], dtype=float)
    dx = float(mesh_spacing[0]) if mesh_spacing is not None else 1.0
    dy = float(mesh_spacing[1]) if mesh_spacing is not None else 1.0
    dz = float(mesh_spacing[2]) if mesh_spacing is not None else 1.0
    total_charge = float(rho.sum() * dx * dy * dz)
    return field_origin, mesh_spacing, total_charge


def frame_number_for_dump(dump: dict, fallback_index: int) -> int:
    metadata = dump["metadata"]
    if "global_step" in metadata:
        return int(parse_numeric(metadata["global_step"]))
    return fallback_index


def render_dump_figure(
    dump: dict,
    overlay: dict[str, object],
    use_grid: bool,
    render_context: dict[str, object] | None = None,
    force_agg: bool = False,
):
    import matplotlib

    if force_agg:
        matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    rho = dump["rho"]
    xs = dump["is"] if use_grid else dump["xs"]
    ys = dump["js"] if use_grid else dump["ys"]
    zs = dump["ks"] if use_grid else (np.asarray(dump["zs"], dtype=float) +
                                      longitudinal_display_delta(dump, use_grid))

    (
        interaction_point,
        ip_element_z_range,
        collision_window_z_range,
        primary_center,
        mirrored_center,
    ) = transformed_overlay(dump, overlay, use_grid)

    if render_context is None:
        cmap, vmin, vmax = color_limits(rho)
        if vmin == vmax:
            vmax = vmin + 1.0
        render_context = compute_global_limits([dump], [overlay], use_grid)
        render_context["cmap"] = cmap
        render_context["vmin"] = vmin
        render_context["vmax"] = vmax
    else:
        cmap = str(render_context["cmap"])
        vmin = float(render_context["vmin"])
        vmax = float(render_context["vmax"])

    ix = len(xs) // 2
    iy = len(ys) // 2
    iz = len(zs) // 2

    x_label = "i" if use_grid else "x [m]"
    y_label = "j" if use_grid else "y [m]"
    z_label = "k" if use_grid else "s [m]"
    label_x_offset = 1.5 if use_grid else 0.003
    raw_primary_center = overlay["primary_center"]
    raw_mirrored_center = overlay["mirrored_center"]
    b1_text = format_bunch_position(dump, raw_primary_center)
    b2_text = format_bunch_position(dump, raw_mirrored_center)
    b1_label = "b1" if b1_text is None else f"b1 {b1_text}"
    b2_label = "b2" if b2_text is None else f"b2 {b2_text}"

    fig = plt.figure(figsize=(13, 4.5), constrained_layout=True)
    gs = fig.add_gridspec(1, 4, width_ratios=[1.0, 1.0, 1.0, 0.06])

    ax_xy = fig.add_subplot(gs[0, 0])
    extent_xy = axis_extent(xs, ys)
    im = ax_xy.imshow(
        rho[:, :, iz].T,
        origin="lower",
        extent=extent_xy,
        aspect="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    ax_xy.set_xlabel(x_label)
    ax_xy.set_ylabel(y_label)
    draw_domain_box(ax_xy, extent_xy)
    if interaction_point is not None:
        add_ip_overlay(ax_xy, float(interaction_point[0]), float(interaction_point[1]), x_label, y_label)
    if primary_center is not None:
        add_bunch_center(
            ax_xy,
            float(primary_center[0]),
            float(primary_center[1]),
            marker="o",
            edgecolor="deepskyblue",
            facecolor="white",
            label="b1",
            x_offset=label_x_offset,
        )
    if mirrored_center is not None:
        add_bunch_center(
            ax_xy,
            float(mirrored_center[0]),
            float(mirrored_center[1]),
            marker="x",
            edgecolor="magenta",
            facecolor="none",
            label="b2",
            x_offset=label_x_offset,
        )
    ax_xy.set_xlim(*render_context["xlim"])
    ax_xy.set_ylim(*render_context["ylim"])

    ax_xz = fig.add_subplot(gs[0, 1])
    extent_xz = axis_extent(xs, zs)
    ax_xz.imshow(
        rho[:, iy, :].T,
        origin="lower",
        extent=extent_xz,
        aspect="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    ax_xz.set_xlabel(x_label)
    ax_xz.set_ylabel(z_label)
    draw_domain_box(ax_xz, extent_xz)
    if interaction_point is not None:
        add_ip_overlay(ax_xz, float(interaction_point[0]), float(interaction_point[2]), x_label, z_label)
    if primary_center is not None:
        add_bunch_center(
            ax_xz,
            float(primary_center[0]),
            float(primary_center[2]),
            marker="o",
            edgecolor="deepskyblue",
            facecolor="white",
            label=b1_label,
            x_offset=label_x_offset,
        )
    if mirrored_center is not None:
        add_bunch_center(
            ax_xz,
            float(mirrored_center[0]),
            float(mirrored_center[2]),
            marker="x",
            edgecolor="magenta",
            facecolor="none",
            label=b2_label,
            x_offset=label_x_offset,
        )
    if ip_element_z_range is not None:
        draw_longitudinal_band(
            ax_xz,
            extent_xz,
            ip_element_z_range,
            edgecolor="darkorange",
            linestyle="--",
            label="BeamBeam",
        )
    if collision_window_z_range is not None and not same_z_range(
        collision_window_z_range, ip_element_z_range
    ):
        draw_longitudinal_band(
            ax_xz,
            extent_xz,
            collision_window_z_range,
            edgecolor="crimson",
            linestyle="-.",
            label="collision window",
        )
    ax_xz.set_xlim(*render_context["xlim"])
    ax_xz.set_ylim(*render_context["zlim"])

    ax_yz = fig.add_subplot(gs[0, 2])
    extent_yz = axis_extent(ys, zs)
    ax_yz.imshow(
        rho[ix, :, :].T,
        origin="lower",
        extent=extent_yz,
        aspect="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    ax_yz.set_xlabel(y_label)
    ax_yz.set_ylabel(z_label)
    draw_domain_box(ax_yz, extent_yz)
    if interaction_point is not None:
        add_ip_overlay(ax_yz, float(interaction_point[1]), float(interaction_point[2]), y_label, z_label)
    if primary_center is not None:
        add_bunch_center(
            ax_yz,
            float(primary_center[1]),
            float(primary_center[2]),
            marker="o",
            edgecolor="deepskyblue",
            facecolor="white",
            label=b1_label,
            x_offset=label_x_offset,
        )
    if mirrored_center is not None:
        add_bunch_center(
            ax_yz,
            float(mirrored_center[1]),
            float(mirrored_center[2]),
            marker="x",
            edgecolor="magenta",
            facecolor="none",
            label=b2_label,
            x_offset=label_x_offset,
        )
    if ip_element_z_range is not None:
        draw_longitudinal_band(
            ax_yz,
            extent_yz,
            ip_element_z_range,
            edgecolor="darkorange",
            linestyle="--",
            label="BeamBeam",
        )
    if collision_window_z_range is not None and not same_z_range(
        collision_window_z_range, ip_element_z_range
    ):
        draw_longitudinal_band(
            ax_yz,
            extent_yz,
            collision_window_z_range,
            edgecolor="crimson",
            linestyle="-.",
            label="collision window",
        )
    ax_yz.set_xlim(*render_context["ylim"])
    ax_yz.set_ylim(*render_context["zlim"])

    cax = fig.add_subplot(gs[0, 3])
    fig.colorbar(im, cax=cax, label="rho [C/m^3]")
    fig.suptitle(build_title(dump, overlay))
    return fig


def write_single_output(
    dump: dict,
    overlay: dict[str, object],
    output: Path,
    use_grid: bool,
    render_context: dict[str, object] | None = None,
) -> None:
    fig = render_dump_figure(dump, overlay, use_grid, render_context=render_context, force_agg=True)
    fig.savefig(output, dpi=150)
    fig.clf()
    print(f"wrote {output}")


def render_gif(
    dumps: list[dict],
    overlays: list[dict[str, object]],
    output: Path,
    use_grid: bool,
    duration_ms: int,
    save_frames: bool,
) -> None:
    frames: list[Image.Image] = []
    render_context = compute_global_limits(dumps, overlays, use_grid)

    for frame_index, (dump, overlay) in enumerate(zip(dumps, overlays), start=1):
        fig = render_dump_figure(dump, overlay, use_grid, render_context=render_context, force_agg=True)
        if save_frames:
            frame_number = frame_number_for_dump(dump, frame_index)
            frame_path = output.with_name(f"frame-{frame_number:05d}.png")
            fig.savefig(frame_path, dpi=150)
        buffer = io.BytesIO()
        fig.savefig(buffer, format="png", dpi=150)
        fig.clf()
        buffer.seek(0)
        frames.append(Image.open(buffer).convert("RGB"))

    frames[0].save(
        output,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
    )
    print(f"wrote {output}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Visualize rho dumps together with the BeamBeam element and collision-window geometry."
    )
    parser.add_argument("filename", help="Path to a scalar RHO dump.")
    parser.add_argument("--output", help="Write a PNG or GIF instead of showing the figure.")
    parser.add_argument("--start", type=int, help="Start step for a movie or single-step render.")
    parser.add_argument("--end", type=int, help="End step for a movie or single-step render.")
    parser.add_argument("--duration", type=int, default=150, help="GIF frame duration in milliseconds.")
    parser.add_argument(
        "--save-frames",
        action="store_true",
        help="When rendering a GIF, also save each frame as frame-0000x.png next to the output file.",
    )
    parser.set_defaults(grid=True)
    parser.add_argument(
        "--grid",
        dest="grid",
        action="store_true",
        help="Plot against grid indices i,j,k. This is now the default.",
    )
    parser.add_argument(
        "--physical",
        dest="grid",
        action="store_false",
        help="Plot against physical coordinates x,y,s in meters instead of i,j,k.",
    )
    parser.add_argument(
        "--no-geometry",
        action="store_true",
        help="Render rho only, without BeamBeam/window overlays or bunch-center markers.",
    )
    parser.add_argument(
        "--all-states",
        action="store_true",
        help="When rendering a step range, include all snapshot kinds instead of filtering to the seed file state.",
    )
    parser.add_argument(
        "--input",
        help="Input file used to infer the BeamBeam element and collision-window geometry.",
    )
    parser.add_argument(
        "--bb-s",
        dest="bb_s",
        type=float,
        help="Absolute s-position of the BeamBeam interaction point.",
    )
    parser.add_argument(
        "--bb-begin-s",
        dest="bb_begin_s",
        type=float,
        help="Absolute s-position where the BeamBeam element begins.",
    )
    parser.add_argument(
        "--bb-end-s",
        dest="bb_end_s",
        type=float,
        help="Absolute s-position where the BeamBeam element ends.",
    )
    parser.add_argument("--window-begin-s", type=float, help="Absolute s-position where the collision window begins.")
    parser.add_argument("--window-end-s", type=float, help="Absolute s-position where the collision window ends.")
    args = parser.parse_args()

    seed = Path(args.filename)
    geometry = {} if args.no_geometry else build_input_geometry(args, seed)

    if (args.start is None) != (args.end is None):
        raise SystemExit("--start and --end must be given together")

    if args.start is None and STEP_FILE_RE.match(seed.name) is not None:
        inside_ip_paths = infer_inside_ip_paths(seed, geometry)
        if inside_ip_paths:
            dumps = [parse_dump(path) for path in inside_ip_paths]
            overlays = [empty_overlay() if args.no_geometry else derive_overlay_geometry(dump, geometry)
                        for dump in dumps]
            dumps, overlays = align_longitudinal_display_frame(dumps, overlays, args.grid)
            overlays = freeze_active_window_geometry(dumps, overlays)
            output = Path(args.output) if args.output else default_output_path(
                inside_ip_paths[0],
                int(STEP_FILE_RE.match(inside_ip_paths[0].name).group(2)),
                int(STEP_FILE_RE.match(inside_ip_paths[-1].name).group(2)),
            )
            if len(dumps) == 1:
                render_context = compute_global_limits(dumps, overlays, args.grid)
                write_single_output(dumps[0], overlays[0], output, args.grid, render_context)
            else:
                render_gif(dumps, overlays, output, args.grid, args.duration, args.save_frames)
            return

    if args.start is None:
        dump = parse_dump(seed)
        overlay = empty_overlay() if args.no_geometry else derive_overlay_geometry(dump, geometry)
        if args.output:
            write_single_output(dump, overlay, Path(args.output), args.grid)
        else:
            fig = render_dump_figure(dump, overlay, args.grid)
            import matplotlib.pyplot as plt

            plt.show()
        return

    paths = resolve_step_paths(seed, args.start, args.end)
    paths = filter_paths_by_unique_global_step(paths)
    dumps = [parse_dump(path) for path in paths]
    overlays = [empty_overlay() if args.no_geometry else derive_overlay_geometry(dump, geometry)
                for dump in dumps]
    dumps, overlays = align_longitudinal_display_frame(dumps, overlays, args.grid)
    overlays = freeze_active_window_geometry(dumps, overlays)

    output = Path(args.output) if args.output else default_output_path(seed, args.start, args.end)
    output = normalize_output_path(output, args.start, args.end)
    if args.start == args.end:
        render_context = compute_global_limits(dumps, overlays, args.grid)
        write_single_output(dumps[0], overlays[0], output, args.grid, render_context)
        return

    render_gif(dumps, overlays, output, args.grid, args.duration, args.save_frames)


if __name__ == "__main__":
    main()
