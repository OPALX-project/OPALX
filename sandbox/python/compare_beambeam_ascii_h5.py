#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


HEADER_RE = re.compile(r"origin=\s*(\([^)]+\))\s+h=\s*(\([^)]+\))")


def parse_header_value(text: str):
    text = text.strip()
    if text.startswith("(") and text.endswith(")"):
        return [float(item.strip()) for item in text[1:-1].split(",") if item.strip()]
    try:
        return float(text)
    except ValueError:
        return text


def read_ascii_dump(path: Path):
    metadata = {}
    column_count = None

    with path.open("r", encoding="utf-8") as stream:
        for raw_line in stream:
            line = raw_line.rstrip("\n")
            if not line.startswith("#"):
                column_count = len(line.split())
                break

            header = line[1:].strip()
            match = HEADER_RE.match(header)
            if match:
                metadata["mesh_origin"] = parse_header_value(match.group(1))
                metadata["mesh_spacing"] = parse_header_value(match.group(2))
                continue
            if "=" in header:
                key, value = header.split("=", 1)
                metadata[key.strip()] = parse_header_value(value.strip())

    if column_count is None:
        raise SystemExit(f"no field rows found in {path}")

    if column_count == 7:
        columns = ["i", "j", "k", "x", "y", "z", "value"]
    elif column_count == 9:
        columns = ["i", "j", "k", "x", "y", "z", "x_value", "y_value", "z_value"]
    else:
        raise SystemExit(f"unsupported field row width {column_count} in {path}")

    frame = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=columns,
        engine="python",
    )

    i_values = sorted(frame["i"].unique())
    j_values = sorted(frame["j"].unique())
    k_values = sorted(frame["k"].unique())
    shape = (len(i_values), len(j_values), len(k_values))
    i_map = {value: index for index, value in enumerate(i_values)}
    j_map = {value: index for index, value in enumerate(j_values)}
    k_map = {value: index for index, value in enumerate(k_values)}

    if column_count == 7:
        field = np.zeros((shape[2], shape[1], shape[0]), dtype=float)
        for row in frame.itertuples(index=False):
            field[k_map[row.k], j_map[row.j], i_map[row.i]] = row.value
    else:
        field = {
            "Ex": np.zeros((shape[2], shape[1], shape[0]), dtype=float),
            "Ey": np.zeros((shape[2], shape[1], shape[0]), dtype=float),
            "Ez": np.zeros((shape[2], shape[1], shape[0]), dtype=float),
        }
        for row in frame.itertuples(index=False):
            iz = k_map[row.k]
            iy = j_map[row.j]
            ix = i_map[row.i]
            field["Ex"][iz, iy, ix] = row.x_value
            field["Ey"][iz, iy, ix] = row.y_value
            field["Ez"][iz, iy, ix] = row.z_value

    return {"metadata": metadata, "field": field}


def decode_h5_value(value):
    if hasattr(value, "item"):
        try:
            return decode_h5_value(value.item())
        except Exception:
            pass
    if hasattr(value, "shape"):
        if value.shape == ():
            return decode_h5_value(value[()])
        return [decode_h5_value(item) for item in value]
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value


def step_sort_key(name: str) -> tuple[int, str]:
    if name.startswith("Step#"):
        suffix = name.split("#", 1)[1]
        try:
            return (int(suffix), name)
        except ValueError:
            pass
    return (10**18, name)


def select_h5_step(h5file, global_step: int, snapshot_kind: str | None):
    for step_name in sorted(h5file.keys(), key=step_sort_key):
        step = h5file[step_name]
        step_value = int(decode_h5_value(step.attrs.get("global_step")))
        if step_value != global_step:
            continue
        if snapshot_kind is not None:
            current_snapshot = str(decode_h5_value(step.attrs.get("snapshot_kind", "")))
            if current_snapshot != snapshot_kind:
                continue
        return step
    raise SystemExit(f"global_step={global_step} not found in {h5file.filename}")


def read_h5_scalar(path: Path, global_step: int, snapshot_kind: str | None):
    with h5py.File(path, "r") as h5file:
        step = select_h5_step(h5file, global_step, snapshot_kind)
        block = step["Block"]
        field_name = sorted(block.keys())[0]
        field_group = block[field_name]
        dataset_name = sorted(field_group.keys(), key=step_sort_key)[0]
        return np.asarray(field_group[dataset_name][()], dtype=float)


def read_h5_vector(path: Path, global_step: int, snapshot_kind: str | None):
    with h5py.File(path, "r") as h5file:
        step = select_h5_step(h5file, global_step, snapshot_kind)
        block = step["Block"]
        field_name = sorted(block.keys())[0]
        field_group = block[field_name]
        dataset_names = sorted(field_group.keys(), key=step_sort_key)
        if len(dataset_names) < 3:
            raise SystemExit(f"vector field in {path} has fewer than three components")
        return {
            "Ex": np.asarray(field_group[dataset_names[0]][()], dtype=float),
            "Ey": np.asarray(field_group[dataset_names[1]][()], dtype=float),
            "Ez": np.asarray(field_group[dataset_names[2]][()], dtype=float),
        }


def metrics(reference, candidate):
    if reference.shape != candidate.shape:
        raise SystemExit(f"shape mismatch: ASCII {reference.shape}, HDF5 {candidate.shape}")
    diff = candidate - reference
    max_abs = float(np.max(np.abs(diff)))
    l2 = float(np.sqrt(np.sum(diff * diff)))
    ref_l2 = float(np.sqrt(np.sum(reference * reference)))
    rel_l2 = l2 / ref_l2 if ref_l2 > 0.0 else (0.0 if l2 == 0.0 else math.inf)
    return max_abs, rel_l2


def print_metrics(name: str, reference, candidate) -> bool:
    max_abs, rel_l2 = metrics(reference, candidate)
    print(f"{name:>4s}: max_abs={max_abs:.6e}  rel_l2={rel_l2:.6e}")
    return rel_l2 <= 1.0e-9 or np.allclose(candidate, reference, rtol=1.0e-9, atol=1.0e-30)


def build_parser():
    parser = argparse.ArgumentParser(description="Compare OPALX BeamBeam ASCII and HDF5 field dumps.")
    parser.add_argument("--rho-ascii", type=Path, required=True)
    parser.add_argument("--phi-ascii", type=Path, required=True)
    parser.add_argument("--ef-ascii", type=Path, required=True)
    parser.add_argument("--rho-h5", type=Path, required=True)
    parser.add_argument("--phi-h5", type=Path, required=True)
    parser.add_argument("--ef-h5", type=Path, required=True)
    parser.add_argument("--global-step", type=int)
    return parser


def main() -> int:
    args = build_parser().parse_args()

    rho_ascii = read_ascii_dump(args.rho_ascii)
    phi_ascii = read_ascii_dump(args.phi_ascii)
    ef_ascii = read_ascii_dump(args.ef_ascii)

    global_step = args.global_step
    if global_step is None:
        global_step = int(rho_ascii["metadata"]["global_step"])
    snapshot_kind = str(rho_ascii["metadata"].get("snapshot_kind", ""))

    rho_h5 = read_h5_scalar(args.rho_h5, global_step, snapshot_kind)
    phi_h5 = read_h5_scalar(args.phi_h5, global_step, snapshot_kind)
    ef_h5 = read_h5_vector(args.ef_h5, global_step, snapshot_kind)

    print(f"global_step={global_step} snapshot_kind={snapshot_kind}")
    ok = True
    ok &= print_metrics("rho", rho_ascii["field"], rho_h5)
    ok &= print_metrics("phi", phi_ascii["field"], phi_h5)
    for component in ("Ex", "Ey", "Ez"):
        ok &= print_metrics(component, ef_ascii["field"][component], ef_h5[component])

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
