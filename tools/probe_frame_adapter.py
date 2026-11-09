#!/usr/bin/env python3
"""Transform OPALX PROBE records into an explicitly defined Cartesian frame.

Positions use r_target = basis @ (r_lab - origin); momenta use p_target =
basis @ p_lab. Basis rows are the target unit axes expressed in laboratory
coordinates. IDs, reference-turn labels, bunch numbers and time are preserved.
This standalone analysis utility is not a native unit-test dependency.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re

import h5py
import numpy as np
import pandas as pd

COLUMNS = ("x", "y", "z", "px", "py", "pz", "id", "turn", "bunchNumber", "time")
INTEGER_COLUMNS = ("id", "turn", "bunchNumber")
FLOAT_COLUMNS = tuple(name for name in COLUMNS if name not in INTEGER_COLUMNS)
HEADER = "# x (m), y (m), z (m), px ( ), py ( ), pz ( ), id, turn ( ), bunchNumber ( ), time (s)"
_INT64_MAX = np.iinfo(np.int64).max


def _validated(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate the ten-column record without rounding integer identities."""
    if tuple(frame.columns) != COLUMNS:
        raise ValueError(f"Expected exactly these probe columns: {', '.join(COLUMNS)}")
    result = frame.copy()
    for name in INTEGER_COLUMNS:
        # Do not pass IDs through float64, which loses integers above 2**53.
        values = []
        for value in result[name]:
            if isinstance(value, (int, np.integer)):
                integer = int(value)
            elif isinstance(value, str) and re.fullmatch(r"[+]?\d+", value):
                integer = int(value)
            else:
                raise ValueError(f"{name} must contain nonnegative integer values")
            if integer < 0 or integer > _INT64_MAX:
                raise ValueError(f"{name} must fit a nonnegative signed 64-bit integer")
            values.append(integer)
        result[name] = pd.Series(values, index=result.index, dtype="int64")
    for name in FLOAT_COLUMNS:
        # Parse decimal tokens directly as float64. pandas.to_numeric can round
        # a 17-digit token by one ULP, defeating an ASCII round trip.
        result[name] = result[name].astype("float64")
        if not np.isfinite(result[name].to_numpy()).all():
            raise ValueError(f"Nonfinite probe {name}")
    return result


def read_probe(path: str | Path) -> pd.DataFrame:
    """Read one station's ten-column ASCII or HDF5 PROBE output.

    HDF5 Step# groups are flush batches, not reference turns. Concatenate all
    batches, then sort stably by bunch/particle/time. No contiguous-ID, unique
    (ID,turn), or one-crossing-per-turn assumption is made. HDF5 q/m datasets
    are outside the ten-column ASCII format and are not exported.
    """
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"Probe file does not exist: {source}")
    if h5py.is_hdf5(source):
        batches = []
        with h5py.File(source, "r") as file:
            groups = sorted((key for key in file if re.fullmatch(r"Step#\d+", key)),
                            key=lambda key: int(key.split("#")[1]))
            if not groups:
                raise ValueError(f"No Step# probe batches in {source}")
            for key in groups:
                group = file[key]
                if not isinstance(group, h5py.Group):
                    raise ValueError(f"{source}:{key} is not an HDF5 group")
                missing = set(COLUMNS) - set(group)
                if missing:
                    raise ValueError(f"{source}:{key} missing probe datasets: {sorted(missing)}")
                arrays = {name: np.asarray(group[name]) for name in COLUMNS}
                if any(array.ndim != 1 for array in arrays.values()):
                    raise ValueError(f"{source}:{key} probe datasets must be one-dimensional")
                if len({len(array) for array in arrays.values()}) != 1:
                    raise ValueError(f"{source}:{key} has unequal probe dataset lengths")
                for name in INTEGER_COLUMNS:
                    if arrays[name].dtype.kind not in "iu":
                        raise ValueError(f"{source}:{key}/{name} must be an integer dataset")
                batches.append(_validated(pd.DataFrame(arrays, columns=COLUMNS)))
        frame = pd.concat(batches, ignore_index=True)
    else:
        # Read tokens as strings first so malformed widths and large IDs remain
        # visible; providing column names before parsing can hide extra columns.
        try:
            frame = pd.read_csv(source, sep=r"\s+", comment="#", header=None,
                                dtype=str, keep_default_na=False)
        except pd.errors.EmptyDataError:
            frame = pd.DataFrame(columns=COLUMNS)
        if len(frame.columns) != len(COLUMNS):
            raise ValueError(f"{source} must contain exactly ten probe columns")
        frame.columns = COLUMNS
        frame = _validated(frame)
    return frame.sort_values(["bunchNumber", "id", "time"], kind="stable").reset_index(drop=True)


def transform_frame(native: pd.DataFrame, origin_m=(0.0, 0.0, 0.0),
                    basis=None) -> pd.DataFrame:
    """Express laboratory positions and momenta in a caller-defined frame.

    origin_m is the target origin in laboratory metres. The three rows of basis
    are orthonormal target axes in laboratory coordinates (identity by default).
    Either handedness is supported; momentum is a polar vector and uses the same
    basis as position, without translation. Scaling, shear and non-Cartesian
    coordinates are rejected. The input and all nonspatial columns are unchanged.
    """
    origin = np.asarray(origin_m, dtype="float64")
    axes = np.eye(3) if basis is None else np.asarray(basis, dtype="float64")
    if origin.shape != (3,) or not np.isfinite(origin).all():
        raise ValueError("origin_m must contain three finite laboratory coordinates in metres")
    if axes.shape != (3, 3) or not np.isfinite(axes).all():
        raise ValueError("basis must be a finite 3-by-3 matrix of target axes")
    # Permit decimal roundoff in supplied unit axes, but reject scale or shear.
    # This is an input-geometry check, not a tracking or comparison tolerance.
    if not np.allclose(axes @ axes.T, np.eye(3), rtol=0.0, atol=1e-12):
        raise ValueError("basis rows must be orthonormal (no scaling or shear)")
    native = _validated(native)
    result = native.copy()
    try:
        with np.errstate(over="raise", invalid="raise"):
            result[["x", "y", "z"]] = (
                native[["x", "y", "z"]].to_numpy() - origin) @ axes.T
            result[["px", "py", "pz"]] = native[["px", "py", "pz"]].to_numpy() @ axes.T
    except FloatingPointError as error:
        raise ValueError("Frame conversion exceeds finite floating-point range") from error
    return _validated(result)


def write_probe_ascii(frame: pd.DataFrame, path: str | Path) -> None:
    """Write one header and ten standard probe columns to a new file.

    Exclusive creation prevents accidental replacement of native or existing
    analysis data. Seventeen-digit floats retain the precision available in the
    source; native ASCII input itself normally has only six significant digits.
    """
    frame = _validated(frame)
    with Path(path).open("x") as output:
        output.write(HEADER + "\n")
        for row in frame.itertuples(index=False, name=None):
            output.write(" ".join(str(int(value)) if name in INTEGER_COLUMNS
                                  else format(float(value), ".17g")
                                  for name, value in zip(COLUMNS, row)) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="one native laboratory PROBE .loss or .h5 file")
    parser.add_argument("output", type=Path, help="new ten-column .loss file in the target frame")
    parser.add_argument("--origin-m", type=float, nargs=3, default=(0.0, 0.0, 0.0),
                        metavar=("X", "Y", "Z"),
                        help="target origin in laboratory metres; default 0 0 0")
    parser.add_argument("--basis", type=float, nargs=9,
                        metavar=("XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ"),
                        help="target unit axes as laboratory components, row-major; default identity")
    args = parser.parse_args()
    try:
        if args.input.resolve() == args.output.resolve():
            raise ValueError("Input and output must be different files")
        basis = None if args.basis is None else np.asarray(args.basis).reshape(3, 3)
        native = read_probe(args.input)
        converted = transform_frame(native, args.origin_m, basis)
        write_probe_ascii(converted, args.output)
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.exit(1, f"probe_frame_adapter: {error}\n")
    print(f"Converted {len(converted)} hits to {args.output}; IDs and reference-turn labels preserved")


if __name__ == "__main__":
    main()
