#!/usr/bin/env python3
"""Run DT and zero-transverse diagnostics for the chicane distribution case."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np

from manufactured_chicane import P0_BETA_GAMMA


ROOT = Path(__file__).resolve().parent
OPALX_EXE = Path("/Users/adelmann/git/opalx/build/src/opalx")
OPAL_EXE = Path("/Users/adelmann/OPAL-2022.1/bin/opal")
OPAL_PREFIX = Path("/Users/adelmann/OPAL-2022.1")

BODY_TEMPLATE = ROOT / "test-chicane-distribution-1.in"
ELEMEDGE_TEMPLATE = ROOT / "test-chicane-distribution-2.in"
OPAL_TEMPLATE = ROOT / "test-chicane-distribution-1_opal.in"
COMPARE_SCRIPT = ROOT / "compare_opal_opalx.py"


def dt_tag(dt: float) -> str:
    return f"dt_{dt:.1e}".replace("+", "").replace("-", "m").replace(".", "p")


def patch_common(text: str, dt: float, n_particles: int) -> str:
    text = text.replace("DT = 1.0e-11,", f"DT = {dt:.16e},")
    text = text.replace("REAL TIME_STEP = 1.0e-11;", f"REAL TIME_STEP = {dt:.16e};")
    text = text.replace("NPARTDIST = 125", f"NPARTDIST = {n_particles}")
    text = text.replace("NALLOC = 125", f"NALLOC = {n_particles}")
    text = text.replace("REAL NUMBER_OF_PARTICLES = 125;", f"REAL NUMBER_OF_PARTICLES = {n_particles};")
    return text


def write_inputs(run_dir: Path, dt: float, n_particles: int) -> tuple[Path, Path, Path]:
    body_input = run_dir / BODY_TEMPLATE.name
    elemedge_input = run_dir / ELEMEDGE_TEMPLATE.name
    opal_input = run_dir / OPAL_TEMPLATE.name

    body_input.write_text(patch_common(BODY_TEMPLATE.read_text(encoding="utf-8"), dt, n_particles), encoding="utf-8")
    elemedge_input.write_text(
        patch_common(ELEMEDGE_TEMPLATE.read_text(encoding="utf-8"), dt, n_particles),
        encoding="utf-8",
    )
    opal_input.write_text(patch_common(OPAL_TEMPLATE.read_text(encoding="utf-8"), dt, n_particles), encoding="utf-8")
    return body_input, elemedge_input, opal_input


def finite_distribution() -> np.ndarray:
    data = np.loadtxt(ROOT / "test-chicane-distribution-1_distribution.txt", skiprows=2)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def zero_transverse_distribution() -> np.ndarray:
    deltas = np.array([-2.0e-3, -1.0e-3, 0.0, 1.0e-3, 2.0e-3])
    z_values = np.array([-1.0e-3, -5.0e-4, 0.0, 5.0e-4, 1.0e-3])
    rows = []
    for delta in deltas:
        for z in z_values:
            rows.append([0.0, 0.0, 0.0, 0.0, z, P0_BETA_GAMMA * (1.0 + delta)])
    return np.asarray(rows, dtype=float)


def write_distribution_files(run_dir: Path, data: np.ndarray) -> int:
    opalx_path = run_dir / "test-chicane-distribution-1_distribution.txt"
    opal_path = run_dir / "test-chicane-distribution-1_opal_distribution.txt"
    with opalx_path.open("w", encoding="utf-8") as stream:
        stream.write(f"{len(data)}\n")
        stream.write("x px y py z pz\n")
        np.savetxt(stream, data, fmt="%.16e")
    with opal_path.open("w", encoding="utf-8") as stream:
        stream.write(f"{len(data)}\n")
        np.savetxt(stream, data, fmt="%.16e")
    return int(len(data))


def run_command(command: list[str], cwd: Path, env: dict[str, str] | None = None) -> None:
    print(f"$ {' '.join(command)}", flush=True)
    result = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (cwd / "run.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        print(result.stdout, flush=True)
        raise RuntimeError(f"command failed with exit code {result.returncode}: {' '.join(command)}")


def run_case(case_dir: Path, case_name: str, dt: float, force: bool) -> dict[str, Path]:
    run_dir = case_dir / dt_tag(dt)
    body_dir = run_dir / "opalx_body"
    elemedge_dir = run_dir / "opalx_elemedge"
    opal_dir = run_dir / "opal_2022"
    compare_body_dir = run_dir / "compare_body"
    compare_elemedge_dir = run_dir / "compare_elemedge"
    for directory in (body_dir, elemedge_dir, opal_dir, compare_body_dir, compare_elemedge_dir):
        directory.mkdir(parents=True, exist_ok=True)

    data = finite_distribution() if case_name == "finite" else zero_transverse_distribution()
    n_particles = write_distribution_files(body_dir, data)
    write_distribution_files(elemedge_dir, data)
    write_distribution_files(opal_dir, data)
    body_input, _, _ = write_inputs(body_dir, dt, n_particles)
    _, elemedge_input, _ = write_inputs(elemedge_dir, dt, n_particles)
    _, _, opal_input = write_inputs(opal_dir, dt, n_particles)

    if force or not (body_dir / "test-chicane-distribution-1.stat").exists():
        run_command([str(OPALX_EXE), body_input.name], body_dir)
    if force or not (elemedge_dir / "test-chicane-distribution-2.stat").exists():
        run_command([str(OPALX_EXE), elemedge_input.name], elemedge_dir)
    if force or not (opal_dir / "test-chicane-distribution-1_opal.stat").exists():
        env = dict(**__import__("os").environ, OPAL_PREFIX=str(OPAL_PREFIX))
        run_command([str(OPAL_EXE), opal_input.name], opal_dir, env=env)

    run_command(
        [
            str(Path("/Users/adelmann/.venv-h6/bin/python")),
            str(COMPARE_SCRIPT),
            "--opalx-stat",
            str(body_dir / "test-chicane-distribution-1.stat"),
            "--opal-stat",
            str(opal_dir / "test-chicane-distribution-1_opal.stat"),
            "--out-dir",
            str(compare_body_dir),
            "--plots",
        ],
        run_dir,
    )
    run_command(
        [
            str(Path("/Users/adelmann/.venv-h6/bin/python")),
            str(COMPARE_SCRIPT),
            "--opalx-stat",
            str(elemedge_dir / "test-chicane-distribution-2.stat"),
            "--opal-stat",
            str(opal_dir / "test-chicane-distribution-1_opal.stat"),
            "--out-dir",
            str(compare_elemedge_dir),
            "--plots",
        ],
        run_dir,
    )
    return {
        "run": run_dir,
        "compare_body": compare_body_dir,
        "compare_elemedge": compare_elemedge_dir,
    }


def read_summary_rows(compare_dir: Path, case_name: str, dt: float, comparison: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with (compare_dir / "stat_summary.csv").open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row["quantity"] in {"rms_x", "rms_y", "rms_s", "mean_x", "mean_s", "ref_x", "ref_z", "By_ref"}:
                row = dict(row)
                row["case"] = case_name
                row["dt_s"] = f"{dt:.16e}"
                row["comparison"] = comparison
                rows.append(row)
    return rows


def read_interval_rows(compare_dir: Path, case_name: str, dt: float, comparison: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with (compare_dir / "by_field_intervals.csv").open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            row = dict(row)
            row["case"] = case_name
            row["dt_s"] = f"{dt:.16e}"
            row["comparison"] = comparison
            rows.append(row)
    return rows


def write_combined_summary(out_dir: Path, metric_rows: list[dict[str, str]], interval_rows: list[dict[str, str]]) -> None:
    metrics_csv = out_dir / "dt_scan_metrics.csv"
    intervals_csv = out_dir / "dt_scan_by_intervals.csv"
    if metric_rows:
        fieldnames = ["case", "dt_s", "comparison", *[key for key in metric_rows[0] if key not in {"case", "dt_s", "comparison"}]]
        with metrics_csv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(metric_rows)
    if interval_rows:
        fieldnames = ["case", "dt_s", "comparison", *[key for key in interval_rows[0] if key not in {"case", "dt_s", "comparison"}]]
        with intervals_csv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(interval_rows)

    def selected(quantity: str, column: str = "rms_diff") -> list[dict[str, str]]:
        return [row for row in metric_rows if row["quantity"] == quantity and column in row]

    lines = [
        "# DT scan summary",
        "",
        f"- Metrics CSV: `{metrics_csv.name}`",
        f"- By interval CSV: `{intervals_csv.name}`",
        "",
        "## Selected RMS differences",
        "",
        "| case | comparison | DT [s] | By_ref rms [T] | rms_s rms [m] | rms_x rms [m] |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    by = {(row["case"], row["comparison"], row["dt_s"]): row for row in selected("By_ref")}
    rs = {(row["case"], row["comparison"], row["dt_s"]): row for row in selected("rms_s")}
    rx = {(row["case"], row["comparison"], row["dt_s"]): row for row in selected("rms_x")}
    for key in sorted(by, key=lambda item: (item[0], item[1], float(item[2]))):
        lines.append(
            "| {} | {} | {:.1e} | {:.6e} | {:.6e} | {:.6e} |".format(
                key[0],
                key[1],
                float(key[2]),
                float(by[key]["rms_diff"]),
                float(rs.get(key, {"rms_diff": "nan"})["rms_diff"]),
                float(rx.get(key, {"rms_diff": "nan"})["rms_diff"]),
            )
        )
    lines.append("")
    (out_dir / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=ROOT / "comparison" / "dt_scan")
    parser.add_argument("--dt", type=float, nargs="+", default=[1.0e-11, 5.0e-12, 2.0e-12, 1.0e-12])
    parser.add_argument("--case", choices=["finite", "zero"], nargs="+", default=["finite", "zero"])
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    metric_rows: list[dict[str, str]] = []
    interval_rows: list[dict[str, str]] = []
    for case_name in args.case:
        case_dir = args.out_dir / case_name
        for dt in args.dt:
            print(f"\n=== {case_name} DT={dt:.3e} ===", flush=True)
            paths = run_case(case_dir, case_name, dt, args.force)
            metric_rows.extend(read_summary_rows(paths["compare_body"], case_name, dt, "body"))
            metric_rows.extend(read_summary_rows(paths["compare_elemedge"], case_name, dt, "elemedge"))
            interval_rows.extend(read_interval_rows(paths["compare_body"], case_name, dt, "body"))
            interval_rows.extend(read_interval_rows(paths["compare_elemedge"], case_name, dt, "elemedge"))

    write_combined_summary(args.out_dir, metric_rows, interval_rows)
    print(f"\nWrote {args.out_dir / 'summary.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
