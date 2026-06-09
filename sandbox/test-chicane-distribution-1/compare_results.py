#!/usr/bin/env python3
"""Compare OPALX output with the manufactured chicane distribution target."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import h5py
import numpy as np

from manufactured_chicane import P0_BETA_GAMMA, R56_REF


def step_index(name: str) -> int:
    match = re.fullmatch(r"Step#(\d+)", name)
    if not match:
        return -1
    return int(match.group(1))


def h5_attr_scalar(group: h5py.Group, name: str) -> float:
    if name not in group.attrs:
        return float("nan")
    value = np.asarray(group.attrs[name], dtype=float)
    if value.size == 0:
        return float("nan")
    return float(value.reshape(-1)[0])


def read_distribution(path: Path) -> np.ndarray:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 3:
        raise ValueError(f"distribution file is too short: {path}")
    n_particles = int(lines[0].strip())
    data = np.loadtxt(path, skiprows=2)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape != (n_particles, 6):
        raise ValueError(
            f"expected {n_particles} rows and 6 columns in {path}, got {data.shape}"
        )
    ids = np.arange(n_particles, dtype=float).reshape(-1, 1)
    delta = (data[:, 5] / P0_BETA_GAMMA - 1.0).reshape(-1, 1)
    return np.hstack((ids, data, delta))


def read_exit_monitor(path: Path, n_particles: int) -> tuple[str, float, dict[str, np.ndarray]]:
    required = {"id", "x", "px", "y", "py", "z", "pz"}
    candidates = []
    with h5py.File(path, "r") as h5:
        for name in sorted(h5.keys(), key=step_index):
            group = h5[name]
            if not required.issubset(group.keys()):
                continue
            ids = np.asarray(group["id"], dtype=int)
            spos = h5_attr_scalar(group, "SPOS")
            candidates.append((len(ids), spos, step_index(name), name))

        if not candidates:
            raise RuntimeError(f"no usable Step# groups found in {path}")

        candidates.sort(key=lambda row: (row[0], np.nan_to_num(row[1], nan=-1.0), row[2]))
        _, spos, _, selected = candidates[-1]
        group = h5[selected]
        data = {key: np.asarray(group[key]) for key in required}

    ids = np.asarray(data["id"], dtype=int)
    if len(ids) != n_particles:
        raise RuntimeError(
            f"selected {selected} contains {len(ids)} particles, expected {n_particles}"
        )
    if len(np.unique(ids)) != n_particles or ids.min() < 0 or ids.max() >= n_particles:
        raise RuntimeError(f"selected {selected} contains invalid particle ids")

    order = np.argsort(ids)
    ordered = {key: np.asarray(value)[order] for key, value in data.items()}
    return selected, spos, ordered


def rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.asarray(values, dtype=float) ** 2)))


def write_plots(out_dir: Path, initial: np.ndarray, monitor: dict[str, np.ndarray]) -> None:
    from PIL import Image, ImageDraw, ImageFont

    out_dir.mkdir(parents=True, exist_ok=True)
    particle_id = initial[:, 0]
    delta = initial[:, 7]
    z_initial = initial[:, 5]
    z_expected = z_initial - R56_REF * delta
    z_exit = np.asarray(monitor["z"], dtype=float)
    slope, intercept = np.polyfit(delta, z_exit - z_initial, 1)

    def fmt_tick(value: float) -> str:
        if value == 0.0:
            return "0"
        if abs(value) < 1.0e-2 or abs(value) >= 1.0e3:
            return f"{value:.2e}"
        return f"{value:.4g}"

    def padded_range(values: np.ndarray) -> tuple[float, float]:
        values = np.asarray(values, dtype=float)
        lo = float(np.min(values))
        hi = float(np.max(values))
        if lo == hi:
            pad = max(abs(lo), 1.0) * 0.05
        else:
            pad = 0.06 * (hi - lo)
        return lo - pad, hi + pad

    def draw_plot(
            path: Path,
            x: np.ndarray,
            y: np.ndarray,
            xlabel: str,
            ylabel: str,
            title: str,
            lines: list[tuple[np.ndarray, np.ndarray, str]] | None = None) -> None:
        width, height = 900, 600
        left, right, top, bottom = 120, 35, 75, 80
        plot_w = width - left - right
        plot_h = height - top - bottom
        font = ImageFont.load_default()
        img = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(img)

        x_all = np.asarray(x, dtype=float)
        y_all = np.asarray(y, dtype=float)
        if lines:
            for lx, ly, _ in lines:
                x_all = np.concatenate((x_all, np.asarray(lx, dtype=float)))
                y_all = np.concatenate((y_all, np.asarray(ly, dtype=float)))
        x_min, x_max = padded_range(x_all)
        y_min, y_max = padded_range(y_all)

        def sx(value: float) -> float:
            return left + (value - x_min) / (x_max - x_min) * plot_w

        def sy(value: float) -> float:
            return top + plot_h - (value - y_min) / (y_max - y_min) * plot_h

        draw.rectangle((left, top, left + plot_w, top + plot_h), outline=(35, 35, 35), width=1)
        for i in range(6):
            tx = x_min + (x_max - x_min) * i / 5.0
            px = sx(tx)
            draw.line((px, top, px, top + plot_h), fill=(225, 225, 225), width=1)
            label = fmt_tick(tx)
            draw.text((px - 25, top + plot_h + 10), label, fill=(35, 35, 35), font=font)

            ty = y_min + (y_max - y_min) * i / 5.0
            py = sy(ty)
            draw.line((left, py, left + plot_w, py), fill=(225, 225, 225), width=1)
            draw.text((8, py - 6), fmt_tick(ty), fill=(35, 35, 35), font=font)

        if lines:
            for lx, ly, color in lines:
                points = [(sx(float(a)), sy(float(b))) for a, b in zip(lx, ly)]
                if len(points) >= 2:
                    draw.line(points, fill=color, width=2)

        for px, py in zip(x, y):
            cx = sx(float(px))
            cy = sy(float(py))
            draw.ellipse((cx - 3, cy - 3, cx + 3, cy + 3), fill=(31, 119, 180))

        draw.text((left, 28), title, fill=(10, 10, 10), font=font)
        draw.text((left + plot_w / 2 - 70, height - 35), xlabel, fill=(10, 10, 10), font=font)
        label_bbox = draw.textbbox((0, 0), ylabel, font=font)
        label_w = label_bbox[2] - label_bbox[0] + 6
        label_h = label_bbox[3] - label_bbox[1] + 6
        label_img = Image.new("RGBA", (label_w, label_h), (255, 255, 255, 0))
        label_draw = ImageDraw.Draw(label_img)
        label_draw.text((3, 3), ylabel, fill=(10, 10, 10), font=font)
        label_img = label_img.rotate(90, expand=True)
        img.paste(label_img, (28, top + plot_h // 2 - label_img.height // 2), label_img)
        img.save(path)

    lo = min(float(np.min(z_expected)), float(np.min(z_exit)))
    hi = max(float(np.max(z_expected)), float(np.max(z_exit)))
    draw_plot(
            out_dir / "manufactured_z_compare.png",
            z_expected,
            z_exit,
            "manufactured z_exit [m]",
            "OPALX monitor z [m]",
            "Manufactured vs OPALX longitudinal output",
            [(np.array([lo, hi]), np.array([lo, hi]), "black")])

    xline = np.array([float(np.min(delta)), float(np.max(delta))])
    draw_plot(
            out_dir / "manufactured_r56.png",
            delta,
            z_exit - z_initial,
            "initial delta",
            "z_exit - z_initial [m]",
            "R56 manufactured slope",
            [(xline, intercept - R56_REF * xline, "black")])

    residual = z_exit - z_expected - intercept
    draw_plot(
            out_dir / "manufactured_z_residual.png",
            particle_id,
            residual,
            "particle id",
            "intercept-corrected z residual [m]",
            "Manufactured longitudinal residual after monitor offset removal",
            [(np.array([float(np.min(particle_id)), float(np.max(particle_id))]), np.array([0.0, 0.0]), "black")])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--distribution",
        type=Path,
        default=Path(__file__).resolve().parent / "test-chicane-distribution-1_distribution.txt",
    )
    parser.add_argument(
        "--monitor",
        type=Path,
        default=Path(__file__).resolve().parent / "test-chicane-distribution-1_exit.h5",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "comparison",
    )
    parser.add_argument("--r56-abs-tol", type=float, default=1.0e-3)
    parser.add_argument("--z-rms-tol", type=float, default=5.0e-4)
    parser.add_argument("--z-max-tol", type=float, default=1.0e-3)
    parser.add_argument("--y-rms-tol", type=float, default=1.0e-8)
    parser.add_argument("--plots", action="store_true")
    args = parser.parse_args()

    initial = read_distribution(args.distribution)
    n_particles = len(initial)
    selected_step, spos, monitor = read_exit_monitor(args.monitor, n_particles)

    delta = initial[:, 7]
    z_initial = initial[:, 5]
    y_initial = initial[:, 3]
    z_expected = z_initial - R56_REF * delta
    z_exit = np.asarray(monitor["z"], dtype=float)
    y_exit = np.asarray(monitor["y"], dtype=float)
    pz_exit = np.asarray(monitor["pz"], dtype=float)

    slope, intercept = np.polyfit(delta, z_exit - z_initial, 1)
    r56_transport = -float(slope)
    z_residual = z_exit - z_expected
    z_residual_centered = z_residual - intercept
    y_residual = y_exit - y_initial
    pz_residual = pz_exit - initial[:, 6]

    checks = {
        "all_particles_at_exit": True,
        "r56_within_tolerance": abs(r56_transport - R56_REF) <= args.r56_abs_tol,
        "z_rms_within_tolerance": rms(z_residual_centered) <= args.z_rms_tol,
        "z_max_within_tolerance": float(np.max(np.abs(z_residual_centered)))
        <= args.z_max_tol,
        "y_rms_within_tolerance": rms(y_residual) <= args.y_rms_tol,
    }
    passed = all(checks.values())

    summary = {
        "passed": passed,
        "selected_step": selected_step,
        "spos_m": float(spos),
        "n_particles": int(n_particles),
        "r56_ref_m": R56_REF,
        "r56_transport_m": r56_transport,
        "r56_difference_m": r56_transport - R56_REF,
        "z_fit_intercept_m": float(intercept),
        "z_residual_rms_m": rms(z_residual),
        "z_residual_max_abs_m": float(np.max(np.abs(z_residual))),
        "z_residual_centered_rms_m": rms(z_residual_centered),
        "z_residual_centered_max_abs_m": float(np.max(np.abs(z_residual_centered))),
        "y_residual_rms_m": rms(y_residual),
        "y_residual_max_abs_m": float(np.max(np.abs(y_residual))),
        "pz_residual_rms": rms(pz_residual),
        "pz_residual_max_abs": float(np.max(np.abs(pz_residual))),
        "x_exit_mean_m": float(np.mean(monitor["x"])),
        "x_exit_rms_m": float(np.std(monitor["x"])),
        "y_exit_mean_m": float(np.mean(y_exit)),
        "y_exit_rms_m": float(np.std(y_exit)),
        "z_exit_mean_m": float(np.mean(z_exit)),
        "z_exit_rms_m": float(np.std(z_exit)),
        "checks": checks,
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "comparison_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    markdown = "\n".join(
        [
            "# test-chicane-distribution-1 comparison",
            "",
            f"Selected monitor step: `{selected_step}`",
            f"Monitor `SPOS`: `{spos:.12e}` m",
            f"Particles at exit: `{n_particles}`",
            "",
            "## Manufactured longitudinal check",
            "",
            f"- Reference `R56`: `{R56_REF:.12e}` m",
            f"- Fitted transport `R56`: `{r56_transport:.12e}` m",
            f"- Difference: `{r56_transport - R56_REF:.12e}` m",
            f"- Fitted monitor offset: `{float(intercept):.12e}` m",
            f"- Raw `z` residual RMS: `{rms(z_residual):.12e}` m",
            f"- Raw `z` residual max abs: `{float(np.max(np.abs(z_residual))):.12e}` m",
            f"- Intercept-corrected `z` residual RMS: `{rms(z_residual_centered):.12e}` m",
            f"- Intercept-corrected `z` residual max abs: `{float(np.max(np.abs(z_residual_centered))):.12e}` m",
            "",
            "## Transverse closure diagnostics",
            "",
            f"- `y` residual RMS: `{rms(y_residual):.12e}` m",
            f"- `y` residual max abs: `{float(np.max(np.abs(y_residual))):.12e}` m",
            f"- final `x` mean: `{float(np.mean(monitor['x'])):.12e}` m",
            f"- final `x` RMS: `{float(np.std(monitor['x'])):.12e}` m",
            "",
            "## Checks",
            "",
            *[f"- `{name}`: `{value}`" for name, value in checks.items()],
            "",
            f"Overall: `{'PASS' if passed else 'FAIL'}`",
            "",
        ]
    )
    (args.out_dir / "comparison_summary.md").write_text(markdown, encoding="utf-8")

    if args.plots:
        write_plots(args.out_dir, initial, monitor)

    print(markdown)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
