#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

import numpy as np


EPS0 = 8.8541878128e-12
C0 = 299792458.0


def gaussian_density(x, y, z, center, sigmas, charge):
    cx, cy, cz = center
    sx, sy, sz = sigmas
    norm = charge / (((2.0 * math.pi) ** 1.5) * sx * sy * sz)
    exponent = (
        -0.5 * ((x - cx) / sx) ** 2
        -0.5 * ((y - cy) / sy) ** 2
        -0.5 * ((z - cz) / sz) ** 2
    )
    return norm * np.exp(exponent)


def bunch_centers(separation, beta, time_s):
    displacement = 0.5 * separation - beta * C0 * time_s
    return (0.0, 0.0, displacement), (0.0, 0.0, -displacement)


def make_grid(sigmas, nsigma, nz_factor):
    sx, sy, sz = sigmas
    nx = max(8, int(math.ceil(nsigma * 2.0 * nz_factor * sx / sx)))
    ny = max(8, int(math.ceil(nsigma * 2.0 * nz_factor * sy / sy)))
    nz = max(8, int(math.ceil(nsigma * 2.0 * nz_factor * sz / min(sx, sy, sz))))

    xs = np.linspace(-nsigma * sx, nsigma * sx, nx)
    ys = np.linspace(-nsigma * sy, nsigma * sy, ny)
    zs = np.linspace(-nsigma * nz_factor * sz, nsigma * nz_factor * sz, nz)
    dx = xs[1] - xs[0] if nx > 1 else 1.0
    dy = ys[1] - ys[0] if ny > 1 else 1.0
    dz = zs[1] - zs[0] if nz > 1 else 1.0
    xv, yv, zv = np.meshgrid(xs, ys, zs, indexing="ij")
    return xv, yv, zv, dx * dy * dz


def total_density(xv, yv, zv, centers, sigmas, charge):
    rho = np.zeros_like(xv)
    for center in centers:
        rho += gaussian_density(xv, yv, zv, center, sigmas, charge)
    return rho


def electric_field_at_point(point, xv, yv, zv, rho, cell_volume):
    px, py, pz = point
    rx = px - xv
    ry = py - yv
    rz = pz - zv
    r2 = rx * rx + ry * ry + rz * rz

    # Avoid the singular self-cell contribution if the probe lies on a grid point.
    r2 = np.where(r2 == 0.0, np.inf, r2)
    inv_r3 = 1.0 / (r2 * np.sqrt(r2))
    prefactor = cell_volume / (4.0 * math.pi * EPS0)

    ex = prefactor * np.sum(rho * rx * inv_r3)
    ey = prefactor * np.sum(rho * ry * inv_r3)
    ez = prefactor * np.sum(rho * rz * inv_r3)
    return np.array([ex, ey, ez])


def axis_scan(axis, scan_values, xv, yv, zv, rho, cell_volume):
    fields = []
    for value in scan_values:
        if axis == "x":
            point = (value, 0.0, 0.0)
        elif axis == "y":
            point = (0.0, value, 0.0)
        else:
            point = (0.0, 0.0, value)
        fields.append(electric_field_at_point(point, xv, yv, zv, rho, cell_volume))
    return np.asarray(fields)


def write_scan(output_path, axis, scan_values, fields):
    header = (
        "# position[m] Ex[V/m] Ey[V/m] Ez[V/m]\n"
        f"# axis={axis}\n"
    )
    with open(output_path, "w", encoding="utf-8") as stream:
        stream.write(header)
        for value, field in zip(scan_values, fields):
            stream.write(
                f"{value:.12e} {field[0]:.12e} {field[1]:.12e} {field[2]:.12e}\n"
            )


def plot_scan(plot_path, axis, scan_values, fields):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return False

    cm_to_inch = 1.0 / 2.54
    fig, ax = plt.subplots(figsize=(12.0 * cm_to_inch, 9.0 * cm_to_inch), dpi=200)
    ax.plot(scan_values, fields[:, 0], label="Ex", linewidth=2.0)
    ax.plot(scan_values, fields[:, 1], label="Ey", linewidth=2.0)
    ax.plot(scan_values, fields[:, 2], label="Ez", linewidth=2.0)
    ax.set_xlabel(f"{axis} [m]")
    ax.set_ylabel("E [V/m]")
    ax.set_title("Field of Two Counterpropagating Gaussian Bunches")
    ax.grid(True, linestyle="--", linewidth=0.7, alpha=0.5)
    ax.legend()
    fig.tight_layout()
    fig.savefig(plot_path, bbox_inches="tight")
    plt.close(fig)
    return True


def auto_time_window_ps(separation, beta):
    crossing_time_s = 0.5 * separation / (max(beta, 1.0e-12) * C0)
    crossing_time_ps = crossing_time_s * 1.0e12
    return -1.5 * crossing_time_ps, 1.5 * crossing_time_ps


def plot_scan_frame(plot_path, axis, scan_values, fields, centers, time_ps):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return False

    axis_index = {"x": 0, "y": 1, "z": 2}[axis]
    component_labels = ["Ex", "Ey", "Ez"]
    component_colors = ["tab:blue", "tab:orange", "tab:green"]

    cm_to_inch = 1.0 / 2.54
    fig, ax = plt.subplots(figsize=(12.0 * cm_to_inch, 9.0 * cm_to_inch), dpi=160)
    for component_index, (label, color) in enumerate(zip(component_labels, component_colors)):
        ax.plot(
            scan_values,
            fields[:, component_index],
            label=label,
            linewidth=2.0,
            color=color,
        )

    center_positions = [center[axis_index] for center in centers]
    for index, center_position in enumerate(center_positions, start=1):
        if scan_values[0] <= center_position <= scan_values[-1]:
            ax.axvline(
                center_position,
                linestyle="--",
                linewidth=1.5,
                color="black",
                alpha=0.6,
                label=f"bunch {index} center" if index == 1 else None,
            )

    ax.set_xlabel(f"{axis} [m]")
    ax.set_ylabel("E [V/m]")
    ax.set_title(f"Two Counterpropagating Gaussian Bunches, t = {time_ps:.3f} ps")
    ax.grid(True, linestyle="--", linewidth=0.7, alpha=0.5)
    ax.legend()
    fig.tight_layout()
    fig.savefig(plot_path, bbox_inches="tight")
    plt.close(fig)
    return True


def make_movie(frame_paths, movie_output, duration_ms):
    try:
        from PIL import Image
    except ModuleNotFoundError:
        return False

    if not frame_paths:
        return False

    images = [Image.open(path) for path in frame_paths]
    images[0].save(
        movie_output,
        save_all=True,
        append_images=images[1:],
        duration=duration_ms,
        loop=0,
    )
    for image in images:
        image.close()
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Direct 3D Coulomb integration for two head-on Gaussian bunches."
    )
    parser.add_argument("--sigma-x", type=float, default=1.0e-3)
    parser.add_argument("--sigma-y", type=float, default=1.0e-3)
    parser.add_argument("--sigma-z", type=float, default=3.0e-3)
    parser.add_argument("--charge", type=float, default=1.0e-9,
                        help="Charge per bunch [C]")
    parser.add_argument("--separation", type=float, default=6.0e-3,
                        help="Initial bunch-center separation [m]")
    parser.add_argument("--beta", type=float, default=1.0)
    parser.add_argument("--time-ps", type=float, default=0.0,
                        help="Time offset [ps]")
    parser.add_argument("--nsigma", type=float, default=4.0,
                        help="Transverse box half-width in sigma")
    parser.add_argument("--nz-factor", type=float, default=2.0,
                        help="Longitudinal box factor relative to nsigma*sigma_z")
    parser.add_argument("--scan-axis", choices=["x", "y", "z"], default="x")
    parser.add_argument("--scan-min", type=float, default=-5.0e-3)
    parser.add_argument("--scan-max", type=float, default=5.0e-3)
    parser.add_argument("--scan-points", type=int, default=81)
    parser.add_argument("--output", type=Path, default=Path("ip-1/data/two_gaussian_field_scan.dat"))
    parser.add_argument("--plot-output", type=Path, default=Path("ip-1/data/two_gaussian_field_scan.png"))
    parser.add_argument("--make-movie", action="store_true",
                        help="Generate a time sequence and assemble a GIF")
    parser.add_argument("--time-start-ps", type=float, default=None)
    parser.add_argument("--time-end-ps", type=float, default=None)
    parser.add_argument("--time-steps", type=int, default=41)
    parser.add_argument("--frames-dir", type=Path, default=Path("ip-1/data/two_gaussian_frames"))
    parser.add_argument("--movie-output", type=Path, default=Path("ip-1/data/two_gaussian_field_scan.gif"))
    parser.add_argument("--movie-duration-ms", type=int, default=120)
    args = parser.parse_args()

    sigmas = (args.sigma_x, args.sigma_y, args.sigma_z)
    time_s = args.time_ps * 1.0e-12
    centers = bunch_centers(args.separation, args.beta, time_s)

    xv, yv, zv, cell_volume = make_grid(sigmas, args.nsigma, args.nz_factor)
    rho = total_density(xv, yv, zv, centers, sigmas, args.charge)

    scan_values = np.linspace(args.scan_min, args.scan_max, args.scan_points)
    fields = axis_scan(args.scan_axis, scan_values, xv, yv, zv, rho, cell_volume)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_scan(args.output, args.scan_axis, scan_values, fields)
    args.plot_output.parent.mkdir(parents=True, exist_ok=True)
    plotted = plot_scan(args.plot_output, args.scan_axis, scan_values, fields)

    print(f"centers={centers}")
    print(f"grid_shape={xv.shape} cell_volume={cell_volume:.6e}")
    print(f"wrote {args.output}")
    if plotted:
        print(f"wrote {args.plot_output}")
    else:
        print("matplotlib not available, skipped plot output")

    if not args.make_movie:
        return

    time_start_ps = args.time_start_ps
    time_end_ps = args.time_end_ps
    if time_start_ps is None or time_end_ps is None:
        auto_start_ps, auto_end_ps = auto_time_window_ps(args.separation, args.beta)
        if time_start_ps is None:
            time_start_ps = auto_start_ps
        if time_end_ps is None:
            time_end_ps = auto_end_ps

    args.frames_dir.mkdir(parents=True, exist_ok=True)
    frame_paths = []
    plotted_frames = True
    for frame_index, frame_time_ps in enumerate(
        np.linspace(time_start_ps, time_end_ps, args.time_steps)
    ):
        frame_centers = bunch_centers(args.separation, args.beta, frame_time_ps * 1.0e-12)
        frame_rho = total_density(xv, yv, zv, frame_centers, sigmas, args.charge)
        frame_fields = axis_scan(args.scan_axis, scan_values, xv, yv, zv, frame_rho, cell_volume)
        frame_path = args.frames_dir / f"frame-{frame_index:04d}.png"
        if not plot_scan_frame(frame_path, args.scan_axis, scan_values, frame_fields, frame_centers, frame_time_ps):
            plotted_frames = False
            break
        frame_paths.append(frame_path)

    if not plotted_frames:
        print("matplotlib not available, skipped movie frame generation")
        return

    args.movie_output.parent.mkdir(parents=True, exist_ok=True)
    if make_movie(frame_paths, args.movie_output, args.movie_duration_ms):
        print(f"wrote {args.movie_output}")
    else:
        print("Pillow not available, skipped GIF assembly")


if __name__ == "__main__":
    main()
