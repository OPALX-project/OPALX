#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
import io
import json
import mimetypes
import os
import sys
import tempfile
import threading
import traceback
import webbrowser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import parse_qs, urlencode, urlparse


SCRIPT_DIR = Path(__file__).resolve().parent
CHECK_COLLWIN_PATH = SCRIPT_DIR / "checkCollWin.py"
MANUFACTURED_PATH = SCRIPT_DIR / "beam-beam-manufactured-solution.py"
BOOSTED_WITNESS_PATH = SCRIPT_DIR / "boosted_gaussian_witness.py"
READ_H5_PATH = SCRIPT_DIR / "read_beambeam_h5.py"
GUI_STATE_PATH = SCRIPT_DIR / ".beambeam_analysis_state.json"


DEFAULT_GUI_STATE = {
    "collwin": {
        "filename": "data/sandbox/BeamBeam-2-RHO_scalar-collwin_vis-000010.dat",
        "input": "sandbox/BeamBeam-2.in",
        "output": "",
        "start": "",
        "end": "",
        "grid": True,
        "no_geometry": False,
    },
    "analytic": {
        "sigma": "1e-3",
        "charge": "-1e-15",
        "half_separation": "2e-3",
        "nx": "64",
        "ny": "64",
        "nz": "64",
        "output": "",
    },
    "compare": {
        "h5_path": "data/sandbox/beambeam-1-beambeam_diagnostics.h5",
        "step": "24",
        "start": "",
        "end": "",
        "sigma": "1e-3",
        "state": "beambeam tracking",
        "component": "rho+phi",
        "output": "",
    },
}


def load_module(module_name: str, path: Path):
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load module {module_name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def merge_state(defaults: dict, loaded: dict) -> dict:
    merged = {}
    for key, value in defaults.items():
        if isinstance(value, dict):
            merged[key] = merge_state(value, loaded.get(key, {}) if isinstance(loaded.get(key, {}), dict) else {})
        else:
            merged[key] = loaded.get(key, value)
    return merged


def load_gui_state() -> dict:
    if not GUI_STATE_PATH.exists():
        return json.loads(json.dumps(DEFAULT_GUI_STATE))
    try:
        loaded = json.loads(GUI_STATE_PATH.read_text())
    except Exception:
        return json.loads(json.dumps(DEFAULT_GUI_STATE))
    return merge_state(DEFAULT_GUI_STATE, loaded if isinstance(loaded, dict) else {})


def save_gui_state(state: dict) -> None:
    GUI_STATE_PATH.write_text(json.dumps(state, indent=2, sort_keys=True))


def load_check_collwin():
    return load_module("checkCollWin_module", CHECK_COLLWIN_PATH)


def load_manufactured_solution():
    return load_module("beam_beam_manufactured_solution_module", MANUFACTURED_PATH)


def load_boosted_witness():
    return load_module("boosted_gaussian_witness_module", BOOSTED_WITNESS_PATH)


def load_read_beambeam_h5():
    return load_module("read_beambeam_h5_module", READ_H5_PATH)


def load_pillow_image():
    try:
        from PIL import Image  # type: ignore
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "Pillow is not installed. Activate the local environment first:\n"
            "  source .venv-h5/bin/activate"
        ) from exc
    return Image


def configure_matplotlib_backend(force_agg: bool) -> None:
    if not force_agg:
        return
    import matplotlib

    matplotlib.use("Agg")


def normalize_compare_state(state: str | None) -> str | None:
    if state is None:
        return None
    normalized = state.strip()
    if not normalized:
        return None
    mapping = {
        "all states": None,
        "all": None,
        "beambeam tracking": "active_beambeam",
        "normal tracking": "normal_tracking",
        "active_beambeam": "active_beambeam",
        "normal_tracking": "normal_tracking",
        "before_beambeam_mesh_enlarge": "before_beambeam_mesh_enlarge",
    }
    return mapping.get(normalized, normalized)


def display_compare_state(raw_state: str) -> str:
    mapping = {
        "active_beambeam": "beambeam tracking",
        "normal_tracking": "normal tracking",
        "before_beambeam_mesh_enlarge": "normal tracking",
    }
    return mapping.get(raw_state, raw_state)


def normalize_compare_component(component: str | None) -> tuple[str, ...]:
    if component is None:
        return ("rho", "phi")
    normalized = component.strip().lower()
    mapping = {
        "rho": ("rho",),
        "phi": ("phi",),
        "ez": ("Ez",),
        "ex-x-axis": ("Ex-x-axis",),
        "ey-y-axis": ("Ey-y-axis",),
        "rho-z-axis": ("rho-z-axis",),
        "phi-z-axis": ("phi-z-axis",),
        "ez-z-axis": ("Ez-z-axis",),
        "rho+phi": ("rho", "phi"),
        "rho+phi+ez": ("rho", "phi", "Ez"),
        "rho,phi": ("rho", "phi"),
        "rho,phi,ez": ("rho", "phi", "Ez"),
        "all": ("rho", "phi", "Ez"),
        "both": ("rho", "phi"),
    }
    return mapping.get(normalized, ("rho", "phi"))


def build_overview_table_text(path: Path, start: int | None = None, end: int | None = None) -> str:
    reader = load_read_beambeam_h5()
    h5py = reader.load_h5py()
    with h5py.File(path, "r") as h5file:
        step_names, rows = reader.collect_step_rows(h5file)
        step_names, rows = reader.filter_steps(step_names, rows, start, end)

    columns = [
        "step",
        "state",
        "t",
        "s",
        "<s>",
        "IP",
        "extend (lab)",
        "Q",
        "shape",
        "origin_z (beam)",
        "spacing",
    ]
    widths = {
        col: max(len(col), *(len(str(row[col])) for row in rows)) if rows else len(col)
        for col in columns
    }

    lines = [f"file: {path}", f"steps: {len(step_names)}", ""]
    header = "  ".join(f"{col:<{widths[col]}}" for col in columns)
    separator = "  ".join("-" * widths[col] for col in columns)
    lines.append(header)
    lines.append(separator)
    for row in rows:
        lines.append("  ".join(f"{str(row[col]):<{widths[col]}}" for col in columns))
    return "\n".join(lines)


def run_collision_window(
    filename: Path,
    *,
    output: Path | None = None,
    start: int | None = None,
    end: int | None = None,
    grid: bool = True,
    no_geometry: bool = False,
    input_file: Path | None = None,
    duration: int = 150,
    save_frames: bool = False,
    force_agg: bool = False,
) -> None:
    collwin = load_check_collwin()
    geometry_args = SimpleNamespace(
        input=str(input_file) if input_file is not None else None,
        bb_s=None,
        bb_begin_s=None,
        bb_end_s=None,
        window_begin_s=None,
        window_end_s=None,
    )

    geometry = {} if no_geometry else collwin.build_input_geometry(geometry_args, filename)

    if start is None and end is None and collwin.STEP_FILE_RE.match(filename.name) is not None:
        inside_ip_paths = collwin.infer_inside_ip_paths(filename, geometry)
        if inside_ip_paths:
            dumps = [collwin.parse_dump(path) for path in inside_ip_paths]
            overlays = [
                collwin.empty_overlay() if no_geometry else collwin.derive_overlay_geometry(dump, geometry)
                for dump in dumps
            ]
            dumps, overlays = collwin.align_longitudinal_display_frame(dumps, overlays, grid)
            overlays = collwin.freeze_active_window_geometry(dumps, overlays)
            resolved_output = output or collwin.default_output_path(
                inside_ip_paths[0],
                int(collwin.STEP_FILE_RE.match(inside_ip_paths[0].name).group(2)),
                int(collwin.STEP_FILE_RE.match(inside_ip_paths[-1].name).group(2)),
            )
            if len(dumps) == 1:
                render_context = collwin.compute_global_limits(dumps, overlays, grid)
                collwin.write_single_output(dumps[0], overlays[0], resolved_output, grid, render_context)
            else:
                collwin.render_gif(dumps, overlays, resolved_output, grid, duration, save_frames)
            return

    if start is None and end is None:
        dump = collwin.parse_dump(filename)
        overlay = collwin.empty_overlay() if no_geometry else collwin.derive_overlay_geometry(dump, geometry)
        if output is not None:
            collwin.write_single_output(dump, overlay, output, grid)
            return
        collwin.render_dump_figure(dump, overlay, grid, force_agg=force_agg)
        import matplotlib.pyplot as plt
        plt.show()
        return

    if (start is None) != (end is None):
        raise SystemExit("--start and --end must be given together")

    paths = collwin.resolve_step_paths(filename, start, end)
    dumps = [collwin.parse_dump(path) for path in paths]
    overlays = [
        collwin.empty_overlay() if no_geometry else collwin.derive_overlay_geometry(dump, geometry)
        for dump in dumps
    ]
    dumps, overlays = collwin.align_longitudinal_display_frame(dumps, overlays, grid)
    overlays = collwin.freeze_active_window_geometry(dumps, overlays)

    resolved_output = output or collwin.default_output_path(filename, start, end)
    resolved_output = collwin.normalize_output_path(resolved_output, start, end)
    if start == end:
        render_context = collwin.compute_global_limits(dumps, overlays, grid)
        collwin.write_single_output(dumps[0], overlays[0], resolved_output, grid, render_context)
        return

    collwin.render_gif(dumps, overlays, resolved_output, grid, duration, save_frames)


def run_manufactured_analytic(
    *,
    sigma: float,
    charge: float,
    half_separation: float,
    nx: int = 64,
    ny: int = 64,
    nz: int = 64,
    x_span: float | None = None,
    y_span: float | None = None,
    z_span: float | None = None,
    output: Path | None = None,
    force_agg: bool = False,
) -> None:
    configure_matplotlib_backend(force_agg or output is not None)
    manufactured = load_manufactured_solution()
    np = manufactured.load_numpy()
    plt = manufactured.load_matplotlib()
    args = SimpleNamespace(
        sigma=sigma,
        charge=charge,
        half_separation=half_separation,
        nx=nx,
        ny=ny,
        nz=nz,
        x_span=x_span,
        y_span=y_span,
        z_span=z_span,
    )
    origin, spacing, shape = manufactured.default_grid(args)
    x_grid, y_grid, z_grid = manufactured.build_grid(np, origin, spacing, shape)
    centers = [
        (0.0, 0.0, -half_separation),
        (0.0, 0.0, +half_separation),
    ]
    analytic = manufactured.gaussian_pair_fields(np, x_grid, y_grid, z_grid, sigma, charge, centers)
    title = f"Manufactured solution | $Q$={charge:.2e} C, $\\sigma$={sigma:.2e} m"
    manufactured.plot_analytic_solution(plt, analytic, origin, spacing, title, output)


def run_pair_in_fields(
    *,
    sigma_m: float,
    charge_C: float,
    source_kinetic_MeV: float,
    half_separation_m: float,
    pair_energy_keV: float,
    pair_duration_ps: float,
    pair_dt_ps: float,
    pair_fine_dt_ps: float | None,
    pair_fine_duration_ps: float,
    pair_output_dt_ps: float,
    pair_direction: tuple[float, float, float],
    pair_origin_m: tuple[float, float, float],
    field_span_mm: float,
    field_grid: int,
    beambeam_radius_m: float,
    beambeam_length_m: float,
    output_prefix: Path,
    write_pair_latex_table: Path,
) -> int:
    boosted = load_boosted_witness()
    args = SimpleNamespace(
        sigma_m=sigma_m,
        charge_C=charge_C,
        source_kinetic_MeV=source_kinetic_MeV,
        half_separation_m=half_separation_m,
        pair_energy_keV=pair_energy_keV,
        pair_duration_ps=pair_duration_ps,
        pair_dt_ps=pair_dt_ps,
        pair_fine_dt_ps=pair_fine_dt_ps,
        pair_fine_duration_ps=pair_fine_duration_ps,
        pair_output_dt_ps=pair_output_dt_ps,
        pair_direction=pair_direction,
        pair_origin_m=pair_origin_m,
        field_span_mm=field_span_mm,
        field_grid=field_grid,
        beambeam_radius_m=beambeam_radius_m,
        beambeam_length_m=beambeam_length_m,
        output_prefix=output_prefix,
        write_pair_latex_table=write_pair_latex_table,
    )
    return boosted.run_pair_demo(args)


def run_manufactured_compare(
    *,
    h5_path: Path,
    step: int,
    sigma: float,
    state: str | None = None,
    comparison: str | None = None,
    output: Path | None = None,
    force_agg: bool = False,
) -> None:
    configure_matplotlib_backend(force_agg or output is not None)
    manufactured = load_manufactured_solution()
    np = manufactured.load_numpy()
    plt = manufactured.load_matplotlib()

    setup, analytic, origin, spacing, ip_beam_z, title = prepare_manufactured_compare(
        manufactured,
        np,
        h5_path,
        step,
        sigma,
        state,
    )
    render_manufactured_compare_output(
        manufactured,
        np,
        plt,
        analytic,
        setup["opalx"],
        origin,
        spacing,
        ip_beam_z,
        title,
        output,
        comparison,
    )
    components = normalize_compare_component(comparison)
    if components in (("rho", "phi"), ("rho", "phi", "Ez")):
        manufactured.plot_rho_z_axis_diagnostic(
            np,
            plt,
            analytic,
            setup["opalx"],
            origin,
            spacing,
            title,
            manufactured.derived_output_path(output, "-rho-z-axis"),
        )
        manufactured.plot_phi_z_axis_diagnostic(
            np,
            plt,
            analytic,
            setup["opalx"],
            origin,
            spacing,
            title,
            manufactured.derived_output_path(output, "-phi-z-axis"),
        )
        manufactured.plot_ez_z_axis_diagnostic(
            np,
            plt,
            analytic,
            setup["opalx"],
            origin,
            spacing,
            title,
            manufactured.derived_output_path(output, "-ez-z-axis"),
        )
        ip_beam_z = setup["interaction_point_s"] - setup["path_length_s"]
        manufactured.plot_ex_x_axis_diagnostic(
            np,
            plt,
            analytic,
            setup["opalx"],
            origin,
            spacing,
            title,
            manufactured.derived_output_path(output, "-ex-x-axis"),
            ip_beam_z,
        )
        manufactured.plot_ey_y_axis_diagnostic(
            np,
            plt,
            analytic,
            setup["opalx"],
            origin,
            spacing,
            title,
            manufactured.derived_output_path(output, "-ey-y-axis"),
            ip_beam_z,
        )


def prepare_manufactured_compare(manufactured, np, h5_path: Path, step: int, sigma: float, state: str | None):
    setup = manufactured.manufactured_setup_from_h5(np, h5_path, step, normalize_compare_state(state))
    origin = setup["origin"]
    spacing = setup["spacing"]
    shape = setup["shape"]
    charge = setup["charge"]
    centers = setup["centers"]

    x_grid, y_grid, z_grid = manufactured.build_grid(np, origin, spacing, shape)
    analytic = manufactured.gaussian_pair_fields(np, x_grid, y_grid, z_grid, sigma, charge, centers)
    ip_beam_z = setup["interaction_point_s"] - setup["path_length_s"]

    title = (
        f"Manufactured solution vs OPALX | step {step} | {display_compare_state(setup['snapshot_kind'])} | "
        f"$Q$={charge:.2e} C, $\\sigma$={sigma:.2e} m"
    )
    return setup, analytic, origin, spacing, ip_beam_z, title


def render_manufactured_compare_output(
    manufactured,
    np,
    plt,
    analytic,
    opalx,
    origin,
    spacing,
    ip_beam_z: float,
    title: str,
    output: Path | None,
    comparison: str | None,
    axis_limits=None,
    scale_overrides=None,
) -> None:
    components = normalize_compare_component(comparison)
    if components == ("rho-z-axis",):
        manufactured.plot_rho_z_axis_diagnostic(
            np,
            plt,
            analytic,
            opalx,
            origin,
            spacing,
            title,
            output,
            axis_limits=axis_limits,
        )
        return
    if components == ("phi-z-axis",):
        manufactured.plot_phi_z_axis_diagnostic(
            np,
            plt,
            analytic,
            opalx,
            origin,
            spacing,
            title,
            output,
            axis_limits=axis_limits,
        )
        return
    if components == ("Ex-x-axis",):
        manufactured.plot_ex_x_axis_diagnostic(
            np,
            plt,
            analytic,
            opalx,
            origin,
            spacing,
            title,
            output,
            ip_beam_z,
            axis_limits=axis_limits,
        )
        return
    if components == ("Ey-y-axis",):
        manufactured.plot_ey_y_axis_diagnostic(
            np,
            plt,
            analytic,
            opalx,
            origin,
            spacing,
            title,
            output,
            ip_beam_z,
            axis_limits=axis_limits,
        )
        return
    if components == ("Ez-z-axis",):
        manufactured.plot_ez_z_axis_diagnostic(
            np,
            plt,
            analytic,
            opalx,
            origin,
            spacing,
            title,
            output,
            axis_limits=axis_limits,
        )
        return
    manufactured.plot_comparison(
        np,
        plt,
        analytic,
        opalx,
        origin,
        spacing,
        title,
        output,
        components=components,
        extent_override=None if axis_limits is None else axis_limits["extent"],
        scale_overrides=scale_overrides,
    )


def prepare_ez_along_bunch_data(
    *,
    h5_path: Path,
    step: int,
    state: str | None = None,
):
    manufactured = load_manufactured_solution()
    np = manufactured.load_numpy()
    h5py = manufactured.load_h5py()

    with h5py.File(h5_path, "r") as h5file:
        step_name, _snapshot_kind = manufactured.select_h5_step(
            h5file,
            step,
            normalize_compare_state(state),
        )
        step_group = h5file[step_name]
        ez_field, origin, spacing = manufactured.read_h5_scalar_field(np, step_group, "Ez")
        particle_mean_r = manufactured.decode_value(step_group.attrs.get("particle_mean_r"))
        interaction_point_s = manufactured.decode_value(step_group.attrs.get("interaction_point_s"))
        path_length_s = manufactured.decode_value(step_group.attrs.get("path_length_s"))
        snapshot_kind = str(manufactured.decode_value(step_group.attrs.get("snapshot_kind", "")))

    if not (isinstance(particle_mean_r, list) and len(particle_mean_r) >= 3):
        raise SystemExit("particle_mean_r missing or invalid in HDF5 step.")

    nx = ez_field.shape[2]
    ny = ez_field.shape[1]
    nz = ez_field.shape[0]
    x_centers = origin[0] + (np.arange(nx) + 0.5) * spacing[0]
    y_centers = origin[1] + (np.arange(ny) + 0.5) * spacing[1]
    z_centers = origin[2] + (np.arange(nz) + 0.5) * spacing[2]

    x_mean = float(particle_mean_r[0])
    y_mean = float(particle_mean_r[1])
    z_mean = float(particle_mean_r[2])

    ix = int(np.argmin(np.abs(x_centers - x_mean)))
    iy = int(np.argmin(np.abs(y_centers - y_mean)))
    ez_line = ez_field[:, iy, ix]

    ip_beam_z = None
    if isinstance(interaction_point_s, float) and isinstance(path_length_s, float):
        ip_beam_z = interaction_point_s - path_length_s

    title = f"step {step}"
    return {
        "np": np,
        "z_centers": z_centers,
        "ez_line": ez_line,
        "z_mean": z_mean,
        "ip_beam_z": ip_beam_z,
        "title": title,
    }


def plot_ez_along_bunch(
    *,
    plt,
    z_centers,
    ez_line,
    z_mean: float,
    ip_beam_z: float | None,
    title: str,
    output: Path | None = None,
    axis_limits=None,
) -> None:
    peak_indices = []
    for index in range(1, len(ez_line) - 1):
        value = float(ez_line[index])
        if value <= 0.0:
            continue
        if value >= float(ez_line[index - 1]) and value > float(ez_line[index + 1]):
            peak_indices.append(index)

    if len(peak_indices) >= 2:
        peak_indices = sorted(
            peak_indices,
            key=lambda index: float(ez_line[index]),
            reverse=True,
        )[:2]
        peak_indices.sort(key=lambda index: float(z_centers[index]))
    else:
        peak_indices = []

    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    ax.plot(z_centers, ez_line, linewidth=2.0, color="tab:blue", label="OPALX $E_z$")
    if ip_beam_z is not None:
        ax.axvline(ip_beam_z, color="darkorange", linestyle=":", linewidth=1.2, label="IP")

    if len(peak_indices) == 2:
        left_index, right_index = peak_indices
        left_z = float(z_centers[left_index])
        right_z = float(z_centers[right_index])
        left_peak = float(ez_line[left_index])
        right_peak = float(ez_line[right_index])
        span_y = 0.92 * min(left_peak, right_peak)
        ax.hlines(span_y, left_z, right_z, colors="crimson", linewidth=2.0)
        ax.vlines(
            [left_z, right_z],
            span_y * 0.985,
            span_y * 1.015,
            colors="crimson",
            linewidth=1.4,
        )
        ax.text(
            0.5 * (left_z + right_z),
            span_y,
            f"{right_z - left_z:.6g} m",
            color="crimson",
            ha="center",
            va="bottom",
            bbox={"boxstyle": "round,pad=0.18", "fc": "white", "ec": "crimson", "alpha": 0.85},
        )

    ax.set_xlabel("z [m]")
    ax.set_ylabel(r"$E_z$ [V/m]")
    ax.set_title(title)
    ax.grid(True, color="0.75", linewidth=0.9, alpha=0.9)
    if ip_beam_z is not None:
        ax.legend()
    if axis_limits is not None:
        ax.set_xlim(*axis_limits["x"])
        ax.set_ylim(*axis_limits["y"])

    if output is not None:
        fig.savefig(output, dpi=150)
        plt.close(fig)
    else:
        plt.show()


def run_ez_along_bunch_plot(
    *,
    h5_path: Path,
    step: int,
    state: str | None = None,
    output: Path | None = None,
    force_agg: bool = False,
) -> None:
    configure_matplotlib_backend(force_agg or output is not None)
    manufactured = load_manufactured_solution()
    plt = manufactured.load_matplotlib()
    data = prepare_ez_along_bunch_data(h5_path=h5_path, step=step, state=state)
    plot_ez_along_bunch(
        plt=plt,
        z_centers=data["z_centers"],
        ez_line=data["ez_line"],
        z_mean=data["z_mean"],
        ip_beam_z=data["ip_beam_z"],
        title=data["title"],
        output=output,
    )


def run_ez_along_bunch_movie(
    *,
    h5_path: Path,
    start: int,
    end: int,
    state: str | None = None,
    output: Path | None = None,
    duration_ms: int = 250,
    force_agg: bool = False,
) -> None:
    if output is None:
        raise SystemExit("movie output path is required")

    configure_matplotlib_backend(True if force_agg or output is not None else False)
    manufactured = load_manufactured_solution()
    plt = manufactured.load_matplotlib()
    image_cls = load_pillow_image()

    lo = min(start, end)
    hi = max(start, end)
    prepared_steps = []
    xmins = []
    xmaxs = []
    ymins = []
    ymaxs = []

    for step in range(lo, hi + 1):
        data = prepare_ez_along_bunch_data(h5_path=h5_path, step=step, state=state)
        prepared_steps.append((step, data))
        xmins.append(float(data["np"].min(data["z_centers"])))
        xmaxs.append(float(data["np"].max(data["z_centers"])))
        ymins.append(float(data["np"].min(data["ez_line"])))
        ymaxs.append(float(data["np"].max(data["ez_line"])))

    def padded(lo_value: float, hi_value: float) -> tuple[float, float]:
        span = hi_value - lo_value
        pad = 0.05 * max(abs(lo_value), abs(hi_value), span, 1.0e-30)
        return lo_value - pad, hi_value + pad

    axis_limits = {
        "x": (min(xmins), max(xmaxs)),
        "y": padded(min(ymins), max(ymaxs)),
    }

    frames = []
    with tempfile.TemporaryDirectory(prefix="beambeam-ez-along-bunch-") as tmpdir:
        tmpdir_path = Path(tmpdir)
        for step, data in prepared_steps:
            frame_path = tmpdir_path / f"ez-{step:06d}.png"
            plot_ez_along_bunch(
                plt=plt,
                z_centers=data["z_centers"],
                ez_line=data["ez_line"],
                z_mean=data["z_mean"],
                ip_beam_z=data["ip_beam_z"],
                title=data["title"],
                output=frame_path,
                axis_limits=axis_limits,
            )
            frames.append(image_cls.open(frame_path).convert("RGBA"))

        if not frames:
            raise SystemExit("no Ez-along-bunch frames generated")

        frames[0].save(
            output,
            save_all=True,
            append_images=frames[1:],
            duration=duration_ms,
            loop=0,
        )


def run_manufactured_compare_movie(
    *,
    h5_path: Path,
    start: int,
    end: int,
    sigma: float,
    state: str | None = None,
    comparison: str | None = None,
    output: Path | None = None,
    duration_ms: int = 250,
    force_agg: bool = False,
) -> None:
    if output is None:
        raise SystemExit("movie output path is required")

    configure_matplotlib_backend(True if force_agg or output is not None else False)
    manufactured = load_manufactured_solution()
    np = manufactured.load_numpy()
    plt = manufactured.load_matplotlib()
    image_cls = load_pillow_image()

    lo = min(start, end)
    hi = max(start, end)
    frames = []
    components = normalize_compare_component(comparison)
    if len(components) != 1:
        raise SystemExit("compare movie requires exactly one comparison kind: rho, phi, Ez, rho-z-axis, phi-z-axis, Ez-z-axis, Ex-x-axis, or Ey-y-axis")
    comparison_name = components[0]
    prepared_steps = []
    image_axis_limits = None
    line_axis_limits = None
    scale_overrides = None

    if comparison_name in ("rho", "phi", "Ez"):
        xmins = []
        xmaxs = []
        zmins = []
        zmaxs = []
        field_vmax = 0.0
        diff_dmax = 0.0
        for step in range(lo, hi + 1):
            setup, analytic, origin, spacing, ip_beam_z, title = prepare_manufactured_compare(
                manufactured,
                np,
                h5_path,
                step,
                sigma,
                state,
            )
            analytic_panel, extent = manufactured.central_xz_slice(analytic[comparison_name], origin, spacing)
            opalx_panel, _ = manufactured.central_xz_slice(setup["opalx"][comparison_name], origin, spacing)
            diff_panel = opalx_panel - analytic_panel
            xmins.append(extent[0])
            xmaxs.append(extent[1])
            zmins.append(extent[2])
            zmaxs.append(extent[3])
            field_vmax = max(field_vmax, float(np.max(np.abs(analytic_panel))), float(np.max(np.abs(opalx_panel))))
            diff_dmax = max(diff_dmax, float(np.max(np.abs(diff_panel))))
            prepared_steps.append((step, setup, analytic, origin, spacing, ip_beam_z, title))
        image_axis_limits = {"extent": (min(xmins), max(xmaxs), min(zmins), max(zmaxs))}
        scale_overrides = {
            comparison_name: {
                "vmax": field_vmax,
                "dmax": diff_dmax,
            }
        }
    else:
        xmins = []
        xmaxs = []
        main_ymins = []
        main_ymaxs = []
        diff_ymins = []
        diff_ymaxs = []
        for step in range(lo, hi + 1):
            setup, analytic, origin, spacing, ip_beam_z, title = prepare_manufactured_compare(
                manufactured,
                np,
                h5_path,
                step,
                sigma,
                state,
            )
            if comparison_name == "rho-z-axis":
                z_coords, analytic_line = manufactured.central_z_axis_line(np, analytic["rho"], origin, spacing)
                _, opalx_line = manufactured.central_z_axis_line(np, setup["opalx"]["rho"], origin, spacing)
            elif comparison_name == "phi-z-axis":
                z_coords, analytic_line = manufactured.central_z_axis_line(np, analytic["phi"], origin, spacing)
                _, opalx_line = manufactured.central_z_axis_line(np, setup["opalx"]["phi"], origin, spacing)
            elif comparison_name == "Ez-z-axis":
                z_coords, analytic_line = manufactured.central_z_axis_line(np, analytic["Ez"], origin, spacing)
                _, opalx_line = manufactured.central_z_axis_line(np, setup["opalx"]["Ez"], origin, spacing)
            elif comparison_name == "Ex-x-axis":
                z_coords, analytic_line = manufactured.transverse_x_axis_line(
                    np, analytic["Ex"], origin, spacing, ip_beam_z
                )
                _, opalx_line = manufactured.transverse_x_axis_line(
                    np, setup["opalx"]["Ex"], origin, spacing, ip_beam_z
                )
            else:
                z_coords, analytic_line = manufactured.transverse_y_axis_line(
                    np, analytic["Ey"], origin, spacing, ip_beam_z
                )
                _, opalx_line = manufactured.transverse_y_axis_line(
                    np, setup["opalx"]["Ey"], origin, spacing, ip_beam_z
                )
            diff_line = opalx_line - analytic_line
            xmins.append(float(np.min(z_coords)))
            xmaxs.append(float(np.max(z_coords)))
            main_ymins.append(float(min(np.min(analytic_line), np.min(opalx_line))))
            main_ymaxs.append(float(max(np.max(analytic_line), np.max(opalx_line))))
            diff_ymins.append(float(np.min(diff_line)))
            diff_ymaxs.append(float(np.max(diff_line)))
            prepared_steps.append((step, setup, analytic, origin, spacing, ip_beam_z, title))

        def padded(lo_value: float, hi_value: float) -> tuple[float, float]:
            span = hi_value - lo_value
            pad = 0.05 * max(abs(lo_value), abs(hi_value), span, 1.0e-30)
            return lo_value - pad, hi_value + pad

        line_axis_limits = {
            "x": (min(xmins), max(xmaxs)),
            "main_y": padded(min(main_ymins), max(main_ymaxs)),
            "diff_y": padded(min(diff_ymins), max(diff_ymaxs)),
        }

    with tempfile.TemporaryDirectory(prefix="beambeam-compare-") as tmpdir:
        tmpdir_path = Path(tmpdir)
        for step, setup, analytic, origin, spacing, ip_beam_z, title in prepared_steps:
            frame_path = tmpdir_path / f"compare-{step:06d}.png"
            render_manufactured_compare_output(
                manufactured,
                np,
                plt,
                analytic,
                setup["opalx"],
                origin,
                spacing,
                ip_beam_z,
                title,
                frame_path,
                comparison,
                axis_limits=image_axis_limits if image_axis_limits is not None else line_axis_limits,
                scale_overrides=scale_overrides,
            )
            frames.append(image_cls.open(frame_path).convert("P"))

    if not frames:
        raise SystemExit("no comparison frames generated")

    frames[0].save(
        output,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
    )


def create_tk_gui() -> None:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
    from PIL import Image, ImageSequence, ImageTk

    gui_state = load_gui_state()
    root = tk.Tk()
    root.title("BeamBeam Analysis")
    root.geometry("860x520")

    notebook = ttk.Notebook(root)
    notebook.pack(fill="both", expand=True, padx=10, pady=(10, 4))

    footer = ttk.Frame(root, padding=(10, 4, 10, 10))
    footer.pack(fill="x")
    footer.columnconfigure(0, weight=1)

    progress_status = tk.StringVar(value="Idle")
    ttk.Label(footer, textvariable=progress_status).grid(row=0, column=0, sticky="w")
    progressbar = ttk.Progressbar(footer, mode="indeterminate", length=220)
    progressbar.grid(row=0, column=1, padx=(10, 10), sticky="e")
    ttk.Button(footer, text="Exit", command=root.destroy).grid(row=0, column=2, sticky="e")

    def browse_file(variable, filetypes):
        value = filedialog.askopenfilename(filetypes=filetypes)
        if value:
            variable.set(value)

    def browse_save(variable, default_ext):
        value = filedialog.asksaveasfilename(defaultextension=default_ext)
        if value:
            variable.set(value)

    def temp_output_path(stem: str, suffix: str) -> Path:
        return Path(tempfile.gettempdir()) / f"{stem}{suffix}"

    def show_render_window(
        title: str,
        paths: list[Path],
        frame_labels_by_path: dict[str, list[int]] | None = None,
    ) -> None:
        viewer = tk.Toplevel(root)
        viewer.title(title)
        viewer.geometry("1200x900")

        notebook_widget = ttk.Notebook(viewer)
        notebook_widget.pack(fill="both", expand=True, padx=8, pady=8)

        for path in paths:
            if not path.exists():
                continue

            frame = ttk.Frame(notebook_widget, padding=8)
            notebook_widget.add(frame, text=path.name)
            ttk.Label(frame, text=str(path)).pack(anchor="w", pady=(0, 8))

            image = Image.open(path)
            max_size = (1100, 760)

            if path.suffix.lower() == ".gif":
                tk_frames = []
                for frame_image in ImageSequence.Iterator(image):
                    converted = frame_image.convert("RGBA")
                    converted.thumbnail(max_size, Image.Resampling.LANCZOS)
                    tk_frames.append(ImageTk.PhotoImage(converted))

                if not tk_frames:
                    continue

                step_labels = None
                if frame_labels_by_path is not None:
                    step_labels = frame_labels_by_path.get(str(path))
                duration_ms = int(image.info.get("duration", 180))
                controls = ttk.Frame(frame)
                controls.pack(fill="x", pady=(0, 8))

                image_label = ttk.Label(frame)
                image_label.pack(fill="both", expand=True)

                current_index = {"value": 0}
                playing = {"value": True}
                after_id = {"value": None}
                suppress_slider_callback = {"value": False}

                step_text = tk.StringVar()

                def update_step_text(index: int) -> None:
                    if step_labels is not None and index < len(step_labels):
                        step_text.set(f"step {step_labels[index]}")
                    else:
                        step_text.set(f"frame {index}")

                def show_frame(index: int, *, update_slider: bool = True) -> None:
                    frame_index = max(0, min(index, len(tk_frames) - 1))
                    current_index["value"] = frame_index
                    image_label.configure(image=tk_frames[frame_index])
                    image_label.image = tk_frames[frame_index]
                    if update_slider:
                        suppress_slider_callback["value"] = True
                        try:
                            slider.set(frame_index)
                        finally:
                            suppress_slider_callback["value"] = False
                    update_step_text(frame_index)

                def advance() -> None:
                    if not playing["value"]:
                        after_id["value"] = None
                        return
                    next_index = (current_index["value"] + 1) % len(tk_frames)
                    show_frame(next_index)
                    after_id["value"] = viewer.after(duration_ms, advance)

                def toggle_play() -> None:
                    playing["value"] = not playing["value"]
                    play_button.configure(text="Stop" if playing["value"] else "Start")
                    if playing["value"] and after_id["value"] is None:
                        after_id["value"] = viewer.after(duration_ms, advance)

                def on_slider(value: str) -> None:
                    if suppress_slider_callback["value"]:
                        return
                    if after_id["value"] is not None:
                        viewer.after_cancel(after_id["value"])
                        after_id["value"] = None
                    show_frame(int(float(value)), update_slider=False)
                    if playing["value"]:
                        after_id["value"] = viewer.after(duration_ms, advance)

                play_button = ttk.Button(controls, text="Stop", command=toggle_play)
                play_button.pack(side="left")
                ttk.Label(controls, textvariable=step_text).pack(side="left", padx=(10, 0))
                slider = tk.Scale(
                    controls,
                    from_=0,
                    to=len(tk_frames) - 1,
                    orient="horizontal",
                    showvalue=False,
                    command=on_slider,
                    length=520,
                )
                slider.pack(side="right", fill="x", expand=True)

                image_label._frames = tk_frames
                show_frame(0)
                after_id["value"] = viewer.after(duration_ms, advance)
                continue

            image_label = ttk.Label(frame)
            image_label.pack(fill="both", expand=True)
            image.thumbnail(max_size, Image.Resampling.LANCZOS)
            photo = ImageTk.PhotoImage(image)
            image_label.configure(image=photo)
            image_label.image = photo

    def show_text_window(title: str, text: str) -> None:
        viewer = tk.Toplevel(root)
        viewer.title(title)
        viewer.geometry("1300x700")

        container = ttk.Frame(viewer, padding=8)
        container.pack(fill="both", expand=True)

        text_widget = tk.Text(container, wrap="none", font=("Menlo", 11))
        y_scroll = ttk.Scrollbar(container, orient="vertical", command=text_widget.yview)
        x_scroll = ttk.Scrollbar(container, orient="horizontal", command=text_widget.xview)
        text_widget.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)

        text_widget.grid(row=0, column=0, sticky="nsew")
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll.grid(row=1, column=0, sticky="ew")
        container.columnconfigure(0, weight=1)
        container.rowconfigure(0, weight=1)

        text_widget.insert("1.0", text)
        text_widget.configure(state="disabled")

    def add_path_row(parent, row, label, variable, browse_types, save=False, default_ext=".png"):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=6, pady=4)
        ttk.Entry(parent, textvariable=variable, width=70).grid(
            row=row, column=1, sticky="ew", padx=6, pady=4
        )
        command = (
            (lambda: browse_save(variable, default_ext))
            if save
            else (lambda: browse_file(variable, browse_types))
        )
        ttk.Button(parent, text="Browse", command=command).grid(row=row, column=2, padx=6, pady=4)

    def handle_errors(action, on_success=None):
        progress_status.set("Reading and preprocessing...")
        progressbar.start(12)
        root.config(cursor="watch")
        root.update_idletasks()

        def finish_success(result):
            progress_status.set("Completed")
            progressbar.stop()
            root.config(cursor="")
            root.update_idletasks()
            if on_success is not None:
                on_success(result)

        def finish_failure(exc: Exception):
            progress_status.set("Failed")
            progressbar.stop()
            root.config(cursor="")
            root.update_idletasks()
            messagebox.showerror("BeamBeam Analysis", str(exc))

        def worker():
            try:
                result = action()
            except Exception as exc:
                root.after(0, lambda exc=exc: finish_failure(exc))
                return
            root.after(0, lambda result=result: finish_success(result))

        threading.Thread(target=worker, daemon=True).start()

    collwin_frame = ttk.Frame(notebook, padding=10)
    collwin_frame.columnconfigure(1, weight=1)
    notebook.add(collwin_frame, text="Collision Window")

    collwin_state = gui_state["collwin"]
    collwin_file = tk.StringVar(value=str(collwin_state["filename"]))
    collwin_input = tk.StringVar(value=str(collwin_state["input"]))
    collwin_output = tk.StringVar(value=str(collwin_state["output"]))
    collwin_start = tk.StringVar(value=str(collwin_state["start"]))
    collwin_end = tk.StringVar(value=str(collwin_state["end"]))
    collwin_grid = tk.BooleanVar(value=bool(collwin_state["grid"]))
    collwin_no_geometry = tk.BooleanVar(value=bool(collwin_state["no_geometry"]))

    add_path_row(
        collwin_frame,
        0,
        "RHO dump",
        collwin_file,
        [("Scalar dumps", "*.dat"), ("All files", "*")],
    )
    add_path_row(
        collwin_frame,
        1,
        "Input file",
        collwin_input,
        [("Input files", "*.in"), ("All files", "*")],
    )
    add_path_row(
        collwin_frame,
        2,
        "Output",
        collwin_output,
        [("Images", "*.png *.gif"), ("All files", "*")],
        save=True,
        default_ext=".png",
    )
    ttk.Label(collwin_frame, text="Start step").grid(row=3, column=0, sticky="w", padx=6, pady=4)
    ttk.Entry(collwin_frame, textvariable=collwin_start, width=12).grid(
        row=3, column=1, sticky="w", padx=6, pady=4
    )
    ttk.Label(collwin_frame, text="End step").grid(row=4, column=0, sticky="w", padx=6, pady=4)
    ttk.Entry(collwin_frame, textvariable=collwin_end, width=12).grid(
        row=4, column=1, sticky="w", padx=6, pady=4
    )
    ttk.Checkbutton(collwin_frame, text="Plot in grid coordinates", variable=collwin_grid).grid(
        row=5, column=0, columnspan=2, sticky="w", padx=6, pady=4
    )
    ttk.Checkbutton(collwin_frame, text="Hide BeamBeam geometry overlays", variable=collwin_no_geometry).grid(
        row=6, column=0, columnspan=2, sticky="w", padx=6, pady=4
    )

    def run_collwin_from_gui():
        filename = Path(collwin_file.get())
        if not filename.exists():
            raise FileNotFoundError(f"missing file: {filename}")
        input_file = Path(collwin_input.get()) if collwin_input.get().strip() else None
        output = Path(collwin_output.get()) if collwin_output.get().strip() else None
        start = int(collwin_start.get()) if collwin_start.get().strip() else None
        end = int(collwin_end.get()) if collwin_end.get().strip() else None
        display_output = output
        if display_output is None:
            suffix = ".png" if start is None or start == end else ".gif"
            display_output = temp_output_path("beambeam-collwin-view", suffix)
        frame_labels = None
        if start is not None and end is not None and display_output.suffix.lower() == ".gif":
            lo = min(start, end)
            hi = max(start, end)
            frame_labels = {str(display_output): list(range(lo, hi + 1))}
        state_snapshot = {
            "filename": collwin_file.get(),
            "input": collwin_input.get(),
            "output": collwin_output.get(),
            "start": collwin_start.get(),
            "end": collwin_end.get(),
            "grid": collwin_grid.get(),
            "no_geometry": collwin_no_geometry.get(),
        }

        def job():
            run_collision_window(
                filename,
                output=display_output,
                start=start,
                end=end,
                grid=state_snapshot["grid"],
                no_geometry=state_snapshot["no_geometry"],
                input_file=input_file,
                force_agg=True,
            )
            gui_state["collwin"] = state_snapshot
            save_gui_state(gui_state)
            return ("Collision Window", [display_output], frame_labels)

        handle_errors(job, lambda result: show_render_window(*result))

    ttk.Button(
        collwin_frame,
        text="Run Collision Window View",
        command=run_collwin_from_gui,
    ).grid(row=7, column=0, columnspan=3, sticky="ew", padx=6, pady=12)

    analytic_frame = ttk.Frame(notebook, padding=10)
    analytic_frame.columnconfigure(1, weight=1)
    notebook.add(analytic_frame, text="Manufactured Analytic")

    analytic_state = gui_state["analytic"]
    sigma_var = tk.StringVar(value=str(analytic_state["sigma"]))
    charge_var = tk.StringVar(value=str(analytic_state["charge"]))
    separation_var = tk.StringVar(value=str(analytic_state["half_separation"]))
    nx_var = tk.StringVar(value=str(analytic_state["nx"]))
    ny_var = tk.StringVar(value=str(analytic_state["ny"]))
    nz_var = tk.StringVar(value=str(analytic_state["nz"]))
    analytic_output = tk.StringVar(value=str(analytic_state["output"]))

    for row, (label, variable) in enumerate(
        (
            ("sigma", sigma_var),
            ("charge", charge_var),
            ("half separation", separation_var),
            ("nx", nx_var),
            ("ny", ny_var),
            ("nz", nz_var),
        )
    ):
        ttk.Label(analytic_frame, text=label).grid(row=row, column=0, sticky="w", padx=6, pady=4)
        ttk.Entry(analytic_frame, textvariable=variable, width=16).grid(
            row=row, column=1, sticky="w", padx=6, pady=4
        )

    add_path_row(
        analytic_frame,
        6,
        "Output",
        analytic_output,
        [("Images", "*.png"), ("All files", "*")],
        save=True,
        default_ext=".png",
    )

    def run_analytic_from_gui():
        output = Path(analytic_output.get()) if analytic_output.get().strip() else None
        display_output = output or temp_output_path("beambeam-manufactured-analytic", ".png")
        state_snapshot = {
            "sigma": sigma_var.get(),
            "charge": charge_var.get(),
            "half_separation": separation_var.get(),
            "nx": nx_var.get(),
            "ny": ny_var.get(),
            "nz": nz_var.get(),
            "output": analytic_output.get(),
        }

        def job():
            run_manufactured_analytic(
                sigma=float(state_snapshot["sigma"]),
                charge=float(state_snapshot["charge"]),
                half_separation=float(state_snapshot["half_separation"]),
                nx=int(state_snapshot["nx"]),
                ny=int(state_snapshot["ny"]),
                nz=int(state_snapshot["nz"]),
                output=display_output,
                force_agg=True,
            )
            gui_state["analytic"] = state_snapshot
            save_gui_state(gui_state)
            return ("Manufactured Analytic", [display_output])

        handle_errors(job, lambda result: show_render_window(*result))

    ttk.Button(
        analytic_frame,
        text="Run Manufactured Analytic",
        command=run_analytic_from_gui,
    ).grid(row=7, column=0, columnspan=3, sticky="ew", padx=6, pady=12)

    compare_frame = ttk.Frame(notebook, padding=10)
    compare_frame.columnconfigure(1, weight=1)
    notebook.add(compare_frame, text="Manufactured Compare")

    compare_tab_state = gui_state["compare"]
    compare_h5 = tk.StringVar(value=str(compare_tab_state["h5_path"]))
    compare_step = tk.StringVar(value=str(compare_tab_state["step"]))
    compare_start = tk.StringVar(value=str(compare_tab_state["start"]))
    compare_end = tk.StringVar(value=str(compare_tab_state["end"]))
    compare_sigma = tk.StringVar(value=str(compare_tab_state["sigma"]))
    compare_state = tk.StringVar(value=str(compare_tab_state["state"]))
    compare_component = tk.StringVar(value=str(compare_tab_state["component"]))
    compare_output = tk.StringVar(value=str(compare_tab_state["output"]))

    add_path_row(
        compare_frame,
        0,
        "BeamBeam HDF5",
        compare_h5,
        [("HDF5 files", "*.h5 *.hdf5"), ("All files", "*")],
    )
    ttk.Label(compare_frame, text="step").grid(row=1, column=0, sticky="w", padx=6, pady=4)
    ttk.Entry(compare_frame, textvariable=compare_step, width=16).grid(
        row=1, column=1, sticky="w", padx=6, pady=4
    )
    ttk.Label(compare_frame, text="start").grid(row=2, column=0, sticky="w", padx=6, pady=4)
    ttk.Entry(compare_frame, textvariable=compare_start, width=16).grid(
        row=2, column=1, sticky="w", padx=6, pady=4
    )
    ttk.Label(compare_frame, text="end").grid(row=3, column=0, sticky="w", padx=6, pady=4)
    ttk.Entry(compare_frame, textvariable=compare_end, width=16).grid(
        row=3, column=1, sticky="w", padx=6, pady=4
    )
    ttk.Label(compare_frame, text="sigma").grid(row=4, column=0, sticky="w", padx=6, pady=4)
    ttk.Entry(compare_frame, textvariable=compare_sigma, width=16).grid(
        row=4, column=1, sticky="w", padx=6, pady=4
    )
    ttk.Label(compare_frame, text="state").grid(row=5, column=0, sticky="w", padx=6, pady=4)
    state_row = ttk.Frame(compare_frame)
    state_row.grid(row=5, column=1, sticky="w", padx=6, pady=4)
    ttk.Combobox(
        state_row,
        textvariable=compare_state,
        values=("", "beambeam tracking", "normal tracking", "all states"),
        width=22,
    ).grid(row=0, column=0, sticky="w")
    ttk.Button(
        state_row,
        text="All States",
        command=lambda: compare_state.set("all states"),
    ).grid(row=0, column=1, sticky="w", padx=(6, 0))
    ttk.Label(compare_frame, text="comparison").grid(row=6, column=0, sticky="w", padx=6, pady=4)
    ttk.Combobox(
        compare_frame,
        textvariable=compare_component,
        values=("rho+phi", "rho+phi+Ez", "rho", "phi", "Ez", "rho-z-axis", "phi-z-axis", "Ez-z-axis", "Ex-x-axis", "Ey-y-axis"),
        width=28,
        state="readonly",
    ).grid(row=6, column=1, sticky="w", padx=6, pady=4)
    add_path_row(
        compare_frame,
        7,
        "Output",
        compare_output,
        [("Images", "*.png"), ("All files", "*")],
        save=True,
        default_ext=".png",
    )

    def run_compare_from_gui():
        h5_path = Path(compare_h5.get())
        if not h5_path.exists():
            raise FileNotFoundError(f"missing file: {h5_path}")
        output = Path(compare_output.get()) if compare_output.get().strip() else None
        state = compare_state.get().strip() or None
        start = compare_start.get().strip()
        end = compare_end.get().strip()
        state_snapshot = {
            "h5_path": compare_h5.get(),
            "step": compare_step.get(),
            "start": compare_start.get(),
            "end": compare_end.get(),
            "sigma": compare_sigma.get(),
            "state": compare_state.get(),
            "component": compare_component.get(),
            "output": compare_output.get(),
        }
        if start or end:
            if not (start and end):
                raise ValueError("both start and end are required for a compare movie")
            display_output = output or temp_output_path("beambeam-compare-movie", ".gif")
            lo = min(int(start), int(end))
            hi = max(int(start), int(end))
            frame_labels = {str(display_output): list(range(lo, hi + 1))}

            def job():
                run_manufactured_compare_movie(
                    h5_path=h5_path,
                    start=int(start),
                    end=int(end),
                    sigma=float(state_snapshot["sigma"]),
                    state=state,
                    comparison=state_snapshot["component"],
                    output=display_output,
                    force_agg=True,
                )
                gui_state["compare"] = state_snapshot
                save_gui_state(gui_state)
                return ("Manufactured Compare Movie", [display_output], frame_labels)

            handle_errors(job, lambda result: show_render_window(*result))
            return
        display_output = output or temp_output_path("beambeam-compare", ".png")

        def job():
            run_manufactured_compare(
                h5_path=h5_path,
                step=int(state_snapshot["step"]),
                sigma=float(state_snapshot["sigma"]),
                state=state,
                comparison=state_snapshot["component"],
                output=display_output,
                force_agg=True,
            )
            gui_state["compare"] = state_snapshot
            save_gui_state(gui_state)
            preview_paths = [display_output]
            components = normalize_compare_component(state_snapshot["component"])
            if components in (("rho", "phi"), ("rho", "phi", "Ez")):
                preview_paths.extend(
                    [
                        display_output.with_name(display_output.stem + "-rho-z-axis" + display_output.suffix),
                        display_output.with_name(display_output.stem + "-phi-z-axis" + display_output.suffix),
                        display_output.with_name(display_output.stem + "-ez-z-axis" + display_output.suffix),
                        display_output.with_name(display_output.stem + "-ex-x-axis" + display_output.suffix),
                        display_output.with_name(display_output.stem + "-ey-y-axis" + display_output.suffix),
                    ]
                )
            return (
                "Manufactured Compare",
                preview_paths,
            )

        handle_errors(job, lambda result: show_render_window(*result))

    def show_compare_overview_from_gui():
        h5_path = Path(compare_h5.get())
        if not h5_path.exists():
            raise FileNotFoundError(f"missing file: {h5_path}")
        overview_text = build_overview_table_text(h5_path, None, None)
        show_text_window("BeamBeam Overview Table", overview_text)

    def show_ez_along_bunch_from_gui():
        h5_path = Path(compare_h5.get())
        if not h5_path.exists():
            raise FileNotFoundError(f"missing file: {h5_path}")
        state_value = compare_state.get().strip() or None
        start_value = compare_start.get().strip()
        end_value = compare_end.get().strip()
        step_value = compare_step.get().strip()

        if start_value or end_value:
            if not (start_value and end_value):
                raise ValueError("both start and end are required for Ez along bunch range output")
            start_int = int(start_value)
            end_int = int(end_value)
            is_movie = start_int != end_int
            selected_step = start_int
        else:
            if not step_value:
                raise ValueError("either step or both start and end are required for Ez along bunch")
            selected_step = int(step_value)
            start_int = selected_step
            end_int = selected_step
            is_movie = False

        display_output = temp_output_path(
            "beambeam-ez-along-bunch",
            ".gif" if is_movie else ".png",
        )

        def job():
            if is_movie:
                run_ez_along_bunch_movie(
                    h5_path=h5_path,
                    start=start_int,
                    end=end_int,
                    state=state_value,
                    output=display_output,
                    force_agg=True,
                )
                return (
                    "Ez Along Bunch",
                    [display_output],
                    {str(display_output): list(range(min(start_int, end_int), max(start_int, end_int) + 1))},
                )
            run_ez_along_bunch_plot(
                h5_path=h5_path,
                step=selected_step,
                state=state_value,
                output=display_output,
                force_agg=True,
            )
            return ("Ez Along Bunch", [display_output], None)

        handle_errors(job, lambda result: show_render_window(*result))

    button_row = ttk.Frame(compare_frame)
    button_row.grid(row=8, column=0, columnspan=3, sticky="ew", padx=6, pady=12)
    button_row.columnconfigure(0, weight=1)
    button_row.columnconfigure(1, weight=1)
    button_row.columnconfigure(2, weight=1)

    ttk.Button(
        button_row,
        text="Run Manufactured Compare",
        command=run_compare_from_gui,
    ).grid(row=0, column=0, sticky="ew", padx=(0, 6))
    ttk.Button(
        button_row,
        text="Show Overview Table",
        command=show_compare_overview_from_gui,
    ).grid(row=0, column=1, sticky="ew", padx=(6, 0))
    ttk.Button(
        button_row,
        text="Show Ez Along Bunch",
        command=show_ez_along_bunch_from_gui,
    ).grid(row=0, column=2, sticky="ew", padx=(6, 0))

    root.mainloop()


def html_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def build_web_gui_page(
    status: str = "",
    values: dict[str, str] | None = None,
    previews: list[tuple[str, Path]] | None = None,
) -> str:
    gui_state = load_gui_state()
    web_defaults = {
        "collwin_filename": str(gui_state["collwin"]["filename"]),
        "collwin_input": str(gui_state["collwin"]["input"]),
        "collwin_output": str(gui_state["collwin"]["output"]),
        "collwin_start": str(gui_state["collwin"]["start"]),
        "collwin_end": str(gui_state["collwin"]["end"]),
        "collwin_physical": "" if bool(gui_state["collwin"]["grid"]) else "1",
        "collwin_no_geometry": "1" if bool(gui_state["collwin"]["no_geometry"]) else "",
        "analytic_sigma": str(gui_state["analytic"]["sigma"]),
        "analytic_charge": str(gui_state["analytic"]["charge"]),
        "analytic_half_separation": str(gui_state["analytic"]["half_separation"]),
        "analytic_nx": str(gui_state["analytic"]["nx"]),
        "analytic_ny": str(gui_state["analytic"]["ny"]),
        "analytic_nz": str(gui_state["analytic"]["nz"]),
        "analytic_output": str(gui_state["analytic"]["output"]),
        "compare_h5_path": str(gui_state["compare"]["h5_path"]),
        "compare_step": str(gui_state["compare"]["step"]),
        "compare_start": str(gui_state["compare"]["start"]),
        "compare_end": str(gui_state["compare"]["end"]),
        "compare_sigma": str(gui_state["compare"]["sigma"]),
        "compare_state": str(gui_state["compare"]["state"]),
        "compare_component": str(gui_state["compare"]["component"]),
        "compare_output": str(gui_state["compare"]["output"]),
    }
    values = {**web_defaults, **(values or {})}
    previews = previews or []

    def get(name: str, default: str = "") -> str:
        return html_escape(values.get(name, default))

    status_html = ""
    if status:
        status_html = f"<pre class='status'>{html_escape(status)}</pre>"

    preview_html = ""
    if previews:
        cards = []
        for label, path in previews:
            url = "/preview?" + urlencode({"path": str(path)})
            cards.append(
                f"<div class='preview-card'><h3>{html_escape(label)}</h3>"
                f"<a href='{html_escape(url)}' target='_blank'>{html_escape(str(path))}</a>"
                f"<img src='{html_escape(url)}' alt='{html_escape(label)}'></div>"
            )
        preview_html = f"<div class='preview-grid'>{''.join(cards)}</div>"

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>BeamBeam Analysis</title>
  <style>
    :root {{
      --bg: #f3efe3;
      --panel: #fffdf7;
      --ink: #1f1f1f;
      --accent: #0d5c63;
      --line: #d4c9b5;
    }}
    body {{
      margin: 0;
      padding: 24px;
      background: linear-gradient(180deg, #ece6d7, var(--bg));
      color: var(--ink);
      font-family: Georgia, "Iowan Old Style", serif;
    }}
    h1 {{
      margin: 0 0 12px 0;
      font-size: 28px;
    }}
    p {{
      margin: 0 0 18px 0;
      max-width: 72ch;
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 18px;
    }}
    .panel {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 12px;
      padding: 16px;
      box-shadow: 0 12px 30px rgba(0, 0, 0, 0.05);
    }}
    .panel h2 {{
      margin: 0 0 10px 0;
      font-size: 20px;
    }}
    label {{
      display: block;
      margin-top: 10px;
      font-size: 14px;
    }}
    input[type="text"], input[type="number"] {{
      width: 100%;
      box-sizing: border-box;
      margin-top: 4px;
      padding: 8px 10px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: white;
      font-family: Menlo, Consolas, monospace;
      font-size: 13px;
    }}
    .inline-row {{
      display: flex;
      gap: 8px;
      align-items: end;
    }}
    .inline-row input[type="text"], .inline-row input[type="number"] {{
      flex: 1;
    }}
    .check {{
      display: flex;
      gap: 8px;
      align-items: center;
      margin-top: 10px;
    }}
    button {{
      margin-top: 14px;
      padding: 10px 14px;
      border: none;
      border-radius: 999px;
      background: var(--accent);
      color: white;
      font-weight: 600;
      cursor: pointer;
    }}
    .status {{
      white-space: pre-wrap;
      background: #101418;
      color: #f3f7fa;
      padding: 14px;
      border-radius: 10px;
      overflow-x: auto;
      margin: 0 0 18px 0;
    }}
    .preview-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
      gap: 18px;
      margin: 0 0 18px 0;
    }}
    .preview-card {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 12px;
      padding: 12px;
      box-shadow: 0 12px 30px rgba(0, 0, 0, 0.05);
    }}
    .preview-card h3 {{
      margin: 0 0 8px 0;
      font-size: 16px;
    }}
    .preview-card a {{
      display: block;
      margin-bottom: 8px;
      color: var(--accent);
      text-decoration: none;
      font-family: Menlo, Consolas, monospace;
      font-size: 12px;
      word-break: break-all;
    }}
    .preview-card img {{
      width: 100%;
      height: auto;
      border-radius: 8px;
      border: 1px solid var(--line);
      background: white;
    }}
  </style>
</head>
<body>
  <h1>BeamBeam Analysis</h1>
  <p>Local browser fallback for systems where Tk cannot open. Enter the paths you want to analyze, then run the selected mode.</p>
  {status_html}
  {preview_html}
  <div class="grid">
    <form class="panel" method="post" action="/run">
      <input type="hidden" name="mode" value="collwin">
      <h2>Collision Window</h2>
      <label>RHO dump
        <input type="text" name="collwin_filename" value="{get('collwin_filename')}">
      </label>
      <label>Input file
        <input type="text" name="collwin_input" value="{get('collwin_input')}">
      </label>
      <label>Output
        <input type="text" name="collwin_output" value="{get('collwin_output')}">
      </label>
      <label>Start step
        <input type="number" name="collwin_start" value="{get('collwin_start')}">
      </label>
      <label>End step
        <input type="number" name="collwin_end" value="{get('collwin_end')}">
      </label>
      <label class="check"><input type="checkbox" name="collwin_physical" {"checked" if get("collwin_physical") else ""}>Plot in physical coordinates</label>
      <label class="check"><input type="checkbox" name="collwin_no_geometry" {"checked" if get("collwin_no_geometry") else ""}>Hide BeamBeam geometry overlays</label>
      <button type="submit">Run Collision Window</button>
    </form>

    <form class="panel" method="post" action="/run">
      <input type="hidden" name="mode" value="analytic">
      <h2>Manufactured Analytic</h2>
      <label>sigma
        <input type="text" name="analytic_sigma" value="{get('analytic_sigma')}">
      </label>
      <label>charge
        <input type="text" name="analytic_charge" value="{get('analytic_charge')}">
      </label>
      <label>half separation
        <input type="text" name="analytic_half_separation" value="{get('analytic_half_separation')}">
      </label>
      <label>nx
        <input type="number" name="analytic_nx" value="{get('analytic_nx')}">
      </label>
      <label>ny
        <input type="number" name="analytic_ny" value="{get('analytic_ny')}">
      </label>
      <label>nz
        <input type="number" name="analytic_nz" value="{get('analytic_nz')}">
      </label>
      <label>Output
        <input type="text" name="analytic_output" value="{get('analytic_output')}">
      </label>
      <button type="submit">Run Manufactured Analytic</button>
    </form>

    <form class="panel" method="post" action="/run">
      <input type="hidden" name="mode" value="compare">
      <h2>Manufactured Compare</h2>
      <label>BeamBeam HDF5
        <input type="text" name="compare_h5_path" value="{get('compare_h5_path')}">
      </label>
      <label>step
        <input type="number" name="compare_step" value="{get('compare_step')}">
      </label>
      <label>start
        <input type="number" name="compare_start" value="{get('compare_start')}">
      </label>
      <label>end
        <input type="number" name="compare_end" value="{get('compare_end')}">
      </label>
      <label>sigma
        <input type="text" name="compare_sigma" value="{get('compare_sigma')}">
      </label>
      <label>state
        <div class="inline-row">
          <input type="text" name="compare_state" id="compare_state" value="{get('compare_state')}">
          <button type="button" onclick="document.getElementById('compare_state').value='all states';">All States</button>
        </div>
      </label>
      <label>comparison
        <input type="text" name="compare_component" value="{get('compare_component')}">
      </label>
      <label>Output
        <input type="text" name="compare_output" value="{get('compare_output')}">
      </label>
      <button type="submit">Run Manufactured Compare / Movie</button>
    </form>
  </div>
</body>
</html>
"""


def run_web_gui(reason: BaseException | None = None) -> None:
    def temp_output_path(stem: str, suffix: str) -> Path:
        return Path(tempfile.gettempdir()) / f"{stem}{suffix}"

    def persist_web_form(form: dict[str, str]) -> None:
        gui_state = load_gui_state()
        mode = form.get("mode", "")
        if mode == "collwin":
            gui_state["collwin"] = {
                "filename": form.get("collwin_filename", ""),
                "input": form.get("collwin_input", ""),
                "output": form.get("collwin_output", ""),
                "start": form.get("collwin_start", ""),
                "end": form.get("collwin_end", ""),
                "grid": not bool(form.get("collwin_physical")),
                "no_geometry": bool(form.get("collwin_no_geometry")),
            }
        elif mode == "analytic":
            gui_state["analytic"] = {
                "sigma": form.get("analytic_sigma", ""),
                "charge": form.get("analytic_charge", ""),
                "half_separation": form.get("analytic_half_separation", ""),
                "nx": form.get("analytic_nx", ""),
                "ny": form.get("analytic_ny", ""),
                "nz": form.get("analytic_nz", ""),
                "output": form.get("analytic_output", ""),
            }
        elif mode == "compare":
            gui_state["compare"] = {
                "h5_path": form.get("compare_h5_path", ""),
                "step": form.get("compare_step", ""),
                "start": form.get("compare_start", ""),
                "end": form.get("compare_end", ""),
                "sigma": form.get("compare_sigma", ""),
                "state": form.get("compare_state", ""),
                "component": form.get("compare_component", ""),
                "output": form.get("compare_output", ""),
            }
        save_gui_state(gui_state)

    class Handler(BaseHTTPRequestHandler):
        def _send_html(self, html: str, status: HTTPStatus = HTTPStatus.OK) -> None:
            payload = html.encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)

        def do_GET(self) -> None:
            parsed = urlparse(self.path)
            if parsed.path == "/":
                self._send_html(build_web_gui_page())
                return
            if parsed.path == "/preview":
                query = parse_qs(parsed.query)
                raw_path = query.get("path", [""])[-1]
                path = Path(raw_path)
                if not path.exists() or not path.is_file():
                    self.send_error(HTTPStatus.NOT_FOUND)
                    return
                payload = path.read_bytes()
                mime_type = mimetypes.guess_type(str(path))[0] or "application/octet-stream"
                self.send_response(HTTPStatus.OK)
                self.send_header("Content-Type", mime_type)
                self.send_header("Content-Length", str(len(payload)))
                self.end_headers()
                self.wfile.write(payload)
                return
            self.send_error(HTTPStatus.NOT_FOUND)

        def do_POST(self) -> None:
            if self.path != "/run":
                self.send_error(HTTPStatus.NOT_FOUND)
                return

            length = int(self.headers.get("Content-Length", "0"))
            raw_body = self.rfile.read(length).decode("utf-8")
            form = {key: values[-1] for key, values in parse_qs(raw_body).items()}
            mode = form.get("mode", "")
            previews: list[tuple[str, Path]] = []

            try:
                if mode == "collwin":
                    output = Path(form["collwin_output"]) if form.get("collwin_output") else temp_output_path(
                        "beambeam-collwin", ".png"
                    )
                    run_collision_window(
                        Path(form["collwin_filename"]),
                        output=output,
                        start=int(form["collwin_start"]) if form.get("collwin_start") else None,
                        end=int(form["collwin_end"]) if form.get("collwin_end") else None,
                        grid=not bool(form.get("collwin_physical")),
                        no_geometry=bool(form.get("collwin_no_geometry")),
                        input_file=Path(form["collwin_input"]) if form.get("collwin_input") else None,
                        force_agg=True,
                    )
                    status = f"Completed collision-window render.\nOutput: {output}"
                    previews.append(("Collision window", output))
                elif mode == "analytic":
                    output = Path(form["analytic_output"]) if form.get("analytic_output") else temp_output_path(
                        "beambeam-analytic", ".png"
                    )
                    run_manufactured_analytic(
                        sigma=float(form["analytic_sigma"]),
                        charge=float(form["analytic_charge"]),
                        half_separation=float(form["analytic_half_separation"]),
                        nx=int(form["analytic_nx"]),
                        ny=int(form["analytic_ny"]),
                        nz=int(form["analytic_nz"]),
                        output=output,
                        force_agg=True,
                    )
                    status = f"Completed manufactured analytic render.\nOutput: {output}"
                    previews.append(("Manufactured analytic", output))
                elif mode == "compare":
                    state = form.get("compare_state") or None
                    if form.get("compare_start") or form.get("compare_end"):
                        if not (form.get("compare_start") and form.get("compare_end")):
                            raise ValueError("both start and end are required for a compare movie")
                        component = form.get("compare_component") or "rho+phi"
                        output = Path(form["compare_output"]) if form.get("compare_output") else temp_output_path(
                            "beambeam-compare-movie", ".gif"
                        )
                        run_manufactured_compare_movie(
                            h5_path=Path(form["compare_h5_path"]),
                            start=int(form["compare_start"]),
                            end=int(form["compare_end"]),
                            sigma=float(form["compare_sigma"]),
                            state=state,
                            comparison=component,
                            output=output,
                            force_agg=True,
                        )
                        status = f"Completed manufactured compare movie.\nOutput: {output}"
                        previews.append(("Manufactured compare movie", output))
                    else:
                        component = form.get("compare_component") or "rho+phi"
                        output = Path(form["compare_output"]) if form.get("compare_output") else temp_output_path(
                            "beambeam-compare", ".png"
                        )
                        run_manufactured_compare(
                            h5_path=Path(form["compare_h5_path"]),
                            step=int(form["compare_step"]),
                            sigma=float(form["compare_sigma"]),
                            state=state,
                            comparison=component,
                            output=output,
                            force_agg=True,
                        )
                        component_tuple = normalize_compare_component(component)
                        status = f"Completed manufactured compare.\nOutput: {output}"
                        previews.append(("Manufactured compare", output))
                        if component_tuple in (("rho", "phi"), ("rho", "phi", "Ez")):
                            status += (
                                "\nDerived outputs: "
                                f"{output.with_name(output.stem + '-rho-z-axis' + output.suffix)}, "
                                f"{output.with_name(output.stem + '-phi-z-axis' + output.suffix)}, "
                                f"{output.with_name(output.stem + '-ez-z-axis' + output.suffix)}, "
                                f"{output.with_name(output.stem + '-ex-x-axis' + output.suffix)}, "
                                f"{output.with_name(output.stem + '-ey-y-axis' + output.suffix)}"
                            )
                            previews.extend(
                                [
                                    ("rho(z)", output.with_name(output.stem + "-rho-z-axis" + output.suffix)),
                                    ("phi(z)", output.with_name(output.stem + "-phi-z-axis" + output.suffix)),
                                    ("Ez(z)", output.with_name(output.stem + "-ez-z-axis" + output.suffix)),
                                    ("Ex(x)", output.with_name(output.stem + "-ex-x-axis" + output.suffix)),
                                    ("Ey(y)", output.with_name(output.stem + "-ey-y-axis" + output.suffix)),
                                ]
                            )
                else:
                    raise ValueError(f"unknown mode: {mode}")
                persist_web_form(form)
            except Exception:
                status = traceback.format_exc()

            self._send_html(build_web_gui_page(status=status, values=form, previews=previews))

        def log_message(self, fmt: str, *args) -> None:
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    host, port = server.server_address
    url = f"http://{host}:{port}/"
    if reason is None:
        print(f"Tk GUI unavailable for {sys.executable}; falling back to local browser UI.")
    else:
        print(
            f"Tk GUI unavailable for {sys.executable}: "
            f"{type(reason).__name__}: {reason}; falling back to local browser UI."
        )
    print(f"Open: {url}")
    try:
        webbrowser.open(url)
    except Exception:
        pass
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


def create_gui() -> None:
    if sys.platform == "darwin":
        try:
            import tkinter as tk

            if float(tk.TkVersion) < 8.6:
                raise RuntimeError(f"Tcl/Tk {tk.TkVersion} is too old for a reliable macOS GUI")
            create_tk_gui()
            return
        except Exception as exc:
            run_web_gui(exc)
            return

    try:
        create_tk_gui()
    except Exception as exc:
        run_web_gui(exc)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combined BeamBeam analysis front-end for collision-window and manufactured-solution tools."
    )
    subparsers = parser.add_subparsers(dest="command")

    subparsers.add_parser("gui", help="Open the small Tk GUI.")

    collwin_parser = subparsers.add_parser("collwin", help="Run the collision-window viewer.")
    collwin_parser.add_argument("filename", type=Path)
    collwin_parser.add_argument("--output", type=Path)
    collwin_parser.add_argument("--start", type=int)
    collwin_parser.add_argument("--end", type=int)
    collwin_parser.add_argument("--physical", action="store_true")
    collwin_parser.add_argument("--no-geometry", action="store_true")
    collwin_parser.add_argument("--input", type=Path)

    analytic_parser = subparsers.add_parser("analytic", help="Run the manufactured analytic view.")
    analytic_parser.add_argument("--sigma", type=float, default=1.0e-3)
    analytic_parser.add_argument("--charge", type=float, default=-1.0e-15)
    analytic_parser.add_argument("--half-separation", type=float, default=2.0e-3)
    analytic_parser.add_argument("--nx", type=int, default=64)
    analytic_parser.add_argument("--ny", type=int, default=64)
    analytic_parser.add_argument("--nz", type=int, default=64)
    analytic_parser.add_argument("--output", type=Path)

    pair_parser = subparsers.add_parser(
        "pair",
        help="Integrate a back-to-back e-/e+ pair in the manufactured boosted-Gaussian fields.",
    )
    pair_parser.add_argument("--sigma-m", type=float, default=0.6e-3)
    pair_parser.add_argument("--charge-C", type=float, default=-2.0027207925e-9)
    pair_parser.add_argument("--source-kinetic-MeV", type=float, default=245.0)
    pair_parser.add_argument("--half-separation-m", type=float, default=0.0)
    pair_parser.add_argument("--pair-energy-keV", type=float, default=313.0)
    pair_parser.add_argument("--pair-duration-ps", type=float, default=1050.0)
    pair_parser.add_argument("--pair-dt-ps", type=float, default=0.1)
    pair_parser.add_argument("--pair-fine-dt-ps", type=float, default=None)
    pair_parser.add_argument("--pair-fine-duration-ps", type=float, default=5.0)
    pair_parser.add_argument("--pair-output-dt-ps", type=float, default=0.1)
    pair_parser.add_argument(
        "--pair-direction",
        type=float,
        nargs=3,
        default=(1.0, 0.0, 0.0),
        metavar=("UX", "UY", "UZ"),
    )
    pair_parser.add_argument(
        "--pair-origin-m",
        type=float,
        nargs=3,
        default=(0.0, 0.0, 0.0),
        metavar=("X", "Y", "Z"),
    )
    pair_parser.add_argument("--field-span-mm", type=float, default=4.0)
    pair_parser.add_argument("--field-grid", type=int, default=121)
    pair_parser.add_argument("--beambeam-radius-m", type=float, default=0.15)
    pair_parser.add_argument("--beambeam-length-m", type=float, default=0.32)
    pair_parser.add_argument("--output-prefix", type=Path, default=Path("sandbox/data/boosted_gaussian_witness"))
    pair_parser.add_argument(
        "--write-pair-latex-table",
        type=Path,
        default=Path("sandbox/note/boosted_gaussian_witness_pair_kinematics_table.tex"),
    )

    compare_parser = subparsers.add_parser("compare", help="Run the manufactured-vs-OPALX compare view.")
    compare_parser.add_argument("h5_path", type=Path)
    compare_parser.add_argument("--step", type=int)
    compare_parser.add_argument("--start", type=int)
    compare_parser.add_argument("--end", type=int)
    compare_parser.add_argument("--sigma", type=float, default=1.0e-3)
    compare_parser.add_argument(
        "--comparison",
        default="rho+phi",
        help="Comparison kind: rho, phi, Ez, rho-z-axis, phi-z-axis, Ez-z-axis, Ex-x-axis, Ey-y-axis, rho+phi, or rho+phi+Ez. Movies require exactly one kind.",
    )
    compare_parser.add_argument(
        "--state",
        help="State selector, e.g. 'beambeam tracking' or 'normal tracking'. Raw HDF names are also accepted.",
    )
    compare_parser.add_argument("--output", type=Path)

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.command in (None, "gui"):
        create_gui()
        return 0

    if args.command == "collwin":
        run_collision_window(
            args.filename,
            output=args.output,
            start=args.start,
            end=args.end,
            grid=not args.physical,
            no_geometry=args.no_geometry,
            input_file=args.input,
        )
        return 0

    if args.command == "analytic":
        run_manufactured_analytic(
            sigma=args.sigma,
            charge=args.charge,
            half_separation=args.half_separation,
            nx=args.nx,
            ny=args.ny,
            nz=args.nz,
            output=args.output,
        )
        return 0

    if args.command == "pair":
        run_pair_in_fields(
            sigma_m=args.sigma_m,
            charge_C=args.charge_C,
            source_kinetic_MeV=args.source_kinetic_MeV,
            half_separation_m=args.half_separation_m,
            pair_energy_keV=args.pair_energy_keV,
            pair_duration_ps=args.pair_duration_ps,
            pair_dt_ps=args.pair_dt_ps,
            pair_fine_dt_ps=args.pair_fine_dt_ps,
            pair_fine_duration_ps=args.pair_fine_duration_ps,
            pair_output_dt_ps=args.pair_output_dt_ps,
            pair_direction=tuple(args.pair_direction),
            pair_origin_m=tuple(args.pair_origin_m),
            field_span_mm=args.field_span_mm,
            field_grid=args.field_grid,
            beambeam_radius_m=args.beambeam_radius_m,
            beambeam_length_m=args.beambeam_length_m,
            output_prefix=args.output_prefix,
            write_pair_latex_table=args.write_pair_latex_table,
        )
        return 0

    if args.command == "compare":
        if args.step is not None and (args.start is not None or args.end is not None):
            raise SystemExit("use either --step or --start/--end for compare")
        if args.step is None and (args.start is None or args.end is None):
            raise SystemExit("compare requires either --step or both --start and --end")
        if args.start is not None and args.end is not None:
            if args.output is None:
                raise SystemExit("compare movie requires --output")
            run_manufactured_compare_movie(
                h5_path=args.h5_path,
                start=args.start,
                end=args.end,
                sigma=args.sigma,
                state=args.state,
                comparison=args.comparison,
                output=args.output,
            )
            return 0
        run_manufactured_compare(
            h5_path=args.h5_path,
            step=args.step,
            sigma=args.sigma,
            state=args.state,
            comparison=args.comparison,
            output=args.output,
        )
        return 0

    raise SystemExit(f"unknown command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
