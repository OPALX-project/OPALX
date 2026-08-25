from __future__ import annotations

import argparse
import glob
import math
import re
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


HEADER_RE = re.compile(r"^#\s*([^=]+)=(.*)$")
VECTOR_METADATA_RE = re.compile(r"([A-Za-z0-9_]+)\s*=\s*(\([^)]*\)|\[[^\]]*\])")
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = SCRIPT_DIR / "data"
VECTOR_METADATA_KEYS = (
    "mesh_origin",
    "mesh_spacing",
    "field_domain_rmin",
    "field_domain_rmax",
    "bunch_bounds_rmin",
    "bunch_bounds_rmax",
)


def resolve_paths(pattern: str) -> list[Path]:
    # pathlib.Path.glob() rejects absolute patterns; use glob so both relative
    # and absolute dump patterns work from any launch directory.
    return sorted(Path(path) for path in glob.glob(pattern))


def load_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    font_candidates = [
        "/System/Library/Fonts/Menlo.ttc",
        "/System/Library/Fonts/Courier.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
    ]
    for path in font_candidates:
        p = Path(path)
        if p.exists():
            return ImageFont.truetype(str(p), size)
    return ImageFont.load_default()


def parse_vector(raw: str) -> tuple[float, ...]:
    matches = re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", raw)
    return tuple(float(part) for part in matches)


def extract_metadata_z_bounds(
    metadata_vectors: dict[str, tuple[float, ...]],
    rmin_key: str,
    rmax_key: str,
) -> tuple[float | None, float | None]:
    rmin = metadata_vectors.get(rmin_key, ())
    rmax = metadata_vectors.get(rmax_key, ())

    zmin = rmin[2] if len(rmin) > 2 else None
    zmax = rmax[2] if len(rmax) > 2 else None
    return zmin, zmax


def parse_metadata_line(
    line: str,
    metadata: dict[str, str],
    metadata_vectors: dict[str, tuple[float, ...]],
) -> None:
    body = line[1:].strip()

    for key, raw_value in VECTOR_METADATA_RE.findall(body):
        value = raw_value.strip()
        metadata[key] = value
        if key in VECTOR_METADATA_KEYS:
            metadata_vectors[key] = parse_vector(value)

    match = HEADER_RE.match(line)
    if match:
        key = match.group(1).strip()
        value = match.group(2).strip()
        if key not in metadata:
            metadata[key] = value


def parse_dump(path: Path) -> dict:
    metadata: dict[str, str] = {}
    metadata_vectors: dict[str, tuple[float, ...]] = {}
    rows: list[tuple[float, float, float, float]] = []

    with path.open() as stream:
        for line in stream:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                parse_metadata_line(line, metadata, metadata_vectors)
                continue

            parts = line.split()
            if len(parts) < 7:
                continue

            x = float(parts[3])
            y = float(parts[4])
            z = float(parts[5])
            rho = float(parts[6])
            rows.append((x, y, z, rho))

    if not rows:
        raise ValueError(f"no rho samples found in {path}")

    xs = sorted({row[0] for row in rows})
    ys = sorted({row[1] for row in rows})
    zs = sorted({row[2] for row in rows})

    dz = zs[1] - zs[0] if len(zs) > 1 else 1.0
    sigma = {(x, y): 0.0 for y in ys for x in xs}
    for x, y, _z, rho in rows:
        sigma[(x, y)] += rho * dz

    return {
        "path": path,
        "metadata": metadata,
        "metadata_vectors": metadata_vectors,
        "xs": xs,
        "ys": ys,
        "sigma": sigma,
        "dz": dz,
        "max_abs_sigma": max(abs(value) for value in sigma.values()),
    }


def color_for_value(value: float, vmax: float) -> tuple[int, int, int]:
    if vmax <= 0.0:
        return (255, 255, 255)

    soft = max(vmax / 20.0, 1.0e-30)
    u = math.asinh(value / soft) / math.asinh(vmax / soft)
    u = max(-1.0, min(1.0, u))

    if u >= 0.0:
        mix = u
        red = 255
        green = int(round(255 * (1.0 - mix)))
        blue = int(round(255 * (1.0 - mix)))
    else:
        mix = -u
        red = int(round(255 * (1.0 - mix)))
        green = int(round(255 * (1.0 - mix)))
        blue = 255

    return (red, green, blue)


def draw_colorbar(
    draw: ImageDraw.ImageDraw,
    x0: int,
    y0: int,
    width: int,
    height: int,
    vmax: float,
    font: ImageFont.FreeTypeFont | ImageFont.ImageFont,
) -> None:
    for row in range(height):
        frac = 1.0 - row / max(height - 1, 1)
        value = (2.0 * frac - 1.0) * vmax
        color = color_for_value(value, vmax)
        draw.line((x0, y0 + row, x0 + width, y0 + row), fill=color)

    draw.rectangle((x0, y0, x0 + width, y0 + height), outline="black", width=1)
    draw.text((x0 + width + 8, y0 - 8), f"+{vmax:.3e}", fill="black", font=font)
    draw.text((x0 + width + 8, y0 + height // 2 - 8), "0", fill="black", font=font)
    draw.text((x0 + width + 8, y0 + height - 16), f"-{vmax:.3e}", fill="black", font=font)


def render_frame(
    frame: dict,
    vmax: float,
    cell_size: int,
    grid_width: int,
    grid_height: int,
    font: ImageFont.FreeTypeFont | ImageFont.ImageFont,
    small_font: ImageFont.FreeTypeFont | ImageFont.ImageFont,
) -> Image.Image:
    margin = 24
    title_height = 90
    legend_width = 130

    width = margin * 3 + grid_width * cell_size + legend_width
    height = margin * 2 + title_height + grid_height * cell_size

    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    metadata = frame["metadata"]
    metadata_vectors = frame["metadata_vectors"]
    phase = metadata.get("phase", "unknown")
    step = metadata.get("global_step", "?")
    bunch_zmin, bunch_zmax = extract_metadata_z_bounds(
        metadata_vectors, "bunch_bounds_rmin", "bunch_bounds_rmax"
    )

    draw.text((margin, margin), "Rho XY Projection", fill="black", font=font)
    draw.text(
        (margin, margin + 28),
        f"integrated quantity=sigma_xy [C/m^2]",
        fill="black",
        font=small_font,
    )
    #draw.text((margin, margin + 50), f"rmin={field_rmin}", fill="black", font=small_font)
    #draw.text((margin, margin + 66), f"rmax={field_rmax}", fill="black", font=small_font)
    zmin_label = f"{bunch_zmin:.6e}" if bunch_zmin is not None else "n/a"
    zmax_label = f"{bunch_zmax:.6e}" if bunch_zmax is not None else "n/a"
    draw.text(
        (margin, margin + 50),
        f"z: {zmin_label} .. {zmax_label}",
        fill="black",
        font=small_font,
    )
    #draw.text((margin, margin + 66), f"bunch rmax={bunch_rmax}", fill="black", font=small_font)

    grid_x0 = margin
    grid_y0 = margin + title_height

    xs = frame["xs"]
    ys = frame["ys"]
    sigma = frame["sigma"]

    for row, y in enumerate(reversed(ys)):
        for col, x in enumerate(xs):
            value = sigma[(x, y)]
            color = color_for_value(value, vmax)
            x0 = grid_x0 + col * cell_size
            y0 = grid_y0 + row * cell_size
            draw.rectangle((x0, y0, x0 + cell_size, y0 + cell_size), fill=color)

    draw.rectangle(
        (grid_x0, grid_y0, grid_x0 + len(xs) * cell_size, grid_y0 + len(ys) * cell_size),
        outline="black",
        width=1,
    )

    draw.text(
        (grid_x0, grid_y0 - 5 +len(ys) * cell_size + 8),
        f"x: {1000*xs[0]:.3} .. {1000*xs[-1]:.3} mm",
        fill="black",
        font=small_font,
    )
    draw.text(
        (grid_x0 + 230, grid_y0 - 5 + len(ys) * cell_size + 8),
        f"y: {1000*ys[0]:.3} .. {1000*ys[-1]:.3} mm",
        fill="black",
        font=small_font,
    )

    legend_x0 = grid_x0 + grid_width * cell_size + margin
    legend_y0 = grid_y0
    draw_colorbar(draw, legend_x0, legend_y0, 24, grid_height * cell_size, vmax, small_font)

    return img


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Make a GIF from OPALX rho debug dumps using an (x,y) projection."
    )
    parser.add_argument(
        "--pattern",
        default=str(DEFAULT_DATA_DIR / "*rho_debug*.dat"),
        help="Glob for rho debug dump files.",
    )
    parser.add_argument(
        "--frames-dir",
        default=str(DEFAULT_DATA_DIR / "rho_xy_frames"),
        help="Directory for rendered PNG frames.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_DATA_DIR / "rho_xy.gif"),
        help="Output GIF filename.",
    )
    parser.add_argument(
        "--duration",
        type=int,
        default=120,
        help="GIF frame duration in milliseconds.",
    )
    args = parser.parse_args()

    paths = resolve_paths(args.pattern)
    if not paths:
        raise SystemExit(f"no files matched pattern: {args.pattern}")

    frames = [parse_dump(path) for path in paths]
    vmax = max(frame["max_abs_sigma"] for frame in frames)

    max_nx = max(len(frame["xs"]) for frame in frames)
    max_ny = max(len(frame["ys"]) for frame in frames)
    cell_size = max(4, min(18, 720 // max(max_nx, max_ny, 1)))

    font = load_font(24)
    small_font = load_font(16)

    frames_dir = Path(args.frames_dir)
    frames_dir.mkdir(parents=True, exist_ok=True)

    rendered = []
    for index, frame in enumerate(frames):
        image = render_frame(frame, vmax, cell_size, max_nx, max_ny, font, small_font)
        image_path = frames_dir / f"frame_{index:04d}.png"
        image.save(image_path)
        rendered.append(image)

    rendered[0].save(
        args.output,
        save_all=True,
        append_images=rendered[1:],
        duration=args.duration,
        loop=0,
    )
    print(f"wrote {len(rendered)} frames to {args.output}")


if __name__ == "__main__":
    main()
