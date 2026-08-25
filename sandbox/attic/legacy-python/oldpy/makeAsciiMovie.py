from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = SCRIPT_DIR / "data"


def load_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    font_candidates = [
        "/System/Library/Fonts/Menlo.ttc",
        "/System/Library/Fonts/Courier.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
    ]
    for path in font_candidates:
        font_path = Path(path)
        if font_path.exists():
            return ImageFont.truetype(str(font_path), size)
    return ImageFont.load_default()


def parse_frames(path: Path) -> list[dict[str, object]]:
    lines = path.read_text().splitlines()
    frames: list[dict[str, object]] = []
    index = 0

    while index < len(lines):
        line = lines[index].strip()
        if not line:
            index += 1
            continue
        if not line.startswith("FRAME "):
            raise ValueError(f"unexpected line {index + 1}: {lines[index]!r}")

        frame_id = int(line.split()[1])
        frame_lines: list[str] = []
        index += 1

        while index < len(lines) and lines[index].strip():
            frame_lines.append(lines[index].rstrip("\n"))
            index += 1

        if len(frame_lines) < 2:
            raise ValueError(f"incomplete frame starting at line {index + 1}")

        frames.append({"id": frame_id, "lines": frame_lines})

        while index < len(lines) and not lines[index].strip():
            index += 1

    return frames


def render_frame(
    frame: dict[str, object],
    font: ImageFont.FreeTypeFont | ImageFont.ImageFont,
    small_font: ImageFont.FreeTypeFont | ImageFont.ImageFont,
    cell_width: int,
    cell_height: int,
    canvas_width: int,
    canvas_height: int,
    max_line_count: int,
) -> Image.Image:
    frame_lines = [str(line) for line in frame["lines"]]
    frame_id = int(frame["id"])

    margin_x = 28
    margin_y = 22
    title_height = 30

    image = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(image)

    draw.text((margin_x, margin_y - 4), f"Collision ASCII Frame {frame_id}", fill="black", font=small_font)
    for line_index, line in enumerate(frame_lines):
        draw.text(
            (margin_x, margin_y + title_height + line_index * cell_height),
            line,
            fill="black",
            font=font,
        )

    return image


def main() -> None:
    parser = argparse.ArgumentParser(description="Render collision ASCII frames to PNGs and a GIF.")
    parser.add_argument(
        "--input",
        default=str(DEFAULT_DATA_DIR / "collision_ascii_frames.txt"),
        help="Input collision ASCII frame dump.",
    )
    parser.add_argument(
        "--frames-dir",
        default=str(DEFAULT_DATA_DIR / "collision_ascii_png"),
        help="Directory for rendered PNG frames.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_DATA_DIR / "collision_ascii.gif"),
        help="Output GIF filename.",
    )
    parser.add_argument(
        "--duration",
        type=int,
        default=90,
        help="GIF frame duration in milliseconds.",
    )
    parser.add_argument(
        "--font-size",
        type=int,
        default=22,
        help="Monospace font size.",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    frames = parse_frames(input_path)
    if not frames:
        raise SystemExit(f"no frames found in {input_path}")

    font = load_font(args.font_size)
    small_font = load_font(max(12, args.font_size - 6))
    bbox = font.getbbox("M")
    cell_width = bbox[2] - bbox[0]
    cell_height = bbox[3] - bbox[1] + 8
    margin_x = 28
    margin_y = 22
    title_height = 30
    max_line_width = max(len(str(line)) for frame in frames for line in frame["lines"])
    max_line_count = max(len(frame["lines"]) for frame in frames)
    canvas_width = margin_x * 2 + max_line_width * cell_width
    canvas_height = margin_y * 2 + title_height + max_line_count * cell_height

    frames_dir = Path(args.frames_dir)
    frames_dir.mkdir(parents=True, exist_ok=True)
    for frame_path in frames_dir.glob("frame_*.png"):
        frame_path.unlink()

    rendered: list[Image.Image] = []
    for frame in frames:
        image = render_frame(
            frame,
            font,
            small_font,
            cell_width,
            cell_height,
            canvas_width,
            canvas_height,
            max_line_count,
        )
        frame_path = frames_dir / f"frame_{int(frame['id']):04d}.png"
        image.save(frame_path)
        rendered.append(image)

    rendered[0].save(
        args.output,
        save_all=True,
        append_images=rendered[1:],
        duration=args.duration,
        loop=0,
        disposal=2,
        optimize=False,
    )
    print(f"wrote {len(rendered)} frames to {args.output}")


if __name__ == "__main__":
    main()
