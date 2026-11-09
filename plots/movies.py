#!/usr/bin/env python3
"""Plot and animate three-layer PEQUOD snapshots from a serial or MPI run."""

import argparse
from pathlib import Path
import re
import struct
from subprocess import CalledProcessError
from tempfile import TemporaryDirectory

import matplotlib
import numpy as np
from matplotlib.animation import FFMpegWriter

matplotlib.use("Agg")
from matplotlib import pyplot as plt


SNAPSHOT_NAME = re.compile(r"(qg3_\d+_s)\.(\d+)(?:\.(\d+))?")
FIELD_NAMES = ("psi1", "psi2", "psi3", "zeta1", "zeta2", "zeta3")
ROW_LABELS = ("Streamfunction", "Relative potential vorticity")


def find_snapshots(run_dir):
    snapshots = {}
    prefixes = set()
    for path in run_dir.iterdir():
        match = SNAPSHOT_NAME.fullmatch(path.name)
        if match and path.is_file():
            prefix, index, _ = match.groups()
            prefixes.add(prefix)
            snapshots.setdefault(int(index), []).append(path)

    if not snapshots:
        raise ValueError(f"No three-layer snapshot headers found in {run_dir}")
    if len(prefixes) != 1:
        raise ValueError(f"Found multiple snapshot grids in {run_dir}")
    return [sorted(headers) for _, headers in sorted(snapshots.items())]


def read_record(stream, count):
    # The Ubuntu build writes little-endian float64 values between 4-byte sizes.
    size = count * 8
    marker = struct.pack("<I", size)
    if stream.read(4) != marker:
        raise ValueError(f"{stream.name}: expected a record with {count} values")
    data = stream.read(size)
    if len(data) != size or stream.read(4) != marker:
        raise ValueError(f"{stream.name}: incomplete or invalid record")
    return np.frombuffer(data, dtype="<f8")


def read_tile(data_path, descriptions, tile):
    ni, nj = tile.shape[1:]
    with data_path.open("rb") as stream:
        times = (read_record(stream, 1)[0], read_record(stream, 1)[0])
        for field, name, description in zip(
            tile, FIELD_NAMES, descriptions, strict=True
        ):
            if description.split() == [name, "ii", "jj"]:
                levels = 1
            elif description.split() == [name, "ii", "jj", "leapfrog"]:
                levels = 2
            else:
                raise ValueError(f"{data_path}: unexpected field description {description!r}")
            values = read_record(stream, ni * nj * levels)
            # Centred runs store the previous level first and the current level last.
            field[:] = values[-ni * nj :].reshape((ni, nj), order="F")
        if stream.read(1):
            raise ValueError(f"{stream.name}: unexpected data after the snapshot")
    return times


def read_snapshot(headers):
    grid = None
    times = None
    seen_ranks = set()
    for header in headers:
        lines = header.read_text().splitlines()
        if len(lines) != 12 or lines[4:6] != ["model_time_day", "total_time_day"]:
            raise ValueError(f"{header}: unexpected snapshot header")
        nprocs, _, _, rank, _, _ = map(int, lines[0].split())
        ni, nj, layers = map(int, lines[1].split())
        i0, j0, local_ni, local_nj = map(int, lines[2].split())
        spacing_cm, _ = map(float, lines[3].split())
        spacing_km = spacing_cm / 100_000
        if layers != 3 or ni < 1 or nj < 1 or not 0 < spacing_km < np.inf:
            raise ValueError(f"{header}: invalid three-layer grid metadata")
        if len(headers) != nprocs:
            raise ValueError(f"{header}: expected {nprocs} rank headers, found {len(headers)}")
        rank_suffix = SNAPSHOT_NAME.fullmatch(header.name)[3]
        expected_rank = 0 if rank_suffix is None else int(rank_suffix)
        if (
            (rank_suffix is None) != (nprocs == 1)
            or rank != expected_rank
            or not 0 <= rank < nprocs
            or rank in seen_ranks
        ):
            raise ValueError(f"{header}: unexpected or duplicate rank {rank}")
        seen_ranks.add(rank)

        tile_grid = (ni, nj, spacing_km, nprocs)
        if grid is None:
            grid = tile_grid
            fields = np.empty((len(FIELD_NAMES), ni, nj))
            filled = np.zeros((ni, nj), dtype=bool)
        elif grid != tile_grid:
            raise ValueError(f"{header}: global grid metadata differs between ranks")

        i1, j1 = i0 + local_ni, j0 + local_nj
        if not (0 <= i0 < i1 <= ni and 0 <= j0 < j1 <= nj):
            raise ValueError(f"{header}: tile lies outside the global grid")
        rows, columns = slice(i0, i1), slice(j0, j1)
        covered = filled[rows, columns]
        if covered.any():
            raise ValueError(f"{header}: tile overlaps another rank")

        if rank_suffix is None:
            data_path = header.with_name(f"{header.name}.d")
        else:
            data_path = header.with_name(f"{header.stem}.d.{rank_suffix}")
        tile_times = read_tile(data_path, lines[6:], fields[:, rows, columns])
        if not np.isfinite(tile_times).all() or (times is not None and times != tile_times):
            raise ValueError(f"{data_path}: invalid or inconsistent snapshot times")
        times = tile_times
        covered[:] = True

    if not filled.all():
        raise ValueError(f"{headers[0]}: rank tiles do not cover the full grid")
    return times[1], fields.reshape(2, 3, ni, nj), spacing_km


def make_outputs(run_dir, snapshots, fps):
    limits = np.zeros(2)
    grid = None
    for snapshot in snapshots:
        day, fields, spacing_km = read_snapshot(snapshot)
        frame_grid = (fields.shape[2:], spacing_km)
        if grid is not None and frame_grid != grid:
            raise ValueError(f"{snapshot[0]}: grid changed between snapshots")
        grid = frame_grid
        limits = np.maximum(limits, np.max(np.abs(fields), axis=(1, 2, 3)))
        if not np.isfinite(limits).all():
            raise ValueError(f"{snapshot[0]}: contains non-finite field values")
    limits[limits == 0] = 1

    (ni, nj), spacing_km = grid
    west, east = 0, (nj - 1) * spacing_km
    south, north = 0, (ni - 1) * spacing_km
    half_cell = spacing_km / 2
    extent = (west - half_cell, east + half_cell, south - half_cell, north + half_cell)

    fig, axes = plt.subplots(2, 3, figsize=(12, 8), sharex=True, sharey=True,
                            layout="constrained")
    images = []
    for row, (label, limit) in enumerate(zip(ROW_LABELS, limits, strict=True)):
        for col, ax in enumerate(axes[row]):
            image = ax.imshow(
                fields[row, col], cmap="RdBu_r", vmin=-limit, vmax=limit,
                origin="upper", extent=extent, interpolation="nearest",
            )
            images.append(image)
            ax.set(xlim=(west, east), ylim=(south, north))
            if row == 0:
                ax.set_title(f"Layer {col + 1}")
        fig.colorbar(images[-1], ax=axes[row], shrink=0.83,
                     label=f"{label} (nondimensional)")
    fig.supxlabel("Eastward distance (km)")
    fig.supylabel("Northward distance (km)")
    title = fig.suptitle(f"PEQUOD | Model day {day:.2f}")
    fig.canvas.draw()
    fig.set_layout_engine("none")

    writer = FFMpegWriter(
        fps=fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart"]
    )
    # Keep existing outputs intact if reading a frame or encoding the video fails.
    try:
        with TemporaryDirectory(dir=run_dir) as temp:
            png = Path(temp) / "fields.png"
            mp4 = Path(temp) / "fields.mp4"
            fig.savefig(png, dpi=150)
            with writer.saving(fig, mp4, dpi=100):
                for snapshot in snapshots:
                    day, fields, _ = read_snapshot(snapshot)
                    title.set_text(f"PEQUOD | Model day {day:.2f}")
                    for image, field in zip(images, fields.reshape(6, ni, nj), strict=True):
                        image.set_data(field)
                    writer.grab_frame()
            for output in (png, mp4):
                output.replace(run_dir / output.name)
    finally:
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "run_dir", type=Path, help="Path to a run-serial.XXXXXX or run-parallel.XXXXXX folder"
    )
    parser.add_argument("--fps", type=int, default=20, help="Frames per second (default: 20)")
    args = parser.parse_args()
    if args.fps < 1:
        parser.error("--fps must be positive")
    if not FFMpegWriter.isAvailable():
        parser.error("FFmpeg is required to write the MP4")

    run_dir = args.run_dir.expanduser().resolve()
    try:
        snapshots = find_snapshots(run_dir)
        print(f"Rendering {len(snapshots)} snapshots from {run_dir}", flush=True)
        make_outputs(run_dir, snapshots, args.fps)
    except (OSError, ValueError, CalledProcessError) as exc:
        parser.exit(1, f"{parser.prog}: {exc}\n")
    print(f"Saved {run_dir / 'fields.png'}")
    print(f"Saved {run_dir / 'fields.mp4'}")


if __name__ == "__main__":
    main()
