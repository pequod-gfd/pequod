#!/usr/bin/env python3
"""Plot and animate three-layer PEQUOD snapshots from a serial or MPI run."""

import argparse
from pathlib import Path
import re
import runpy
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


def field_levels(header, descriptions):
    levels = []
    for name, description in zip(FIELD_NAMES, descriptions, strict=True):
        dimensions = description.split()
        if dimensions == [name, "ii", "jj"]:
            levels.append(1)
        elif dimensions == [name, "ii", "jj", "leapfrog"]:
            levels.append(2)
        else:
            raise ValueError(f"{header}: unexpected field description {description!r}")
    return levels


def read_tile(data_path, levels, tile, offset):
    ni, nj = tile.shape[1:]
    count = ni * nj
    with data_path.open("rb") as stream:
        stream.seek(offset)
        times = (read_record(stream, 1)[0], read_record(stream, 1)[0])
        for field, level in zip(tile, levels, strict=True):
            values = read_record(stream, count * level)
            # Centred runs store the previous level first and the current level last.
            field[:] = values[-count:].reshape((ni, nj), order="F")
    return times


def read_dump(headers):
    grid = None
    frame_count = None
    seen_ranks = set()
    tiles = []
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
        levels = field_levels(header, lines[6:])
        # Two scalar time records, then one record for each field.
        frame_size = 32 + sum(8 + 8 * local_ni * local_nj * level for level in levels)
        count, remainder = divmod(data_path.stat().st_size, frame_size)
        if not count or remainder:
            raise ValueError(f"{data_path}: incomplete snapshot data")
        if frame_count is not None and frame_count != count:
            raise ValueError(f"{data_path}: ranks have different numbers of snapshots")
        frame_count = count
        tiles.append((data_path, rows, columns, levels, frame_size))
        covered[:] = True

    if not filled.all():
        raise ValueError(f"{headers[0]}: rank tiles do not cover the full grid")
    for index in range(frame_count):
        fields = np.empty((len(FIELD_NAMES), ni, nj))
        times = None
        for data_path, rows, columns, levels, frame_size in tiles:
            tile_times = read_tile(data_path, levels, fields[:, rows, columns], index * frame_size)
            if not np.isfinite(tile_times).all() or (times is not None and times != tile_times):
                raise ValueError(f"{data_path}: invalid or inconsistent snapshot times")
            times = tile_times
        yield times[1], fields.reshape(2, 3, ni, nj), spacing_km


def read_frames(snapshots):
    for headers in snapshots:
        for day, fields, spacing_km in read_dump(headers):
            yield headers[0], day, fields, spacing_km


def read_cape_mask(run_dir, shape, spacing_km):
    path = run_dir / "cape_geometry.dat"
    if not path.is_file():
        return False
    geometry_reader = runpy.run_path(str(Path(__file__).resolve().parents[1] / "scripts/qg_irregular"))
    geometry = geometry_reader["Geometry"].read(path)
    grid = geometry.grid
    if ((grid.ni, grid.nj) != shape or grid.units != "cm"
            or not np.isclose(grid.h / 100_000, spacing_km, rtol=1e-12, atol=0)):
        raise ValueError(f"{path}: geometry grid does not match the snapshots")
    return geometry.nodes[:, :, 0] == geometry_reader["SOLID"]


def make_outputs(run_dir, snapshots, fps):
    limits = np.zeros(2)
    grid = None
    for count, (source, day, fields, spacing_km) in enumerate(read_frames(snapshots), 1):
        frame_grid = (fields.shape[2:], spacing_km)
        if grid is not None and frame_grid != grid:
            raise ValueError(f"{source}: grid changed between snapshots")
        grid = frame_grid
        limits = np.maximum(limits, np.max(np.abs(fields), axis=(1, 2, 3)))
        if not np.isfinite(limits).all():
            raise ValueError(f"{source}: contains non-finite field values")
    limits[limits == 0] = 1
    print(f"Rendering {count} snapshots", flush=True)

    (ni, nj), spacing_km = grid
    land = read_cape_mask(run_dir, (ni, nj), spacing_km)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("0.8")
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
                np.ma.array(fields[row, col], mask=land), cmap=cmap, vmin=-limit, vmax=limit,
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
                for _, day, fields, _ in read_frames(snapshots):
                    title.set_text(f"PEQUOD | Model day {day:.2f}")
                    for image, field in zip(images, fields.reshape(6, ni, nj), strict=True):
                        image.set_data(np.ma.array(field, mask=land))
                    writer.grab_frame()
            for output in (png, mp4):
                output.replace(run_dir / output.name)
    finally:
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "run_dir", type=Path, help="Path to a rectangular or cape qg3 run folder"
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
        print(f"Reading {len(snapshots)} snapshot files or rank groups from {run_dir}", flush=True)
        make_outputs(run_dir, snapshots, args.fps)
    except (OSError, ValueError, CalledProcessError) as exc:
        parser.exit(1, f"{parser.prog}: {exc}\n")
    print(f"Saved {run_dir / 'fields.png'}")
    print(f"Saved {run_dir / 'fields.mp4'}")


if __name__ == "__main__":
    main()
