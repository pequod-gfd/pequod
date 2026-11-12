"""Inspect and merge experimental PEQUOD CMM half-stage checkpoint shards."""

import argparse
from dataclasses import dataclass
import hashlib
import math
import os
from pathlib import Path
import re
import runpy
import struct
import tempfile

import numpy as np


MAGIC = b"PEQUOD_CMM_CHECKPOINT_1\n"
HEADER = struct.Struct("<64s64s64s4iq3d2q")
CELL_FIELDS = ("pv_previous", "pv_half", "advection_rate", "beta_history",
               "volume_source_previous", "volume_source_current")
FACE_FIELDS = ("pv_previous", "pv_next", "transport_previous", "transport_current", "flux_next")
COMPONENT_FIELDS = ("circulation", "torque_previous", "torque_current", "streamfunction")


@dataclass(frozen=True)
class Stage:
    geometry: str
    operator: str
    closure: str
    ni: int
    nj: int
    layers: int
    islands: int
    step: int
    time: float
    dt_previous: float
    dt_current: float

    def validate(self):
        for digest in (self.geometry, self.operator, self.closure):
            if re.fullmatch(r"[0-9a-f]{64}", digest) is None:
                raise ValueError("expected lowercase SHA256 checkpoint identities")
        if min(self.ni, self.nj) < 3 or self.layers < 1 or self.islands < 0 or self.step < 1:
            raise ValueError("invalid checkpoint dimensions or step")
        if not all(map(math.isfinite, (self.time, self.dt_previous, self.dt_current))):
            raise ValueError("nonfinite checkpoint time metadata")
        if min(self.dt_previous, self.dt_current) <= 0:
            raise ValueError("checkpoint steps must be positive")


@dataclass
class Checkpoint:
    stage: Stage
    mass: np.ndarray
    components: np.ndarray
    cells: np.ndarray
    faces: np.ndarray


def records(fields, layers):
    return np.dtype([("id", "<i8"), ("values", "<f8", (len(fields), layers))])


def check_records(table, fields, layers, limit):
    if table.ndim != 1 or table.dtype != records(fields, layers):
        raise ValueError("checkpoint record layout mismatch")
    ids = table["id"]
    if np.any((ids < 1) | (ids > limit)) or np.any(ids[1:] <= ids[:-1]):
        raise ValueError("checkpoint IDs must be strictly increasing and in range")
    if not np.all(np.isfinite(table["values"])):
        raise ValueError("nonfinite checkpoint record")


def validate(checkpoint):
    s = checkpoint.stage
    s.validate()
    if checkpoint.mass.shape != (s.layers-1,) or checkpoint.components.shape != (
            s.islands+1, len(COMPONENT_FIELDS), s.layers):
        raise ValueError("checkpoint replicated state has incompatible shape")
    if not np.all(np.isfinite(checkpoint.mass)) or not np.all(np.isfinite(checkpoint.components)):
        raise ValueError("nonfinite checkpoint component or mass state")
    check_records(checkpoint.cells, CELL_FIELDS, s.layers, s.ni*s.nj)
    check_records(checkpoint.faces, FACE_FIELDS, s.layers, 2*s.ni*s.nj)


def write_checkpoint(path, checkpoint):
    """Publish one complete shard atomically after encoding a portable stream."""
    validate(checkpoint)
    path = Path(path)
    s = checkpoint.stage
    header = HEADER.pack(s.geometry.encode(), s.operator.encode(), s.closure.encode(), s.ni, s.nj,
                         s.layers, s.islands, s.step, s.time, s.dt_previous, s.dt_current,
                         len(checkpoint.cells), len(checkpoint.faces))
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}-", delete=False) as stream:
        temporary = Path(stream.name)
        try:
            digest = hashlib.sha256()
            for block in (MAGIC, header, np.asarray(checkpoint.mass, dtype="<f8"),
                          np.asarray(checkpoint.components, dtype="<f8"), checkpoint.cells, checkpoint.faces):
                data = block if isinstance(block, bytes) else block.tobytes(order="C")
                stream.write(data)
                digest.update(data)
            stream.write(digest.digest())
            stream.flush()
            os.fsync(stream.fileno())
            stream.close()
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)


def read_checkpoint(path, *, expected=None):
    path = Path(path)
    with path.open("rb") as stream:
        magic = stream.read(len(MAGIC))
        if magic != MAGIC:
            raise ValueError("unknown checkpoint format or version")
        header = stream.read(HEADER.size)
        if len(header) != HEADER.size:
            raise ValueError("truncated checkpoint header")
        geometry, operator, closure, *args, nc, nf = HEADER.unpack(header)
        s = Stage(geometry.decode("ascii"), operator.decode("ascii"), closure.decode("ascii"), *args)
        s.validate()
        if expected is not None and s != expected:
            raise ValueError("checkpoint stage or identity mismatch")
        cells = records(CELL_FIELDS, s.layers)
        faces = records(FACE_FIELDS, s.layers)
        shape = s.islands+1, len(COMPONENT_FIELDS), s.layers
        repeated = s.layers-1+math.prod(shape)
        size = repeated*8+nc*cells.itemsize+nf*faces.itemsize
        if not (0 <= nc <= s.ni*s.nj and 0 <= nf <= 2*s.ni*s.nj):
            raise ValueError("invalid checkpoint record counts")
        if os.fstat(stream.fileno()).st_size != len(MAGIC)+HEADER.size+size+32:
            raise ValueError("checkpoint length mismatch")
        data = stream.read(size)
        digest = hashlib.sha256(magic+header)
        digest.update(data)
        if digest.digest() != stream.read(32):
            raise ValueError("checkpoint SHA256 mismatch")
    mass_and_components = np.frombuffer(data, "<f8", count=repeated)
    mass = mass_and_components[:s.layers-1]
    components = mass_and_components[s.layers-1:].reshape(shape)
    result = Checkpoint(s, mass, components, np.frombuffer(data, cells, nc, repeated*8),
                        np.frombuffer(data, faces, nf, repeated*8+nc*cells.itemsize))
    validate(result)
    return result


def geometry_ids(path):
    """Read the canonical geometry and enumerate its wet cells and shared faces."""
    api = runpy.run_path(str(Path(__file__).with_name("qg_irregular")))
    geom = api["Geometry"].read(path)
    wedge = geom.nodes[:, :, 2]
    cell_ids = np.arange(1, wedge.size+1, dtype=np.int64).reshape(wedge.shape)
    cells = cell_ids[wedge != 0]
    south = 2*cell_ids[:-1][(wedge[:-1] & ((1 << 5) | (1 << 6))) != 0]-1
    east = 2*cell_ids[:, :-1][(wedge[:, :-1] & ((1 << 7) | 1)) != 0]
    return geom, cells, np.sort(np.r_[south, east])


def merge_checkpoints(checkpoints, cell_ids, face_ids):
    """Require exactly the supplied wet-cell and face sets from the geometry."""
    if not checkpoints:
        raise ValueError("at least one checkpoint shard is required")
    first, *others = checkpoints
    validate(first)
    for shard in others:
        validate(shard)
        if shard.stage != first.stage or not np.array_equal(shard.mass, first.mass) or not np.array_equal(
                shard.components, first.components):
            raise ValueError("checkpoint shards disagree on the replicated stage")
    merged = Checkpoint(first.stage, first.mass, first.components,
                        np.sort(np.concatenate([p.cells for p in checkpoints]), order="id"),
                        np.sort(np.concatenate([p.faces for p in checkpoints]), order="id"))
    validate(merged)
    if not np.array_equal(merged.cells["id"], cell_ids) or not np.array_equal(merged.faces["id"], face_ids):
        raise ValueError("checkpoint records do not match the required geometry IDs")
    return merged


def repartition(checkpoint, cell_ids, face_ids):
    """Select IDs for one new partition, including an empty dry rank."""
    selected = []
    for table, wanted in ((checkpoint.cells, cell_ids), (checkpoint.faces, face_ids)):
        wanted = np.asarray(wanted)
        if wanted.ndim != 1 or not np.issubdtype(wanted.dtype, np.integer):
            raise ValueError("expected a vector of integer checkpoint IDs")
        if np.any(wanted[1:] <= wanted[:-1]):
            raise ValueError("requested IDs must be strictly increasing")
        rows = np.searchsorted(table["id"], wanted)
        if np.any(rows >= len(table)) or not np.array_equal(table["id"][rows], wanted):
            raise ValueError("requested checkpoint IDs are absent")
        selected.append(table[rows])
    return Checkpoint(checkpoint.stage, checkpoint.mass, checkpoint.components, *selected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    inspect = sub.add_parser("inspect")
    inspect.add_argument("checkpoint", type=Path)
    merge = sub.add_parser("merge")
    merge.add_argument("--geometry", type=Path, required=True, help="version-1 irregular geometry .dat")
    merge.add_argument("--output", type=Path, required=True)
    merge.add_argument("shards", type=Path, nargs="+")
    args = parser.parse_args()
    try:
        if args.command == "inspect":
            chk = read_checkpoint(args.checkpoint)
            print(chk.stage)
            print(f"{len(chk.cells)} cell records, {len(chk.faces)} face records")
        else:
            if args.output.resolve() in {path.resolve() for path in (*args.shards, args.geometry)}:
                raise ValueError("merged output must be a separate file")
            geom, cells, faces = geometry_ids(args.geometry)
            chk = merge_checkpoints([read_checkpoint(path) for path in args.shards], cells, faces)
            if (chk.stage.geometry, chk.stage.ni, chk.stage.nj, chk.stage.islands) != (
                    geom.digest(), geom.grid.ni, geom.grid.nj, len(geom.loops)-1):
                raise ValueError("checkpoint geometry identity mismatch")
            write_checkpoint(args.output, chk)
    except (OSError, ValueError, TypeError, KeyError, struct.error) as error:
        parser.exit(1, f"cmm_checkpoint: {error}\n")


if __name__ == "__main__":
    main()
