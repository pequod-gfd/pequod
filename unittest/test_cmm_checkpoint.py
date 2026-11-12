"""Validate the experimental CMM stream, exact coverage, and repartition."""
from dataclasses import replace
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from cmm_checkpoint import (CELL_FIELDS, COMPONENT_FIELDS, FACE_FIELDS, Checkpoint, Stage, geometry_ids,
                            merge_checkpoints, read_checkpoint, records, repartition, write_checkpoint)


def table(ids, fields, layers):
    result = np.zeros(len(ids), records(fields, layers))
    result["id"] = ids
    result["values"] = ids[:, None, None]*0.125+np.arange(len(fields)*layers).reshape(-1, layers)*0.25
    return result


def reject(function, text):
    try:
        function()
    except ValueError as error:
        assert text in str(error), str(error)
    else:
        raise AssertionError(f"checkpoint case accepted without {text!r}")


def check(geometry, work):
    geom, cells, faces = geometry_ids(geometry)
    layers = 3
    islands = len(geom.loops)-1
    assert islands == 2
    stage = Stage(geom.digest(), "1"*64, "2"*64, geom.grid.ni, geom.grid.nj, layers,
                  islands, 9, 1.375, 0.375, 0.125)
    mass = np.array([0.25, -0.75])
    components = np.arange((islands+1)*len(COMPONENT_FIELDS)*layers).reshape(islands+1, -1, layers)/8
    original = Checkpoint(stage, mass, components, table(cells, CELL_FIELDS, layers),
                          table(faces, FACE_FIELDS, layers))
    files = [work / f"shard-{i}.bin" for i in range(5)]
    for rank, path in enumerate(files):
        write_checkpoint(path, repartition(original, cells[rank::5], faces[rank::5]))
    shards = [read_checkpoint(path, expected=stage) for path in files[::-1]]
    merged = merge_checkpoints(shards, cells, faces)
    for attribute in ("mass", "components", "cells", "faces"):
        np.testing.assert_array_equal(getattr(merged, attribute), getattr(original, attribute))
    dry = np.array([], dtype=np.int64)
    assert len(repartition(merged, dry, dry).cells) == 0
    for rank in range(3):
        part = repartition(merged, cells[rank::3], faces[rank::3])
        np.testing.assert_array_equal(part.faces, original.faces[rank::3])
        np.testing.assert_array_equal(part.cells, original.cells[rank::3])

    # Exercise the CLI against the geometry, then verify byte reproducibility.
    script = Path(__file__).resolve().parents[1] / "scripts/cmm_checkpoint.py"
    out = work / "merged.bin"
    subprocess.run([sys.executable, script, "merge", "--geometry", geometry, "--output", out, *files], check=True)
    direct = work / "direct.bin"
    write_checkpoint(direct, original)
    assert direct.read_bytes() == out.read_bytes()

    reject(lambda: read_checkpoint(out, expected=replace(stage, closure="3"*64)), "stage or identity")
    reject(lambda: merge_checkpoints(shards[:-1], cells, faces), "required geometry IDs")
    reject(lambda: merge_checkpoints(shards+[shards[0]], cells, faces), "strictly increasing")
    mixed = [replace(shards[0], stage=replace(stage, dt_current=0.5)), *shards[1:]]
    reject(lambda: merge_checkpoints(mixed, cells, faces), "replicated stage")
    mixed = [replace(shards[0], mass=mass*2), *shards[1:]]
    reject(lambda: merge_checkpoints(mixed, cells, faces), "replicated stage")
    data = bytearray(out.read_bytes())
    data[-40] ^= 1
    corrupt = work / "corrupt.bin"
    corrupt.write_bytes(data)
    reject(lambda: read_checkpoint(corrupt), "SHA256 mismatch")
    corrupt.write_bytes(data[:-1])
    reject(lambda: read_checkpoint(corrupt), "length mismatch")
    data[0] ^= 1
    corrupt.write_bytes(data)
    reject(lambda: read_checkpoint(corrupt), "unknown checkpoint format")
    print("CMM checkpoint shard, coverage, identity and repartition PASS")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("requires geometry file and scratch parent")
    with tempfile.TemporaryDirectory(prefix="cmm-checkpoint-", dir=sys.argv[2]) as folder:
        check(Path(sys.argv[1]), Path(folder))
