"""Geometry oracles for qg_irregular. Run with python3 unittest/test_irregular_geometry.py."""
import hashlib
from fractions import Fraction
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "qg_irregular"
API = runpy.run_path(str(SCRIPT))
Grid, Geometry, GeometryError = (API[name] for name in ("Grid", "Geometry", "GeometryError"))
RAYS = ((1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1))


def expand(vertices):
    """Fixture helper only. Subdivide straight edges of an explicit polygon."""
    loop = []
    for (x, y), (u, v) in zip(vertices, vertices[1:] + vertices[:1]):
        steps = max(abs(u - x), abs(v - y))
        assert steps and (u - x) % steps == (v - y) % steps == 0
        loop.extend((x + k * (u - x) // steps, y + k * (v - y) // steps) for k in range(steps))
    return loop


def to_indices(loop, n):
    return tuple((n - y, x + 1) for x, y in loop)


def inside(point, loop):
    """Independent winding-number oracle, in x east and y north."""
    x, y = point
    winding = 0
    for (a, b), (c, d) in zip(loop, loop[1:] + loop[:1]):
        cross = (c - a) * (y - b) - (d - b) * (x - a)
        if b <= y < d and cross > 0:
            winding += 1
        if d <= y < b and cross < 0:
            winding -= 1
    return winding != 0


def raster_oracle(loops, n):
    """Test at rational sector centroids. No polygon can cross a sector interior."""
    result = np.zeros((n, n), dtype=int)
    polygons = [list((j - 1, n - i) for i, j in loop) for loop in loops]
    offsets = [(Fraction(x + u, 6), Fraction(y + v, 6))
               for (x, y), (u, v) in zip(RAYS, RAYS[1:] + RAYS[:1])]
    for i in range(n):
        for j in range(n):
            for bit, (dx, dy) in enumerate(offsets):
                point = j + dx, n - 1 - i + dy
                if inside(point, polygons[0]) and not any(inside(point, loop) for loop in polygons[1:]):
                    result[i, j] |= 1 << bit
    return result


def transforms(loop, n):
    for mirrored in (False, True):
        points = [(n - 1 - x if mirrored else x, y) for x, y in loop]
        for _ in range(4):
            yield to_indices(points, n)
            points = [(y, n - 1 - x) for x, y in points]


class GeometryTests(unittest.TestCase):
    def test_parser(self):
        predicate = API["parse_predicate"]("0 <= x <= 1 AND NOT (x^2 + y**2 > 1) OR abs(y) == 2")
        np.testing.assert_array_equal(predicate(np.array([0.0, 1.0, 0.0, 3.0, -1.0]),
                                                np.array([0.0, 0.0, 1.1, -2.0, 0.0])),
                                      [True, True, False, True, False])
        for source in ("open('file')", "x.real > 0", "x[0] > 0", "x + 1", "sqrt(-1) < x", "1e999 < x",
                       "1 / 0 > x", "sin(x, y) == 0", "x if y else 1"):
            with self.subTest(source=source), self.assertRaises((ValueError, ArithmeticError)):
                API["parse_predicate"](source)(0.0, 0.0)
        with self.assertRaisesRegex(GeometryError, "equal coordinate spacing"):
            Grid(5, 6, (0., 1., 0., 1.), "normalized")
        translated = Grid(5, 7, (12., 24., -7., 1.), "cm")
        np.testing.assert_array_equal(translated.coordinates(np.array([0, 6]), np.array([0, 4])),
                                      [[12, 24], [1, -7]])

    def test_fig8_and_exact_sectors(self):
        n = 9
        grid = Grid(n, n, (0., 8., 0., 8.), "normalized")
        paths = [[(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (3, 2)],
                 [(0, 0), (1, 1), (2, 2), (3, 2)],
                 [(0, 0), (1, 0), (2, 1), (3, 2)]]
        for path in paths:
            loop = expand(path + [(7, 2), (7, 7), (0, 7)])
            for transformed in transforms(loop, n):
                geometry = Geometry.from_loops(grid, [transformed])
                self.assertEqual(len(geometry.loops[0]), len(loop))
                np.testing.assert_array_equal(geometry.nodes[:, :, 2], raster_oracle(geometry.loops, n))
        invalid = to_indices([(0, 0), (3, 2), (4, 2), (4, 3)], n)
        with self.assertRaisesRegex(GeometryError, "nonadjacent"):
            Geometry.from_loops(grid, [invalid])

    def test_all_sectors_and_orientations(self):
        n = 11
        grid = Grid(n, n, (0., 10., 0., 10.), "normalized")
        for count in range(1, 8):
            loop = expand([(5, 5)] + [(5 + 4*x, 5 + 4*y) for x, y in RAYS[:count + 1]])
            for transformed in transforms(loop, n):
                geometry = Geometry.from_loops(grid, [transformed])
                self.assertEqual(int(geometry.nodes[5, 5, 2]).bit_count(), count)
                np.testing.assert_array_equal(geometry.nodes[:, :, 2], raster_oracle(geometry.loops, n))

    def test_invalid_topology(self):
        n = 17
        grid = Grid(n, n, (0., 16., 0., 16.), "normalized")
        box = to_indices(expand([(0, 0), (16, 0), (16, 16), (0, 16)]), n)
        crossing = to_indices([(2, 2), (3, 3), (2, 3), (3, 2)], n)
        with self.assertRaisesRegex(GeometryError, "crossing diagonals"):
            Geometry.from_loops(grid, [crossing])
        touching = to_indices(expand([(0, 8), (2, 6), (4, 8), (2, 10)]), n)
        with self.assertRaisesRegex(GeometryError, "shared boundary node"):
            Geometry.from_loops(grid, [box, touching])
        gap = to_indices(expand([(1, 8), (3, 6), (5, 8), (3, 10)]), n)
        with self.assertRaisesRegex(GeometryError, "one-spacing gap"):
            Geometry.from_loops(grid, [box, gap])
        dumbbell = to_indices(expand([(1, 1), (5, 1), (5, 3), (8, 3), (8, 1), (13, 1),
                                      (13, 6), (8, 6), (8, 4), (5, 4), (5, 6), (1, 6)]), n)
        with self.assertRaisesRegex(GeometryError, "disconnected strict-fluid graph"):
            Geometry.from_loops(grid, [dumbbell])
        predicate = API["parse_predicate"]("NOT ((x - 8.5)^2 + (y - 8.5)^2 < 0.15^2)")
        with self.assertRaisesRegex(GeometryError, "subgrid closed contour"):
            API["sampled_loops"](grid, predicate, 1e-8, 4)
        oscillatory = API["parse_predicate"]("y > 8.25 + 0.7*cos(9*x)")
        with self.assertRaisesRegex(GeometryError, "multiple sampled crossings"):
            API["sampled_loops"](grid, oscillatory, 1e-8, 4)

    def test_roundtrip_corruption_and_plot(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            out = directory / "two-islands.dat"
            plot = directory / "two-islands.png"
            predicate = ("(x > 0.093 AND y > 0.093) AND NOT (abs(x-0.34375)+abs(y-0.4375) <= 0.125 "
                         "OR abs(x-0.75)+abs(y-0.65625) <= 0.125)")
            command = [sys.executable, str(SCRIPT), "33", "33", "--bounds", "0", "1", "0", "1",
                       "--coordinate-units", "normalized", "--fluid", predicate,
                       "--output", str(out), "--plot", str(plot)]
            subprocess.run(command, check=True, capture_output=True, text=True)
            snapshot = out.read_bytes()
            self.assertEqual(len(Geometry.read(out).loops), 3)
            image = plot.read_bytes()
            self.assertTrue(image.startswith(b"\x89PNG"))
            subprocess.run(command, check=True, capture_output=True, text=True)
            self.assertEqual(hashlib.sha256(snapshot).digest(), hashlib.sha256(out.read_bytes()).digest())
            failed = subprocess.run(command + ["--max-snap-distance", "1e-5"], capture_output=True, text=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("exceeds limit", failed.stderr)
            self.assertEqual(snapshot, out.read_bytes())
            lines = snapshot.decode("ascii").splitlines()
            for before, after, expected in (("PEQUOD_IRREGULAR 1", "PEQUOD_IRREGULAR 2", "invalid record"),
                                            ("1 1 0 -1 0", "1 1 1 -1 255", "flags disagree"),
                                            ("END", "END\nstray", "unexpected data")):
                with self.subTest(before=before):
                    out.write_text("\n".join(after if line == before else line for line in lines) + "\n")
                    with self.assertRaisesRegex(GeometryError, expected):
                        Geometry.read(out)
            digest = next(line for line in lines if line.startswith("HASH SHA256 "))
            out.write_bytes(snapshot.replace(digest.encode(), b"HASH SHA256 " + b"a"*64))
            with self.assertRaisesRegex(GeometryError, "SHA256 mismatch"):
                Geometry.read(out)


if __name__ == "__main__":
    unittest.main()
