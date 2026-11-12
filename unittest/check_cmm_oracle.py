"""Run prescribed-Dirichlet Fortran fixtures against independent SciPy sparse LU."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


GENERATOR = Path(__file__).resolve().parents[1] / "scripts" / "qg_irregular"
CARDINAL = [(-1, 0), (0, -1), (1, 0), (0, 1)]
SHIFTS = [0, 0.3, 0.3]
ALGEBRAIC_FIXTURES = (
    ("rectangle", 17, "True", (1, 4)),
    ("triangle", 33, "x+y < 1", (1, 2, 4, 12)),
    ("bay", 33, "(x <= 0.5 AND y >= 0.3) OR (x > 0.5 AND y > 0.1)", (1, 4)),
    ("islands", 33, "(x > 0.093 AND y > 0.093) AND NOT (abs(x-0.34375)+abs(y-0.4375) <= 0.125 "
                    "OR abs(x-0.75)+abs(y-0.65625) <= 0.125)", (1, 2, 4, 12)),
    ("near-islands", 33, "NOT (abs(x-0.34375)+abs(y-0.5) <= 0.125 "
                         "OR abs(x-0.65625)+abs(y-0.5) <= 0.125)", (1, 4)),
)


def generate(folder, n, predicate, name):
    dat = folder / f"{name}-{n}.dat"
    subprocess.run([sys.executable, str(GENERATOR), str(n), str(n),
                    "--bounds", "0", "1", "0", "1", "--coordinate-units", "normalized",
                    "--fluid", predicate, "--output", str(dat), "--plot", str(dat.with_suffix(".png"))],
                   check=True, stdout=subprocess.DEVNULL)
    return dat


def fields(dat, path, n):
    # Only the flags and ordered indices are file inputs to the independent
    # solve. Assemble row numbers and all couplings here.
    records = np.loadtxt(dat, skiprows=8, max_rows=n*n, dtype=int)
    kinds = records[:, 2].reshape((n, n))
    components = records[:, 3].reshape((n, n))
    interior_wall = np.argwhere(kinds[1:-1, 1:-1] == 2) + 1
    m = len(interior_wall)
    data = np.fromfile(path, dtype=np.float64)
    bulk = 9*n*n
    return (kinds, components, interior_wall, data[:bulk].reshape((n, n, 3, 3), order="F"),
            data[bulk:].reshape((m, m, 3), order="F"))


def operator(nodes, shift, rhs=None, components=None, g=None):
    row_of = {tuple(point): row for row, point in enumerate(nodes)}
    rows, columns, values = [], [], []
    b = np.empty(len(nodes)) if rhs is not None else None
    for row, (i, j) in enumerate(nodes):
        rows.append(row)
        columns.append(row)
        values.append(-4 - shift)
        if b is not None:
            b[row] = rhs[i, j]
        for di, dj in CARDINAL:
            neighbor = i + di, j + dj
            if neighbor in row_of:
                rows.append(row)
                columns.append(row_of[neighbor])
                values.append(1)
            elif b is not None:
                component = components[neighbor]
                assert component >= 0, (i, j, neighbor)
                b[row] -= g[component]
    return coo_matrix((values, (rows, columns)), shape=(len(nodes), len(nodes))).tocsc(), b


def check_oracle(dat, output, n):
    kinds, components, walls, values, capacitance = fields(dat, output, n)
    fluid = np.argwhere(kinds == 1)
    fy, fx = fluid.T
    rows, cols = np.indices(kinds.shape) + 1
    ncomponents = np.max(components) + 1
    reference_error = residual = wall_error = matrix_error = 0.0
    for mode, shift in enumerate(SHIFTS, 1):
        rhs = ((23*rows+17*cols+13*rows*cols+5*mode) % 43-21)/31
        g = 0.25*(1 + 2*np.arange(ncomponents) - mode)
        matrix, b = operator(fluid, shift, rhs, components, g)
        expected = splu(matrix).solve(b)
        for trial in (0, 1):
            phi = values[:, :, mode-1, trial]
            reference_error = max(reference_error, np.max(abs(phi[fy, fx] - expected)))
            residual = max(residual, np.max(abs(matrix @ phi[fy, fx] - b)))
            wall_error = max(wall_error, np.max(abs(phi[kinds == 2] - g[components[kinds == 2]])))

        if len(walls):
            # Full homogeneous background Green block by independent sparse LU.
            mesh = [(i, j) for i in range(1, n-1) for j in range(1, n-1)]
            full_matrix, _ = operator(mesh, shift)
            indices = (walls[:, 0]-1)*(n-2) + walls[:, 1]-1
            impulses = np.zeros((len(mesh), len(walls)))
            impulses[indices, np.arange(len(walls))] = 1
            c = -splu(full_matrix).solve(impulses)[indices]
            matrix_error = max(matrix_error, np.max(abs(c - capacitance[:, :, mode-1])))
            assert np.min(np.linalg.eigvalsh(capacitance[:, :, mode-1])) > 0
    assert reference_error < 1e-10
    assert max(residual, wall_error, matrix_error) < 1e-11
    print(f"  oracle {reference_error:.3e}, PDE {residual:.3e}, wall {wall_error:.3e}, C {matrix_error:.3e}")
    return values[kinds != 0]


def run_fixture(executable, ranks, dat, output, n, case):
    mpirun = os.environ.get("MPIEXEC", "mpirun")
    result = subprocess.run([mpirun, "-np", str(ranks), str(executable), str(dat),
                             str(output), str(n), case], text=True, capture_output=True)
    if result.returncode:
        print(result.stdout, result.stderr, file=sys.stderr)
        result.check_returncode()
    print(executable.parent.parent.name, case, dat.stem, ranks, result.stdout.strip(), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--petsc", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    for name, n, predicate, rank_counts in ALGEBRAIC_FIXTURES:
        dat = generate(args.work, n, predicate, name)
        reference = None
        for executable in (args.native, args.petsc):
            for ranks in rank_counts:
                output = args.work / f"{executable.parent.parent.name}-{name}-{ranks}.bin"
                run_fixture(executable, ranks, dat, output, n, "algebraic")
                values = check_oracle(dat, output, n)
                if reference is None:
                    reference = values
                else:
                    error = np.max(abs(values-reference))
                    assert error < 1e-10
                    print(f"  agreement with native serial {error:.3e}")

    errors = []
    for n in (17, 33, 65):
        dat = generate(args.work, n, "x+y < 1", "triangle")
        output = args.work / f"smooth-{n}.bin"
        run_fixture(args.native, 4, dat, output, n, "smooth")
        kinds, _, _, values, _ = fields(dat, output, n)
        y, x = np.indices(kinds.shape) / (n-1)
        y = 1-y
        exact = (x*y*(1-x-y))**2
        errors.append(float(np.max(abs(values[:, :, :, 0][kinds == 1] - exact[kinds == 1, None]))))
    rates = np.log2(np.array(errors[:-1])/errors[1:])
    print("smooth streamfunction max errors:", errors, "orders:", rates)
    assert np.all(rates >= 1.8)


if __name__ == "__main__":
    main()
