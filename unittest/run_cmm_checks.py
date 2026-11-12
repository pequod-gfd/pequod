"""Build and run the CMM development checks from a fresh source snapshot."""

import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile

import contourpy
import matplotlib
import numpy as np
import scipy

from check_cmm_oracle import ALGEBRAIC_FIXTURES, check_oracle, fields


ROOT = Path(__file__).resolve().parents[1]
FLAGS = ("-O0 -g -fcheck=all -fbacktrace -finit-integer=-666 -finit-real=nan "
         "-fdefault-real-8 -fdefault-double-8 -fallow-argument-mismatch "
         "-ffree-line-length-none -Wall -std=gnu -DDDEBUG")
ARG_FIXTURES = ("test_weak_diffusion", "test_cmm_profile", "test_cmm_dirichlet", "test_irregular_mass", "test_irregular_mpi_geometry",
                "test_irregular_straight_wall", "test_irregular_operators", "test_irregular_vertex_order",
                "test_irregular_smooth_corners", "test_irregular_island_constraints", "test_irregular_transport",
                "test_irregular_forcing", "test_cape_cabaret")
VERTEX_FIXTURES = {
    "rectangle": "True",
    "diamond": "abs(x-0.5)+abs(y-0.5) < 0.375",
    "triangle": "x+y < 1",
    "concave": "x < 0.5 OR y < 0.5",
}
SQUARE_ISLANDS = ("NOT ((x >= 0.1875 AND x <= 0.3125 AND y >= 0.4375 AND y <= 0.5625) "
                  "OR (x >= 0.6875 AND x <= 0.8125 AND y >= 0.4375 AND y <= 0.5625))")
GEOMETRIES = (
    ("rectangle", 33, "True"),
    ("dry", 33, "abs(x-0.1875)+abs(y-0.1875) < 0.125"),
    ("diamond", 33, "abs(x-0.5)+abs(y-0.5) < 0.375"),
    ("triangle", 17, "x+y < 1"),
    ("triangle", 65, "x+y < 1"),
)


class CheckFailed(RuntimeError):
    pass


class Runner:
    def __init__(self, work, args):
        self.work = work
        self.args = args
        self.env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", MPLCONFIGDIR=str(work / "mpl"),
                        OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
        for key in ("MAKEFLAGS", "MFLAGS", "PYTHONOPTIMIZE"):
            self.env.pop(key, None)
        self.env.setdefault("OMPI_MCA_btl", "self,sm")
        self.env.setdefault("PRTE_MCA_rmaps_default_mapping_policy", ":oversubscribe")
        self.events = []
        self.metadata = {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
                         "matplotlib": matplotlib.__version__, "contourpy": contourpy.__version__,
                         "flags": FLAGS, "timeout_seconds": args.timeout,
                         "local_mpi_environment": {k: self.env[k] for k in
                            ("OMPI_MCA_btl", "PRTE_MCA_rmaps_default_mapping_policy")}}
        (work / "logs").mkdir()
        self.source = work / "source"
        self.source.mkdir()
        inputs = [*ROOT.glob("*.F90"), *ROOT.glob("*.h"), *ROOT.glob("makefile*"),
                  *(ROOT / "unittest").glob("*.F90"), *(ROOT / "unittest").glob("*.py"),
                  ROOT / "unittest" / "makefile", ROOT / "scripts" / "qg_irregular",
                  ROOT / "scripts" / "cmm_checkpoint.py"]
        self.metadata["source_sha256"] = {}
        for source in sorted(inputs):
            relative = source.relative_to(ROOT)
            target = self.source / relative
            target.parent.mkdir(exist_ok=True, parents=True)
            data = source.read_bytes()
            target.write_bytes(data)
            self.metadata["source_sha256"][str(relative)] = hashlib.sha256(data).hexdigest()

    def event(self, name, status, detail=""):
        self.events.append(dict(case=name, status=status, detail=detail))
        print(f"{status:6} {name}" + (f" — {detail}" if detail else ""), flush=True)

    def run(self, name, command, *, expected=0, contains=None, timeout=None, cwd=None, report_skip=False):
        command = list(map(str, command))
        log = self.work / "logs" / f"{name}.log"
        with log.open("w") as stream:
            stream.write("$ " + shlex.join(command) + "\n")
            stream.flush()
            with subprocess.Popen(command, cwd=cwd or self.source, env=self.env, stdout=stream,
                                  stderr=subprocess.STDOUT, start_new_session=True) as child:
                try:
                    status = child.wait(timeout=timeout or self.args.timeout)
                except subprocess.TimeoutExpired as error:
                    os.killpg(child.pid, signal.SIGTERM)
                    try:
                        child.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        os.killpg(child.pid, signal.SIGKILL)
                        child.wait()
                    self.event(name, "FAIL", f"timeout ({log})")
                    raise CheckFailed(f"{name} timed out") from error
        output = log.read_text(errors="replace").split("\n", 1)[1]
        if report_skip and status == 0 and output.startswith("SKIP "):
            self.event(name, "SKIP", output.splitlines()[0].removeprefix("SKIP "))
            return None
        if (status != expected if expected is not None else status == 0) or (contains and contains not in output):
            self.event(name, "FAIL", f"exit {status} ({log})")
            raise CheckFailed(f"{name}: {output[-3000:]}")
        self.event(name, "PASS", str(log))
        return output.strip()

    def build(self, backend):
        directory = self.work / backend
        shutil.copytree(self.source, directory)
        config = directory / "config.h"
        text = config.read_text()
        if text.count("#define PETSC") != 1:
            raise CheckFailed("expected one PETSC switch in config.h")
        cflags = libs = ""
        if backend == "native":
            config.write_text(text.replace("#define PETSC", "#undef PETSC"))
        else:
            self.metadata["PETSc"] = self.run("petsc-version", ["pkg-config", "--modversion", "PETSc"])
            cflags = self.run("petsc-cflags", ["pkg-config", "--cflags", "PETSc"])
            libs = self.run("petsc-libs", ["pkg-config", "--libs", "PETSc"])
        self.metadata.setdefault("effective_config", {})[backend] = config.read_text()
        self.run(f"build-{backend}", ["make", "-C", directory, f"-j{self.args.jobs}", "qg3", "unittest-bin",
                 f"FC={self.args.mpif90}", f"LD={self.args.mpif90}", f"FCFLAGS={FLAGS} {cflags}",
                 "LDFLAGS=-O0 -g", f"LIBS={libs} -llapack -lblas"], timeout=300)
        return directory / "unittest"

    def mpi(self, name, ranks, executable, *args, **kwargs):
        return self.run(name, [self.args.mpiexec, "-np", ranks, executable, *args], **kwargs)

    def generate(self, name, n, predicate, *, depth=4):
        folder = self.work / "geometry"
        folder.mkdir(exist_ok=True)
        dat = folder / f"{name}-{n}.dat"
        self.run(f"geometry-{name}-{n}", [sys.executable, self.source / "scripts/qg_irregular", n, n, "--bounds", 0, 1, 0, 1,
                 "--coordinate-units", "normalized", "--fluid", predicate, "--max-depth", depth,
                 "--output", dat, "--plot", dat.with_suffix(".png")])
        return dat

    def report(self):
        (self.work / "report.json").write_text(json.dumps(
            {"metadata": self.metadata, "checks": self.events}, indent=2) + "\n")
        counts = {status: sum(row["status"] == status for row in self.events)
                  for status in ("PASS", "FAIL", "XFAIL", "SKIP")}
        print(f"CMM development checks {counts}. Report {self.work / 'report.json'}", flush=True)


def check(runner):
    work, args = runner.work, runner.args
    runner.metadata["compiler"] = runner.run("compiler-version", [args.mpif90, "--version"]).splitlines()[0]
    runner.metadata["mpi"] = runner.run("mpi-version", [args.mpiexec, "--version"]).splitlines()[0]
    backends = ("native", "petsc") if args.backend == "both" else (args.backend,)
    if args.backend != "both":
        runner.event("other backend", "SKIP", f"explicit --backend {args.backend}")

    binaries = {backend: runner.build(backend) for backend in backends}
    specs = {(name, n): predicate for name, n, predicate in
             (*[entry[:3] for entry in ALGEBRAIC_FIXTURES], *GEOMETRIES)}
    specs.update({(name, n): predicate for name, predicate in VERTEX_FIXTURES.items() for n in (17, 33, 65)})
    specs.update({("square-islands", n): SQUARE_ISLANDS for n in (33, 65, 129, 257)})
    paths = {(name, n): runner.generate(name, n, predicate, depth=2 if n == 257 else 4)
             for (name, n), predicate in specs.items()}
    runner.metadata["geometry_sha256"] = {
        f"{name}-{n}": path.read_text().splitlines()[6].split()[-1] for (name, n), path in paths.items()}
    triangle = paths["triangle", 33]
    corrupt = work / "corrupt.dat"
    correct = triangle.read_text()
    corrupt.write_text(correct.replace(runner.metadata["geometry_sha256"]["triangle-33"], "0" * 64, 1))
    if args.inject_bad_fixture:
        triangle.write_bytes(corrupt.read_bytes())

    runner.run("python-geometry", [sys.executable, runner.source / "unittest/test_irregular_geometry.py", "-v"],
               contains="Ran 7 tests")
    runner.run("checkpoint-shards", [sys.executable, runner.source / "unittest/test_cmm_checkpoint.py",
                                    paths["square-islands", 33], work], contains="repartition PASS")
    for backend, binary in binaries.items():
        for source in sorted((runner.source / "unittest").glob("test_*.F90")):
            if source.stem in ARG_FIXTURES or source.stem.startswith(("test_cmm_", "test_irregular_")):
                continue
            if "helmholtz_petsc" in source.stem and backend == "native":
                runner.event(f"{backend}-legacy-{source.stem}", "SKIP", "PETSc fixture disabled in this build")
                continue
            match = re.search(r"^!!NP (\d+)$", source.read_text(), re.M)
            ranks = int(match[1]) if match else 1
            result = runner.mpi(f"{backend}-legacy-{source.stem}", ranks, binary / source.stem)
            lines = [line.split()[-1].upper() for line in result.splitlines() if line.split()]
            if "FAIL" in lines or "PASS" not in lines:
                runner.event(f"{backend}-legacy-{source.stem}-assertions", "FAIL",
                             "explicit failed assertion or missing PASS output")
        if backend == "petsc":
            for ranks in (1, 4):
                name = f"petsc-scalapack-{ranks}"
                runner.mpi(name, ranks, binary / "test_cmm_scalapack", report_skip=True,
                           contains="ScaLAPACK repeated Cholesky solve error:")
        reader = binary / "test_irregular_geometry_reader"
        runner.run(f"{backend}-sha-vectors", [reader], contains=": PASS")
        runner.run(f"{backend}-degree-four-moments", [binary / "test_irregular_moments"], contains="PASS")
        runner.run(f"{backend}-smooth-corner-scalar", [binary / "test_irregular_smooth_corners", "--scalar"],
                   contains="168 oriented")
        for (name, n), path in paths.items():
            runner.run(f"{backend}-reader-{name}-{n}", [reader, path, n, n, 1],
                       contains=runner.metadata["geometry_sha256"][f"{name}-{n}"])
        runner.run(f"{backend}-reader-corrupt", [reader, corrupt, 33, 33, 1],
                   expected=1, contains="SHA256 mismatch")
        for program in ARG_FIXTURES:
            runner.run(f"{backend}-{program}-missing-arguments", [binary / program],
                       expected=2, contains="requires")
        runner.run(f"{backend}-straight-scalar", [binary / "test_irregular_straight_wall", "--scalar"],
                   contains=": PASS")
        for fixture in ("triangle", "dry"):
            for ranks in (1, 2, 4, 12):
                runner.mpi(f"{backend}-ownership-{fixture}-{ranks}", ranks,
                           binary / "test_irregular_mpi_geometry", paths[fixture, 33], fixture,
                           contains="PASS")
        for fixture in ("triangle", "islands", "dry"):
            for ranks in (1, 4, 12):
                runner.mpi(f"{backend}-flux-balance-{fixture}-{ranks}", ranks, binary / "test_irregular_operators",
                           paths[fixture, 33], contains="max cubic face defect:")
        for fixture in ("triangle", "concave", "islands", "dry"):
            for ranks in (1, 4, 12):
                runner.mpi(f"{backend}-cut-transport-{fixture}-{ranks}", ranks, binary / "test_irregular_transport",
                           paths[fixture, 33], contains="cut-transport errors, CFL, area, time order:")
        cape_errors = []
        for n, variant, ranks in ((33, "affine", 1), (33, "affine", 4),
                                  (33, "smooth", 4), (65, "smooth", 4), (129, "smooth", 4)):
            case = f"{backend}-cape-cabaret-{n}-{variant}-{ranks}"
            folder = work / case
            folder.mkdir()
            output = runner.mpi(case, ranks, binary / "test_cape_cabaret", n, variant, cwd=folder,
                                contains="cape CABARET errors:")
            if variant == "smooth":
                fields = output.split("cape CABARET errors:", 1)[1].split()
                cape_errors.append([float(value) for value in fields[:2]])
        orders = np.log2(np.asarray(cape_errors[:-1]) / np.asarray(cape_errors[1:]))
        if not np.all(orders > 1.7):
            raise CheckFailed(f"{backend} cape CABARET smooth convergence regressed: {orders}")
        runner.event(f"{backend}-cape-cabaret-order", "PASS", str(orders.tolist()))
        for fixture in ("triangle", "islands", "dry"):
            for ranks in (1, 4, 12):
                runner.mpi(f"{backend}-forcing-{fixture}-{ranks}", ranks, binary / "test_irregular_forcing",
                           paths[fixture, 33], contains="irregular wind and drag errors:")
        for ranks in (1, 4):
            case = f"{backend}-forcing-cape-{ranks}"
            folder = work / case
            folder.mkdir()
            runner.mpi(case, ranks, binary / "test_irregular_forcing", "--cape", cwd=folder,
                       contains="cape wind errors:")
        for ranks in (1, 4, 12):
            runner.mpi(f"{backend}-smooth-corner-{ranks}", ranks, binary / "test_irregular_smooth_corners",
                       paths["concave", 33], contains="MPI smooth corner max error:")
        for fixture in ("triangle", "rectangle", "islands"):
            ranks_to_test = (1, 2, 4) if fixture == "triangle" else (1, 4)
            for ranks in ranks_to_test:
                output = runner.mpi(f"{backend}-mass-{fixture}-{ranks}", ranks, binary / "test_irregular_mass",
                                    paths[fixture, 33], fixture, contains="potential total:")
                if fixture == "triangle" and ranks > 1:
                    if "Collective invalid inputs (inactive cache, shape, NaN, disagreement) PASS" not in output:
                        raise CheckFailed("mass fixture did not exercise collective invalid inputs")
        for fixture in ("rectangle", "diamond", "islands", "triangle"):
            for ranks in (1, 4, 12):
                runner.mpi(f"{backend}-wall-{fixture}-{ranks}", ranks, binary / "test_irregular_straight_wall",
                           paths[fixture, 33], fixture, contains="PASS")
        vertex = binary / "test_irregular_vertex_order"
        for fixture in VERTEX_FIXTURES:
            errors = []
            for n in (17, 33, 65):
                out = runner.mpi(f"{backend}-vertex-{fixture}-{n}", 4, vertex,
                                 paths[fixture, n], n, fixture, "direct", contains="vertex-errors")
                errors.append(read_vertex_errors(out)[:2])
            orders = np.log2(np.asarray(errors[:-1])/errors[1:])
            assert np.all(orders >= 1.8), (fixture, orders)
            runner.event(f"{backend}-vertex-{fixture}-order", "PASS", str(orders.tolist()))
        for fixture in ("cubic", "concave", "islands", "near-islands"):
            reference = None
            for ranks in (1, 4, 12):
                mode = "inverted" if fixture in ("cubic", "concave", "islands") else "direct"
                out = runner.mpi(f"{backend}-vertex-{fixture}-{mode}-{ranks}", ranks, vertex,
                                 paths["triangle" if fixture == "cubic" else fixture, 33], 33, fixture, mode,
                                 contains="vertex-errors")
                errors = read_vertex_errors(out)
                if reference is None:
                    reference = errors
                assert np.max(abs(reference-errors)) < 1e-9, (fixture, ranks, reference, errors)
        island_errors = []
        for n in (33, 65, 129, 257):
            output = runner.mpi(f"{backend}-weak-islands-{n}", 4, binary / "test_irregular_island_constraints",
                                paths["square-islands", n], n, contains="reused influence factor PASS", timeout=max(args.timeout, 450))
            line, = [line for line in output.splitlines() if line.startswith("island-errors ")]
            island_errors.append(list(map(float, line.split()[1:])))
        orders = np.log2(np.array(island_errors[:-1])[:, :2]/np.array(island_errors[1:])[:, :2])
        for coarse, fine, order in zip((33, 65, 129), (65, 129, 257), orders):
            status = "PASS" if np.all(order >= 1.8) else "FAIL"
            runner.event(f"{backend}-weak-island-convergence-{coarse}-{fine}", status,
                         f"orders {order.tolist()}; all grid errors {island_errors}")

        physical_errors = []
        slip_errors = []
        for n in (33, 65, 129):
            result = runner.mpi(f"{backend}-physical-island-means-{n}", 4, binary / "test_irregular_island_constraints",
                                paths["square-islands", n], n, "averages", contains="independent-island-slip")
            line, = [line for line in result.splitlines() if line.startswith("island-errors ")]
            physical_errors.append(list(map(float, line.split()[1:3])))
            line, = [line for line in result.splitlines() if line.startswith("independent-island-slip ")]
            slip_errors.append(float(line.split()[1]))
        errors = np.asarray(physical_errors)
        orders = np.log2(errors[:-1]/errors[1:])
        if np.all(orders >= 1.8) and all(a > b for a,b in zip(slip_errors[:-1],slip_errors[1:])):
            status = "PASS"
        else:
            status = "FAIL"
        runner.event(f"{backend}-physical-island-refinement", status,
                     f"psi and constant orders {orders.tolist()}, independent slip residuals {slip_errors}")

    for backend, binary in binaries.items():
        output = work / f"{backend}-reordered.bin"
        runner.mpi(f"{backend}-reordered-inversion", 4, binary / "test_cmm_dirichlet",
                   triangle, output, 33, "triangle", "reordered",
                   contains="aggregate explicit two-solve, Green, dense setup bytes:")
        check_oracle(triangle, output, 33)
    for backend, binary in binaries.items():
        rank_counts = (1, 4, 12) if backend == "petsc" else (4,)
        for ranks in rank_counts:
            output = work / f"{backend}-distributed-{ranks}.bin"
            result = runner.mpi(f"{backend}-distributed-CMM-{ranks}", ranks, binary / "test_cmm_dirichlet",
                                paths["triangle", 65], output, 65, "algebraic", "distributed-reordered",
                                report_skip=True, contains="aggregate explicit two-solve, Green, dense setup bytes:")
            if result is not None:
                check_oracle(paths["triangle", 65], output, 65)
    for name, n, _, rank_counts in ALGEBRAIC_FIXTURES:
        reference = None
        for backend, binary in binaries.items():
            for ranks in rank_counts:
                label = f"{backend}-CMM-{name}-{ranks}"
                output = work / f"{label}.bin"
                runner.mpi(label, ranks, binary / "test_cmm_dirichlet",
                           paths[name, n], output, n, "algebraic", contains="rcond:")
                values = check_oracle(paths[name, n], output, n)
                if reference is None:
                    reference = values
                discrepancy = np.max(abs(reference-values))
                assert discrepancy < 1e-10
                runner.event(f"{label}-oracle", "PASS", f"rank/backend difference {discrepancy:.3e}")

    errors = {}
    for backend, binary in binaries.items():
        maxima = []
        for n in (17, 33, 65):
            label = f"{backend}-smooth-triangle-{n}"
            output = work / f"{label}.bin"
            runner.mpi(label, 4, binary / "test_cmm_dirichlet",
                       paths["triangle", n], output, n, "smooth", contains="rcond:")
            kinds, _, _, values, _ = fields(paths["triangle", n], output, n)
            y, x = np.indices(kinds.shape)/(n-1)
            exact = (x*(1-y)*(y-x))**2
            fluid = kinds == 1
            maxima.append(float(np.max(abs(values[:, :, :, 0][fluid]-exact[fluid, None]))))
        orders = np.log2(np.array(maxima[:-1])/maxima[1:])
        assert np.all(orders >= 1.8)
        errors[backend] = dict(errors=maxima, orders=orders.tolist())
        runner.event(f"{backend}-smooth-convergence", "PASS", str(errors[backend]))
    runner.metadata["elliptic_convergence"] = errors

    if "petsc" in binaries:
        runner.mpi("petsc-reject-inaccurate-background", 4, binaries["petsc"] / "test_cmm_dirichlet",
                   triangle, work / "inaccurate.bin", 33, "algebraic", "-ksp_type", "preonly",
                   "-pc_type", "jacobi", expected=None, contains="rectangle solve backward error")
    else:
        runner.event("petsc-reject-inaccurate-background", "SKIP", "PETSc backend was excluded")

    report = work / "corner-experiments.json"
    runner.run("corner-research-fixtures", [sys.executable, runner.source / "unittest/check_corner_reconstruction.py",
               "--discrete", "--report", report], contains="No-slip singular patches")
    for vertex in json.loads(report.read_text())["vertex_order"]:
        orientation = "diagonal" if vertex["diagonal"] else "axial"
        orders = np.asarray(vertex["orders"])[:, 0]
        if np.all(orders >= 1.8):
            raise CheckFailed(f"R2 {orientation}: unexpected pass. Replace the known-failure designation with runtime acceptance.")
        if not np.all((orders > 0.9) & (orders < 1.1)):
            raise CheckFailed(f"R2 {orientation}: changed failure signature {orders}")
        runner.event(f"R2-four-centre-comparison-{orientation}", "PASS",
                     f"reproduced the retired comparator's order {orders}; cubic path tested above")


def read_vertex_errors(output):
    line, = [line for line in output.splitlines() if line.startswith("vertex-errors ")]
    return np.array(list(map(float, line.split()[1:-1])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("native", "petsc", "both"), default="both")
    parser.add_argument("--work-parent", type=Path, help="directory in which to retain a fresh run")
    parser.add_argument("--mpif90", default=os.environ.get("MPIF90", "mpif90"))
    parser.add_argument("--mpiexec", default=os.environ.get("MPIEXEC", "mpirun"))
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=60, help="per-case deadline in seconds")
    parser.add_argument("--inject-bad-fixture", action="store_true", help="exercise failure propagation")
    args = parser.parse_args()
    if not __debug__ or args.timeout <= 0 or args.jobs < 1:
        parser.error("assertions, positive timeout and at least one build job are required")
    work = Path(tempfile.mkdtemp(prefix="pequod-cmm-check.", dir=args.work_parent)).resolve()
    runner = Runner(work, args)
    try:
        check(runner)
    except (AssertionError, OSError, ValueError, CheckFailed) as error:
        if not runner.events or runner.events[-1]["status"] != "FAIL":
            runner.event("orchestrator", "FAIL", repr(error))
        print(error, file=sys.stderr)
        return 1
    finally:
        runner.report()
    return int(any(row["status"] == "FAIL" for row in runner.events))


if __name__ == "__main__":
    sys.exit(main())
