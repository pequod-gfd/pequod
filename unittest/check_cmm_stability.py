"""Independent energy and spectral checks for the coupled no-slip source."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, eigvals, solve, svdvals
from scipy.sparse import coo_matrix, diags

from check_cmm_weak_corner import Reference


def measure(n, case):
    slope = {"225": -1, "315": 1}.get(case, 0)
    ref = Reference(n, concave=case in ("225", "270", "315"), islands=case=="islands", notch_slope=slope)
    strict = np.flatnonzero(ref.strict)
    rows = strict.tolist()
    columns = list(range(len(strict)))
    for k,(lo,hi) in enumerate(ref.boxes):
        wall = np.flatnonzero(((ref.indices>=lo)&(ref.indices<=hi)).all(axis=1))
        rows.extend(wall)
        columns.extend([len(strict)+k]*len(wall))
    inject = coo_matrix((np.ones(len(rows)),(rows,columns)),shape=(len(ref.area),len(strict)+len(ref.boxes))).tocsr()
    stiffness = ref.stiffness@inject
    energy = (-inject.T@stiffness).toarray()
    dissipation = (inject.T@ref.diffusion@diags(1/(ref.h**2*ref.area))@stiffness).toarray()
    symmetric = (dissipation+dissipation.T)/2
    # Map the energy inner product into the Euclidean norm. The actual
    # linear generator is -E^-1 K. Its spectrum is evaluated independently
    # from the symmetric Rayleigh bounds.
    e,u = eigh(energy)
    assert np.min(e)>0, "the Dirichlet outer gauge must remove the constant nullspace"
    transform = (u/np.sqrt(e))@u.T
    scaled = transform@dissipation@transform
    eig = eigvals(-scaled)
    coercivity = eigh(symmetric,energy,eigvals_only=True,subset_by_index=[0,0])[0]
    norm = svdvals(scaled)[0]
    euler_limit = 2*coercivity/norm**2 if coercivity>0 else None

    # Compare the actual worst single-step amplification with the energy
    # proof on three dt values, including one outside an explicit CFL.
    identity = np.eye(len(e))
    gain = {}
    for dt in (0.25*ref.h**2,2*ref.h**2,0.1):
        cn = solve(identity+dt*scaled/2,identity-dt*scaled/2)
        gain[str(dt)] = float(svdvals(cn)[0])
    return dict(n=n,case=case,coercivity=float(coercivity),operator_norm=float(norm),
                rightmost_eigenvalue=float(np.max(eig.real)),
                sufficient_euler_dt=None if euler_limit is None else float(euler_limit),cn_gains=gain)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grids", type=int, nargs="+", default=(17,33))
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for n in args.grids:
        for case in ("rectangle","225","270","315","islands"):
            row = measure(n,case)
            rows.append(row)
            print(row,flush=True)
    args.report.write_text(json.dumps(rows,indent=2)+"\n")


if __name__ == "__main__":
    main()
