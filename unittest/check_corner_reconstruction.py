"""Corner coefficient and manufactured-wedge experiments for M6.

This fixture develops local operators on an ideal sector. It does not install
the experimental coefficients as a wall closure in the time-stepping code.
"""
import argparse
import json
from math import factorial
from pathlib import Path

import numpy as np
from numpy.polynomial.polynomial import polyder, polyval2d


RAYS = np.array([(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)])
DEGREE = 4
POWERS = np.array([(a, total-a) for total in range(1, DEGREE+1) for a in range(total+1)])
TOTAL = POWERS.sum(axis=1)
CONDITION_LIMIT = 1e6
GAUSS_X, GAUSS_W = np.polynomial.legendre.leggauss(3)


def geometry(start, sectors):
    tangents = RAYS[[start, (start+sectors) % 8]]
    normals = tangents[:, ::-1] * [[-1, 1], [1, -1]]
    normals = normals / np.linalg.norm(normals, axis=1)[:, None]
    return tangents, normals


def samples(start, sectors, radius):
    xy = np.array([(x, y) for x in range(-radius, radius+1) for y in range(-radius, radius+1)
                   if x or y])
    theta = (np.arctan2(xy[:, 1], xy[:, 0]) - start*np.pi/4) % (2*np.pi)
    return xy[(theta > 1e-12) & (theta < sectors*np.pi/4-1e-12)]


def derivative_row(xy, dx=0, dy=0):
    x, y = xy
    return np.array([x**(a-dx) * y**(b-dy) / (factorial(a-dx)*factorial(b-dy))
                     if a >= dx and b >= dy else 0 for a, b in POWERS], dtype=float)


def normal_row(xy, normal):
    nx, ny = normal
    return nx*derivative_row(xy, 1, 0) + ny*derivative_row(xy, 0, 1)


def normal2_row(xy, normal):
    nx, ny = normal
    return (nx*nx*derivative_row(xy, 2, 0) + 2*nx*ny*derivative_row(xy, 1, 1)
            + ny*ny*derivative_row(xy, 0, 2))


def torque_row(xy, normal):
    nx, ny = normal
    return nx*(derivative_row(xy, 3, 0)+derivative_row(xy, 1, 2)) + ny*(
        derivative_row(xy, 2, 1)+derivative_row(xy, 0, 3))


def constraints(tangents, normals, slip, alpha):
    rows = []
    for t, n in zip(tangents, normals):
        for degree in range(1, DEGREE+1):
            rows.append(np.where(TOTAL == degree, derivative_row(t), 0))
        for k in range(DEGREE if slip == "no-slip" else DEGREE-1):
            first = np.where(TOTAL == k+1, normal_row(t, n), 0)
            second = np.where(TOTAL == k+2, normal2_row(t, n), 0)
            rows.append(first if slip == "no-slip" else second-alpha*first)
    rows = np.array(rows)
    norms = np.linalg.norm(rows, axis=1)
    return rows[norms > 0] / norms[norms > 0, None]


def functionals(tangents, normals):
    output = [derivative_row((0, 0), 2, 0)+derivative_row((0, 0), 0, 2)]
    for derivative in (normal_row, torque_row):
        for t, n in zip(tangents, normals):
            # Three-point Gaussian quadrature integrates each derivative
            # polynomial exactly over the half-segment adjacent to the corner.
            length = np.linalg.norm(t)/2
            output.append(sum(weight*derivative((node+1)*t/4, -n)
                              for node, weight in zip(GAUSS_X, GAUSS_W))*length/2)
    return np.array(output)


def fit_corner(start, sectors, slip, alpha):
    tangents, normals = geometry(start, sectors)
    constraint = constraints(tangents, normals, slip, alpha)
    _, singular, vectors = np.linalg.svd(constraint, full_matrices=True)
    rank = np.count_nonzero(singular > 64*np.finfo(float).eps*max(constraint.shape)*singular[0])
    basis = vectors[rank:].T
    assert np.max(np.abs(constraint @ basis), initial=0) < 1e-12
    for radius in range(1, 5):
        xy = samples(start, sectors, radius)
        if not len(xy):
            continue
        design = np.stack(list(map(derivative_row, xy)))
        if not basis.shape[1]:
            return xy, radius, basis, np.zeros((len(POWERS), len(xy))), 1.0
        weight = 1 / (1 + np.sum(xy**2, axis=1))
        weighted = weight[:, None]*design
        u, s, vh = np.linalg.svd(weighted @ basis, full_matrices=False)
        rank_floor = 64*np.finfo(float).eps*max(weighted.shape)*np.linalg.norm(weighted)
        if len(s) == basis.shape[1] and s[-1] > max(rank_floor, s[0]/CONDITION_LIMIT):
            return xy, radius, basis, ((basis @ vh.T / s) @ u.T)*weight, s[0]/s[-1]
    raise AssertionError(f"insufficient admissible donor rank at sector {sectors}, ray {start}, {slip}")


def polynomial(coefficients):
    table = np.zeros((DEGREE+1, DEGREE+1))
    for (a, b), coefficient in zip(POWERS, coefficients):
        table[a, b] = coefficient / (factorial(a)*factorial(b))
    return table


def independent_functionals(coefficients, tangents, normals, slip, alpha):
    table = polynomial(coefficients)
    gradient = [polyder(table, axis=k) for k in range(2)]
    laplacian = np.zeros_like(table)
    laplacian[:-2, :] += polyder(table, m=2, axis=0)
    laplacian[:, :-2] += polyder(table, m=2, axis=1)
    lap_gradient = [polyder(laplacian, axis=k) for k in range(2)]
    output = [laplacian[0, 0]]
    for derivatives in (gradient, lap_gradient):
        for t, n in zip(tangents, normals):
            # Different nodes from the coefficient generator. Simpson 3/8 is
            # exact for the cubic-or-lower derivative along each ray.
            values = [sum(-n[k]*polyval2d(*(fraction*t/2), derivatives[k]) for k in range(2))
                      for fraction in (0, 1/3, 2/3, 1)]
            output.append(np.dot([1, 3, 3, 1], values)*np.linalg.norm(t)/16)
    # Check the ray laws independently using the polynomial derivatives.
    for t, n in zip(tangents, normals):
        first = gradient[0]*n[0]
        second = gradient[1]*n[1]
        for r in (0.0, 0.25, 1.1):
            x, y = r*t
            value = polyval2d(x, y, table)
            assert abs(value) < 2e-11
            velocity = polyval2d(x, y, first)+polyval2d(x, y, second)
            if slip == "no-slip":
                assert abs(velocity) < 2e-11
            else:
                normal2 = sum(n[i]*n[j]*polyval2d(x, y, polyder(gradient[i], axis=j))
                              for i in range(2) for j in range(2))
                # A p=4 expansion enforces degrees 0..2 of the relation.
                # Remove the allowed cubic term in alpha*d_n before testing.
                top = table.copy()
                for a, b in POWERS[TOTAL < DEGREE]:
                    top[a, b] = 0
                cubic = sum(n[i]*polyval2d(x, y, polyder(top, axis=i)) for i in range(2))
                assert abs(normal2-alpha*(velocity-cubic)) < 2e-11
    return output


def smooth_experiment():
    report = []
    for slip, alpha in (("free", 0), ("partial", 0.7), ("no-slip", 0)):
        for sectors in range(1, 8):
            for start in range(8):
                xy, radius, basis, fit, condition = fit_corner(start, sectors, slip, alpha)
                tangents, normals = geometry(start, sectors)
                rows = functionals(tangents, normals) @ fit
                residual = 0.0
                for coeff in basis.T:
                    known = independent_functionals(coeff, tangents, normals, slip, alpha)
                    data = polyval2d(xy[:, 0], xy[:, 1], polynomial(coeff))
                    # Shift by a nonzero C to exercise gauge invariance.
                    actual = rows @ ((0.375+data)-0.375)
                    residual = max(residual, np.max(abs(actual-known)))
                    assert residual < 1e-10
                    if sectors != 4:
                        assert abs(known[0]) < 1e-12
                report.append(dict(slip=slip, sectors=sectors, start=start, radius=radius,
                                   dofs=basis.shape[1], condition=condition, residual=residual,
                                   sites=xy.tolist(), weights=rows.tolist(),
                                   wall_value_weights=(-rows.sum(axis=1)).tolist()))
    by_key = {(item["slip"], item["sectors"], item["start"]): item for item in report}
    for item in report:
        for reflected in (False, True):
            start, sectors = item["start"], item["sectors"]
            start = (-(start+sectors) if reflected else start+2) % 8
            other = by_key[(item["slip"], sectors, start)]
            assert other["radius"] == item["radius"]
            a = {tuple(p): w for p, w in zip(other["sites"], np.array(other["weights"]).T)}
            ordering = [0, 2, 1, 4, 3] if reflected else [0, 1, 2, 3, 4]
            for (x, y), w in zip(item["sites"], np.array(item["weights"]).T):
                assert np.allclose(a[(x,-y) if reflected else (-y,x)], w[ordering], atol=1e-10, rtol=0)
    return report


def singular_experiment():
    """Free-slip harmonic wedge modes, whose integrated J is known exactly."""
    report = []
    for sectors in (5, 6, 7):
        exponent = 4 / sectors
        for start in range(8):
            xy = samples(start, sectors, 3)
            r = np.linalg.norm(xy, axis=1)
            theta = (np.arctan2(xy[:, 1], xy[:, 0])-start*np.pi/4) % (2*np.pi)
            tangents, _ = geometry(start, sectors)
            ends = np.linalg.norm(tangents, axis=1)/2
            powers = exponent * np.arange(1, 5)
            design = r[:, None]**powers * np.sin(theta[:, None]*powers)
            integrals = np.vstack([-ends[0]**powers, ends[1]**powers * [-1, 1, -1, 1]])
            for retained in (2, 3):
                assert np.linalg.matrix_rank(design[:, :retained]) == retained
                rows = integrals[:, :retained] @ np.linalg.pinv(design[:, :retained])
                assert np.max(abs(rows @ design[:, :retained] - integrals[:, :retained])) < 2e-14
                # Compare with a manufactured field including the first
                # omitted harmonic, keeping the physical amplitudes fixed.
                errors = []
                for h in (0.2, 0.1, 0.05):
                    coefficients = (h**powers * [1, 0.3, 0.1, 0.07])[:retained+1]
                    exact = integrals[:, :retained+1] @ coefficients
                    errors.append(float(np.max(abs(rows @ (design[:, :retained+1] @ coefficients)-exact))))
                rate = float(np.log2(errors[-2]/errors[-1]))
                assert abs(rate - (retained+1)*exponent) < 1e-8
                report.append(dict(sectors=sectors, start=start, exponent=exponent, retained=retained,
                                   rate=rate, errors=errors, sites=xy.tolist(), weights=rows.tolist()))
    return report


def biharmonic_experiment():
    """No-slip wedge: test a finite part against a volume/open-face identity."""
    quadrature, weights = np.polynomial.legendre.leggauss(24)
    result = []
    for sectors in (5, 6, 7):
        beta = sectors*np.pi/4
        low, high = 1e-8, np.pi/beta
        for _ in range(55):
            mid = (low+high)/2
            if np.sin(beta*mid)+mid*np.sin(beta) > 0:
                low = mid
            else:
                high = mid
        lam = 1+(low+high)/2
        p = lam-2
        a, b = np.cos(p*beta/2), -np.cos(lam*beta/2)
        amplitude = b*(lam**2-p**2)
        assert abs(a*np.cos(lam*beta/2)+b*np.cos(p*beta/2)) < 1e-15
        assert abs(a*lam*np.sin(lam*beta/2)+b*p*np.sin(p*beta/2)) < 1e-14

        for start in range(8):
            centre = start*np.pi/4 + beta/2
            tangents, normals = geometry(start, sectors)
            ends = np.linalg.norm(tangents, axis=1)/2
            vertices = RAYS[(start+np.arange(sectors+1)) % 8]/2
            smooth_cross = normals[0] @ normals[1]

            def smooth_gradient(points):
                u, v = (points @ normals.T).T
                return (4*u[:, None]*normals[0] + 4*v[:, None]*normals[1] +
                        8*smooth_cross*(v[:, None]*normals[0]+u[:, None]*normals[1]))

            def singular_gradient(points):
                r = np.linalg.norm(points, axis=1)
                phi = (np.arctan2(points[:, 1], points[:, 0])-centre+np.pi) % (2*np.pi)-np.pi
                radial = points/r[:, None]
                angular = radial[:, ::-1]*[-1, 1]
                return amplitude*p*r[:, None]**(p-1)*(
                    np.cos(p*phi)[:, None]*radial-np.sin(p*phi)[:, None]*angular)

            flux = np.zeros(2)
            for v0, v1 in zip(vertices[:-1], vertices[1:]):
                side = v1-v0
                outward = side[::-1]*[1, -1]
                points = v0 + (quadrature[:, None]+1)*side/2
                flux += [weights @ (gradient(points) @ outward)/2
                         for gradient in (singular_gradient, smooth_gradient)]
            volume = np.array([0, (8+16*smooth_cross**2)*sectors/8])
            weak = volume-flux

            # Directly integrate the truncated wall rays and the inner circular
            # cap. The two wall integrals alone diverge as the cap shrinks.
            raw_torques, residuals = [], []
            for epsilon in (1/16, 1/64, 1/256):
                raw = 0.0
                for end, ray, normal in zip(ends, tangents, normals):
                    lo, hi = np.log(epsilon), np.log(end)
                    r = np.exp((lo+hi)/2+(hi-lo)*quadrature/2)
                    points = r[:, None]*ray/np.linalg.norm(ray)
                    raw += (weights*r) @ (singular_gradient(points) @ -normal)*(hi-lo)/2
                angles = centre + beta/2*quadrature
                radial = np.column_stack([np.cos(angles), np.sin(angles)])
                inner = weights @ np.sum(singular_gradient(epsilon*radial)*(-radial), axis=1)*epsilon*beta/2
                raw_torques.append(float(raw))
                residuals.append(float(abs(raw+inner-weak[0])))
            assert abs(raw_torques[-1]) > abs(raw_torques[0])
            assert max(residuals) < 1e-11

            smooth_wall = sum(end*np.dot(smooth_gradient((ray/4)[None, :])[0],-normal)
                              for end, ray, normal in zip(ends, tangents, normals))
            assert abs(smooth_wall-weak[1]) < 1e-13
            xy = samples(start, sectors, 3)
            phi = (np.arctan2(xy[:, 1], xy[:, 0])-centre+np.pi) % (2*np.pi)-np.pi
            singular = np.linalg.norm(xy, axis=1)**lam*(a*np.cos(lam*phi)+b*np.cos(p*phi))
            smooth = np.prod((xy @ normals.T)**2, axis=1)
            design = np.column_stack([singular, smooth])
            assert np.linalg.matrix_rank(design) == 2
            row = weak @ np.linalg.pinv(design)
            assert np.max(abs(row @ design-weak)) < 1e-12
            result.append(dict(sectors=sectors, start=start, exponent=float(lam), raw_torques=raw_torques,
                               identity_residuals=residuals, weak_torque=weak.tolist(),
                               sites=xy.tolist(), weights=row.tolist()))
    return result


def finite_slip_experiment():
    """An exact finite-slip singular field with integrable biharmonic source."""
    nodes, weights = np.polynomial.legendre.leggauss(24)
    alpha = 0.7
    report = []
    for sectors in (5, 6, 7):
        beta = sectors*np.pi/4
        lam = np.pi/beta
        for start in range(8):
            angle = start*np.pi/4
            tangents, normals = geometry(start, sectors)
            a = np.linalg.solve(normals, np.full(2,alpha/2))
            a2 = a@a
            rotation = np.exp(-1j*angle)
            ac = complex(*a)*rotation
            vertices = RAYS[(start+np.arange(sectors+1)) % 8]/2

            def fields(r, theta):
                directions = np.column_stack([np.cos(angle+theta),np.sin(angle+theta)])
                linear = r*(directions@a)
                # Keep the branch angle across the entire reentrant sector.
                value = r**lam*np.sin(lam*theta)
                first = lam*r**(lam-1)*np.exp(1j*(lam-1)*theta)
                second = lam*(lam-1)*r**(lam-2)*np.exp(1j*(lam-2)*theta)
                omega = np.exp(linear)*(2*np.imag(ac*first)+a2*value)
                return linear, value, first, second, omega

            def flux(r, theta, normal):
                linear, _, first, second, omega = fields(r,theta)
                nc = (normal[...,0]+1j*normal[...,1])*rotation
                return omega*(normal@a)+np.exp(linear)*np.imag(nc*(2*ac*second+a2*first))

            def radial_integral(linear, power):
                # Integral t**(power-1)*exp(t*linear), t in [0,1].
                return sum(linear**m/(factorial(m)*(power+m)) for m in range(24))

            def volume_radial(r, theta):
                linear, value, first, second, _ = fields(r,theta)
                return (4*np.imag(ac**2*second)*radial_integral(linear,lam) +
                        4*a2*np.imag(ac*first)*radial_integral(linear,lam+1) +
                        a2*a2*value*radial_integral(linear,lam+2))

            face_integral = 0.0
            volume = 0.0
            for left, right in zip(vertices[:-1],vertices[1:]):
                side = right-left
                points = left+(nodes[:,None]+1)*side/2
                r = np.linalg.norm(points,axis=1)
                theta = (np.arctan2(points[:,1],points[:,0])-angle) % (2*np.pi)
                face_integral += weights @ flux(r,theta,side[::-1]*[1,-1])/2
                cross = left[0]*right[1]-left[1]*right[0]
                volume += cross*(weights @ volume_radial(r,theta))/2
            weak = volume-face_integral

            raw_torques, residuals = [], []
            for epsilon in (1/16,1/64,1/256):
                raw = 0.0
                for end, normal, theta_wall in zip(np.linalg.norm(tangents,axis=1)/2,normals,(0,beta)):
                    lo, hi = np.log(epsilon), np.log(end)
                    r = np.exp((lo+hi)/2+(hi-lo)*nodes/2)
                    theta = np.full_like(r,theta_wall)
                    linear, value, first, _, omega = fields(r,theta)
                    velocity = np.exp(linear)*(np.imag(complex(*normal)*rotation*first)+(a@normal)*value)
                    assert np.max(abs(omega-alpha*velocity)) < 1e-12
                    raw += (weights*r) @ flux(r,theta,-normal)*(hi-lo)/2
                theta = beta*(nodes+1)/2
                radial = np.column_stack([np.cos(angle+theta),np.sin(angle+theta)])
                radius = np.full_like(theta,epsilon)
                cap = (weights @ flux(radius,theta,-radial))*epsilon*beta/2
                small_volume = (weights @ volume_radial(radius,theta))*epsilon**2*beta/2
                raw_torques.append(float(raw))
                residuals.append(float(abs(raw+cap+small_volume-weak)))
            assert abs(raw_torques[-1]) > abs(raw_torques[0])
            assert max(residuals) < 1e-11
            report.append(dict(sectors=sectors,start=start,exponent=float(lam),raw_torques=raw_torques,
                               identity_residuals=residuals,weak_torque=float(weak)))
    return report


def discrete_wedge_experiment():
    """Separate analytic sample reproduction from actual five-point inversion."""
    from scipy.sparse.linalg import splu
    from check_cmm_oracle import operator

    exponent = 2/3
    inner, outer = 0.16, 0.32
    stencil = samples(2, 6, 3)
    stencil_radius = np.linalg.norm(stencil, axis=1)
    stencil_theta = (np.arctan2(stencil[:, 1], stencil[:, 0])-np.pi/2) % (2*np.pi)
    stencil_powers = exponent*np.arange(1, 4)
    stencil_design = stencil_radius[:, None]**stencil_powers*np.sin(stencil_theta[:, None]*stencil_powers)
    fit = np.linalg.pinv(stencil_design)
    calibration = []
    for radius in (32, 64, 128):
        row, column = np.indices((2*radius+1, 2*radius+1))
        x, y = column-radius, radius-row
        theta = (np.arctan2(y, x)-np.pi/2) % (2*np.pi)
        exact = np.hypot(x, y)**exponent*np.sin(exponent*theta)
        boundary = (abs(x) == radius) | (abs(y) == radius)
        strict = (~boundary) & ((x < 0) | (y < 0))
        box_values = np.where(boundary & ((x < 0) | (y < 0)), exact, 0)
        coords = np.argwhere(strict)
        rhs = np.array([-sum(box_values[i+di, j+dj] for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)))
                        for i, j in coords])
        matrix, _ = operator(coords, 0)
        lattice = np.zeros_like(exact)
        lattice[strict] = splu(matrix).solve(rhs)
        calibration.append(float(fit[0] @ lattice[radius-stencil[:, 1],radius+stencil[:, 0]]))
    report = []
    for n in (33, 65, 129, 257):
        h = 1/(n-1)
        row, column = np.indices((n, n))
        x, y = column*h-0.5, 0.5-row*h
        r = np.hypot(x, y)
        theta = (np.arctan2(y, x)-np.pi/2) % (2*np.pi)
        strict = ((row > 0) & (row < n-1) & (column > 0) & (column < n-1) &
                  ((x < 0) | (y < 0)))
        coords = np.argwhere(strict)
        singular = r**exponent*np.sin(exponent*theta)
        ramp = (r > inner) & (r < outer)
        t = (r[ramp]-inner)/(outer-inner)
        derivative = (-30*t**2+60*t**3-30*t**4)/(outer-inner)
        second = (-60*t+180*t**2-120*t**3)/(outer-inner)**2
        rhs = np.zeros((n, n))
        rhs[ramp] = h*h*singular[ramp]*(second+(1+2*exponent)*derivative/r[ramp])
        cutoff = np.zeros((n, n))
        cutoff[r <= inner] = 1
        cutoff[ramp] = 1-10*t**3+15*t**4-6*t**5
        matrix, _ = operator(coords, 0)
        actual = np.zeros((n,n))
        actual[strict] = splu(matrix).solve(rhs[strict])
        assert np.max(abs(matrix @ actual[strict]-rhs[strict])) < 5e-15
        mid = (n-1)//2
        patch_data = actual[mid-stencil[:, 1],mid+stencil[:, 0]]
        estimate = (fit @ patch_data)[0]/h**exponent
        annulus = strict & (r >= 0.08) & (r <= 0.14)
        far_estimate = (singular[annulus] @ actual[annulus])/np.sum(singular[annulus]**2)
        report.append(dict(grid=n, solution_error=float(np.max(abs(actual[strict]-(cutoff*singular)[strict]))),
                           radius_three_relative_amplitude_error=float(abs(estimate-1)),
                           calibrated_relative_amplitude_error=float(abs(estimate/calibration[-1]-1)),
                           annulus_relative_amplitude_error=float(abs(far_estimate-1))))
    assert report[-1]["annulus_relative_amplitude_error"] < report[0]["annulus_relative_amplitude_error"]
    assert report[-1]["calibrated_relative_amplitude_error"] < report[-1]["radius_three_relative_amplitude_error"]
    return dict(grids=report, lattice_amplitude_normalization=calibration)


def vertex_order_experiment():
    """Measure the velocity error from switching to an exact wall vertex."""
    degree3 = [(a, d-a) for d in range(4) for a in range(d+1)]
    target = np.r_[1, np.zeros(9)]
    result = []
    for diagonal in (False, True):
        normal = np.array([-1,1])/np.sqrt(2) if diagonal else np.array([0,1])
        tangent = normal[::-1]*[1,-1]

        def exact(points):
            distance = points @ normal
            along = points @ tangent
            return distance**2*(1+0.2*along+0.3*distance**2)

        errors = []
        for h in (1/64,1/128,1/256):
            ys = np.r_[0.5 if diagonal else 0,np.arange(1 if diagonal else 0,round(0.2/h))+0.5]
            vertices = np.column_stack([np.full(len(ys),0.5),ys])
            raw = [0.0]
            corrected = [0.0]
            for vertex in vertices[1:]:
                cell = np.floor(vertex)
                regular = cell+np.array([[0,0],[0,1],[1,0],[1,1]])
                average = np.mean(exact(regular*h))
                distance = (vertex@normal)*h
                raw.append(average)
                blend = min(1.0,distance/0.08)
                if blend == 1.0:
                    corrected.append(average)
                    continue
                donors = cell+np.array([(i,j) for i in range(-2,4) for j in range(-2,4)])
                donors = donors[donors@normal >= -1e-14]
                x,y = (donors-vertex).T
                moments = np.array([x**a*y**b for a,b in degree3])
                weights, _, rank, _ = np.linalg.lstsq(moments,target,rcond=None)
                assert rank == len(degree3) and np.max(abs(moments@weights-target)) < 1e-12
                cubic = weights@exact(donors*h)
                corrected.append((1-blend)*cubic+blend*average)
            reference = exact(vertices*h)
            face_length = h*np.diff(ys)
            errors.append([np.max(abs(np.diff(value-reference)/face_length))
                           for value in (np.array(raw),np.array(corrected))])
        orders = np.log2(np.array(errors[:-1])/errors[1:])
        assert np.all((orders[:,0] > 0.9) & (orders[:,0] < 1.1))
        assert np.all(orders[:,1] > 1.8)
        result.append(dict(diagonal=diagonal,errors=np.array(errors).tolist(),orders=orders.tolist()))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--discrete", action="store_true", help="also run the independent SciPy discrete wedge solve")
    args = parser.parse_args()
    smooth = smooth_experiment()
    singular = singular_experiment()
    biharmonic = biharmonic_experiment()
    finite_slip = finite_slip_experiment()
    discrete = discrete_wedge_experiment() if args.discrete else {}
    vertices = vertex_order_experiment()
    args.report.write_text(json.dumps({"smooth": smooth, "free_slip_wedge": singular,
                                       "no_slip_wedge": biharmonic, "finite_slip_wedge": finite_slip,
                                       "discrete_wedge": discrete, "vertex_order": vertices}, indent=2) + "\n")
    print(f"Smooth patches: {len(smooth)}, max condition {max(x['condition'] for x in smooth):.3g}, "
          f"max moment residual {max(x['residual'] for x in smooth):.3e}")
    for retained in (2, 3):
        rates = [x["rate"] for x in singular if x["retained"] == retained]
        print(f"Free-slip singular patches with {retained} modes: {len(rates)}, "
              f"integrated error orders {min(rates):.3f}–{max(rates):.3f}")
    print(f"No-slip singular patches: {len(biharmonic)}, max truncated-contour identity residual "
          f"{max(max(x['identity_residuals']) for x in biharmonic):.3e}")
    print(f"Finite-slip singular patches: {len(finite_slip)}, max truncated-contour identity residual "
          f"{max(max(x['identity_residuals']) for x in finite_slip):.3e}")
    for row in discrete.get("grids", []):
        print("Discrete 270-degree wedge:", row)
    if discrete:
        print("Radius 32/64/128 lattice calibrations:", discrete["lattice_amplitude_normalization"])
    for item in vertices:
        print("Straight-wall vertex velocity orders, original and blended cubic:", item)


if __name__ == "__main__":
    main()
