"""Independent monolithic steady Stokes reference for the C1 weak wall equations."""

import argparse
from dataclasses import dataclass
import json
from math import comb
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
from scipy.sparse import bmat, coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import splu


RAYS = np.array([(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]) / 2
POWERS = np.array([(a, k-a) for k in range(4) for a in range(k, -1, -1)])


def monomials(points):
    return np.prod(points[..., None, :]**POWERS, axis=-1)


@dataclass
class Reference:
    n: int
    concave: bool
    islands: bool = False
    notch_slope: int = 0

    def __post_init__(self):
        self.h = 1/(self.n-1)
        self.centre = (self.n-1)//2
        self.boxes = []
        if self.islands:
            radius = (self.n-1)/16
            for x in (0.25,0.75):
                centre = (self.n-1)*np.array([x,0.5])
                self.boxes.append((centre-radius,centre+radius))
        self.indices = np.array([(x, y) for y in range(self.n) for x in range(self.n)])
        self.offsets, self.weights = self.quadrature(6)
        shape = self.indices[:, None, :] + self.offsets[None, :, :]
        inside = self.wet(shape)
        area = inside @ self.weights
        keep = area > 0
        self.indices = self.indices[keep]
        self.area = area[keep]
        self.strict = np.isclose(self.area, 1, atol=1e-14, rtol=0)
        self.number = {tuple(q): i for i, q in enumerate(self.indices)}
        self.masked_weights = inside[keep]*self.weights/self.area[:, None]
        self.moments = self.masked_weights @ monomials(self.offsets)
        self.face_data = list(self.faces())
        self.stiffness = self.assemble(False)
        self.diffusion = self.assemble(True)

    @staticmethod
    def quadrature(order):
        z, w = np.polynomial.legendre.leggauss(order)
        z, w = (z+1)/2, w/2
        u, v = np.meshgrid(z, z, indexing="ij")
        weight = (u*np.outer(w, w)).ravel()/4
        offsets = [u[..., None]*(p+v[..., None]*(q-p))
                   for p, q in zip(RAYS, np.roll(RAYS,-1,axis=0))]
        return np.reshape(offsets, (-1,2)), np.tile(weight, 8)

    def wet(self, points):
        wet = ((points > 0) & (points < self.n-1)).all(axis=-1)
        if self.concave:
            p = points-self.centre
            wet &= (p[...,1] < 0) | (p[...,0]-self.notch_slope*p[...,1] < 0)
        for lo,hi in self.boxes:
            wet &= ~((points >= lo) & (points <= hi)).all(axis=-1)
        return wet

    def visible(self, origin, candidate):
        for lo,hi in self.boxes:
            a,b = 0.0,1.0
            for x,dx,left,right in zip(origin,candidate-origin,lo,hi):
                if dx == 0:
                    if not left < x < right:
                        a,b = 1.0,0.0
                        break
                else:
                    entry,exit = sorted(((left-x)/dx,(right-x)/dx))
                    a,b = max(a,entry),min(b,exit)
            if a < b:
                return False
        if not self.concave:
            return True
        # Candidate centres and face midpoints lie in the box. Test for an
        # interior intersection with the closed NE notch by splitting at
        # both coordinate crossings. No sampled coastline geometry is used.
        p = origin-self.centre
        q = candidate-origin
        start = np.array([p[1],p[0]-self.notch_slope*p[1]])
        delta = np.array([q[1],q[0]-self.notch_slope*q[1]])
        t = [0.0, 1.0]
        for v, d in zip(start, delta):
            if d and 0 < -v/d < 1:
                t.append(-v/d)
        t.sort()
        return all(np.min(start+(a+b)/2*delta) < 0 for a, b in zip(t[:-1],t[1:]))

    def faces(self):
        for first, p in enumerate(self.indices):
            for normal, tangent in ((np.array([1, 0]), np.array([0, 1])), (np.array([0, 1]), np.array([1, 0]))):
                second = self.number.get(tuple(p+normal))
                if second is None:
                    continue
                midpoint = p+normal/2
                wet = self.wet(midpoint + np.array([[-0.25], [0.25]])*tangent)
                if wet.any():
                    yield first, second, midpoint, normal, tangent, wet

    def assemble(self, cubic):
        size = len(self.indices)
        rows, columns, values = [], [], []
        for i, j, midpoint, normal, tangent, wet in self.face_data:
            length = wet.sum()/2
            if cubic:
                donors, coeff = self.face_row(midpoint, normal, tangent, wet)
            else:
                donors, coeff = [i,j], [-length,length]
            rows.extend([i]*len(donors)+[j]*len(donors))
            columns.extend(donors+donors)
            values.extend(coeff)
            values.extend(-np.asarray(coeff))
        return coo_matrix((values,(rows,columns)),shape=(size,size)).tocsr()

    def face_row(self, origin, normal, tangent, wet):
        points = origin + np.array([[-0.5 if wet[0] else 0], [0.5 if wet[1] else 0]])*tangent
        origin = points.mean(axis=0)
        line = (points-origin).mean(axis=0)
        displacement = (points[1]-points[0])/(2*np.sqrt(3))
        target = np.zeros(len(POWERS))
        for k, (a,b) in enumerate(POWERS):
            if a:
                target[k] += a*normal[0]*np.mean((line[0]+np.array([-1,1])*displacement[0])**(a-1) *
                                                 (line[1]+np.array([-1,1])*displacement[1])**b)
            if b:
                target[k] += b*normal[1]*np.mean((line[0]+np.array([-1,1])*displacement[0])**a *
                                                 (line[1]+np.array([-1,1])*displacement[1])**(b-1))
        target *= wet.sum()/2
        return self.fit_row(origin,target)

    def wall_trace(self):
        rows,columns,values = [],[],[]
        for i in np.flatnonzero(~self.strict):
            p = self.indices[i]
            for tangent in np.array([(1,0),(0,1),(-1,0),(0,-1)]):
                midpoint = p+tangent/4
                normal = tangent[::-1]*[-1,1]
                wet = self.wet(midpoint+np.array([[-0.125],[0.125]])*normal)
                if wet[0] == wet[1]:
                    continue
                offset = np.array([[-1],[1]])*tangent/(4*np.sqrt(3))
                target = 0.5*monomials(offset).mean(axis=0)
                donors,coefficients = self.fit_row(midpoint,target)
                rows.extend([i]*len(donors))
                columns.extend(donors)
                values.extend(coefficients)
        n = len(self.area)
        return coo_matrix((values,(rows,columns)),shape=(n,n)).tocsr()

    def fit_row(self, origin, target):
        centre = np.floor(origin).astype(int)
        for radius in range(2,5):
            donors = []
            for x in range(centre[0]-radius,centre[0]+radius+1):
                for y in range(centre[1]-radius,centre[1]+radius+1):
                    index = self.number.get((x,y))
                    if index is not None and self.visible(origin, np.array([x,y])):
                        donors.append(index)
            if len(donors) < len(POWERS):
                continue
            offsets = self.indices[donors]-origin
            design = np.zeros((len(donors),len(POWERS)))
            for k,(a,b) in enumerate(POWERS):
                for j,(r,s) in enumerate(POWERS):
                    if r <= a and s <= b:
                        design[:,k] += (comb(a,r)*comb(b,s)*offsets[:,0]**(a-r)*offsets[:,1]**(b-s) *
                                        self.moments[donors,j])
            distance = 1/(1+np.sum(offsets**2,axis=1))
            weights, _, rank, _ = np.linalg.lstsq(design.T*distance,target,rcond=1e-6)
            if rank == len(POWERS):
                assert np.max(abs(design.T@(weights*distance)-target)) < 2e-11
                return donors, weights*distance
        raise ValueError(f"insufficient reference donors at {origin}")

    def means(self, function):
        points = self.h*(self.indices[:,None,:]+self.offsets[None,:,:])
        return np.sum(function(points)*self.masked_weights, axis=1)


def smooth_solution(concave):
    p = Polynomial([0, 0, 1, -2, 1])
    if concave:
        p *= Polynomial([0.25, -1, 1])

    def fields(q):
        x, y = q[...,0], q[...,1]
        psi = p(x)*p(y)
        vorticity = p.deriv(2)(x)*p(y)+p(x)*p.deriv(2)(y)
        forcing = p.deriv(4)(x)*p(y)+2*p.deriv(2)(x)*p.deriv(2)(y)+p(x)*p.deriv(4)(y)
        return psi,vorticity,forcing

    return fields


def singular_solution(beta, notch_slope):
    lo,hi = 1e-8, np.pi/beta
    for _ in range(55):
        mid = (lo+hi)/2
        if np.sin(beta*mid)+mid*np.sin(beta) > 0:
            lo = mid
        else:
            hi = mid
    lam = 1+(lo+hi)/2
    lower = lam-2
    a,b = np.cos(lower*beta/2), -np.cos(lam*beta/2)
    c = b*(lam**2-lower**2)
    inner, outer = 0.23,0.4
    chi = Polynomial([1,0,0,0,0,-126,420,-540,315,-70])

    def polar(q):
        point = q-0.5
        r = np.linalg.norm(point,axis=-1)
        angle = (np.arctan2(point[...,1],point[...,0])-(2*np.pi-beta)) % (2*np.pi)-beta/2
        return r,angle

    def fields(q):
        r,angle = polar(q)
        f = a*np.cos(lam*angle)+b*np.cos(lower*angle)
        safe_r = np.where(r > 0,r,1)
        singular = r**lam*f
        w = c*safe_r**lower*np.cos(lower*angle)
        t = np.clip((r-inner)/(outer-inner),0,1)
        derivatives = [chi(t)] + [chi.deriv(k)(t)/(outer-inner)**k*((t>0)&(t<1)) for k in range(1,5)]
        ch, cr, crr, crrr, crrrr = derivatives
        d = crr+(2*lam+1)*cr/safe_r
        dr = crrr+(2*lam+1)*(crr/safe_r-cr/safe_r**2)
        drr = crrrr+(2*lam+1)*(crrr/safe_r-2*crr/safe_r**2+2*cr/safe_r**3)
        rhs = w*(crr+(2*lower+1)*cr/safe_r+d)+singular*(drr+(2*lam+1)*dr/safe_r)
        return ch*singular, ch*w+singular*d, rhs

    def wall_torque(radius):
        # True finite-part wall torque in the smaller square patch from the
        # integrated derivative of analytic vorticity on the four open sides.
        z, w = np.polynomial.legendre.leggauss(30)
        flux = 0.0
        boundary = np.array([[-1,-1],[1,-1],[1,1],[-1,1]])*radius
        for start,end in zip(boundary,np.roll(boundary,-1,axis=0)):
            delta = end-start
            p = [start[1],start[0]-notch_slope*start[1]]
            d = [delta[1],delta[0]-notch_slope*delta[1]]
            points = sorted([0.0,1.0]+[-v/u for v,u in zip(p,d) if u and 0 < -v/u < 1])
            for a,b in zip(points[:-1],points[1:]):
                middle = start+(a+b)/2*delta
                if middle[1] >= 0 and middle[0]-notch_slope*middle[1] >= 0:
                    continue
                step = (b-a)*delta
                normal = [step[1],-step[0]]
                q = start+a*delta+(1+z[:,None])/2*step
                r, angle = polar(q+0.5)
                er = q/r[:,None]
                et = er[:,::-1]*[-1,1]
                derivative = (c*lower*r[:,None]**(lower-1) *
                              (np.cos(lower*angle)[:,None]*er-np.sin(lower*angle)[:,None]*et))
                flux += w @ (derivative@normal)/2
        return -flux

    return fields, wall_torque


def island_solution(aligned=False):
    """Two clamped smooth island flows with unequal constants and torques."""
    z = Polynomial([0,1])
    smoothstep = 1-630*(z**4*(1-z)**4).integ()
    width = 0.1
    half_side = 0.0625
    support = 0.25 if aligned else 0.075
    taper = [smoothstep.deriv(k)/width**k for k in range(5)]
    aligned_step = 140*(z**3*(1-z)**3).integ()
    tangent = (1-z*z)**5
    wall = z*z*(1+z)*(1-z*z)**5
    tangent_deriv = [tangent.deriv(k)/half_side**k for k in range(5)]
    wall_deriv = [wall.deriv(k)*support**(2-k) for k in range(5)]
    levels = np.array([0.6,-0.35])
    amplitudes = np.array([2.0,-1.25])
    exact_torque = -amplitudes*6/support*half_side*(tangent.integ()(1)-tangent.integ()(-1))

    def plateau(x):
        t = (abs(x)-half_side-0.0125)/width
        return [np.where(t<=0,float(k==0),np.where(t<1,p(np.clip(t,0,1)),0.0)) for k,p in enumerate(taper)]

    def aligned_plateau(x, flat_half, edge_half):
        width = edge_half-flat_half
        t = (edge_half-abs(x))/width
        return [np.where(t>=1,float(k==0),np.where(t>0,
                  aligned_step.deriv(k)(np.clip(t,0,1))/width**k,0.0)) for k in range(5)]

    def evaluate(points):
        x,y = points[...,0],points[...,1]
        result = [np.zeros_like(x) for _ in range(3)]
        py = aligned_plateau(y-0.5,half_side,0.5) if aligned else plateau(y-0.5)
        dy = (y-(0.5+half_side))/support
        wy = [np.where((dy>=0)&(dy<1),p(np.clip(dy,0,1)),0.0) for p in wall_deriv]
        for centre,level,amplitude in zip((0.25,0.75),levels,amplitudes):
            px = aligned_plateau(x-centre,half_side,0.25) if aligned else plateau(x-centre)
            dx = (x-centre)/half_side
            tx = [np.where(abs(dx)<1,p(np.clip(dx,-1,1)),0.0) for p in tangent_deriv]
            for coefficient,g,h in ((level,px,py),(amplitude,tx,wy)):
                result[0] += coefficient*g[0]*h[0]
                result[1] += coefficient*(g[2]*h[0]+g[0]*h[2])
                result[2] += coefficient*(g[4]*h[0]+2*g[2]*h[2]+g[0]*h[4])
        return result

    return evaluate,levels,exact_torque


def solve_slip(n, free):
    ref = Reference(n, concave=False)
    alpha = 0.0 if free else 3.0
    strict = np.flatnonzero(ref.strict)
    h = ref.h
    area = ref.area

    def basis(x):
        b = 0.5*alpha
        exponential = np.exp(b*x*(1-x)+1j*np.pi*x)
        log_first = b*(1-2*x)+1j*np.pi
        log_second = -2*b
        return (exponential.imag, (exponential*(log_first**2+log_second)).imag,
                (exponential*(log_first**4+6*log_first**2*log_second+3*log_second**2)).imag)

    def forcing(q):
        gx,gxx,gxxxx = basis(q[...,0])
        gy,gyy,gyyyy = basis(q[...,1])
        return gxxxx*gy+2*gxx*gyy+gx*gyyyy

    trace = ref.wall_trace()
    # Multiply Gamma=-ell*int(omega) by alpha on wall rows. For
    # alpha=0 the same equation is the zero wall vorticity integral.
    scale = np.where(ref.strict,1,alpha)
    mass = diags(scale*h*h*area)+h*trace
    stiffness = -diags(scale)@ref.stiffness[:,strict]
    matrix = bmat([[mass,stiffness],[ref.diffusion[strict,:],None]],format="csc")
    rhs = h*h*area*ref.means(forcing)
    solution,residual,pivot = scaled_solve(matrix,np.r_[np.zeros(len(area)),rhs[strict]])
    psi = np.zeros_like(area)
    psi[strict] = solution[len(area):]
    actual = basis(h*ref.indices[:,0])[0]*basis(h*ref.indices[:,1])[0]
    x = range(n//4,3*n//4)
    values = [psi[[ref.number[(i,j)] for i in x]] for j in (1,2,3,4)]
    inward = (18*values[0]-9*values[1]+2*values[2])/(6*h)
    laplace = (-104*values[0]+114*values[1]-56*values[2]+11*values[3])/(12*h**2)
    # Third-order normal slope and second-order normal second derivative
    # from interior psi alone. Along the wall the tangential second derivative is zero.
    return dict(n=n,psi_error=float(np.max(abs(psi-actual))),
                independent_slip=float(np.max(abs(laplace-alpha*inward))),
                linear_residual=residual,min_abs_scaled_pivot=pivot)


def solve_islands(n, transient=False, aligned=False):
    ref = Reference(n, concave=False, islands=True)
    exact,levels,torque = island_solution(aligned=aligned)
    area = ref.area
    h = ref.h
    strict = np.flatnonzero(ref.strict)
    on_island = np.array([((ref.indices>=lo)&(ref.indices<=hi)).all(axis=1) for lo,hi in ref.boxes])
    assert np.count_nonzero(on_island) == n-1
    trace = csr_matrix(on_island.T.astype(float))
    component = csr_matrix(on_island.astype(float))
    rhs = h*h*area*ref.means(lambda p: exact(p)[2])
    mass = diags(h*h*area)
    strict_stiffness = -ref.stiffness[:,strict]
    trace_stiffness = -ref.stiffness@trace
    open_laplace = ref.diffusion[strict,:]
    component_open = component@ref.diffusion
    # d = Af-D_open omega. Stationary pressure compatibility requires
    # sum(d)+T_parent=0. The independent analytic T_parent is -torque.
    op = bmat([[mass,strict_stiffness,trace_stiffness],
               [open_laplace,None,csr_matrix((len(strict),2))],
               [component_open,csr_matrix((2,len(strict))),csr_matrix((2,2))]], format="csc")
    b = np.r_[np.zeros(len(area)),rhs[strict],component@rhs-torque]
    solution,residual,pivot = scaled_solve(op,b)
    if transient:
        return transient_islands(ref, exact, torque, trace, component, mass, strict_stiffness,
                                 trace_stiffness, rhs, solution[:len(area)])
    recovered = solution[-2:]
    psi = trace@recovered
    psi[strict] = solution[len(area):-2]
    error = np.max(abs(psi-exact(h*ref.indices)[0]))

    # Solve independently with prescribed exact constants, leaving the
    # torque rows out of the linear system. Compare derived wall torque
    # against analytic line integration, not the enforced residual above.
    fixed = bmat([[mass,strict_stiffness],[open_laplace,None]], format="csc")
    predicted,_,_ = scaled_solve(fixed,np.r_[-trace_stiffness@levels,rhs[strict]])
    predicted_torque = component@(rhs-ref.diffusion@predicted[:len(area)])
    slip = []
    for k,(lo,hi) in enumerate(ref.boxes):
        x = round((lo[0]+hi[0])/2)
        sample = [psi[ref.number[x,round(hi[1])+j]] for j in (1,2,3)]
        slip.append(abs(np.dot([-11,18,-9,2],[recovered[k],*sample])/(6*h)))
    return dict(n=n,psi_error=float(error),constant_error=float(np.max(abs(recovered-levels))),
                prescribed_torque_error=float(np.max(abs(predicted_torque-torque))), independent_slip=max(slip),
                linear_residual=residual,min_abs_scaled_pivot=pivot)


def transient_islands(ref, exact, torque, trace, component, mass, strict_stiffness,
                      trace_stiffness, steady_source, initial):
    """Coupled Crank–Nicolson diffusion with stationary no-slip torque rows."""
    strict = np.flatnonzero(ref.strict)
    boundary = ~ref.strict
    num_cells = len(ref.area)
    omega = ref.means(lambda p: exact(p)[1])
    source0 = -ref.h**2*ref.area*omega-steady_source
    # The manufactured parent force has independently specified circulation
    # on each island and curl from the analytic vorticity equation.
    def run(steps):
        state = initial.copy()
        dt_small,dt_large = np.array([0.75,1.25])*0.0625/steps
        factors = {}
        for dt in (dt_small,dt_large):
            left = mass-0.5*dt*ref.diffusion
            right = mass+0.5*dt*ref.diffusion
            matrix = bmat([[mass,strict_stiffness,trace_stiffness],
                           [left[strict],None,csr_matrix((len(strict),2))],
                           [component@left,csr_matrix((2,len(strict))),csr_matrix((2,2))]], format="csc")
            row_scale = 1/np.asarray(abs(matrix).max(axis=1).toarray()).ravel()
            factors[dt] = right,splu((diags(row_scale)@matrix).tocsc()),row_scale

        time = 0.0
        torque_defect = 0.0
        curl_defect = 0.0
        for k in range(steps):
            dt = (dt_small,dt_large)[k%2]
            right,factor,row_scale = factors[dt]
            avg = 0.5*(np.exp(-time)+np.exp(-time-dt))
            source = avg*source0
            parent = -avg*torque
            known = right@state+dt*source
            b = np.r_[np.zeros(num_cells),known[strict],component@known-dt*parent]
            solution = factor.solve(row_scale*b)
            new = solution[:num_cells]
            # Diagnostic wall flux follows from the actual cell mean update.
            # Independent tests below use analytic forcing and compare psi.
            balance = mass@(new-state)/dt-source-ref.diffusion@((state+new)/2)
            torque_defect = max(torque_defect,np.max(abs(component@balance+parent)))
            curl_defect = max(curl_defect,abs(np.sum(balance[boundary])+np.sum(source)))
            state = new
            time += dt
        psi = trace@solution[-2:]
        psi[strict] = solution[num_cells:-2]
        return psi,torque_defect,curl_defect

    reference,_,_ = run(128)
    results = [run(steps) for steps in (8,16,32)]
    temporal_errors = [float(np.max(abs(psi-reference))) for psi,_,_ in results]
    orders = np.log2(np.asarray(temporal_errors[:-1])/temporal_errors[1:])
    psi_error = np.max(abs(reference-np.exp(-0.0625)*exact(ref.h*ref.indices)[0]))
    torque_defect = max(row[1] for row in results)
    curl_defect = max(row[2] for row in results)
    assert max(torque_defect,curl_defect) < 1e-8
    assert np.all(orders>1.8), orders
    return dict(n=ref.n,psi_error=float(psi_error),temporal_errors=temporal_errors,
                temporal_orders=orders.tolist(),torque_defect=float(torque_defect),curl_defect=float(curl_defect))


def scaled_solve(op,b):
    row_scale = 1/np.asarray(abs(op).max(axis=1).toarray()).ravel()
    factor = splu((diags(row_scale)@op).tocsc())
    solution = factor.solve(row_scale*b)
    residual = np.max(abs(op@solution-b))
    assert residual < 1e-10
    return solution,float(residual),float(np.min(abs(factor.U.diagonal())))


def solve(case, n):
    if case in ("islands","island-transient","islands-aligned"):
        return solve_islands(n,transient=case == "island-transient",aligned=case == "islands-aligned")
    if case in ("finite-slip", "free-slip"):
        return solve_slip(n, case == "free-slip")
    slope = {"singular225":-1,"singular315":1}.get(case,0)
    ref = Reference(n, concave=case != "rectangle", notch_slope=slope)
    if case.startswith("singular"):
        beta = {"singular225":1.25,"singular":1.5,"singular315":1.75}[case]*np.pi
        exact,local_torque = singular_solution(beta,slope)
    else:
        exact,local_torque = smooth_solution(ref.concave),None
    strict = np.flatnonzero(ref.strict)
    strict_count = len(strict)
    area = ref.area
    h = ref.h
    psi_exact, _, _ = exact(ref.indices*h)
    rhs = h*h*area*ref.means(lambda p: exact(p)[2])
    # Unknowns are all cell-mean omega and strict psi. The two row
    # blocks are A h² omega - B psi = 0, D_open omega = A h² f on F.
    # Homogeneous no-slip sets the weak Poisson boundary remainder to zero.
    op = bmat([[diags(h*h*area), -ref.stiffness[:,strict]],
               [ref.diffusion[strict,:], None]], format="csc")
    b = np.r_[np.zeros(len(area)),rhs[strict]]
    solution,residual,pivot = scaled_solve(op,b)
    omega, psi_f = solution[:-strict_count], solution[-strict_count:]
    psi = np.zeros_like(area)
    psi[strict] = psi_f
    wall_residual = rhs-ref.diffusion@omega
    assert abs(np.sum(rhs)-np.sum(wall_residual[~ref.strict])) < 1e-10
    i1 = [ref.number[(x,1)] for x in range(n//4,3*n//4)]
    i2 = [ref.number[(x,2)] for x in range(n//4,3*n//4)]
    i3 = [ref.number[(x,3)] for x in range(n//4,3*n//4)]
    slip = np.max(abs((18*psi[i1]-9*psi[i2]+2*psi[i3])/(6*h)))
    if ref.concave:
        # Inward normal at the eastward reentrant ray, away from its tip.
        # Use independently differenced interior samples, including within
        # the nonzero part of the singular manufactured streamfunction.
        lo = ref.centre+max(1,round(0.07/h))
        hi = ref.centre+round(0.17/h)
        samples = [np.array([ref.number[(x,ref.centre-k)] for x in range(lo,hi+1)])
                   for k in (1,2,3)]
        slip = max(slip,np.max(abs((18*psi[samples[0]]-9*psi[samples[1]]+2*psi[samples[2]])/(6*h))))
    result = dict(n=n, psi_error=float(np.max(abs(psi-psi_exact))), independent_slip=float(slip),
                  linear_residual=residual, min_abs_scaled_pivot=pivot)
    if local_torque is not None:
        m = int(0.125/h)
        patch = (np.abs(ref.indices-ref.centre) <= m).all(axis=1)
        measured = np.sum(wall_residual[patch & ~ref.strict])
        true = local_torque((m+0.5)*h)
        result.update(corner_torque=float(measured), exact_corner_torque=float(true),
                      relative_corner_torque_error=float(abs(measured-true)/abs(true)))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--case", choices=("rectangle","concave","singular","singular225","singular315",
                                 "islands","island-transient","islands-aligned","finite-slip","free-slip"), nargs="+",
                        default=["rectangle","concave","singular"])
    parser.add_argument("--grids", type=int, nargs="+", default=[17,33,65])
    args = parser.parse_args()
    result = {}
    for case in args.case:
        result[case] = []
        for n in args.grids:
            row = solve(case,n)
            result[case].append(row)
            print(case,row,flush=True)
    args.report.write_text(json.dumps(result,indent=2)+"\n")


if __name__ == "__main__":
    main()
