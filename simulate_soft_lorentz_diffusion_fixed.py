
import numpy as np

r0 = 1.0
sigma = 0.01
w = 0.05
L = 2.0*r0 + w

a1 = np.array([L, 0.0])
a2 = np.array([L/2.0, np.sqrt(3.0)*L/2.0])
A = np.column_stack((a1, a2))
Ainv = np.linalg.inv(A)

def V_fermi(r, r0=r0, sigma=sigma):
    return 1.0 / (1.0 + np.exp((r - r0)/sigma))

def dVdr_fermi(V, sigma=sigma):
    return -(1.0/sigma) * V * (1.0 - V)

def lattice_points(nrange=2):
    pts = []
    for i in range(-nrange, nrange+1):
        for j in range(-nrange, nrange+1):
            pts.append(i*a1 + j*a2)
    return np.array(pts, dtype=float)

LP = lattice_points(2)

def potential_at(r):
    r = np.atleast_2d(r)
    Vtot = np.zeros(r.shape[0])
    for R in LP:
        d = r - R
        dist = np.linalg.norm(d, axis=1)
        Vtot += V_fermi(dist)
    return Vtot

def force_at(r):
    r = np.atleast_2d(r)
    F = np.zeros_like(r)
    for R in LP:
        d = r - R
        dist = np.linalg.norm(d, axis=1) + 1e-14
        V = V_fermi(dist)
        dVdr = dVdr_fermi(V)
        coeff = -dVdr / dist
        F += (coeff[:,None] * d)
    return F

def wrap_to_cell(r):
    uv = (Ainv @ r.T).T
    n = np.floor(uv)
    uv_wrapped = uv - n
    r_wrapped = (A @ uv_wrapped.T).T
    shifts = n[...,0][:,None]*a1 + n[...,1][:,None]*a2
    return r_wrapped, shifts

def sample_positions(N, E, max_tries=10000, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    u = rng.random((N,)); v = rng.random((N,))
    r = u[:,None]*a1 + v[:,None]*a2
    V = potential_at(r)
    tries = 0
    mask = (E - V) <= 0.0
    while np.any(mask):
        M = np.sum(mask)
        u = rng.random((M,)); v = rng.random((M,))
        r_new = u[:,None]*a1 + v[:,None]*a2
        V_new = potential_at(r_new)
        r[mask] = r_new; V[mask] = V_new
        mask = (E - V) <= 0.0
        tries += 1
        if tries > max_tries:
            raise RuntimeError("Could not place enough initial positions for E={}".format(E))
    return r, V

def sample_velocities(N, E, V_at_r, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    K = E - V_at_r
    vmag = np.sqrt(2.0 * K)
    phi = 2.0*np.pi*rng.random((N,))
    return np.column_stack((vmag*np.cos(phi), vmag*np.sin(phi)))

def simulate_msd(E, N=200, dt=1e-3, steps=20000, sample_stride=20, rng=None):
    if rng is None:
        rng = np.random.default_rng(2025)
    r, V = sample_positions(N, E, rng=rng)
    v = sample_velocities(N, E, V, rng=rng)
    F = force_at(r)

    # FIX: keep cumulative lattice shift and absolute positions
    L_cum = np.zeros_like(r)
    r_abs = r.copy()
    r_abs0 = r_abs.copy()

    times = [0.0]; msd = [0.0]
    for k in range(1, steps+1):
        v_half = v + 0.5*dt*F
        r_new = r + dt*v_half
        r_wrapped, shifts = wrap_to_cell(r_new)
        L_cum += shifts
        r = r_wrapped
        F_new = force_at(r)
        v = v_half + 0.5*dt*F_new
        F = F_new

        if (k % sample_stride)==0:
            r_abs = r + L_cum
            t = k*dt
            disp = r_abs - r_abs0
            msd_t = np.mean(np.sum(disp*disp, axis=1))
            times.append(t); msd.append(msd_t)
    return np.array(times), np.array(msd)

def estimate_D_from_msd(times, msd, frac_window=0.5):
    n = len(times)
    start = int(n*(1.0 - frac_window))
    if start < 1: start = 1
    t_fit = times[start:]; m_fit = msd[start:]
    A = np.vstack([t_fit, np.ones_like(t_fit)]).T
    slope, intercept = np.linalg.lstsq(A, m_fit, rcond=None)[0]
    D = slope/4.0
    resid = m_fit - (slope*t_fit + intercept)
    s_err = np.sqrt(np.mean(resid**2))
    denom = np.sum((t_fit - np.mean(t_fit))**2) + 1e-12
    slope_err = s_err / np.sqrt(denom)
    D_err = slope_err / 4.0
    return D, D_err, slope, intercept
