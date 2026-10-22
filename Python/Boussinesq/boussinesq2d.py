"""
boussinesq2d.py — a small 2-D (x–z) Boussinesq / anelastic model in pure NumPy.

Formulation
-----------
Vorticity–streamfunction form.  With a reference density rho0(z) (constant for
Boussinesq) the velocity is derived from a *mass* streamfunction psi:

    rho0 u' = -dpsi/dz ,     rho0 w = +dpsi/dx          (so that div(rho0 u) = 0)

The total horizontal wind is U(z) + u', where U(z) is an optional prescribed
background wind.  The prognostic variables are the perturbation vorticity

    eta = du'/dz - dw/dx

and the buoyancy b (= g*theta'/theta0, or -g*rho'/rho0 in Boussinesq).  The
equations solved are

    d eta/dt = -(1/rho0) div(rho0 u eta) + eta w (dln rho0/dz)
               - w (U_zz - U_z dln rho0/dz) - db/dx
               + nu  Lap(eta) - sigma(x,z) eta

    d b/dt   = -(1/rho0) div(rho0 u b) - w N^2(z)
               + kappa Lap(b) - sigma(x,z) b

    d/dz( (1/rho0) dpsi/dz ) + (1/rho0) d^2psi/dx^2 = -eta

where u = (U + u', w) is the full advecting velocity, nu/kappa are optional
Laplacian diffusivities, and sigma(x,z) is a Rayleigh "sponge" damping rate that
is zero in the interior box |x|<=S, 0<=z<=D and ramps linearly to sigma_max at
the edges of the domain.  Set rho0 = None (or a constant) for Boussinesq.

Boundary conditions
-------------------
  x : periodic on [-L, L).
  z = 0 : w = f(x, t)   (user-supplied callable, x-mean removed if nonzero)
          implemented as psi(x,0,t) = rho0(0) * int f dx   (spectral integral)
  z = H : w = 0          (psi = 0)
  eta = b = 0 on z = 0 and z = H (inflowing fluid is unperturbed; the sponge
          should normally cover the top anyway).

Numerics
--------
  * Colocated grid: x_i = -L + i*dx (i = 0..Nx-1), z_k = k*dz (k = 0..Nz).
    Fields are arrays of shape (Nz+1, Nx)  ->  index [k, i].
  * Second-order centered differences.  Advection uses the skew-symmetric
    (half flux-form, half advective-form) discretisation, which conserves both
    the first and second moments of eta and b on this grid, so there is no
    implicit numerical dissipation — use nu/kappa and the sponge for that.
  * Elliptic solve: FFT in x, tridiagonal (Thomas) solve in z, vectorised over
    wavenumbers.  The x-eigenvalue is that of the discrete 3-point second
    difference so the solve is consistent with the finite differences.
  * Time stepping: 3-stage strong-stability-preserving Runge–Kutta (SSPRK3).

Typical use
-----------
    from boussinesq2d import Boussinesq2D
    m = Boussinesq2D(Nx=256, Nz=128, L=50e3, H=20e3, N2=1e-4, U=10.0,
                     w_bottom=lambda x, t: ...,
                     sponge=dict(S=40e3, D=12e3, sigma_max=0.005))
    dt = 0.5 * m.suggest_dt()
    out = m.run(t_end=3600*4, dt=dt, output_every=600)   # dict of arrays
    # ...or drive it yourself:
    for n in range(1000):
        m.step(dt)
        if n % 50 == 0:
            u, w = m.u, m.w   # diagnostics are refreshed after every step
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Callable, Optional, Union

import numpy as np

ArrayLike = Union[float, np.ndarray, Callable[[np.ndarray], np.ndarray]]


# ----------------------------------------------------------------------------
# Grid
# ----------------------------------------------------------------------------
@dataclass
class Grid:
    Nx: int
    Nz: int
    L: float
    H: float

    def __post_init__(self):
        self.dx = 2.0 * self.L / self.Nx
        self.dz = self.H / self.Nz
        self.x = -self.L + self.dx * np.arange(self.Nx)          # periodic
        self.z = self.dz * np.arange(self.Nz + 1)                # includes 0 and H
        self.X, self.Z = np.meshgrid(self.x, self.z)             # shape (Nz+1, Nx)
        # Wavenumbers for rfft along x
        self.kx = 2.0 * np.pi * np.fft.rfftfreq(self.Nx, d=self.dx)
        # Eigenvalue of the discrete second difference  (q[i+1]-2q[i]+q[i-1])/dx^2
        self.lam = -(2.0 - 2.0 * np.cos(self.kx * self.dx)) / self.dx**2


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def _as_profile(val: Optional[ArrayLike], z: np.ndarray, default: float) -> np.ndarray:
    """Turn a scalar / array / callable into a 1-D profile on the z nodes."""
    if val is None:
        return np.full_like(z, default, dtype=float)
    if callable(val):
        return np.asarray(val(z), dtype=float) * np.ones_like(z)
    val = np.asarray(val, dtype=float)
    if val.ndim == 0:
        return np.full_like(z, float(val))
    if val.shape != z.shape:
        raise ValueError(f"profile has shape {val.shape}, expected {z.shape}")
    return val


# ----------------------------------------------------------------------------
# Model
# ----------------------------------------------------------------------------
class Boussinesq2D:
    def __init__(
        self,
        Nx: int,
        Nz: int,
        L: float,
        H: float,
        N2: ArrayLike = 1.0e-4,
        U: Optional[ArrayLike] = None,
        rho0: Optional[ArrayLike] = None,
        w_bottom: Optional[Callable[[np.ndarray, float], np.ndarray]] = None,
        nu: float = 0.0,
        kappa: float = 0.0,
        sponge: Optional[dict] = None,
    ):
        """
        Parameters
        ----------
        Nx, Nz : number of x points (periodic) and z *intervals* (Nz+1 levels).
        L, H   : x in [-L, L), z in [0, H].
        N2     : buoyancy frequency squared; scalar, array on z nodes, or f(z).
        U      : background wind U(z); scalar, array on z nodes, or f(z).  None -> 0.
        rho0   : reference density; None or scalar -> Boussinesq.  Array/f(z) -> anelastic.
        w_bottom : f(x, t) -> array of shape (Nx,), the vertical velocity at z = 0.
                   None -> w = 0 at the bottom.
        nu, kappa : Laplacian diffusivities for eta and b (m^2/s).  0 = none.
        sponge : dict(S=..., D=..., sigma_max=..., combine='max'|'sum').
                 Damping is 0 for |x|<=S, z<=D and increases linearly to
                 sigma_max at |x|=L and at z=H.  None -> no sponge.
        """
        self.g = Grid(Nx, Nz, L, H)
        g = self.g
        self.t = 0.0
        self.nstep = 0

        # --- background profiles (on z nodes) -------------------------------
        self.N2 = _as_profile(N2, g.z, 0.0)
        self.U = _as_profile(U, g.z, 0.0)
        self.rho0 = _as_profile(rho0, g.z, 1.0)
        if np.any(self.rho0 <= 0):
            raise ValueError("rho0 must be positive")
        self.anelastic = not np.allclose(self.rho0, self.rho0[0])

        # d ln rho0 / dz, U_z, U_zz  (centered inside, one-sided at the ends)
        self.dlnrho = np.gradient(np.log(self.rho0), g.dz)
        self.Uz = np.gradient(self.U, g.dz)
        self.Uzz = np.gradient(self.Uz, g.dz)
        # the combination that multiplies -w in the vorticity equation
        self._Ucoef = (self.Uzz - self.Uz * self.dlnrho)[1:-1, None]

        # --- bottom boundary condition --------------------------------------
        self.w_bottom = w_bottom
        self._warned_mean = False

        # --- diffusion -----------------------------------------------------
        self.nu = float(nu)
        self.kappa = float(kappa)

        # --- sponge --------------------------------------------------------
        self.sigma = np.zeros_like(g.X)
        self.sponge = sponge
        if sponge is not None:
            S = sponge.get("S", L)
            D = sponge.get("D", H)
            smax = sponge.get("sigma_max", 0.0)
            combine = sponge.get("combine", "max")
            rx = np.clip((np.abs(g.X) - S) / max(L - S, 1e-300), 0.0, 1.0)
            rz = np.clip((g.Z - D) / max(H - D, 1e-300), 0.0, 1.0)
            if combine == "max":
                self.sigma = smax * np.maximum(rx, rz)
            elif combine == "sum":
                self.sigma = smax * (rx + rz)
            else:
                raise ValueError("sponge['combine'] must be 'max' or 'sum'")

        # --- elliptic solver set-up ----------------------------------------
        self._setup_elliptic()

        # --- prognostic state ----------------------------------------------
        self.eta = np.zeros_like(g.X)
        self.b = np.zeros_like(g.X)
        # diagnostics (filled by _diagnose)
        self.psi = np.zeros_like(g.X)
        self.u = np.zeros_like(g.X)        # perturbation u'
        self.w = np.zeros_like(g.X)
        self._diagnose()

    # ------------------------------------------------------------------
    # Elliptic solver:  d/dz( a dpsi/dz ) + a lam psi = rhs,   a = 1/rho0
    # Dirichlet psi(0) = psi_bot(x), psi(H) = 0.  Unknowns k = 1..Nz-1.
    # ------------------------------------------------------------------
    def _setup_elliptic(self):
        g = self.g
        a = 1.0 / self.rho0
        # a on half levels (geometric mean of rho0 is exact for exponential rho0)
        a_half = 1.0 / np.sqrt(self.rho0[:-1] * self.rho0[1:])     # k+1/2, k = 0..Nz-1
        dz2 = g.dz**2
        n = g.Nz - 1                                               # interior rows
        lo = a_half[0:n] / dz2            # coefficient of psi[k-1], k = 1..Nz-1
        up = a_half[1:n + 1] / dz2        # coefficient of psi[k+1]
        diag = -(lo + up)[:, None] + a[1:-1, None] * g.lam[None, :]   # (n, M)
        self._ell_lo = lo
        self._ell_up = up
        # Thomas forward elimination of the (m-independent) off-diagonals
        dprime = np.empty_like(diag)
        mult = np.zeros((n, g.lam.size))
        dprime[0] = diag[0]
        for k in range(1, n):
            mult[k] = lo[k] / dprime[k - 1]
            dprime[k] = diag[k] - mult[k] * up[k - 1]
        self._ell_dprime = dprime
        self._ell_mult = mult
        self._ell_n = n

    def _solve_psi(self, eta: np.ndarray, psi_bot: np.ndarray) -> np.ndarray:
        g = self.g
        n = self._ell_n
        rhs = np.fft.rfft(-eta[1:-1], axis=1)                       # (n, M)
        rhs[0] -= self._ell_lo[0] * np.fft.rfft(psi_bot)            # bottom BC
        # (top BC psi = 0 contributes nothing)
        # forward sweep
        y = np.empty_like(rhs)
        y[0] = rhs[0]
        mult = self._ell_mult
        for k in range(1, n):
            y[k] = rhs[k] - mult[k] * y[k - 1]
        # back substitution
        dprime, up = self._ell_dprime, self._ell_up
        sol = np.empty_like(rhs)
        sol[n - 1] = y[n - 1] / dprime[n - 1]
        for k in range(n - 2, -1, -1):
            sol[k] = (y[k] - up[k] * sol[k + 1]) / dprime[k]
        psi = np.empty_like(eta)
        psi[0] = psi_bot
        psi[-1] = 0.0
        psi[1:-1] = np.fft.irfft(sol, n=g.Nx, axis=1)
        return psi

    # ------------------------------------------------------------------
    # Bottom boundary: psi(x,0,t) = rho0(0) * int f dx  (spectral integration)
    # ------------------------------------------------------------------
    def _psi_bottom(self, t: float) -> np.ndarray:
        g = self.g
        if self.w_bottom is None:
            return np.zeros(g.Nx)
        f = np.asarray(self.w_bottom(g.x, t), dtype=float) * np.ones(g.Nx)
        fbar = f.mean()
        if abs(fbar) > 1e-6 * (np.abs(f).max() + 1e-300) and not self._warned_mean:
            warnings.warn(
                "w_bottom has a nonzero x-mean, which is incompatible with "
                "periodic x and w=0 at the top; the mean is being removed."
            )
            self._warned_mean = True
        fh = np.fft.rfft(f - fbar)
        ph = np.zeros_like(fh)
        ph[1:] = fh[1:] / (1j * g.kx[1:])
        # the realised bottom w is the centered difference of psi_bot (2nd-order
        # accurate approximation of f); drop the Nyquist mode for safety.
        if g.Nx % 2 == 0:
            ph[-1] = 0.0
        return self.rho0[0] * np.fft.irfft(ph, n=g.Nx)

    # ------------------------------------------------------------------
    # Finite-difference operators (colocated grid, shape (Nz+1, Nx))
    # ------------------------------------------------------------------
    def _ddx(self, q: np.ndarray) -> np.ndarray:
        return (np.roll(q, -1, axis=1) - np.roll(q, 1, axis=1)) / (2.0 * self.g.dx)

    def _ddz_full(self, q: np.ndarray) -> np.ndarray:
        """Centered inside, 2nd-order one-sided at k=0 and k=Nz."""
        dz = self.g.dz
        out = np.empty_like(q)
        out[1:-1] = (q[2:] - q[:-2]) / (2.0 * dz)
        out[0] = (-3.0 * q[0] + 4.0 * q[1] - q[2]) / (2.0 * dz)
        out[-1] = (3.0 * q[-1] - 4.0 * q[-2] + q[-3]) / (2.0 * dz)
        return out

    def _lap_interior(self, q: np.ndarray) -> np.ndarray:
        g = self.g
        qxx = (np.roll(q, -1, axis=1) - 2.0 * q + np.roll(q, 1, axis=1)) / g.dx**2
        qzz = (q[2:] - 2.0 * q[1:-1] + q[:-2]) / g.dz**2
        return qxx[1:-1] + qzz

    def _advect(self, q: np.ndarray, Fx: np.ndarray, Fz: np.ndarray) -> np.ndarray:
        """
        (1/rho0) div(rho0 u q) at interior rows, skew-symmetric form:
            1/2 [ d(Fx q)/dx + Fx dq/dx + d(Fz q)/dz + Fz dq/dz ] / rho0
        Fx, Fz are the mass fluxes rho0*u, rho0*w on the full grid.
        """
        dz = self.g.dz
        xpart = 0.5 * (self._ddx(Fx * q) + Fx * self._ddx(q))[1:-1]
        Fzq = Fz * q
        zpart = 0.5 * ((Fzq[2:] - Fzq[:-2]) + Fz[1:-1] * (q[2:] - q[:-2])) / (2.0 * dz)
        return (xpart + zpart) / self.rho0[1:-1, None]

    # ------------------------------------------------------------------
    # Tendencies
    # ------------------------------------------------------------------
    def _velocities(self, eta: np.ndarray, t: float):
        """Return psi, total mass flux (Fx, Fz) and w on the full grid."""
        psi = self._solve_psi(eta, self._psi_bottom(t))
        Fz = self._ddx(psi)                                   # rho0 w
        Fx = -self._ddz_full(psi) + self.rho0[:, None] * self.U[:, None]   # rho0 (U+u')
        w = Fz / self.rho0[:, None]
        return psi, Fx, Fz, w

    def _tendencies(self, eta: np.ndarray, b: np.ndarray, t: float):
        psi, Fx, Fz, w = self._velocities(eta, t)
        wi = w[1:-1]
        sig = self.sigma[1:-1]

        deta = -self._advect(eta, Fx, Fz) - self._ddx(b)[1:-1] - sig * eta[1:-1]
        if self.anelastic:
            deta += eta[1:-1] * wi * self.dlnrho[1:-1, None]
        if np.any(self.U != 0.0):
            deta -= wi * self._Ucoef
        if self.nu > 0.0:
            deta += self.nu * self._lap_interior(eta)

        db = -self._advect(b, Fx, Fz) - wi * self.N2[1:-1, None] - sig * b[1:-1]
        if self.kappa > 0.0:
            db += self.kappa * self._lap_interior(b)
        return deta, db

    # ------------------------------------------------------------------
    # Time stepping (SSPRK3)
    # ------------------------------------------------------------------
    def step(self, dt: float):
        """Advance the model by one time step of length dt."""
        eta0, b0, t = self.eta, self.b, self.t

        d1, e1 = self._tendencies(eta0, b0, t)
        eta1 = eta0.copy(); b1 = b0.copy()
        eta1[1:-1] += dt * d1; b1[1:-1] += dt * e1

        d2, e2 = self._tendencies(eta1, b1, t + dt)
        eta2 = eta1.copy(); b2 = b1.copy()
        eta2[1:-1] += dt * d2; b2[1:-1] += dt * e2
        eta2 = 0.75 * eta0 + 0.25 * eta2
        b2 = 0.75 * b0 + 0.25 * b2

        d3, e3 = self._tendencies(eta2, b2, t + 0.5 * dt)
        eta3 = eta2.copy(); b3 = b2.copy()
        eta3[1:-1] += dt * d3; b3[1:-1] += dt * e3
        self.eta = (eta0 + 2.0 * eta3) / 3.0
        self.b = (b0 + 2.0 * b3) / 3.0

        self.t += dt
        self.nstep += 1
        self._diagnose()

    def _diagnose(self):
        """Refresh psi, u', w from the current state."""
        psi, Fx, Fz, w = self._velocities(self.eta, self.t)
        self.psi = psi
        self.w = w
        self.u = Fx / self.rho0[:, None] - self.U[:, None]     # perturbation only

    def run(self, t_end: float, dt: float, output_every: Optional[float] = None,
            callback: Optional[Callable[["Boussinesq2D"], None]] = None,
            fields=("u", "w", "b", "eta", "psi"), verbose: bool = False) -> dict:
        """
        Integrate until t >= t_end.  Snapshots of `fields` are stored every
        `output_every` seconds (None -> only the final state).  Returns a dict
        with 't' (1-D) and one (nt, Nz+1, Nx) array per field, plus 'x', 'z'.
        `callback(model)` is called after each stored snapshot.
        """
        nsteps = int(np.ceil((t_end - self.t) / dt - 1e-9))
        store_every = (max(1, int(round(output_every / dt)))
                       if output_every is not None else nsteps)
        out = {f: [] for f in fields}
        out["t"] = []

        def _store():
            out["t"].append(self.t)
            for f in fields:
                out[f].append(getattr(self, f).copy())
            if callback is not None:
                callback(self)

        _store()
        for n in range(1, nsteps + 1):
            self.step(dt)
            if not np.isfinite(self.eta).all():
                raise FloatingPointError(f"model blew up at step {self.nstep}, t={self.t:.1f}")
            if n % store_every == 0 or n == nsteps:
                _store()
                if verbose:
                    print(f"t = {self.t/3600:7.3f} h   max|w| = {np.abs(self.w).max():.3e}")
        for k in out:
            out[k] = np.asarray(out[k])
        out["x"] = self.g.x
        out["z"] = self.g.z
        return out

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------
    def suggest_dt(self, cfl: float = 0.5, umax_guess: float = 0.0) -> float:
        """
        A conservative maximum time step from (i) the buoyancy oscillation
        (SSPRK3 is stable for N dt < sqrt(3)), (ii) the advective CFL with the
        background wind plus `umax_guess`, (iii) explicit diffusion, and
        (iv) the sponge (sigma dt < 2).  Multiply by a safety factor yourself.
        """
        g = self.g
        lims = []
        Nmax = np.sqrt(max(self.N2.max(), 0.0))
        if Nmax > 0:
            lims.append(np.sqrt(3.0) / Nmax)
        umax = np.abs(self.U).max() + umax_guess
        if umax > 0:
            lims.append(cfl * g.dx / umax)
        for coef in (self.nu, self.kappa):
            if coef > 0:
                lims.append(0.5 * min(g.dx, g.dz) ** 2 / coef / 2.0)
        if self.sigma.max() > 0:
            lims.append(2.0 / self.sigma.max())
        return min(lims) if lims else np.inf

    def momentum_flux(self) -> np.ndarray:
        """Horizontally averaged vertical flux of horizontal momentum rho0 <u'w'>(z)."""
        return self.rho0 * np.mean(self.u * self.w, axis=1)

    def energy(self) -> tuple:
        """Domain-integrated perturbation kinetic and available potential energy."""
        g = self.g
        wgt = self.rho0[:, None] * np.ones_like(self.u)
        wgt[0] *= 0.5; wgt[-1] *= 0.5
        ke = 0.5 * np.sum(wgt * (self.u**2 + self.w**2)) * g.dx * g.dz
        N2 = np.where(self.N2 > 0, self.N2, np.inf)[:, None]
        pe = 0.5 * np.sum(wgt * self.b**2 / N2) * g.dx * g.dz
        return ke, pe


def to_xarray(out: dict):
    """Convert the dict returned by Boussinesq2D.run into an xarray.Dataset."""
    import xarray as xr   # optional dependency
    data = {}
    for k, v in out.items():
        if k in ("t", "x", "z"):
            continue
        data[k] = (("t", "z", "x"), v)
    return xr.Dataset(data, coords={"t": out["t"], "z": out["z"], "x": out["x"]})
