import numpy as np
from boussinesq2d import Boussinesq2D

# --- 1. Elliptic solver check (anelastic): manufactured solution ------------
L, H = 50e3, 20e3
Hr = 7e3
m = Boussinesq2D(Nx=128, Nz=100, L=L, H=H, rho0=lambda z: np.exp(-z / Hr))
g = m.g
k = 2 * np.pi / L
mz = np.pi / H
psi_true = np.sin(k * g.X) * np.sin(2 * mz * g.Z) * 1e4          # psi=0 at z=0,H
rho = m.rho0[:, None]
# eta = -[ d/dz( (1/rho) psi_z ) + (1/rho) psi_xx ]
psi_z = np.sin(k * g.X) * 2 * mz * np.cos(2 * mz * g.Z) * 1e4
psi_zz = -np.sin(k * g.X) * (2 * mz) ** 2 * np.sin(2 * mz * g.Z) * 1e4
psi_xx = -k**2 * psi_true
eta = -((psi_zz + psi_z / Hr) / rho + psi_xx / rho)
psi_num = m._solve_psi(eta, np.zeros(g.Nx))
err = np.abs(psi_num - psi_true).max() / np.abs(psi_true).max()
print(f"elliptic solver rel. error (Nz=100): {err:.3e}")
m2 = Boussinesq2D(Nx=256, Nz=200, L=L, H=H, rho0=lambda z: np.exp(-z / Hr))
g2 = m2.g
psi_true2 = np.sin(k * g2.X) * np.sin(2 * mz * g2.Z) * 1e4
rho2 = m2.rho0[:, None]
psi_z2 = np.sin(k * g2.X) * 2 * mz * np.cos(2 * mz * g2.Z) * 1e4
psi_zz2 = -np.sin(k * g2.X) * (2 * mz) ** 2 * np.sin(2 * mz * g2.Z) * 1e4
eta2 = -((psi_zz2 + psi_z2 / Hr) / rho2 - k**2 * psi_true2 / rho2)
err2 = np.abs(m2._solve_psi(eta2, np.zeros(g2.Nx)) - psi_true2).max() / np.abs(psi_true2).max()
print(f"elliptic solver rel. error (Nz=200): {err2:.3e}   ratio = {err/err2:.2f} (expect ~4)")

# --- 2. Bottom BC check: realised w(x,0) vs requested f ---------------------
f = lambda x, t: 0.3 * np.exp(-(x / 5e3) ** 2) * np.cos(2 * np.pi * x / 20e3) - 0.0
m3 = Boussinesq2D(Nx=256, Nz=50, L=L, H=H, w_bottom=f)
fx = f(g2.x, 0.0); fx -= fx.mean()
print(f"bottom w error (rel.): {np.abs(m3.w[0] - fx).max() / np.abs(fx).max():.3e}")

# --- 3. Conservation with no forcing/diffusion/sponge (Boussinesq) ----------
m4 = Boussinesq2D(Nx=96, Nz=64, L=10e3, H=10e3, N2=1e-4)
g4 = m4.g
rng = np.random.default_rng(0)
blob = np.exp(-((g4.X) ** 2 + (g4.Z - 5e3) ** 2) / (1.5e3) ** 2)
m4.b[1:-1] = 0.05 * blob[1:-1]
m4.eta[1:-1] = 1e-3 * np.sin(2 * np.pi * g4.X / 10e3)[1:-1] * blob[1:-1]
m4._diagnose()
dt = 0.3 * m4.suggest_dt(umax_guess=2.0)
E0 = sum(m4.energy())
bsum0 = m4.b.sum()
for n in range(400):
    m4.step(dt)
E1 = sum(m4.energy())
print(f"energy drift over 400 steps: {(E1 - E0) / E0:+.3e}   sum(b) drift: {(m4.b.sum()-bsum0):+.3e}")
print(f"max|w| = {np.abs(m4.w).max():.3f} m/s, dt = {dt:.1f} s")
