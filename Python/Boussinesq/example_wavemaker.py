"""
Example 2: an oscillating, localised bottom "piston" in an anelastic
atmosphere with no mean wind.

    w(x, 0, t) = w0 exp(-x^2/sx^2) sin(omega t),   omega = N / sqrt(2)

A monochromatic source radiates energy along beams at an angle theta from the
vertical given by cos(theta) = omega/N, i.e. 45 degrees here (the familiar
"St Andrew's cross").  Because rho0 decreases with height the wave amplitude
grows like rho0^(-1/2) until the beams are absorbed in the sponge.

This also shows how to drive the model step by step from Python and collect
your own diagnostics along the way.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from boussinesq2d import Boussinesq2D

N = 0.01
omega = N / np.sqrt(2.0)
w0, sx = 0.05, 2.0e3
Hrho = 7.0e3
L, H = 40.0e3, 30.0e3
S, D = 30.0e3, 20.0e3


def w_bottom(x, t):
    return w0 * np.exp(-(x / sx) ** 2) * np.sin(omega * t)


m = Boussinesq2D(Nx=400, Nz=200, L=L, H=H, N2=N**2,
                 rho0=lambda z: np.exp(-z / Hrho),
                 w_bottom=w_bottom,
                 sponge=dict(S=S, D=D, sigma_max=0.2 * N))

dt = 10.0
nsteps = int(3 * 3600 / dt)
snap_times = [1800.0, 3600.0, 7200.0, 10800.0]
snaps = {}
ke_hist = []

# --- hand-driven time loop ---------------------------------------------------
for n in range(nsteps):
    m.step(dt)
    ke_hist.append((m.t, m.energy()[0]))
    for ts in snap_times:
        if abs(m.t - ts) < 0.5 * dt:
            snaps[ts] = m.w.copy()

# --- figure: w scaled by sqrt(rho0) so the beams look uniform ----------------
g = m.g
scale = np.sqrt(m.rho0)[:, None]
fig, axs = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
lev = np.linspace(-0.02, 0.02, 21)
for ax, ts in zip(axs.flat, snap_times):
    cf = ax.contourf(g.x / 1e3, g.z / 1e3, snaps[ts] * scale, lev, cmap="RdBu_r", extend="both")
    ax.axhline(D / 1e3, color="k", ls="--", lw=0.8)
    for s in (-S, S):
        ax.axvline(s / 1e3, color="k", ls="--", lw=0.8)
    # expected beam direction: 45 degrees
    for sgn in (-1, 1):
        ax.plot([0, sgn * 25], [0, 25], color="gray", ls=":", lw=1)
    ax.set_title(f"t = {ts/3600:.1f} h")
    ax.set_xlabel("x (km)"); ax.set_ylabel("z (km)")
fig.colorbar(cf, ax=axs, shrink=0.6, label=r"$w\,\sqrt{\rho_0/\rho_0(0)}$  (m/s)")
fig.suptitle("Oscillating piston, omega = N/sqrt(2): beams at 45 deg, absorbed in the sponge")
fig.savefig("wavemaker.png", dpi=130)
print("wrote wavemaker.png")

t_ke = np.array(ke_hist)
print(f"KE at end: {t_ke[-1, 1]:.3e}  (max over run {t_ke[:, 1].max():.3e}) -> "
      "levels off once the sponge absorbs as much as the piston radiates")
