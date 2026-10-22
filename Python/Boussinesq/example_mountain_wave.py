"""
Example 1: hydrostatic mountain wave over a Witch-of-Agnesi ridge, compared
with linear theory.

The mountain enters only through the linearised lower boundary condition
    w(x, 0, t) = U * dh/dx ,     h(x) = h0 a^2 / (a^2 + x^2),
which is exactly the kind of bottom forcing the model is built around.

Parameters: U = 10 m/s, N = 0.01 1/s, h0 = 100 m, a = 10 km
  -> N a / U = 10   (hydrostatic),   N h0 / U = 0.1   (linear)
  -> vertical wavelength 2 pi U / N = 6.28 km

Two reference solutions are computed:
  (a) the classical infinite-domain hydrostatic Queney solution, and
  (b) the linear steady solution for the *same periodic forcing* the model
      actually sees (the 240 km box cannot hold the longest Agnesi waves,
      which carry ~10 % of the infinite-domain momentum flux).
(b) is the fair test of the numerics; (a) shows the effect of the finite box.
"""
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from boussinesq2d import Boussinesq2D

U0, N, h0, a = 10.0, 0.01, 100.0, 10.0e3
L, H = 120.0e3, 25.0e3
S, D = 90.0e3, 12.0e3                 # sponge-free box |x|<S, z<D
Nx, Nz = 480, 200                     # dx = 500 m, dz = 125 m
tau_ramp = 3600.0                     # ramp the forcing in over 1 h
t_end = 12 * 3600.0


def w_bottom(x, t):
    ramp = min(t / tau_ramp, 1.0)
    dhdx = -2.0 * h0 * a**2 * x / (a**2 + x**2) ** 2
    return ramp * U0 * dhdx


m = Boussinesq2D(Nx=Nx, Nz=Nz, L=L, H=H, N2=N**2, U=U0, w_bottom=w_bottom,
                 sponge=dict(S=S, D=D, sigma_max=0.2 * N))
dt = 20.0
print(f"dt = {dt} s  (suggest_dt = {m.suggest_dt():.1f} s)")

t0 = time.time()
out = m.run(t_end=t_end, dt=dt, output_every=3600.0, verbose=True)
print(f"wall time: {time.time() - t0:.1f} s for {m.nstep} steps")

# ---- linear steady solutions ----------------------------------------------
g = m.g
mz = N / U0
# (a) classical infinite-domain hydrostatic Queney solution
w_queney = U0 * h0 * a * ((g.X**2 - a**2) * np.sin(mz * g.Z)
                          - 2.0 * a * g.X * np.cos(mz * g.Z)) / (a**2 + g.X**2) ** 2
# (b) linear, non-hydrostatic, steady solution for the same periodic forcing
#     w_hat(k,z) = w_hat(k,0) exp(i m z),  m = sqrt(N^2/U^2 - k^2)  (k < N/U)
#                                          m = i sqrt(k^2 - N^2/U^2) (evanescent)
wh0 = np.fft.rfft(m.w[0])                      # realised bottom w (ramp finished)
k = g.kx
mm = np.where(k < N / U0, np.sqrt(np.maximum(N**2 / U0**2 - k**2, 0.0)) + 0j,
              1j * np.sqrt(np.maximum(k**2 - N**2 / U0**2, 0.0)))
wh = wh0[None, :] * np.exp(1j * mm[None, :] * g.z[:, None])
w_lin = np.fft.irfft(wh, n=g.Nx, axis=1)
with np.errstate(divide="ignore", invalid="ignore"):
    uh = np.where(k > 0, -mm[None, :] / k[None, :] * wh, 0.0)   # continuity
u_lin = np.fft.irfft(uh, n=g.Nx, axis=1)
mf_lin = np.mean(u_lin * w_lin, axis=1)

# compare inside the sponge-free box only
box = (np.abs(g.X) <= S) & (g.Z <= D)
rms = lambda f: np.sqrt(np.mean(f[box] ** 2))
print(f"RMS error of w in the clean box: {100*rms(m.w - w_lin)/rms(w_lin):.2f} % "
      f"(vs periodic linear),  {100*rms(m.w - w_queney)/rms(w_queney):.2f} % (vs Queney)")

mf = m.momentum_flux()
mf_agnesi = -(np.pi / 4.0) * N * U0 * h0**2 / (2.0 * L)   # infinite-domain drag / 2L
kz = int(3e3 / g.dz)
print(f"momentum flux at z=3 km:  model {mf[kz]:.3e}   periodic linear {mf_lin[kz]:.3e}"
      f"   infinite-domain Agnesi {mf_agnesi:.3e}")

# ---- figure ---------------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
lev = np.linspace(-0.1, 0.1, 21)
xk, zk = g.x / 1e3, g.z / 1e3
for ax, fld, ttl in zip(axs[:2], (m.w, w_lin),
                        (f"model, t = {m.t/3600:.0f} h", "linear theory (periodic)")):
    cf = ax.contourf(xk, zk, fld, lev, cmap="RdBu_r", extend="both")
    ax.axhline(D / 1e3, color="k", ls="--", lw=0.8)
    for s in (-S, S):
        ax.axvline(s / 1e3, color="k", ls="--", lw=0.8)
    ax.set_xlabel("x (km)"); ax.set_ylabel("z (km)"); ax.set_title(f"w (m/s): {ttl}")
fig.colorbar(cf, ax=axs[:2], shrink=0.9)
axs[2].plot(mf, zk, label="model")
axs[2].plot(mf_lin, zk, "k--", label="linear, periodic box")
axs[2].axvline(mf_agnesi, color="gray", ls=":", label="linear, infinite domain")
axs[2].axhline(D / 1e3, color="gray", ls=":")
axs[2].set_xlabel(r"$\rho_0\,\overline{u'w'}$  (m$^2$ s$^{-2}$)"); axs[2].set_ylabel("z (km)")
axs[2].set_title("momentum flux"); axs[2].legend(fontsize=8)
fig.savefig("mountain_wave.png", dpi=130)
print("wrote mountain_wave.png")
