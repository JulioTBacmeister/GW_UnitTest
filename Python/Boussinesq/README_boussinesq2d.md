# boussinesq2d — a small 2-D (x–z) Boussinesq / anelastic model in NumPy

Files

| file | what it is |
|---|---|
| `boussinesq2d.py` | the model (one class, no dependencies beyond NumPy) |
| `example_mountain_wave.py` | hydrostatic mountain wave vs linear theory (validation) |
| `example_wavemaker.py` | oscillating piston in an anelastic atmosphere, hand-driven time loop |
| `test_basic.py` | elliptic-solver convergence, bottom-BC accuracy, inviscid conservation |

## What it solves

Vorticity–streamfunction form.  With a reference density ρ₀(z) (constant → Boussinesq)
the velocity comes from a *mass* streamfunction ψ:

    ρ₀ u' = −∂ψ/∂z ,   ρ₀ w = +∂ψ/∂x        so that  ∇·(ρ₀ u) = 0 .

Total horizontal wind is U(z) + u'.  Prognostic variables are the perturbation
vorticity η = ∂u'/∂z − ∂w/∂x and the buoyancy b (= g θ'/θ₀):

    ∂η/∂t = −(1/ρ₀)∇·(ρ₀ u η) + η w (d ln ρ₀/dz) − w (U_zz − U_z d ln ρ₀/dz) − ∂b/∂x
            + ν ∇²η − σ(x,z) η
    ∂b/∂t = −(1/ρ₀)∇·(ρ₀ u b) − w N²(z) + κ ∇²b − σ(x,z) b
    ∂/∂z((1/ρ₀) ∂ψ/∂z) + (1/ρ₀) ∂²ψ/∂x² = −η

Set `rho0=None` (default) for Boussinesq; the extra terms switch off automatically.

## Boundary conditions

* x: periodic on [−L, L).
* z = 0: `w = w_bottom(x, t)` (any callable returning an array over x).  Implemented
  as ψ(x,0,t) = ρ₀(0)·∫f dx (spectral integral).  Because the box is periodic in x and
  closed at the top, the x-mean of f must be zero; a nonzero mean is removed with a warning.
* z = H: w = 0 (ψ = 0).
* η = b = 0 on both z-boundaries — inflowing fluid is unperturbed.  Normally the sponge
  covers the top anyway.

## Sponge

`sponge=dict(S=..., D=..., sigma_max=..., combine='max')`: the Rayleigh damping rate is
zero for |x| ≤ S and z ≤ D and increases **linearly** to `sigma_max` at |x| = L and at
z = H (`combine='max'` takes the larger of the two ramps; `'sum'` adds them).  η and b
are damped toward zero.

Tuning lesson from the validation run: the sponge must be *gentle*.  With a vertical ramp
of ~2 vertical wavelengths, `sigma_max = 0.5 N` reflected enough to put a 12 % error in the
interior w field; `sigma_max = 0.1–0.2 N` brought that to ~3 % (the remainder is
transients and the finite periodic box).  Longer ramps let you use stronger damping.

## Numerics (all second order)

* Colocated grid, fields are arrays of shape `(Nz+1, Nx)`, index `[k, i]`;
  `x_i = −L + i·dx`, `z_k = k·dz`, both boundaries are grid rows.
* Centered differences.  Advection is in skew-symmetric form (½ flux-form + ½ advective),
  which conserves the first and second moments of η and b exactly on this grid — there is
  **no implicit numerical damping**, so for long or nonlinear runs give it a small ν/κ
  (a few ×10⁻² dx²/dt is plenty) or rely on the sponge.
* Elliptic solve: real FFT in x, tridiagonal solve in z vectorised over wavenumbers; the
  x-eigenvalue is that of the discrete 3-point second difference, so the discrete
  divergence of the mass flux is exactly zero.  The factorisation is done once.
* Time stepping: SSP-RK3.  Stability limits: N·dt < √3, advective CFL,
  σ_max·dt < 2, explicit diffusion.  `suggest_dt()` returns the most restrictive
  of these — multiply by your own safety factor (I used ~0.8).

Cost: the mountain-wave case (480 × 201 points, 2160 steps) takes ~70 s on one core.

## Usage

```python
from boussinesq2d import Boussinesq2D, to_xarray

m = Boussinesq2D(Nx=480, Nz=200, L=120e3, H=25e3,
                 N2=1e-4,                      # scalar, array on z nodes, or f(z)
                 U=10.0,                       # background wind, same options
                 rho0=lambda z: np.exp(-z/7e3),  # omit for Boussinesq
                 w_bottom=lambda x, t: ...,    # bottom w, array over x
                 nu=0.0, kappa=0.0,
                 sponge=dict(S=90e3, D=12e3, sigma_max=2e-3))

# one-liner run, snapshots in a dict (or xarray via to_xarray(out))
out = m.run(t_end=12*3600, dt=20.0, output_every=3600, verbose=True)
out['w'].shape          # (nt, Nz+1, Nx);  also 'u', 'b', 'eta', 'psi', 't', 'x', 'z'

# ...or drive it yourself
for n in range(1000):
    m.step(dt)
    m.u, m.w, m.b, m.eta, m.psi      # refreshed after every step
    m.momentum_flux()                # rho0 <u'w'>(z)
    m.energy()                       # (KE, APE), domain integrated
```

`m.u` is the perturbation u' (U is not included).  You can set initial conditions by
writing to `m.eta[1:-1]` and `m.b[1:-1]` and calling `m._diagnose()`.

## Validation (example_mountain_wave.py)

Witch of Agnesi, U = 10 m/s, N = 0.01 s⁻¹, h₀ = 100 m, a = 10 km, forcing
w(x,0) = U·dh/dx.  After 12 h the modelled w agrees with the linear steady solution for the
same periodic forcing to 6.7 % RMS over the whole sponge-free box (3 % below 6 km, where
the slow long-wave components have arrived), and the momentum flux at 3 km is within 1 %
of that solution (−3.14e-3 vs −3.18e-3 m² s⁻²).  The infinite-domain Agnesi value is
−3.27e-3; the 3 % gap is the ~10 % of the Agnesi spectrum at wavelengths longer than the
240 km box — a property of the box, not of the numerics.

## Things you may want to change

* The bottom-vorticity condition (η = 0) is the "wavemaker" choice; for a no-slip wall you
  would diagnose η from ψ instead.
* The sponge ramp is linear as requested; a smooth (e.g. sin²) ramp reflects less.
* Laplacian ν/κ act on all scales; a ∇⁴ hyperdiffusion would be more scale-selective.
