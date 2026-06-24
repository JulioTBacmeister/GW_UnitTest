"""
gw_movmtn_calc_mod.py

Python/NumPy translation of gw_movmtn_calc_mod.F90
Driver module for the moving mountain gravity wave parameterisation.

Indexing convention
-------------------
All level indices are 0-based (Fortran k=1 → Python k=0).
  - Midpoint arrays : shape (..., pver)     [levels 0 .. pver-1]
  - Interface arrays: shape (..., pver+1)   [levels 0 .. pver  ]
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional

from gw_common import (GWBand, Coords1D, pver,
                        gw_drag_prof, calc_taucd, momentum_flux, momentum_fixer,
                        energy_change,
                        west, east, south, north)
from gw_movmtn import MovMtnSourceDesc, gw_movmtn_src

# ---------------------------------------------------------------------------
# Module-level state  (analogous to Fortran module variables)
# ---------------------------------------------------------------------------
band_movmtn: Optional[GWBand] = None
vramp:       Optional[np.ndarray] = None

# ---------------------------------------------------------------------------
# Namelist parameters (set by the caller or read from a namelist file).
# These match the variables in namelist_mod.F90.
# ---------------------------------------------------------------------------
alpha_gw_movmtn:  float = 1.0
movmtn_plaunch:   float = 95000.0   # launch pressure [Pa]
movmtn_psteer:    float = 85000.0   # steering pressure [Pa]
movmtn_source:    int   = 1         # 1 = vorticity, 2 = PBL mom flux
gw_apply_tndmax:  bool  = True


# ===========================================================================
# PhysicsPtend  (lightweight Python equivalent of physics_ptend in F90)
# ===========================================================================

@dataclass
class PhysicsPtend:
    """
    Parameterisation tendency accumulator.

    Attributes
    ----------
    u : ndarray, shape (pcols, pver)   zonal wind tendency [m s⁻²]
    v : ndarray, shape (pcols, pver)   meridional wind tendency [m s⁻²]
    s : ndarray, shape (pcols, pver)   heating tendency [J kg⁻¹ s⁻¹]
    q : ndarray, shape (pcols, pver, pcnst)  constituent tendency [kg kg⁻¹ s⁻¹]
    """
    u: np.ndarray
    v: np.ndarray
    s: np.ndarray
    q: np.ndarray

    @classmethod
    def zeros(cls, pcols: int, pver_: int, pcnst: int) -> "PhysicsPtend":
        return cls(
            u=np.zeros((pcols, pver_)),
            v=np.zeros((pcols, pver_)),
            s=np.zeros((pcols, pver_)),
            q=np.zeros((pcols, pver_, pcnst)),
        )


# ===========================================================================
# set_band_movmtn
# ===========================================================================

def set_band_movmtn(band_movmtn_in: GWBand) -> None:
    """Set the module-level GWBand for moving mountains."""
    global band_movmtn
    band_movmtn = band_movmtn_in
    print(f" Band MovMTN {band_movmtn.ngwv}")


# ===========================================================================
# set_vramp
# ===========================================================================

def set_vramp() -> None:
    """Allocate vramp and set to 1 for all levels."""
    global vramp
    vramp = np.ones(pver)


# ===========================================================================
# report_from_within
# ===========================================================================

def report_from_within() -> None:
    """Diagnostic: print module-level state."""
    if band_movmtn is not None:
        print(f" Inside gw_movmtn_calc_mod: Band MovMtn {band_movmtn.ngwv}")
    else:
        print(" Inside gw_movmtn_calc_mod: band_movmtn not set")


# ===========================================================================
# gw_init_movmtn
# ===========================================================================

def gw_init_movmtn(file_name: str, band: GWBand) -> MovMtnSourceDesc:
    """
    Initialise a MovMtnSourceDesc with dummy (stub) lookup table values.

    This is a Python port of the abridged ``gw_init_movmtn`` in
    gw_movmtn_calc_mod.F90, which fills the table with sentinel values
    (-999) rather than reading from a file.

    Parameters
    ----------
    file_name : str   (unused; kept for signature compatibility)
    band      : GWBand

    Returns
    -------
    desc : MovMtnSourceDesc
    """
    ngwv = band.ngwv
    desc = MovMtnSourceDesc()

    # Heating-depth dimension (Fortran: maxh = 15)
    desc.maxh = 15

    # Mean-wind dimension (Fortran: maxuh = (241-1)/2 = 120)
    desc.maxuh = (241 - 1) // 2   # = 120

    # Allocate and fill heating-depth axis [km → m conversion; stub = 1000 m]
    desc.hd = np.full(desc.maxh, 1000.0)   # Fortran: desc%hd = 1000._r8

    # Allocate and fill wind axis (positive half; stub = 0)
    desc.uh = np.zeros(desc.maxuh)

    # Allocate and fill momentum flux table with -999 (stub / unread)
    desc.mfcc = np.full(
        (desc.maxh, 2 * desc.maxuh + 1, 2 * ngwv + 1),
        -999.0)

    return desc


# ===========================================================================
# gw_movmtn_calc
# ===========================================================================

def gw_movmtn_calc(ncol, lchnk, dt, pref_edge,
                   u, v, t, vort4gw, p, piln, zm, zi,
                   nm, ni, rhoi, kvtt, q, dse,
                   effgw_movmtn_pbl, ptend):
    """
    Top-level driver for the moving mountain GW parameterisation.

    Parameters
    ----------
    ncol : int
        Number of active atmospheric columns.
    lchnk : int
        Chunk identifier.
    dt : float
        Time step [s].
    pref_edge : ndarray, shape (pver+1,)
        Reference pressure at interface levels [Pa].
    u, v : ndarray, shape (ncol, pver)
        Midpoint zonal / meridional winds [m s⁻¹].
    t : ndarray, shape (ncol, pver)
        Midpoint temperatures [K].
    vort4gw : ndarray, shape (ncol, pver)
        3-D relative vorticity [s⁻¹].
    p : Coords1D
        Pressure coordinate object.
    piln : ndarray, shape (ncol, pver+1)
        Log of interface pressures.
    zm : ndarray, shape (ncol, pver)
        Midpoint altitudes above ground [m].
    zi : ndarray, shape (ncol, pver+1)
        Interface altitudes above ground [m].
    nm : ndarray, shape (ncol, pver)
        Midpoint Brunt-Väisälä frequency [rad s⁻¹].
    ni : ndarray, shape (ncol, pver+1)
        Interface Brunt-Väisälä frequency [rad s⁻¹].
    rhoi : ndarray, shape (ncol, pver+1)
        Interface density [kg m⁻³].
    kvtt : ndarray, shape (ncol, pver+1)
        Molecular thermal diffusivity [m² s⁻¹].
    q : ndarray, shape (ncol, pver, pcnst)
        Constituent mixing ratios.
    dse : ndarray, shape (ncol, pver)
        Dry static energy [J kg⁻¹].
    effgw_movmtn_pbl : float
        Tendency efficiency.
    ptend : PhysicsPtend
        Accumulated tendencies — modified in place.

    Returns
    -------
    flx_heat : ndarray, shape (ncol,)
        Net energy flux [J m⁻²].
    """
    global band_movmtn, vramp

    itime = 1   # time index for diagnostic output (stub)

    # ----------------------------------------------------------------
    # Initialise moving mountain source descriptor
    # ----------------------------------------------------------------
    movmtn_desc = gw_init_movmtn("DummyName.nc", band_movmtn)

    # Locate the 950 hPa level (Fortran: do k = 0, pver; pref_edge(k+1) < 95000)
    for k in range(pver + 1):
        # Fortran: pref_edge(k+1) with k from 0 to pver
        #   pref_edge is 1-based size pver+1; pref_edge(k+1) → Python pref_edge[k]
        if pref_edge[k] < 95000.0:
            movmtn_desc.k = k   # 0-based interface level

    movmtn_desc.min_hdepth = 1.0

    # ----------------------------------------------------------------
    # Find steering and launch levels from pref_edge
    # Use -1 as "not found" sentinel (analogous to Fortran uninitialized = 0)
    # ----------------------------------------------------------------
    movmtn_ksteer  = -1
    movmtn_klaunch = -1

    for k in range(pver):
        # Fortran: pref_edge(k+1) and pref_edge(k) — 1-based, k=1..pver
        # Python : pref_edge[k]   and pref_edge[k-1] for the same levels,
        #          but we iterate k over 0..pver-1 so the interfaces are
        #          pref_edge[k] (top) and pref_edge[k+1] (bottom).
        if pref_edge[k + 1] >= movmtn_psteer and pref_edge[k] < movmtn_psteer:
            movmtn_ksteer = k
    for k in range(pver):
        if pref_edge[k + 1] >= movmtn_plaunch and pref_edge[k] < movmtn_plaunch:
            movmtn_klaunch = k

    print(" ooooohhh  ... in gw_movmtn_calc")

    # ----------------------------------------------------------------
    # Initialise accumulated momentum flux and diagnostics
    # ----------------------------------------------------------------
    taummx   = np.zeros((ncol, pver + 1))
    taummy   = np.zeros((ncol, pver + 1))
    tau_diag = np.full((ncol, pver + 1), -9999.0)

    # Zero-out dummy source arrays for the non-PBL path
    upwp_clubb_gw = np.zeros((ncol, pver))
    vpwp_clubb_gw = np.zeros((ncol, pver))
    ttend_clubb   = np.zeros((ncol, pver))
    ttend_dp      = np.zeros((ncol, pver))
    xpwp_clubb    = np.sqrt(upwp_clubb_gw**2 + vpwp_clubb_gw**2)

    # ----------------------------------------------------------------
    # Compute GW source
    # ----------------------------------------------------------------
    src_level, tend_level, tau, ubm, ubi, xv, yv, phase_speeds, hdepth = \
        gw_movmtn_src(
            ncol, lchnk, band_movmtn, movmtn_desc,
            u, v, ttend_dp, ttend_clubb, xpwp_clubb, vort4gw,
            zm, alpha_gw_movmtn, movmtn_source,
            movmtn_ksteer, movmtn_klaunch)

    # Project stress into cardinal directional components
    taucd = calc_taucd(ncol, band_movmtn.ngwv, tend_level,
                       tau, phase_speeds, xv, yv, ubi)

    print(" Now back in gw_movmtn_calc  topi       ")
    print(f" range src_level    {src_level.min():d}  {src_level.max():d}")
    print(f" range tend_level   {tend_level.min():d}  {tend_level.max():d}")
    print(f" range phase_speeds {phase_speeds.min():.4e}  {phase_speeds.max():.4e}")
    print(f" range tau          {tau.min():.4e}  {tau.max():.4e}")
    print(f" band_movmtn%effkwv {band_movmtn.effkwv:.6e}")

    # Stub for ncfile_put_col3d diagnostic output
    # (In Fortran these write to a NetCDF file; here we skip them.)
    # ncfile_put_col3d('TAU_SRC_MOVMVTN', tau[:, band_movmtn.ngwv, :], ...)

    # ----------------------------------------------------------------
    # Propagate waves: compute drag tendencies
    # ----------------------------------------------------------------
    effgw = np.full(ncol, effgw_movmtn_pbl)
    print(f" range effgw       {effgw.min():.4e}  {effgw.max():.4e}")

    (utgw, vtgw, ttgw, qtgw, egwdffi,
     gwut, dttdf, dttke, tau_diag) = gw_drag_prof(
        ncol, band_movmtn, p, src_level, tend_level, dt,
        t, vramp, piln, rhoi, nm, ni, ubm, ubi, xv, yv,
        effgw, phase_speeds, kvtt, q, dse, tau,
        lapply_effgw_in=gw_apply_tndmax,
        tau_diag=tau_diag)

    # Project stress into cardinal directional components (post-drag)
    taucd = calc_taucd(ncol, band_movmtn.ngwv, tend_level,
                       tau, phase_speeds, xv, yv, ubi)

    # ----------------------------------------------------------------
    # Accumulate diffusivity (Fortran: egwdffi_tot was uninitialized —
    # initialised to zero here to avoid undefined behaviour)
    # ----------------------------------------------------------------
    egwdffi_tot = np.zeros((ncol, pver + 1))
    for k in range(pver + 1):
        egwdffi_tot[:, k] += egwdffi[:, k]

    # ----------------------------------------------------------------
    # Add tendencies to ptend
    # ----------------------------------------------------------------
    for k in range(pver):
        ptend.u[:ncol, k] += utgw[:, k]
        ptend.v[:ncol, k] += vtgw[:, k]
        ptend.s[:ncol, k] += ttgw[:, k]

    pcnst = q.shape[2] if q.ndim == 3 else 1
    for m in range(pcnst):
        for k in range(pver):
            ptend.q[:ncol, k, m] += qtgw[:, k, m]

    # ----------------------------------------------------------------
    # Momentum fixer
    # ----------------------------------------------------------------
    um_flux, vm_flux = momentum_flux(tend_level, taucd)
    momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

    # Stress in zonal / meridional directions for diagnostics
    ngwv = band_movmtn.ngwv
    for k in range(pver + 1):
        taummx[:, k] = tau[:, ngwv, k] * xv
        taummy[:, k] = tau[:, ngwv, k] * yv

    # Diagnostic stubs (ncfile_put_col3d omitted)

    # ----------------------------------------------------------------
    # Energy change for CAM energy checker
    # ----------------------------------------------------------------
    de = energy_change(dt, p, u, v,
                       ptend.u[:ncol, :], ptend.v[:ncol, :],
                       ptend.s[:ncol, :])
    flx_heat = np.zeros(ncol)
    flx_heat[:ncol] = de

    print(" wheeeewww   ... survived gw_movmtn_calc")

    return flx_heat
