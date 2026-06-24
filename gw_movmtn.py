"""
gw_movmtn.py

Python/NumPy translation of gw_movmtn.F90
Gravity wave source parameterization from the obstacle (moving mountain) effect
produced by internal atmospheric circulations.

Indexing convention
-------------------
All level indices are 0-based (Fortran k=1 → Python k=0).
  - Midpoint arrays : shape (..., pver)     [levels 0 .. pver-1]
  - Interface arrays: shape (..., pver+1)   [levels 0 .. pver  ]
Wavenumber index l ranges from -ngwv to +ngwv.
In arrays the wavenumber dimension has size 2*ngwv+1;
  l maps to array index  l + ngwv.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional

from gw_common import GWBand, pver, gravit, rair, midpoint_interp

# ---------------------------------------------------------------------------
# Module-level flag  (analogous to use_gw_movmtn_pbl in phys_control.F90)
# Default matches phys_control.F90: use_gw_movmtn_pbl = .true.
# ---------------------------------------------------------------------------
use_gw_movmtn_pbl: bool = True


# ===========================================================================
# MovMtnSourceDesc
# ===========================================================================

@dataclass
class MovMtnSourceDesc:
    """
    Settings and lookup table for the moving mountain GW source.

    Attributes
    ----------
    storm_shift : bool
        Whether wind speeds are shifted relative to storm cells.
    k : int
        Level index where wind speed is used as the source speed (0-based).
    min_hdepth : float
        Heating depths below this value [m] will be ignored.
    maxh : int
        Size of the lookup table heating-depth dimension.
    maxuh : int
        Half-size of the lookup table wind dimension
        (full wind axis = -maxuh .. +maxuh, size 2*maxuh+1).
    hd : ndarray, shape (maxh,)
        Heating depths in the lookup table [m].
    uh : ndarray, shape (maxuh,)
        Wind speeds in the lookup table (positive half only).
    mfcc : ndarray, shape (maxh, 2*maxuh+1, 2*ngwv+1)
        Lookup table of source spectra f(depth, wind, phase_speed).
    """
    storm_shift: bool  = False
    k:           int   = 0
    min_hdepth:  float = 0.0
    maxh:        int   = 0
    maxuh:       int   = 0
    hd:   Optional[np.ndarray] = None
    uh:   Optional[np.ndarray] = None
    mfcc: Optional[np.ndarray] = None


# ===========================================================================
# index_of_nearest
# ===========================================================================

def index_of_nearest(x: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """
    Return 0-based indices of the grid points nearest to each value in x.

    Mirrors ``index_of_nearest`` in gw_movmtn.F90
    (which returns 1-based indices in Fortran; here we return 0-based).

    Parameters
    ----------
    x    : ndarray, shape (n,)
    grid : ndarray, shape (m,)

    Returns
    -------
    idx : ndarray of int, shape (n,)
    """
    n = len(grid)
    interfaces = 0.5 * (grid[:-1] + grid[1:])   # midpoints between grid nodes
    idx = np.zeros(len(x), dtype=int)
    for i in range(n - 1):
        idx = np.where(x > interfaces[i], i + 1, idx)
    return idx


# ===========================================================================
# vorticity_flux_src
# ===========================================================================

def vorticity_flux_src(vorticity: np.ndarray, ncol: int,
                       alpha_gw_movmtn: float):
    """
    Calculate GW flux source from vorticity.

    Fortran pverx argument = pver (midpoint dimension);
    level indices here are already 0-based.

    Parameters
    ----------
    vorticity : ndarray, shape (ncol, pver)
    ncol : int
    alpha_gw_movmtn : float   tuning parameter

    Returns
    -------
    vort_src       : ndarray, shape (ncol,)
    steering_level : ndarray of int, shape (ncol,)  0-based
    launch_level   : ndarray of int, shape (ncol,)  0-based
    """
    scale_factor = 1.e4   # scales vorticity amplitude to u'w' in CLUBB
    nlayers = 10

    # Fortran: steering_level = pverx - 20  (pverx=pver, 1-based)
    # Python  (0-based): pver - 1 - 20 = pver - 21
    steering_level = np.full(ncol, pver - 21, dtype=int)
    launch_level   = steering_level - 10        # pver - 31

    vort_src = np.zeros(ncol)
    for k in range(nlayers):
        # Fortran: vort_src += scale_factor * abs(vorticity(:, pverx-k))
        # pverx-k (1-based) → pver-1-k (0-based)
        vort_src[:ncol] += scale_factor * np.abs(vorticity[:ncol, pver - 1 - k])
    vort_src[:ncol] = alpha_gw_movmtn * vort_src[:ncol] / nlayers

    return vort_src, steering_level, launch_level


# ===========================================================================
# shcu_flux_src
# ===========================================================================

def shcu_flux_src(xpwp_shcu: np.ndarray, ncol: int,
                  alpha_gw_movmtn: float):
    """
    Calculate GW flux source from ShCu/PBL and set steering/launch levels.

    Fortran pverx argument = pver+1 (interface dimension);
    level indices here are 0-based.

    Parameters
    ----------
    xpwp_shcu      : ndarray, shape (ncol, pver+1)
    ncol           : int
    alpha_gw_movmtn: float

    Returns
    -------
    xpwp_src       : ndarray, shape (ncol,)
    steering_level : ndarray of int, shape (ncol,)  0-based
    launch_level   : ndarray of int, shape (ncol,)  0-based
    """
    # Fortran (pverx=pver+1, 1-based): steering_level = (pverx-1) - 5 = pver - 5
    # Python  (0-based): pver - 5 - 1 = pver - 6
    steering_level = np.full(ncol, pver - 6, dtype=int)
    launch_level   = steering_level - 10        # pver - 16

    nlayers = 5
    xpwp_src = np.zeros(ncol)
    for k in range(nlayers):
        # Fortran: xpwp_shcu(:, pverx-k) where pverx=pver+1 (1-based, interface)
        # 0-based interface index: pver - k
        xpwp_src[:ncol] += xpwp_shcu[:ncol, pver - k]
    xpwp_src[:ncol] = alpha_gw_movmtn * xpwp_src[:ncol] / nlayers

    return xpwp_src, steering_level, launch_level


# ===========================================================================
# gw_movmtn_src
# ===========================================================================

def gw_movmtn_src(ncol, lchnk, band, desc,
                  u, v, netdt, netdt_shcu, xpwp_shcu, vorticity,
                  zm, alpha_gw_movmtn, movmtn_source,
                  ksteer_in, klaunch_in):
    """
    Flexible driver for GW source from obstacle effects produced by internal
    circulations.

    Parameters
    ----------
    ncol : int
    lchnk : int
    band : GWBand
    desc : MovMtnSourceDesc
    u, v : ndarray, shape (ncol, pver)
        Midpoint zonal/meridional winds [m s⁻¹].
    netdt : ndarray, shape (ncol, pver)
        Heating rate from deep convection [K s⁻¹].
    netdt_shcu : ndarray, shape (ncol, pver)
        Heating from ShCu/PBL [K s⁻¹].
    xpwp_shcu : ndarray, shape (ncol, pver+1)
        Higher-order momentum flux from ShCu/PBL [m² s⁻²].
    vorticity : ndarray, shape (ncol, pver)
        Relative vorticity [s⁻¹].
    zm : ndarray, shape (ncol, pver)
        Midpoint altitudes above ground [m].
    alpha_gw_movmtn : float
        Tuning parameter controlling proportion of PBL flux emitted as GW.
    movmtn_source : int
        Source type: 1 = vorticity, 2 = PBL momentum fluxes.
    ksteer_in : int
        Steering level override (0-based). Pass -1 to use internally computed level.
    klaunch_in : int
        Launch level override (0-based). Pass -1 to use internally computed level.

    Returns
    -------
    src_level  : ndarray of int, shape (ncol,)    0-based
    tend_level : ndarray of int, shape (ncol,)    0-based
    tau        : ndarray, shape (ncol, 2*ngwv+1, pver+1)
        Wave Reynolds stress [m² s⁻²].
    ubm        : ndarray, shape (ncol, pver)
        Wind projected onto wave direction at midpoints [m s⁻¹].
    ubi        : ndarray, shape (ncol, pver+1)
        Wind projected onto wave direction at interfaces [m s⁻¹].
    xv, yv     : ndarray, shape (ncol,)
        Zonal/meridional components of unit wave propagation vector.
    c          : ndarray, shape (ncol, 2*ngwv+1)
        Phase speeds [m s⁻¹].
    hdepth     : ndarray, shape (ncol,)
        Heating depth [m].
    """
    ngwv   = band.ngwv
    nwaves = 2 * ngwv + 1

    # ---- initialise output arrays ----------------------------------------
    tau    = np.zeros((ncol, nwaves, pver + 1))
    hdepth = np.zeros(ncol)
    c      = np.zeros((ncol, nwaves))
    ubm    = np.zeros((ncol, pver))
    ubi    = np.zeros((ncol, pver + 1))
    xv     = np.zeros(ncol)
    yv     = np.zeros(ncol)

    print(" Debugging gw_movmtn.py ")
    print(f" use_gw_movmtn_pbl {use_gw_movmtn_pbl}")
    print(f" movmtn_source     {movmtn_source}")
    print(f" min max vort      {vorticity.min():.4e}  {vorticity.max():.4e}")

    # ------------------------------------------------------------------
    # Source flux and initial steering / launch levels
    # ------------------------------------------------------------------
    source_type = movmtn_source
    if source_type == 1:
        xpwp_src, Steer_k, Launch_k = vorticity_flux_src(
            vorticity, ncol, alpha_gw_movmtn)
    elif source_type == 2:
        xpwp_src, Steer_k, Launch_k = shcu_flux_src(
            xpwp_shcu, ncol, alpha_gw_movmtn)
    else:
        xpwp_src = np.zeros(ncol)
        Steer_k  = np.zeros(ncol, dtype=int)
        Launch_k = np.zeros(ncol, dtype=int)

    # Override steering / launch levels when valid inputs are supplied.
    # (Fortran uses > 0 for 1-based; Python uses >= 0 for 0-based.)
    if klaunch_in >= 0:
        Launch_k[:ncol] = klaunch_in
    if ksteer_in >= 0:
        Steer_k[:ncol] = ksteer_in

    # ------------------------------------------------------------------
    # Wind and unit vector at the steering level
    # ------------------------------------------------------------------
    usteer = np.array([u[i, Steer_k[i]] for i in range(ncol)])
    vsteer = np.array([v[i, Steer_k[i]] for i in range(ncol)])
    # steer_level kept as float for potential diagnostics
    steer_level = Steer_k.astype(float)

    # Unit vector at steering level (inline get_unit_vector)
    mag_steer  = np.sqrt(usteer**2 + vsteer**2)
    xv_steer   = np.where(mag_steer > 0.0, usteer / mag_steer, 0.0)
    yv_steer   = np.where(mag_steer > 0.0, vsteer / mag_steer, 0.0)

    # Cell retrograde speed — always 0 in current implementation
    Cell_Retro_Speed = np.minimum(np.sqrt(usteer**2 + vsteer**2), 0.0)
    usteer -= xv_steer * Cell_Retro_Speed
    vsteer -= yv_steer * Cell_Retro_Speed
    # After the above, (usteer, vsteer) is the cell speed (= ground-based
    # wave phase speed for moving mountain GW).

    # ------------------------------------------------------------------
    # Heating depth calculation
    # ------------------------------------------------------------------
    # boti: bottom index of heating region (0-based), sentinel = -1 (unset)
    # topi: top    index of heating region (0-based), sentinel = -1 (unset)
    boti = np.full(ncol, -1, dtype=int)
    topi = np.full(ncol, -1, dtype=int)

    if use_gw_movmtn_pbl:
        # Fortran: boti=pver (1-based) → 0-based: pver-1
        boti[:] = pver - 1
        topi[:] = Launch_k[:ncol]
    else:
        # Scan from surface upward; k decreases as we go toward the top
        for k in range(pver - 1, -1, -1):
            for i in range(ncol):
                if boti[i] == -1:
                    if zm[i, k] >= 20000.0:
                        boti[i] = k
                        topi[i] = k
                    else:
                        if netdt[i, k] > 0.0:
                            boti[i] = k
                elif topi[i] == -1:
                    if zm[i, k] >= 20000.0:
                        topi[i] = k
                    else:
                        if not (netdt[i, k] > 0.0):
                            topi[i] = k
            if np.all(topi != -1):
                break
        # Safety: any column where topi was never found → set equal to boti
        topi = np.where(topi == -1, boti, topi)

    # Heating depth [m]
    idx    = np.arange(ncol)
    hdepth = zm[idx, topi] - zm[idx, boti]

    # Index into the lookup table heating-depth axis
    hd_idx = index_of_nearest(hdepth, desc.hd)
    # Use -1 as "too shallow / inactive" sentinel (Fortran uses 0 for 1-based)
    below_min = hdepth < max(desc.min_hdepth, desc.hd[0])
    hd_idx    = np.where(below_min, -1, hd_idx)

    # Maximum heating rate over the heating region
    q0 = np.zeros(ncol)
    if not use_gw_movmtn_pbl:
        k_lo = int(np.min(topi))
        k_hi = int(np.max(boti))
        for k in range(k_lo, k_hi + 1):
            mask = (k >= topi) & (k <= boti)
            q0   = np.where(mask, np.maximum(q0, netdt[:ncol, k]), q0)

    CF = 20.0      # heating rate conversion factor (1 / 5%)
    AL = 1.0e5     # averaging length [m]
    q0 = q0 * CF
    qj = gravit / rair * q0    # unit conversion to m s⁻³

    # Cell speed: CS1 = |usteer|, CS = CS1*(xv_steer + yv_steer)
    CS1 = np.sqrt(usteer**2 + vsteer**2)
    CS  = CS1 * xv_steer + CS1 * yv_steer

    # ------------------------------------------------------------------
    # Wave-relative wind profiles  (U - c_cell)
    # ------------------------------------------------------------------
    uwavef = u[:ncol, :] - usteer[:, np.newaxis]   # (ncol, pver)
    vwavef = v[:ncol, :] - vsteer[:, np.newaxis]

    # Wave-relative wind at the source (topi) level
    udiff = np.array([uwavef[i, topi[i]] for i in range(ncol)])
    vdiff = np.array([vwavef[i, topi[i]] for i in range(ncol)])

    # Unit vector in direction of wavevector (inline get_unit_vector)
    mag2   = np.sqrt(udiff**2 + vdiff**2)
    ubisrc = mag2.copy()
    xv = np.where(mag2 > 0.0, udiff / mag2, 0.0)
    yv = np.where(mag2 > 0.0, vdiff / mag2, 0.0)

    # Project wave-relative midpoint winds onto the wavevector
    ubm = uwavef * xv[:, np.newaxis] + vwavef * yv[:, np.newaxis]

    # Source-level on-crest wind
    ubmsrc = np.array([ubm[i, topi[i]] for i in range(ncol)])

    # Ensure source-level wind is always positive; flip sign of ubm, xv, yv
    sgn = np.sign(ubmsrc)
    sgn = np.where(sgn == 0.0, 1.0, sgn)
    ubm *= sgn[:, np.newaxis]
    xv  *= sgn
    yv  *= sgn

    # Interface wind by averaging adjacent midpoint winds
    # Fortran: ubi(:,1) = ubm(:,1); ubi(:,2:pver) = midpoint_interp(ubm)
    ubi[:, 0]      = ubm[:, 0]                   # top interface
    ubi[:, 1:pver] = midpoint_interp(ubm)        # interior interfaces
    # ubi[:, pver] (bottom) is intentionally left at 0

    # ------------------------------------------------------------------
    # Wind for lookup table: wind at cell-top in frame of moving cell
    # ------------------------------------------------------------------
    uh = np.zeros(ncol)
    for i in range(ncol):
        ut    = ubm[i, topi[i]]
        uh[i] = ut - CS[i]

    # Phase speeds: c = 0 for l=0 (only bin for ngwv=0)
    c[:, ngwv] = 0.0

    # ------------------------------------------------------------------
    # GW source: loop over columns
    # ------------------------------------------------------------------
    for i in range(ncol):
        if not use_gw_movmtn_pbl:
            if hd_idx[i] >= 0:
                # Look up momentum flux spectrum using depth and wind.
                # Note: the Fortran code appears to swap the depth and wind
                # dimensions when accessing mfcc — this translation mirrors
                # that access faithfully.
                uhmm_idx_arr = index_of_nearest(uh[i:i+1], desc.uh)
                uhmm_idx = int(uhmm_idx_arr[0])
                # Fortran: mfcc(uhmm_idx, hd_idx, 0) — dim order (wind, depth, phase)
                # Python (matching Fortran access): mfcc[uhmm_idx, hd_idx+maxuh, ngwv]
                taumm = abs(float(desc.mfcc[uhmm_idx, hd_idx[i] + desc.maxuh, ngwv]))
                taumm = taumm * qj[i] * qj[i] / AL / 1000.0
                # Assign sign based on ground-based phase speed (= CS)
                taumm = -np.copysign(taumm, CS[i])

                # Ground-based phase speed bin
                c0    = np.full(nwaves, CS[i])
                c_idx = index_of_nearest(c0, c[i, :])   # (nwaves,)

                # Assign tau at source interface levels
                # Fortran: tau(i, c_idx(i,:), topi(i):topi(i)+1) = taumm(i)
                k_src = topi[i]
                for lj in range(nwaves):
                    tau[i, c_idx[lj], k_src]     = taumm
                    tau[i, c_idx[lj], k_src + 1] = taumm
        else:
            # PBL branch: assign flux just above the launch level interface
            tau[i, ngwv, topi[i] + 1] = xpwp_src[i]

    print(f" range topi       {topi.min():d}  {topi.max():d}")
    print(f" range xpwp_src   {xpwp_src.min():.4e}  {xpwp_src.max():.4e}")
    print(f" range tau        {tau.min():.4e}  {tau.max():.4e}")

    src_level  = topi.copy()
    tend_level = topi.copy()

    return src_level, tend_level, tau, ubm, ubi, xv, yv, c, hdepth
