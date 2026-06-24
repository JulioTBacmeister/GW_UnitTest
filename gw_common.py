"""
gw_common.py

Python/NumPy translation of gw_common.F90
Gravity wave common parameterization routines.

Indexing convention
-------------------
All level indices are 0-based (Fortran k=1 → Python k=0).
  - Midpoint arrays : shape (..., pver)      [levels 0 .. pver-1]
  - Interface arrays: shape (..., pver+1)    [levels 0 .. pver  ]
Wavenumber index l ranges from -ngwv to +ngwv.
In arrays the wavenumber dimension has size 2*ngwv+1;
  l maps to array index  l + ngwv.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional

# ---------------------------------------------------------------------------
# Cardinal direction indices
# ---------------------------------------------------------------------------
west  = 0
east  = 1
south = 2
north = 3

pi = np.pi

# ---------------------------------------------------------------------------
# Module-level state  (set by gw_common_init)
# ---------------------------------------------------------------------------
pver:   int   = 0
gravit: float = np.finfo(float).max
rair:   float = np.finfo(float).max
qbo_hdepth_scaling: float = 0.0

_tau_0_ubc: bool                  = False
_ktop:      int                   = 0
_rog:       float                 = np.finfo(float).max   # rair / gravit
_alpha:     Optional[np.ndarray]  = None
_prndl:     float                 = 0.0

# Hard-coded numerical limits
_dback   = 0.05
_taumin  = 1.e-10
_tndmax  = 400.0 / 86400.0   # 400 m/s/day → m/s/s
_umcfac  = 0.5
_ubmc2mn = 0.01


# ===========================================================================
# Helper
# ===========================================================================

def midpoint_interp(arr: np.ndarray) -> np.ndarray:
    """
    Average adjacent values along the last axis.

    Mirrors ``midpoint_interp`` in ``gw_utils.F90``.

    Parameters
    ----------
    arr : ndarray, shape (..., N)

    Returns
    -------
    ndarray, shape (..., N-1)
    """
    return 0.5 * (arr[..., :-1] + arr[..., 1:])


# ===========================================================================
# GWBand
# ===========================================================================

@dataclass
class GWBand:
    """
    A band of wavelengths into which gravity waves can be emitted.

    Parameters
    ----------
    ngwv : int
        Half-number of phase speed bins (spectrum size = 2*ngwv+1).
    dc : float
        Phase speed increment [m/s].
    wavelength : float
        Characteristic horizontal wavelength [m].
    fcrit2 : float, optional
        Critical Froude number squared.  Forced to 1.0 to match Fortran.

    Attributes
    ----------
    cref : ndarray, shape (2*ngwv+1,)
        Reference phase speeds [m/s], uniformly spaced from -ngwv*dc to +ngwv*dc.
    kwv : float
        Horizontal wavenumber [1/m].
    effkwv : float
        Effective wavenumber = fcrit2 * kwv.
    """
    ngwv:      int
    dc:        float
    wavelength: float
    fcrit2:    float = 1.0          # forced to 1.0 in the Fortran source

    cref:    np.ndarray = field(init=False, repr=False)
    kwv:     float      = field(init=False)
    effkwv:  float      = field(init=False)

    def __post_init__(self):
        self.fcrit2  = 1.0          # enforce regardless of input
        ls           = np.arange(-self.ngwv, self.ngwv + 1, dtype=float)
        self.cref    = self.dc * ls
        self.kwv     = 2.0 * pi / self.wavelength
        self.effkwv  = self.fcrit2 * self.kwv


# ===========================================================================
# Coords1D
# ===========================================================================

@dataclass
class Coords1D:
    """
    Pressure coordinate object (analogous to the Fortran ``Coords1D`` type).

    Attributes
    ----------
    ifc : ndarray, shape (ncol, pver+1)
        Interface pressures [Pa].
    del_ : ndarray, shape (ncol, pver)
        Layer pressure thickness  del_[:,k] = ifc[:,k+1] - ifc[:,k].
    rdel : ndarray, shape (ncol, pver)
        1 / del_.
    rdst : ndarray, shape (ncol, pver-1)
        Reciprocal of pressure distance between adjacent midpoints:
        rdst[:,k] = 1 / (pmid[:,k+1] - pmid[:,k])  for k = 0 .. pver-2.
    """
    ifc:  np.ndarray
    del_: np.ndarray
    rdel: np.ndarray
    rdst: np.ndarray

    @classmethod
    def from_ifc(cls, ifc: np.ndarray) -> "Coords1D":
        """Construct from interface pressure array of shape (ncol, pver+1)."""
        del_  = np.diff(ifc, axis=1)                      # (ncol, pver)
        rdel  = 1.0 / del_
        pmid  = 0.5 * (ifc[:, :-1] + ifc[:, 1:])         # (ncol, pver)
        rdst  = 1.0 / np.diff(pmid, axis=1)               # (ncol, pver-1)
        return cls(ifc=ifc, del_=del_, rdel=rdel, rdst=rdst)


# ===========================================================================
# Initialization
# ===========================================================================

def gw_common_init(pver_in, tau_0_ubc_in, ktop_in, gravit_in, rair_in,
                   alpha_in, prndl_in, qbo_hdepth_scaling_in):
    """
    Initialise module-level parameters.

    Parameters
    ----------
    pver_in : int
        Number of vertical levels.
    tau_0_ubc_in : bool
        If True, enforce tau = 0 at the upper boundary.
    ktop_in : int
        Top-most level index for GW propagation (0-based).
    gravit_in : float
        Gravitational acceleration [m s⁻²].
    rair_in : float
        Gas constant for dry air [J kg⁻¹ K⁻¹].
    alpha_in : array_like, shape (pver+1,)
        Newtonian cooling coefficients at each interface level.
    prndl_in : float
        Inverse Prandtl number.
    qbo_hdepth_scaling_in : float
        Scaling factor for QBO heating depth.
    """
    global pver, gravit, rair, qbo_hdepth_scaling
    global _tau_0_ubc, _ktop, _rog, _alpha, _prndl

    pver               = int(pver_in)
    _tau_0_ubc         = bool(tau_0_ubc_in)
    _ktop              = int(ktop_in)
    gravit             = float(gravit_in)
    rair               = float(rair_in)
    _alpha             = np.asarray(alpha_in, dtype=float)
    _prndl             = float(prndl_in)
    qbo_hdepth_scaling = float(qbo_hdepth_scaling_in)
    _rog               = rair / gravit


# ===========================================================================
# gw_prof
# ===========================================================================

def gw_prof(ncol, p, cpair, t):
    """
    Compute background-state profiles for the GW parameterisation.

    Parameters
    ----------
    ncol : int
        Number of columns.
    p : Coords1D
        Pressure coordinate object.
    cpair : float
        Specific heat of dry air at constant pressure [J kg⁻¹ K⁻¹].
    t : ndarray, shape (ncol, pver)
        Midpoint temperatures [K].

    Returns
    -------
    rhoi : ndarray, shape (ncol, pver+1)
        Interface densities [kg m⁻³].
    nm : ndarray, shape (ncol, pver)
        Midpoint Brunt-Väisälä frequencies [rad s⁻¹].
    ni : ndarray, shape (ncol, pver+1)
        Interface Brunt-Väisälä frequencies [rad s⁻¹].
    """
    n2min = 5.e-5

    ti   = np.empty((ncol, pver + 1))
    rhoi = np.empty((ncol, pver + 1))
    ni   = np.empty((ncol, pver + 1))

    # Top interface (k=0): isothermal assumption above top level
    ti[:, 0]   = t[:, 0]
    rhoi[:, 0] = p.ifc[:, 0] / (rair * ti[:, 0])
    ni[:, 0]   = np.sqrt(gravit**2 / (cpair * ti[:, 0]))

    # Interior interfaces (k = 1 .. pver-1): centred differences
    # ti at interior interfaces is the average of adjacent midpoint temps.
    ti[:, 1:pver]   = midpoint_interp(t)                        # (ncol, pver-1)
    rhoi[:, 1:pver] = p.ifc[:, 1:pver] / (rair * ti[:, 1:pver])
    # dtdp at interior interfaces; p.rdst has shape (ncol, pver-1)
    dtdp = (t[:, 1:] - t[:, :-1]) * p.rdst                     # (ncol, pver-1)
    n2   = gravit**2 / ti[:, 1:pver] * (1.0 / cpair - rhoi[:, 1:pver] * dtdp)
    ni[:, 1:pver] = np.sqrt(np.maximum(n2min, n2))

    # Bottom interface (k=pver): use bottom midpoint temperature
    ti[:, pver]   = t[:, pver - 1]
    rhoi[:, pver] = p.ifc[:, pver] / (rair * ti[:, pver])
    ni[:, pver]   = ni[:, pver - 1]

    # Midpoint B-V frequencies: average adjacent interface values
    nm = midpoint_interp(ni)                                     # (ncol, pver)

    return rhoi, nm, ni


# ===========================================================================
# Stubs for gw_diffusion routines
# ===========================================================================

class _DummyDecomp:
    """Placeholder for the TriDiagDecomp object returned by gw_ediff."""
    pass


def _gw_ediff_stub(ncol, pver, ngwv, kbot_tend, ktop, tend_level,
                   gwut, ubm, nm, rhoi, dt, prndl, gravit_arg, p, c,
                   vramp, ro_adjust=None):
    """
    Stub for gw_ediff (gw_diffusion.F90).

    Returns zero diffusivity and a dummy decomposition object.
    Replace this with the real implementation when available.
    """
    egwdffi = np.zeros((ncol, pver + 1))
    decomp  = _DummyDecomp()
    return egwdffi, decomp


def _gw_diff_tend_stub(ncol, pver, kbot_tend, ktop, field_in, dt, decomp):
    """
    Stub for gw_diff_tend (gw_diffusion.F90).

    Returns a zero tendency array.
    Replace this with the real implementation when available.
    """
    return np.zeros_like(field_in)


# ===========================================================================
# gw_drag_prof
# ===========================================================================

def gw_drag_prof(ncol, band, p, src_level, tend_level, dt,
                 t, vramp, piln, rhoi, nm, ni, ubm, ubi, xv, yv,
                 effgw, c, kvtt, q, dse, tau,
                 ro_adjust=None, kwvrdg=None, satfac_in=None,
                 lapply_effgw_in=None, lapply_vdiff=None):
    """
    Compute gravity wave drag tendencies.

    tau is modified in place (intent(inout) in Fortran).

    Parameters
    ----------
    ncol : int
        Number of columns.
    band : GWBand
    p : Coords1D
    src_level : ndarray of int, shape (ncol,)
        Level from which GWs propagate upward (0-based).
    tend_level : ndarray of int, shape (ncol,)
        Lowest level where wind tendencies are calculated (0-based).
    dt : float
        Time step [s].
    t : ndarray (ncol, pver)
        Midpoint temperatures [K].
    vramp : ndarray (pver,) or None
        Ramp-down coefficient for diffusion.
    piln : ndarray (ncol, pver+1)
        Log of interface pressures.
    rhoi : ndarray (ncol, pver+1)
        Interface densities [kg m⁻³].
    nm : ndarray (ncol, pver)
        Midpoint Brunt-Väisälä frequency [rad s⁻¹].
    ni : ndarray (ncol, pver+1)
        Interface Brunt-Väisälä frequency [rad s⁻¹].
    ubm : ndarray (ncol, pver)
        Midpoint wind projected onto wave direction [m s⁻¹].
    ubi : ndarray (ncol, pver+1)
        Interface wind projected onto wave direction [m s⁻¹].
    xv : ndarray (ncol,)
        Zonal component of unit wave propagation vector.
    yv : ndarray (ncol,)
        Meridional component of unit wave propagation vector.
    effgw : ndarray (ncol,)
        Tendency efficiency factor.
    c : ndarray (ncol, 2*ngwv+1)
        Phase speeds [m s⁻¹].  Wave l maps to index l+ngwv.
    kvtt : ndarray (ncol, pver+1)
        Molecular thermal diffusivity [m² s⁻¹].
    q : ndarray (ncol, pver, nq)
        Constituent mixing ratios.
    dse : ndarray (ncol, pver)
        Dry static energy [J kg⁻¹].
    tau : ndarray (ncol, 2*ngwv+1, pver+1)
        Wave Reynolds stress [Pa].  Modified in place.
    ro_adjust : ndarray (ncol, 2*ngwv+1, pver+1), optional
        Inertial adjustment factors.
    kwvrdg : ndarray (ncol,), optional
        Diagnosed horizontal wavenumber for ridges [1/m].
    satfac_in : float, optional
        Saturation factor (default 2.0).
    lapply_effgw_in : bool, optional
        Apply efficiency and tendency limiters (default True).
    lapply_vdiff : bool, optional
        Apply vertical diffusion (default True).

    Returns
    -------
    utgw : ndarray (ncol, pver)
        Zonal wind tendency [m s⁻²].
    vtgw : ndarray (ncol, pver)
        Meridional wind tendency [m s⁻²].
    ttgw : ndarray (ncol, pver)
        Temperature tendency [K s⁻¹].
    qtgw : ndarray (ncol, pver, nq)
        Constituent tendency [kg kg⁻¹ s⁻¹].
    egwdffi : ndarray (ncol, pver+1)
        GW effective diffusivity at interfaces.
    gwut : ndarray (ncol, pver, 2*ngwv+1)
        Per-wave wind tendency [m s⁻²].
    dttdf : ndarray (ncol, pver)
        Temperature tendency from diffusion [K s⁻¹].
    dttke : ndarray (ncol, pver)
        Temperature tendency from KE conversion [K s⁻¹].
    tau_diag : ndarray (ncol, pver+1)
        Pre-adjustment Reynolds stress for the l=0 wave (diagnostic).
    """
    ngwv = band.ngwv

    # Optional argument defaults
    satfac                = satfac_in       if satfac_in       is not None else 2.0
    do_vertical_diffusion = lapply_vdiff    if lapply_vdiff    is not None else True
    lapply_effgw          = lapply_effgw_in if lapply_effgw_in is not None else True

    kbot_tend = int(np.max(tend_level))
    kbot_src  = int(np.max(src_level))

    # -----------------------------------------------------------------------
    # Allocate output arrays
    # -----------------------------------------------------------------------
    utgw    = np.zeros((ncol, pver))
    vtgw    = np.zeros((ncol, pver))
    ttgw    = np.zeros((ncol, pver))
    qtgw    = np.zeros_like(q)
    egwdffi = np.zeros((ncol, pver + 1))
    gwut    = np.zeros((ncol, pver, 2 * ngwv + 1))
    dttdf   = np.zeros((ncol, pver))
    dttke   = np.zeros((ncol, pver))

    # -----------------------------------------------------------------------
    # Stress profiles: scan bottom → top
    # (Fortran: do k = kbot_src, ktop, -1)
    # -----------------------------------------------------------------------
    for k in range(kbot_src, _ktop - 1, -1):

        d = _dback + kvtt[:, k]     # total diffusivity (ncol,)

        for l in range(-ngwv, ngwv + 1):
            li = l + ngwv           # array index for wavenumber l

            ubmc     = ubi[:, k] - c[:, li]          # (ncol,)
            mask_src = src_level >= k                 # (ncol,) bool

            # ---- saturation stress -----------------------------------------
            sign_same = (ubmc > 0.0) == (ubi[:, k + 1] > c[:, li])

            if kwvrdg is not None:
                tausat = np.where(
                    mask_src & sign_same,
                    np.abs(kwvrdg * rhoi[:, k] * ubmc**3 / (satfac * ni[:, k])),
                    0.0)
            else:
                tausat = np.where(
                    mask_src & sign_same,
                    np.abs(band.effkwv * rhoi[:, k] * ubmc**3 / (satfac * ni[:, k])),
                    0.0)

            if ro_adjust is not None:
                tausat = np.where(mask_src,
                                  tausat * np.sqrt(ro_adjust[:, li, k]),
                                  tausat)

            # ---- damped stress and new tau ----------------------------------
            ubmc2 = np.maximum(ubmc**2, _ubmc2mn)

            if kwvrdg is not None:
                mi  = (ni[:, k] / (2.0 * kwvrdg * ubmc2) *
                       (_alpha[k] + ni[:, k]**2 / ubmc2 * d))
                # wrk computed but not multiplied into taudmp in the kwvrdg branch
                # (matches the Fortran: taudmp = tau(:,l,k+1) with no exp)
                taudmp = tau[:, li, k + 1]
            else:
                mi     = (ni[:, k] / (2.0 * band.kwv * ubmc2) *
                          (_alpha[k] + ni[:, k]**2 / ubmc2 * d))
                wrk    = -2.0 * mi * _rog * t[:, k] * (piln[:, k + 1] - piln[:, k])
                taudmp = tau[:, li, k + 1] * np.exp(wrk)

            # Zero out values below taumin before taking the minimum
            tausat = np.where(tausat <= _taumin, 0.0, tausat)
            taudmp = np.where(taudmp <= _taumin, 0.0, taudmp)

            tau[:, li, k] = np.where(mask_src,
                                      np.minimum(taudmp, tausat),
                                      tau[:, li, k])

    # Upper boundary condition
    if _tau_0_ubc:
        tau[:, :, _ktop] = 0.0

    # Diagnostic: pre-adjustment tau for l=0 wave
    tau_diag = tau[:, ngwv, :].copy()   # l=0 → index ngwv

    # -----------------------------------------------------------------------
    # Apply efficiency to the completed stress profile
    # (Fortran: do k = ktop, kbot_tend+1)
    # -----------------------------------------------------------------------
    if lapply_effgw:
        for k in range(_ktop, kbot_tend + 2):
            mask = ((k - 1) <= tend_level)                    # (ncol,)
            tau[mask, :, k] *= effgw[mask, np.newaxis]

    # -----------------------------------------------------------------------
    # Compute tendencies from stress divergence: scan top → bottom
    # (Fortran: do k = ktop, kbot_tend)
    # -----------------------------------------------------------------------
    ubt = np.zeros((ncol, pver))

    for k in range(_ktop, kbot_tend + 1):

        ubt[:, k] = 0.0

        for l in range(-ngwv, ngwv + 1):
            li = l + ngwv

            # Stress-divergence wind tendency
            ubtl = gravit * (tau[:, li, k + 1] - tau[:, li, k]) * p.rdel[:, k]

            # Limiter 1: prevent sign reversal of (u - c)
            ubtl = np.minimum(ubtl,
                              _umcfac * np.abs(c[:, li] - ubm[:, k]) / dt)

            if not lapply_effgw:
                ubtl = np.minimum(ubtl, _tndmax)

            mask_tend = k <= tend_level                        # (ncol,)

            # Fortran sign(a,b) = copysign: apply sign of (c-ubm) to |ubtl|
            gwut_kl = np.copysign(ubtl, c[:, li] - ubm[:, k])
            gwut[:, k, li] = np.where(mask_tend, gwut_kl, gwut[:, k, li])
            ubt[:, k]      = np.where(mask_tend,
                                       ubt[:, k] + gwut[:, k, li],
                                       ubt[:, k])

        # Limiter 2: cap total tendency at tndmax
        if lapply_effgw:
            abs_ubt       = np.abs(ubt[:, k])
            ubt_lim_ratio = np.where(abs_ubt > _tndmax,
                                     _tndmax / abs_ubt, 1.0)
            ubt[:, k]     = ubt_lim_ratio * ubt[:, k]
        else:
            ubt_lim_ratio = np.ones(ncol)

        for l in range(-ngwv, ngwv + 1):
            li = l + ngwv

            gwut[:, k, li] *= ubt_lim_ratio

            # Protect against floating-point underflow
            gwut[:, k, li] = np.where(
                np.abs(gwut[:, k, li]) < 1.e-15, 0.0, gwut[:, k, li])

            # Adjust tau on the interface below to reflect bounded tendency
            mask_tend = k <= tend_level
            tau[:, li, k + 1] = np.where(
                mask_tend,
                tau[:, li, k] + np.abs(gwut[:, k, li]) * p.del_[:, k] / gravit,
                tau[:, li, k + 1])

        # Project mean wind tendency onto zonal / meridional components
        mask_tend  = k <= tend_level
        utgw[:, k] = np.where(mask_tend, ubt[:, k] * xv, utgw[:, k])
        vtgw[:, k] = np.where(mask_tend, ubt[:, k] * yv, vtgw[:, k])

        if vramp is not None:
            utgw[:, k] *= vramp[k]
            vtgw[:, k] *= vramp[k]

    # -----------------------------------------------------------------------
    # Legacy path: undo Sean Santos efficiency mods (lapply_effgw = False)
    # -----------------------------------------------------------------------
    if not lapply_effgw:
        for k in range(_ktop, kbot_tend + 2):
            mask = ((k - 1) <= tend_level)
            tau[mask, :, k] *= effgw[mask, np.newaxis]
        for k in range(_ktop, kbot_tend + 1):
            gwut[:, k, :] *= effgw[:, np.newaxis]
            utgw[:, k]    *= effgw
            vtgw[:, k]    *= effgw

    # -----------------------------------------------------------------------
    # Vertical diffusion (calls to gw_diffusion routines)
    # -----------------------------------------------------------------------
    if do_vertical_diffusion:
        egwdffi, decomp = _gw_ediff_stub(
            ncol, pver, ngwv, kbot_tend, _ktop, tend_level,
            gwut, ubm, nm, rhoi, dt, _prndl, gravit, p, c, vramp,
            ro_adjust=ro_adjust)

        if q.ndim == 3:
            for m in range(q.shape[2]):
                qtgw[:, :, m] = _gw_diff_tend_stub(
                    ncol, pver, kbot_tend, _ktop, q[:, :, m], dt, decomp)

        dttdf = _gw_diff_tend_stub(
            ncol, pver, kbot_tend, _ktop, dse, dt, decomp)

    # -----------------------------------------------------------------------
    # KE → heat conversion:  dttke -= (ubm - c) * gwut
    # -----------------------------------------------------------------------
    for l in range(-ngwv, ngwv + 1):
        li = l + ngwv
        dttke[:, _ktop:kbot_tend + 1] -= (
            (ubm[:, _ktop:kbot_tend + 1] - c[:, li:li + 1]) *
            gwut[:, _ktop:kbot_tend + 1, li])

    ttgw = dttke + dttdf

    if vramp is not None:
        for k in range(_ktop, kbot_tend + 1):
            ttgw[:, k] *= vramp[k]

    return utgw, vtgw, ttgw, qtgw, egwdffi, gwut, dttdf, dttke, tau_diag


# ===========================================================================
# calc_taucd
# ===========================================================================

def calc_taucd(ncol, ngwv, tend_level, tau, c, xv, yv, ubi):
    """
    Calculate Reynolds stress for waves propagating in each cardinal direction.

    Parameters
    ----------
    ncol : int
    ngwv : int
    tend_level : ndarray of int, shape (ncol,)   [0-based]
    tau : ndarray (ncol, 2*ngwv+1, pver+1)
    c   : ndarray (ncol, 2*ngwv+1)
    xv  : ndarray (ncol,)
    yv  : ndarray (ncol,)
    ubi : ndarray (ncol, pver+1)

    Returns
    -------
    taucd : ndarray (ncol, pver+1, 4)
        Columns in dimension-2 correspond to [west, east, south, north].
    """
    taucd = np.zeros((ncol, pver + 1, 4))

    # Interface wind at tend_level+1 for each column
    idx      = np.arange(ncol)
    ubi_tend = ubi[idx, tend_level + 1]   # (ncol,)

    for k in range(_ktop, int(np.max(tend_level)) + 2):

        taub = np.zeros(ncol)
        tauf = np.zeros(ncol)
        mask_k = (k - 1) <= tend_level       # (ncol,)

        for l in range(-ngwv, ngwv + 1):
            li    = l + ngwv
            tausg = np.copysign(tau[:, li, k], c[:, li] - ubi[:, k])

            below = c[:, li] < ubi_tend       # waves "behind" the wind
            taub  = np.where(mask_k &  below, taub + tausg, taub)
            tauf  = np.where(mask_k & ~below, tauf + tausg, tauf)

        # Assign to cardinal directions
        xpos = xv > 0.0
        ypos = yv > 0.0

        taucd[:, k, east]  = np.where(mask_k,
                                       np.where(xpos, tauf, taub) * xv,
                                       taucd[:, k, east])
        taucd[:, k, west]  = np.where(mask_k,
                                       np.where(xpos, taub, tauf) * xv,
                                       taucd[:, k, west])
        taucd[:, k, north] = np.where(mask_k,
                                       np.where(ypos, tauf, taub) * yv,
                                       taucd[:, k, north])
        taucd[:, k, south] = np.where(mask_k,
                                       np.where(ypos, taub, tauf) * yv,
                                       taucd[:, k, south])

    return taucd


# ===========================================================================
# momentum_flux
# ===========================================================================

def momentum_flux(tend_level, taucd):
    """
    Calculate the momentum flux from below the GW region.

    Parameters
    ----------
    tend_level : ndarray of int, shape (ncol,)
    taucd : ndarray (ncol, pver+1, 4)

    Returns
    -------
    um_flux : ndarray (ncol,)
    vm_flux : ndarray (ncol,)
    """
    idx      = np.arange(len(tend_level))
    k_bottom = tend_level + 1
    um_flux  = taucd[idx, k_bottom, east]  + taucd[idx, k_bottom, west]
    vm_flux  = taucd[idx, k_bottom, north] + taucd[idx, k_bottom, south]
    return um_flux, vm_flux


# ===========================================================================
# momentum_fixer
# ===========================================================================

def momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw):
    """
    Subtract GW momentum change from wind tendencies below the GW region
    to enforce momentum conservation.

    Parameters
    ----------
    tend_level : ndarray of int, shape (ncol,)
    p : Coords1D
    um_flux : ndarray (ncol,)
    vm_flux : ndarray (ncol,)
    utgw : ndarray (ncol, pver)   modified in place
    vtgw : ndarray (ncol, pver)   modified in place
    """
    idx = np.arange(len(tend_level))
    # Reciprocal total mass from surface to source level (dp/g)
    rdm = gravit / (p.ifc[idx, pver] - p.ifc[idx, tend_level + 1])
    du  = -um_flux * rdm
    dv  = -vm_flux * rdm

    for k in range(int(np.min(tend_level)) + 1, pver):
        mask       = k > tend_level
        utgw[:, k] = np.where(mask, utgw[:, k] + du, utgw[:, k])
        vtgw[:, k] = np.where(mask, vtgw[:, k] + dv, vtgw[:, k])


# ===========================================================================
# energy_change
# ===========================================================================

def energy_change(dt, p, u, v, dudt, dvdt, dsdt):
    """
    Calculate the change in total energy from tendencies.

    Parameters
    ----------
    dt : float
    p : Coords1D
    u, v : ndarray (ncol, pver)   winds at start of time step
    dudt, dvdt : ndarray (ncol, pver)
    dsdt : ndarray (ncol, pver)   heating tendency

    Returns
    -------
    de : ndarray (ncol,)   net energy change per column [J m⁻²]
    """
    de = np.sum(
        p.del_ / gravit * (
            dsdt
            + dudt * (u + dudt * 0.5 * dt)
            + dvdt * (v + dvdt * 0.5 * dt)
        ),
        axis=1)
    return de


# ===========================================================================
# energy_fixer
# ===========================================================================

def energy_fixer(tend_level, p, de, ttgw):
    """
    Subtract net energy gain/loss from heating tendency below the GW region.

    Parameters
    ----------
    tend_level : ndarray of int, shape (ncol,)
    p : Coords1D
    de : ndarray (ncol,)
    ttgw : ndarray (ncol, pver)   modified in place
    """
    idx   = np.arange(len(tend_level))
    de_dm = -de * gravit / (p.ifc[idx, pver] - p.ifc[idx, tend_level + 1])

    for k in range(int(np.min(tend_level)) + 1, pver):
        mask       = k > tend_level
        ttgw[:, k] = np.where(mask, ttgw[:, k] + de_dm, ttgw[:, k])


# ===========================================================================
# coriolis_speed
# ===========================================================================

def coriolis_speed(band, lat):
    """
    Absolute value of the local Coriolis frequency divided by kwv [m/s].

    Parameters
    ----------
    band : GWBand
    lat  : ndarray (ncol,)   latitude in radians

    Returns
    -------
    ndarray (ncol,)
    """
    omega_earth = 2.0 * pi / 86400.0
    return np.abs(np.sin(lat) * 2.0 * omega_earth / band.kwv)


# ===========================================================================
# adjust_inertial
# ===========================================================================

def adjust_inertial(band, tend_level, u_coriolis, c, ubi, tau, ro_adjust):
    """
    Adjust tau and compute ro_adjust for inertial gravity waves.

    Blocks upward wave propagation through levels where the Coriolis effect
    is too strong.  tau is modified in place; ro_adjust is overwritten.

    Parameters
    ----------
    band : GWBand
    tend_level : ndarray of int, shape (ncol,)
    u_coriolis : ndarray (ncol,)
        |f| / kwv  [m/s].
    c : ndarray (ncol, 2*ngwv+1)
    ubi : ndarray (ncol, pver+1)
    tau : ndarray (ncol, 2*ngwv+1, pver+1)   modified in place
    ro_adjust : ndarray (ncol, 2*ngwv+1, pver+1)   overwritten
    """
    ngwv   = band.ngwv
    ncol   = len(tend_level)
    nwaves = 2 * ngwv + 1

    # unblocked_mask[i, li] = True until Coriolis blocks wave propagation
    unblocked_mask = np.ones((ncol, nwaves), dtype=bool)
    ro_adjust[:]   = 0.0

    # Iterate bottom → top through interface levels where tau is defined
    for k in range(int(np.max(tend_level)) + 1, _ktop - 1, -1):

        valid = k <= tend_level + 1              # (ncol,) bool

        # Speed of wind relative to each wave at this level
        speed_diff    = np.abs(ubi[:, k:k+1] - c)        # (ncol, nwaves)
        unblocked_here = speed_diff > u_coriolis[:, np.newaxis]

        # Cumulative AND: once blocked, stays blocked
        unblocked_mask = np.where(
            valid[:, np.newaxis],
            unblocked_mask & unblocked_here,
            unblocked_mask)

        valid_2d   = valid[:, np.newaxis]                # (ncol, 1)
        unblocked  = valid_2d &  unblocked_mask
        blocked    = valid_2d & ~unblocked_mask

        # ro_adjust = 1 - (u_coriolis / (ubi - c))^2  where unblocked
        # Guard against division by zero with a safe denominator
        denom       = np.where(
            unblocked,
            ubi[:, k:k+1] - c,        # (ncol, nwaves)
            1.0)                       # dummy, won't be used
        uc_ratio    = u_coriolis[:, np.newaxis] / denom
        ro_adjust[:, :, k] = np.where(unblocked, 1.0 - uc_ratio**2, ro_adjust[:, :, k])

        tau[:, :, k] = np.where(blocked, 0.0, tau[:, :, k])
