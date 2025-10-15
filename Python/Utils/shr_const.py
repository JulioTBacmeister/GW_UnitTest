# shr_const.py

import numpy as np

class ShrConst:
    # Physical constants
    PI = 3.14159265358979323846  # pi
    CDAY = 86400.0  # sec in calendar day ~ sec
    SDAY = 86164.0  # sec in siderial day ~ sec
    OMEGA = 2.0 * PI / SDAY  # earth rot ~ rad/sec
    REARTH = 6.37122e6  # radius of earth ~ m
    G = 9.80616  # acceleration of gravity ~ m/s^2

    STEBOL = 5.67e-8  # Stefan-Boltzmann constant ~ W/m^2/K^4
    BOLTZ = 1.38065e-23  # Boltzmann's constant ~ J/K/molecule
    AVOGAD = 6.02214e26  # Avogadro's number ~ molecules/kmole
    RGAS = AVOGAD * BOLTZ  # Universal gas constant ~ J/K/kmole
    MWDAIR = 28.966  # molecular weight dry air ~ kg/kmole
    MWWV = 18.016  # molecular weight water vapor
    RDAIR = RGAS / MWDAIR  # Dry air gas constant ~ J/K/kg
    RWV = RGAS / MWWV  # Water vapor gas constant ~ J/K/kg
    ZVIR = (RWV / RDAIR) - 1.0  # RWV/RDAIR - 1.0
    KARMAN = 0.4  # Von Karman constant
    PSTD = 101325.0  # standard pressure ~ pascals
    PDB = 0.0112372  # ratio of 13C/12C in Pee Dee Belemnite (C isotope standard)

    TKTRIP = 273.16  # triple point of fresh water ~ K
    TKFRZ = 273.15  # freezing T of fresh water ~ K
    TKFRZSW = TKFRZ - 1.8  # freezing T of salt water ~ K
    ZSRFLYR = 3.0  # ocn surf layer depth for diurnal SST cal ~ m

    RHODAIR = PSTD / (RDAIR * TKFRZ)  # density of dry air at STP ~ kg/m^3
    RHOFW = 1.000e3  # density of fresh water ~ kg/m^3
    RHOSW = 1.026e3  # density of sea water ~ kg/m^3
    RHOICE = 0.917e3  # density of ice ~ kg/m^3
    CPDAIR = 1.00464e3  # specific heat of dry air ~ J/kg/K
    CPWV = 1.810e3  # specific heat of water vap ~ J/kg/K
    CPVIR = (CPWV / CPDAIR) - 1.0  # CPWV/CPDAIR - 1.0
    CPFW = 4.188e3  # specific heat of fresh h2o ~ J/kg/K
    CPSW = 3.996e3  # specific heat of sea h2o ~ J/kg/K
    CPICE = 2.11727e3  # specific heat of fresh ice ~ J/kg/K
    LATICE = 3.337e5  # latent heat of fusion ~ J/kg
    LATVAP = 2.501e6  # latent heat of evaporation ~ J/kg
    LATSUB = LATICE + LATVAP  # latent heat of sublimation ~ J/kg
    CONDICE = 2.1  # thermal conductivity of ice ~ W/m/K
    KAPPA_LAND_ICE = CONDICE / (RHOICE * CPICE)  # Diffusivity of heat in land ice

    TF0 = 6.22e-2  # The freezing temperature at zero pressure in sub-ice-shelf ocean cavities ~ C
    DTF_DP = -7.43e-8  # The coefficient for the term proportional to the (limited) pressure in the freezing temperature in sub-ice-shelf ocean cavities ~ C Pa^{-1}
    DTF_DS = -5.63e-2  # The coefficient for the term proportional to salinity in the freezing temperature in sub-ice-shelf ocean cavities ~ C PSU^{-1}
    DTF_DPDS = -1.74e-10  # The coefficient for the term proportional to salinity times pressure in the freezing temperature in sub-ice-shelf ocean cavities ~ C PSU^{-1} Pa^{-1}
    OCN_REF_SAL = 34.7  # ocn ref salinity (psu)
    ICE_REF_SAL = 4.0  # ice ref salinity (psu)

    SPVAL = 1.0e30  # special missing value
    SPVAL_TOLMIN = 0.99 * SPVAL  # min spval tolerance
    SPVAL_TOLMAX = 1.01 * SPVAL  # max spval tolerance
    SPVAL_AERODEP = 1.e29  # special aerosol deposition

    # Water Isotope Ratios in Vienna Standard Mean Ocean Water (VSMOW)
    VSMOW_18O = 2005.2e-6  # 18O/16O in VMSOW
    VSMOW_17O = 379.e-6  # 17O/16O in VMSOW
    VSMOW_16O = 0.997628  # 16O/Tot in VMSOW
    VSMOW_D = 155.76e-6  # 2H/1H in VMSOW
    VSMOW_T = 1.85e-6  # 3H/1H in VMSOW
    VSMOW_H = 0.99984426  # 1H/Tot in VMSOW
    RSTD_H2ODEV = 1.0  # Rstd Dev Use

# Usage
if __name__ == "__main__":
    print("PI:", ShrConst.PI)
    print("Acceleration of gravity (G):", ShrConst.G)
    print("Stefan-Boltzmann constant (STEBOL):", ShrConst.STEBOL)
    # Add more prints to test other constants if needed
