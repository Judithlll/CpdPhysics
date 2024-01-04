import cgs
import numpy as np

## A very simple PPD...

rinn = 0.1 *cgs.au
rout = 100 *cgs.au

#this is necessary apparently...
sigmol=2e-15


def Mcp_t (t):
    return cgs.mSun


def M_influx (t0, tEnd):
    return 0.0


def key_disk_properties (rad, t, dold=None):
    """
    simple power-law solutions
    """
    rc = 30*cgs.au
    sigma = 1700 *(rad/cgs.au)**-1.5 *np.exp(-rad/rc)
    mgas = 2.34 *np.ones_like(rad)
    temp = 300 *(rad/cgs.au)**-0.5

    return sigma, temp, mgas


def Omega_K (rad, t, Mcp):
    return 2e-7 *(rad/cgs.au)**-1.5

