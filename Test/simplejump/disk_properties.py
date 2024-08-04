import cgs
import numpy as np
import parameters as pars

## A very simple PPD...

rinn = pars.dgasgrid['rinn']
rout = pars.dgasgrid['rout']

#this is necessary apparently...
sigmol=2e-15


def Mcp_t (t):
    return cgs.Msun

def dot_Mg (t, mode=None):
    return 0


def user_add_var ():
    """
    a list of attributes to be added to the disk class
    """
    return {'alpha':1e-3}

def user_add_fun ():
    return [dot_Md]

def user_add_eval ():
    return [eta, m_gas]


def key_disk_properties (rad, t, dold=None):
    """
    simple power-law solutions
    """
    sigma = 1700 *(rad/cgs.RJ)**-1.5 *np.exp(-rad/rout)
    mgas = m_gas(rad)
    temp = 300 *(rad/cgs.RJ)**-0.5

    return sigma, temp, mgas

def eta (disk):
    return 0.01*np.ones_like(disk.loc)

def m_gas (rad):
    return 2.34*np.ones_like(rad)

def H_d (Hg, St):
    return 1.0*np.ones_like(disk.loc)

def dot_Md (disk):
    return 1e-3*0.1*cgs.MJ/1e6/cgs.yr


