import cgs
import numpy as np
import parameters as pars
## A very simple PPD...

rinn = pars.dgasgrid['rinn'] 
#rout = 800 *cgs.au

#this is necessary apparently...
#sigmol=2e-15

def r_out (t):
    return pars.dgasgrid['rout'] 
def Mcp_t (t):
    return cgs.Msun

def dot_Mg (t, mode=None):
    return 0

def M_influx (t0, tEnd):
    return 0.0

def user_add_var ():
    """
    a list of attributes to be added to the disk class
    """
    #follow Huang et al. 2024
    return {'alpha':pars.dgasprop['alpha']}

def user_add_fun ():
    return [dot_Md]

def user_add_eval ():
    return [eta, m_gas]


def key_disk_properties (rad, t, dold=None):
    """
    simple power-law solutions
    """
    #see Huang 2024

    rc = pars.dgasgrid['rout'] 
    Mdisk = pars.dgasprop['Mcp0'] 
    sigma =  Mdisk/(2*np.pi*rc*(1-np.exp(-1)))*np.exp(-rad/rc)/rad
    mgas = m_gas(rad)

    temp = 200 *(rad/cgs.au)**-0.5

    return sigma, temp, mgas

def eta (disk):
    #eq.20 YS
    hgas = disk.cs/disk.vK
    #hgas = 0.033 *(disk.loc/cgs.au)**0.25
    return 1.6 *hgas**2

def m_gas (rad):
    return 2.34*np.ones_like(rad)

def dot_Md (disk):
    return 1.0

