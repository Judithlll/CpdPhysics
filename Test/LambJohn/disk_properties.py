import cgs
import numpy as np
import parameters as pars 
## A very simple PPD...

rinn = 0.1 *cgs.au
#rout = 800 *cgs.au

fraction = 0.0 
#this is necessary apparently...
#sigmol=2e-15

def r_out (t):
    return pars.dgasgrid['rout']

def Mcp_t (t):
    return cgs.Msun

def dot_Mg (t, mode=None):
    return fraction*cgs.Msun/1e6/cgs.yr


def M_influx (t0, tEnd):
    return 0.0

def user_add_var ():
    """
    a list of attributes to be added to the disk class
    """
    return {'alpha':1e-3}

def user_add_fun ():
    return [dot_Md]

def user_add_eval ():
    return [eta, m_gas]


def be_ta(t):
    beta0 = 500
    return beta0*np.exp(-t/3e6/cgs.yr)

def sigma_gas (rad, t):
    """
    simple power-law solutions
    """
    #
    sigma = be_ta(t) *(rad/cgs.au)**(-1)
    return sigma 

def key_disk_properties (rad, t, dold=None):
    """
    simple power-law solutions
    """
    #
    sigma = sigma_gas(rad, t)
    mgas = m_gas(rad)

    #from schoonenberg 
    Tvisc = 350*(rad/cgs.au)**(-3/4)
    Tirr = 177*(rad/cgs.au)**(-0.5)
    temp = (Tvisc**4 + Tirr**4)**0.25 

    return sigma, temp, mgas

eta0 = 0.0015 
def eta (disk):
    #from L&J 2014
    eta = eta0*(disk.loc/cgs.au)**(1/2)
    return eta 

def m_gas (rad):
    return 2.34*np.ones_like(rad)

#def H_d (Hg, St):
#    return 1.0*np.ones_like(disk.loc)

def dot_Md (time):
    return 1.0+0.01*dot_Mg(time)

