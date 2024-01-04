import numpy as np
import cgs

def get_auxiliary (mcp, loc, sigmaG, temp, mu, sigmol):

    omega = np.sqrt(cgs.gC *mcp/loc**3)      
    cs =  np.sqrt(cgs.kB*temp /(mu*cgs.mp))
    vth = np.sqrt(8/np.pi)*cs 
    Hg = cs/omega
    rhog = sigmaG/np.sqrt(2*np.pi) /Hg
    lmfp = mu*cgs.mp /(sigmol*rhog)

    return omega, cs, vth, Hg, rhog, lmfp

