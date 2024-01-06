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


def sig_mol (mu):
    """
    molecular cross section
    This probably depends (weakly) on the makeup of the gas
    But we leave this to future use
    """

    return 2e-15 *np.ones_like(mu)


class DISK (object):
    """
    Disk object including every disk properties use the methods definded in the disk_properties.py
    """

    def __init__ (self, sigmaG, temp, mu, loc, time, mcp):
        """
        disk classes are initialized ONLY
        by (sigma,temp,mu) on a 1D grid loc at a certain time
        """
        self.loc = loc
        self.sigmaG = sigmaG
        self.temp = temp
        self.mu = mu
        self.time = time
        self.mcp = mcp

    def add_auxiliary (self):
        """
        this are auxiliary disk properties that directly follow
        from the initialization
        """
        self.OmegaK = np.sqrt(cgs.gC *self.mcp/self.loc**3)      
        
        self.cs =  np.sqrt(cgs.kB*self.temp/(self.mu*cgs.mp))
        self.vth = np.sqrt(8/np.pi)*self.cs 
 
        self.vK = self.loc *self.OmegaK
        self.Hg = self.cs/self.OmegaK 
        self.rhog = self.sigmaG/(2*np.pi)**0.5/self.Hg

        self.sigmol = sig_mol (self.mu)
        self.lmfp = self.mu*cgs.mp/(self.sigmol*self.rhog)

        #alpha is too specific, must be added through add_more
        #
        #self.nu = self.alpha*self.cs*self.Hg


    def add_uservar (self, dd):
        """
        this add variables to the disk object
        """
        for key, val in dd.items():
            setattr(self, key, val)


    def add_userfun (self, dfunL):
        """
        this add functions to the disk object
        """
        for fun in dfunL:
            setattr(self, fun.__name__, fun)


    def add_user_eval (self, dfunL):
        """
        this evaluates the functions using the disk class
        (this is a bit complex)
        """
        for fun in dfunL:
            val = fun(self)
            setattr(self, fun.__name__, val)

