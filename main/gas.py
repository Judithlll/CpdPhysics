import parameters as pars
import disk_properties as dp
import numpy as np
import cgs
import sys

class GAS ():
    """
    [23.12.13]: copied the essentials from /Newlagrange
    """

    def __init__(self, gasL, dcomposL, time=0, rinn=-1, rout=None, mode='prescribed', ngrid=200, 
                    solmode='default', ttol=1e-4, dcompos={}, **kwargs):
        """
        initialize disk class. 
        
        Modes:
            prescribed: will look for *disk.py* to find the gas surface 
                            density at (r,t). No grids are used
            grid:       evolves the gas AND vapor surface density *on a 
                        grid*; can accounts for source terms.
        """
        self.mode = mode
        self.oldkeyprop = None


        if mode=='gridstatic':
            #need to initialize the grid
            self.loc = 10**np.linspace(np.log10(rinn),np.log10(rout), ngrid)
            
            #adds the key properties
            #update_grid (self,time)

        else:
            #no structure is made
            pass


    def get_inner_radius (self, time):
        if self.mode=='prescribed':
            return dp.get_rinn (time)
    #TBD: evolve mcp and rcp here?          

    def get_key_disk_properties (self, loc, time):
        """
        depending on the two modes:

        -- prescribed:  take directly for disk properties
        -- grid:        call the for GAS instance

        loc are the locations of the superparticles or planets
        """

        if self.mode=='prescribed':
            sigmaG, temp, muout = dp.key_disk_properties(loc, time, dold=self.oldkeyprop)
            if type(temp)==np.ndarray:
                self.oldkeyprop = {'sigmaG':sigmaG, 'temp':temp, 'muout':muout}
        else:

            if self.mode=='gridevolve':
                self.advance_grid(time)

            elif self.mode=='gridstatic':
                self.update_grid(time)


            else:
                print('No valid mode specified')
                sys.exit()


            #now intrapolate to the particle positions!
            sigmaG = np.interp(loc, self.loc, self.sigmaG)
            temp = np.interp(loc, self.loc, self.temp)
            muout = np.interp(loc, self.loc, self.mu)


        ## return this (perhaps as a class instance??)
        return sigmaG, temp, muout


    def update_grid (self, time):
        """
        this advances the key properties of the grid to the new time.
        we look for a disk function ...
        """
        self.sigmaG, self.temp, self.mu =\
                dp.key_disk_properties (self.loc, time)


    def advance_grid (self, dt, mode='implicit', split=False, source_func=None, **kwargs):
        """
        grid-only // TBD

        advance grid by dt.  Save the state as newSigma (intermediat surfac density)
        """

        #...
        # do not split operator and solve the PDE all together
        newSigmaL, vgas = self.adv_implicit(dt, source_func, **kwargs)

        self.newSigmaL = newSigmaL
        self.new_dt = dt
