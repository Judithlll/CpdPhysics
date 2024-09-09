import parameters as pars
import matplotlib.pyplot as plt
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


        if mode == 'prescribed':
            pass
        else:
            #not sure whether we should do it here or in disk object which has 
            #more parameters
            #need to initialize the grid
            self.locgrid = 10**np.linspace(np.log10(rinn), np.log10(rout), ngrid)
            self.ngrid = ngrid
            #adds the key properties

            #update_grid (self,time)


    def get_inner_radius (self, time):
        if self.mode=='prescribed':
            return dp.get_rinn (time)

    def update_oldkeyprop(self, time):
        if type(self.temp)==np.ndarray:
            self.oldkeyprop = {'sigmaG':self.sigmaG, 'temp':self.temp, 'muout':self.muout, 'time':time}


    def get_key_disk_properties (self, loc, time):
        """
        depending on the two modes:

        -- prescribed:  take directly for disk properties
        -- grid:        call the for GAS instance

        loc are the locations of the superparticles or planets
        """

        if self.mode=='prescribed':
            self.sigmaG, self.temp, self.muout = dp.key_disk_properties(loc, time, dold=self.oldkeyprop)
            self.update_oldkeyprop(time)
        else:

            if self.mode=='gridevolve':
                self.advance_grid(time)

            elif self.mode=='gridstatic':
                sigmaG, temp, muout = dp.key_disk_properties(self.locgrid, time, dold = self.oldkeyprop)
                self.sigmaG = sigmaG 
                self.temp = temp
                self.muout = muout 

                self.update_oldkeyprop(time)
                
                self.update_grid(time)


            else:
                print('No valid mode specified')
                sys.exit()


            #now intrapolate to the particle positions!
            sigmaG = np.interp(loc, self.locgrid, self.sigmaG)
            temp = np.interp(loc, self.locgrid, self.temp)
            muout = np.interp(loc, self.locgrid, self.muout)


        ## return this (perhaps as a class instance??)
        return sigmaG, temp, muout


    def update_grid (self, time):
        """
        this advances the key properties of the grid to the new time.
        we look for a disk function ...
        """
        pass


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
