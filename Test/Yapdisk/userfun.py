import cgs 
import numpy as np 
import subprocess as sp 
import physics 
import matplotlib.pyplot as plt


def do_stuff (system, init=False, final=False):

    if init:
        plt.figure()
        plt.title('Turbulence-Dominated')
        plt.xlabel('Distance')
        plt.ylabel(r'Surface Density $\Sigma (g/cm^2)$')
        #plt.ylim(1e1, 1e4)
        #plt.xlim(1e-1, 4e1)
        plt.yscale('log')
        plt.xscale('log')
        pass

    elif final:
        plt.legend()
        plt.savefig('turbulence-dominated.png')
    else:
        tkeyL = system.minTimes.nameL
        tminarr = system.minTimes.tminarr

        #partices drift/growth/rel.motion
        imin = system.minTimes.dpart['imin']

        sfmt = '{:8d} {:5d} {:10.2e} {:3d} {:2d} {:2d} {:10.2e}'
        line = sfmt.format(system.ntime, len(system.particles.massL), system.deltaT, 
                                            tminarr.argmin(), imin[0],imin[1], system.time/cgs.yr)
        
        if system.ntime%100 ==0:
            plt.plot(system.gas.locgrid/cgs.au, system.gas.sigmaG, label = 't = {:.2e}'.format(system.time/cgs.yr)) 
        #output = sp.run('tail -n1 log/system_evol.log', shell=True)
        print(line)

        #if len(system.messages.msgL)>0: import pdb; pdb.set_trace()



def dm_dt (*args):
    return 0 

def H_d (St, disk):
    return disk.Hg*(1+St/disk.alpha)**(-0.5)

def Stokes_number (**kwargs):
    #here the stokes number should be very complex
    return physics.Stokes_Epstein(**kwargs)
