import cgs 
import numpy as np 
import subprocess as sp 
import physics 


def do_stuff (system, init=False, final=False):

    if init:
        pass
    else:
        tkeyL = system.minTimes.nameL
        tminarr = system.minTimes.tminarr

        #partices drift/growth/rel.motion
        imin = system.minTimes.dpart['imin']

        sfmt = '{:8d} {:5d} {:10.2e} {:3d} {:2d} {:2d} {:10.2e}'
        line = sfmt.format(system.ntime, len(system.particles.massL), system.deltaT, 
                                            tminarr.argmin(), imin[0],imin[1], system.time/cgs.yr)

        #output = sp.run('tail -n1 log/system_evol.log', shell=True)
        print(line)

        #if len(system.messages.msgL)>0: import pdb; pdb.set_trace()



def dm_dt (*args):
    return 0 

def Stokes_number (**kwargs):
    #here the stokes number should be very complex
    pass  
