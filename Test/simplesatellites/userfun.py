import cgs
import numpy as np
import physics
import fileio


def Stokes_number(delv, Rd, vth, lmfp, OmegaK, rhog, rhoint = 1.4):
    v_dg=np.abs(delv)
    Rep=4*Rd*v_dg/vth/lmfp
    #following Shibaike, eq....
    CD=24/Rep*(1+0.27*Rep)**0.43+0.47*(1-np.exp(-0.04*Rep**0.38))
    St=8/3/CD*rhoint*Rd/rhog/v_dg*OmegaK
    return St

def do_stuff(system, init=False, final = False):
    if init: 
        pass 
    elif final:
        fileio.store_class(system, 'system')
        #fileio.store_class(data, 'data')##CWO: this and all stuff below does not seem to be general. Move to do_stuff perhaps
    else:
        #here we print something about planet and their resonance 
        tkeyL = system.minTimes.nameL
        tminarr = system.minTimes.tminarr

        #partices drift/growth/rel.motion
        imin = system.minTimes.dpart['imin']

        #planet and resonance 
        jres = [system.planetL[i].inxt+1 for i in range(1,system.nplanet)]
        prat = [(system.planetL[i+1].loc/system.planetL[i].loc)**(3/2) for i in range(system.nplanet-1) if system.nplanet>=2]
        prat_form = ['{:.4f}'.format(rat) for rat in prat]


        sfmt = '{:8d} {:5d} {:10.2e} {:10.2e} {} {}'
        line = sfmt.format(system.ntime, len(system.particles.massL), 
                           system.deltaT, system.time/cgs.yr, jres, prat_form)


        #output = sp.run('tail -n1 log/system_evol.log', shell=True)
        print(line)

def H_d (St, disk):
    return physics.H_d(disk.Hg,St,disk.alpha)

def del_v (St, disk):
    return 0.0
        
def dm_dt (*args):
    return 0.0

def planet_migration (*args):
    """
    LZX[24.08.04]: Now that we want to also test the resonance, so here we set the 
    migration rate a random value, not 0.0
    """
    return -0.2

def init_planets ():
    """
    this initializes the planets 

    return planet properties: 
        - time it appears
        - location
        - mass
        - composition

    [23.12.05]:it's probably better to format this as dictionaries (TBD:!)
    """
    #composition should follow what you have defined in parameters
    #and be normalized

    #[23.12.05]:let's worry about composition issues later...
    #fcomp = np.ones_like(pars.composL, dtype=float)
    #fcomp = fcomp /sum(fcomp)

    #return lists for the N-planets we have 
    timeL = [1e6*cgs.yr, 2e6*cgs.yr, 3e6*cgs.yr, 4e6*cgs.yr] 
    #some things wrong with the initial location is set to the out edge
    #about particles number
    locationL = [33*cgs.au, 50*cgs.au, 78*cgs.au, 100*cgs.au] 
    massL = [3e26, 3e26, 3e26, 3e26] 
    compoL = np.array([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])

    return timeL, locationL, massL, compoL
