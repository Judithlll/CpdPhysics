import cgs
import numpy as np
import subprocess as sp
import physics
import pylab as pl
import parameters as pars
import matplotlib.pyplot as plt

def init_compos_Z (material):
    """
    set the initial composition of the disk Z

    this is a 2D array, such that Z(i,j) gives
    the composition of particle i for species j
    """

    #repackage
    if material=='silicates':
        def f_compos (rad):
            Zsil = 0.0134
            return Zsil

    return f_compos

def PIM(*args):
    return 1e99

def M_critical(*args):
    return 0.0

def add_planet(system, mode):
    """
    The add planet here should have some options, but both of them can
    apply the add_planet under System, which add planet from a candidate 
    list, so here we should modify the candidate list.

    The options are:
    1. given: the add from the given information
    2. random: the add from the planetesimal generating rate
    """
    #First modify the candidate 
    ##LZX[24.08.27]:maybe check the planetesimal formation models
    
    #Then add the planet to the planet list
    if mode == 'capture':
        #maybe also some newly captured seeds
        #pnew = core.PLANET(time, mcp0, rho)
        #system.planet_candidate.append(pnew) 
        
        if len(system.planet_candidate)>0:
            for planet in system.planet_candidate:
                if planet.starttime <= system.time:
                    system.add_planet(planet)
                    system.planet_candidate.remove(planet)
    elif mode == 'insitu':
        pass

    #Finally we should sort the planets according to the location.

    return system

def Stokes_number(delv, Rd, vth, lmfp, Omega_K, rho_g, rhoint):
    return np.pi*Rd*rhoint/2/rho_g

def init_compos (compos):
    dcompos = {}
    if compos=='silicates':
        dcompos['Z_init'] = init_compos_Z (compos)

    return dcompos


def g_r (rad, gas, d=1.5):
    Z_0 = init_compos_Z ('silicates')
    Z0 = Z_0(rad)

    sig0, temp, mmw = gas.get_key_disk_properties(rad,0.0)
    #sig0, mmw = disk.sigma_gas_ini (rad)
    #sig0 = Z0*disk.sigma_gas_ini (rad)
    return rad**(d+1) *sig0*Z0


def sigma_rt (rad, vdr, time, gas, d=1.5):
    """
    Youdin & Shu solution
    """
    ri = rad *(1 -(d-1)*vdr*time/rad)**(-1/(d-1))
    sig = rad**(-d-1) *g_r(ri,gas)
    return sig



def do_stuff (system, init=False, final=False):
    global iplot, tplot, fg, ax

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

        #sfmt = '{:10.3e} {:10.3e} {:10.3e}'        
        #line = sfmt.format(system.time, system.particles.v_r[-1], system.particles.locL[-1]/cgs.au)

    if final: 
        plt.xlim(0,200)
        plt.plot(system.particles.locL/cgs.au, system.particles.sfd, 
                  label= r'$\alpha$ = {:.0e}  St={:.0e}'.format(pars.dgasprop['alpha'], pars.fixed_St))
        for p in system.planetL:
            plt.axvline(p.loc/cgs.au, c='r')

        filename = r'final_sfd+{:.0e}+{:.0e}'.format(pars.dgasprop['alpha'], pars.fixed_St)
        plt.legend()
        plt.savefig(filename+'.jpg')
        plt.close()

        #need  to save something to input into RadMC3D
        np.savez('saved_data/'+filename+'.npz', loc=system.particles.locL, sfd=system.particles.sfd, 
                                  hg=system.particles.Hg, mstar = pars.dgasprop['Mcp0'], 
                                  alpha = pars.dgasprop['alpha'], St = pars.fixed_St)
        




def del_v (St, disk):
    return 0.0

def epsilon_PA (planetLoc,planetMass,cross_p):
    """
    Get the pebble accretion rate
    #"""

    St = cross_p.St
    eta = cross_p.eta

    #prefer to call this qp
    mus = planetMass/cross_p.mcp
    Hp = physics.H_d(cross_p.Hg, St, cross_p.alpha)
    # Hp=cross_p.Hg*(1+cross_p.St/ cross_p.alpha*(1+2*cross_p.St)/(1+cross_p.St))
    hp = Hp/planetLoc

    #expression from eq.18 S+19
    delv_o_vK = 0.52*(mus*St)**(1/3)+eta/(1+5.7*(mus/eta**3*cross_p.St))
    
    P_eff = 1/np.sqrt( (0.32*np.sqrt(mus*delv_o_vK/ St/ eta**2))**(-2)
                      +(0.39*mus/ eta/hp)**(-2)) #Liu & Ormel 2018

    return P_eff


def dm_dt (particles):
    return 0 


#def Stokes_number(v_r, Rd, v_th, lmfp, Omega_K, rho_g, rhoint):
#    return Rd*rhoint /(rho_g*v_th) *Omega_K

#def Stokes_number (**kwargs):
#    return physics.Stokes_Epstein(**kwargs)

#def Stokes_number (**kwargs):
#    return physics.Stokes_general(**kwargs)


def planet_migration(*args):
    return 0.0

def init_planets (): 
    timeL = [0.5*cgs.yr, 1*cgs.yr, 2*cgs.yr, 3*cgs.yr]
    locationL = [20.0*cgs.au, 26.2*cgs.au, 34.3*cgs.au, 75.5*cgs.au]
    massL = [16.1*cgs.Mea, 19.8*cgs.Mea, 24.2*cgs.Mea, 30.2*cgs.Mea] #tbd 
    compoL =np.array([[1.0], [1.0], [1.0], [1.0]])

    return timeL, locationL, massL, compoL

def plot_sfd(*args):
    return 
