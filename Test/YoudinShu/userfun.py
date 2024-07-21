import cgs
import numpy as np
import subprocess as sp
import physics
import pylab as pl

def init_compos_Z (material):
    """
    set the initial composition of the disk Z

    this is a 2D array, such that Z(i,j) gives
    the composition of particle i for species j
    """

    #repackage
    rc = 200*cgs.au
    if material=='silicates':
        def f_compos (rad):
            Zsil = 0.007 *np.exp(-(rad/rc)**2)
            return Zsil

    return f_compos


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
        #print(line)

        #if len(system.messages.msgL)>0: import pdb; pdb.set_trace()

        if system.time>=tplot[iplot]:
            colL = ['k', '0.8', '0.6', '0.4', 'r', 'm', 'b', 'g', (0.8,0.4,0.)]
            print('should plot stuff')
            iplot += 1

            ax.loglog(system.particles.locL/cgs.au, system.particles.sfd, '.', 
                        ms=2, lw=0.5, c=colL[iplot])

            sigana = sigma_rt (system.particles.locL, -system.particles.v_r, system.time, system.gas)

            ax.loglog(system.particles.locL/cgs.au, sigana, c='k', lw=0.5)
            ax1.semilogx(system.particles.locL/cgs.au, system.particles.sfd/sigana-1, c=colL[iplot], lw=0.5)
            ax2.semilogx(system.particles.locL/cgs.au, system.particles.sfd/sigana-1, c=colL[iplot], lw=0.5)

            for aa in [ax,ax1]:
                aa.set_xlim(0.4, 250)


            ax.set_ylim(1e-3, 2e3)
            ax1.set_ylim(-0.1, 0.1)
            ax2.set_ylim(-0.01, 0.01)

            ax.set_xlabel('distance')
            ax.set_ylabel('surface density')

            fg.savefig('testYS-plot.png', dpi=180)


def del_v (St, disk):
    return 0.0


def H_d (St, disk):
    return disk.Hg


def dm_dt (*args):
    return 0


#def Stokes_number(v_r, Rd, v_th, lmfp, Omega_K, rho_g, rhoint):
#    return Rd*rhoint /(rho_g*v_th) *Omega_K

#def Stokes_number (**kwargs):
#    return physics.Stokes_Epstein(**kwargs)

#def Stokes_number (**kwargs):
#    return physics.Stokes_general(**kwargs)


tplot = np.array([0,1e4,2e4,5e4,7e4,1e5,1.5e5,2e5,1e99]) *cgs.yr
iplot = 0

fg, (ax2, ax1, ax) = pl.subplots(3,1, figsize=(4,5),
                gridspec_kw={'height_ratios':[1,1,3]}, sharex=True
                )

