import cgs
import numpy as np

def do_stuff (system, init=False, final=False):

    if init:
        pass
    else:
        tkeyL = system.minTimes.nameL
        tminarr = system.minTimes.tminarr

        sfmt = '{:5d} {:5d} {:10.2e} {:3d} {:10.2e}'
        line = sfmt.format(system.ntime, len(system.particles.massL), system.deltaT, 
                                            tminarr.argmin(), system.time/cgs.yr)
        print(line)


def del_v (St, disk):
    return 0.0


def H_d (St, disk):
    return disk.Hg


def dm_dt (*args):
    return 0


def Stokes_number(v_r, Rd, v_th, lmfp, Omega_K, rho_g, rhoint):
    return Rd*rhoint /(rho_g*Rd) *Omega_K
