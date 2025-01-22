import numpy as np
import cgs
import core 
import scipy.integrate as sciint
import scipy.optimize as sciop
import disk_properties as dp
import os, sys
import parameters as pars
import shutil
import userfun
import physics
import subprocess as sp

def clear_dir(directory):
    """
    mainly used in plot a series of figures. If the directory is not empty, 
    then clear the directory
    """
    if os.path.exists(directory):
        #[24.04.21]cwo: I have no idea what is going on here...
        #               is shutil really necessary?
        #               also put this in fileio.py
        for filename in os.listdir(directory):
            file_path = os.path.join(directory, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    #the directory does not exist yet. Create
    else:
        sp.run(['mkdir', directory])


def get_res_chain(res_setL):
    newsetL = []
    for ss in res_setL:
        for k,nwss in enumerate(newsetL):
            if nwss.intersection(ss):
                newsetL[k] = nwss.union(ss) #merges
                ss = None
                break

        if ss is not None:
            newsetL.append(ss)

    return newsetL

def locate_chain(res_chainL, planetnumber):
    
    numset = {planetnumber}
    for chain in res_chainL:
        if numset.issubset(chain):
            return chain

def sample_equal (fx, x0, xf, Nsample=100, ttol=1e-5):
    """
    [23.12.30]: function copied from /NewLagrance by CWO

     samples the CDF of fx(x) equally over [x0,xf] 
    by Nsample points

    NOTE: I use brentq with a for-loop, which is probably VERY SLOW
    """

    Fnorm, err = sciint.quad(fx, x0, xf)

    def y_root (x,yaim=0,xold=x0,yold=0):
        """
        y:  y-value in range (0,1)
        x:  root
        """
        Fy, err = sciint.quad(fx, xold, x)
        return Fy/Fnorm +yold -yaim


    xL = [x0]; xold = x0; xE = xf; yold=0.
    for k in range(Nsample):
        y = (k+0.5) /Nsample
        xr, info = sciop.brentq(y_root, xold, xE, args=(y,xold,yold), full_output=True)
        xL.append(xr)

        #[22.12.13]some strange bug, perhaps related to non-continuity?!
        yck = y_root(xr,0,x0)
        if abs(y-yck)>ttol:
            xr, info = sciop.brentq(y_root, xold, xE, args=(y,x0,0), full_output=True)

        xE = xf
        xold = xr*1.0
        yold = y*1.0

    return xL

def load_dict_from_file (fname):
    dout = {}
    with open(fname,'r') as f:
        line = f.readline()
        while line:
            if line[0]!='#' and len(line)>=5:
                key, sep, val, *dumL = line.split()

                #boolean
                if val=='True':
                    dout[key] = True
                elif val=='False':
                    dout[key] = False

                elif val.isalnum():
                    dout[key] = val

                else:
                    try:
                        dout[key] = int(val)
                    except:
                        dout[key] = float(val)

            line = f.readline()
    return dout


def Stokes_number (disk, size, rhoint, Sto=0.001, errorX=1e-4, nmax=100):
    """
    Get Stokes number and drift velocity by fixed-point iteration
    
    The Stokes number rely on the size of particles, also the size rely on 
    Stokes_number in the pre-fractor CD (and more clearly the Re), so here 
    we need iteration. 

    maybe there can be an analitically solution, but its a little complex
    """
    #put in physics.py

    if Sto is None:
        St = 1e-4
    elif type(disk.eta)==np.ndarray and type(Sto)==np.ndarray and len(disk.eta)==len(Sto):
        St = Sto
    else:
        St = 1e-4

    niter = 1
    while niter<nmax:
        niter += 1
        vr = physics.radial_v(St, disk.eta, disk.vK)

        #[25.01.20]cwo: made this a bit more general
        #call userfunction for Stokes number

        if pars.dragmodel=='Epstein':
            Stn = physics.Stokes_Epstein (size, rhoint, disk.vth, disk.rhog, disk.OmegaK)
            Stn *= np.sqrt(8/np.pi) #difference b/w sound speed and thermal velocity

            #this is how Youdin & Shu do it..
            vr = -2*St *disk.eta *disk.vK

        #[25.01.01]cwo Stokes number is fixed??
        #just for ism final project
        elif pars.dragmodel=='fixed_St':
            St = np.ones_like(loc)*pars.fixed_St

        #user defines
        else:
            Stn = userfun.Stokes_number(delv=vr, Rd=size, vth=disk.vth, lmfp=disk.lmfp, OmegaK=disk.OmegaK, 
                                    rhog=disk.rhog, rhoint=rhoint)

        #better to do relative error 
        error = abs(St/Stn-1)

        if error.max()<errorX:
            break
        else:
            St = Stn

    #print('# iterations:', niter)
    #if niter>50 and type(eta)==np.ndarray: import pdb; pdb.set_trace()

    return Stn, vr


def get_stokes_number(disk,t,sPLmtot,rhoint):
    mphy = sPLmtot
    Rd=(mphy/(rhoint*4/3*np.pi))**(1/3)

    eta=disk.eta
    v_K=disk.vK
    v_th=disk.vth
    lmfp=disk.lmfp
    rho_g=disk.rhog
    Omega_K=disk.OmegaK

    St,v_r = St_iterate(eta,v_K,v_th,lmfp,rho_g,Omega_K,Rd) 
    return St


def sfd_spline (msup, loc, **kwargs):
    """
    also fails LJ test
    """
    from scipy.interpolate import PchipInterpolator

    locmid = np.sqrt(loc[1:]*loc[:-1])
    locmidext = np.concatenate((
        [loc[0]*np.sqrt(loc[0]/loc[1])],
        locmid,
        [np.sqrt(loc[-1]/loc[-2])*loc[-1]]))

    mcum = np.concatenate(([0],np.cumsum(msup)))

    # Create a monotonic piecewise spline
    pchip = PchipInterpolator(np.log(locmidext), mcum)
    px = pchip.derivative()

    return px(np.log(loc)) /(2*np.pi*loc**2)


def sfd_fixedbin (msup, loc, pgrid, specloc=[], noff=1):
    """
    We simply fix the bins and count...
    """

    rinn = pgrid[0]
    rout = pgrid[-1]

    fac = (rout/rinn)**(1/(len(pgrid)-1))

    sfd_avg = np.zeros_like(loc)
    for koff in range(noff):
        grid = pgrid *fac**(koff/noff)

        if koff!=0:
            grid = np.concatenate(([rinn],grid[:-1],[rout]))

        rbin = np.sqrt(grid[1:]*grid[:-1]) #the midpoints
        wbin = grid[1:] -grid[:-1] #bin spacing
        nbin = len(rbin)

        Abin = 2*np.pi*rbin*wbin
        
        ig = np.searchsorted(grid, loc) -1  #points to the bin index

        ig = np.maximum(0, ig) #hack...


        mbin = np.bincount(ig, weights=msup, minlength=nbin)    #total mass in bin 0, 1, ...
        mbin /= Abin                            #transform to surface density
        sfd = mbin[ig]                            #corresponding mass for particles
        sfd_avg += sfd /noff

    return sfd


def sfd_special (msup, loc, specloc):
    """
    like sfd_simple, but accounting for special locations
    """

    #sort the special locations
    specL = sorted(list(specloc)+[np.inf])

    loc0 = 0
    wdelL = []
    for k, loc1 in enumerate(specL):
        ii = (loc>loc0) *(loc<loc1)

        if sum(ii)==0:
            continue
        if sum(ii)==1:
            print('[functions.sfd_special]BUG: 1 particle -- cannot calculate sfd?')
            sys.exit()

        locmid = np.sqrt(loc[ii][1:]*loc[ii][:-1])
        wdel = np.concatenate(([2*(locmid[0]-loc[ii][0])], 
                                np.diff(locmid),
                               [2*(loc[ii][-1]-locmid[-1])]))

        wdelL.append(wdel)
        loc0 = loc1

    wdel = np.concatenate(wdelL)
    sfd = msup /(2*np.pi *loc) /wdel
    return sfd


def sfd_simple (msup, loc, specloc):
    """
    calculate the surface density at i based on the distance
    between i-1 and i+1
    """
    # import matplotlib.pyplot as plt

    # wdel = 0.5*np.log(loc[2:]/loc[:-2])
    #
    # #[22.08.19]added this hack to prevent the 
    # #traffic jam at the inner boundary
    # w0 = max(np.log(loc[1]/loc[0]), 0.5*np.log(loc[2]/loc[0]) )
    #
    # # extrapolate end points
    # wdel = np.concatenate(([w0], wdel, [np.log(loc[-1]/loc[-2])]))

    #
    locmid = np.sqrt(loc[1:]*loc[:-1])
    wdel = np.diff(locmid)
    wdel = np.concatenate(([2*(locmid[0]-loc[0])], wdel, [2*(loc[-1]-locmid[-1])]))

    ## around the speloc, we also take 2*(locmid-locspe_particle) for wdel.
    # specidx = np.searchsorted(loc, specloc)
    # for idx in specidx:
    #     wdel[idx] = 2*(locmid[idx]- loc[idx]) 
    #     wdel[idx-1] = 2*(loc[idx-1]-locmid[idx-2])

    sfd = msup /(2*np.pi *loc) /wdel
    return sfd


def sfd_chris (msup, loc):
    """
    sfd determined from back
    """
    wdel = np.diff(loc)
    wdel = np.concatenate((wdel, [wdel[-1]]))
    sfd = msup /(2*np.pi *loc) /wdel
    return sfd


def dotMgTscale(radL,deltaTL):
    """
    calculate the time scale of mass flux
    """
    v=np.diff(radL*cgs.RJ,axis=1)/deltaTL
    vTimesr=v*radL[:,0:-1]*cgs.RJ
    first=np.diff(vTimesr,axis=1)/deltaTL[:-1]
    second=v[:,:-1]*np.diff(vTimesr,axis=1)/np.diff(radL,axis=1)[:,:-1]/cgs.RJ
    Tscale=1/((first-second)/v[:,:-1]/radL[:,:-2]/cgs.RJ)
    return Tscale

def determine_type (val):
    """
    this determines the type of val

    hist
        [23.12.03]:copied from /NewLagrange by Chris
    """

    if val=='None':
        return None
    #boolean value
    elif val=='True':
        return True
    elif val=='False':
        return False

    #integer
    elif val.isdigit():
        return int(val)

    else:
        #float/string value
        try:
            out = float(val)
        except:
            out = val

        return out






