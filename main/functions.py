import numpy as np
import cgs
import core 
import scipy.integrate as sciint
import scipy.optimize as sciop
import disk_properties as dp
import os
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



def St_iterate (eta, vK, vth, lmfp, rhog,
                OmegaK, Rd, rhoint=1.4, Sto=0.001, errorX=1e-4, nmax=100):
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
    elif type(eta)==np.ndarray and type(Sto)==np.ndarray and len(eta)==len(Sto):
        St = Sto
    else:
        St = 1e-4

    niter = 1
    while niter<nmax:
        niter += 1
        vr = physics.radial_v(St, eta, vK)

        #call userfunction for Stokes number
        Stn = userfun.Stokes_number(delv=vr, Rd=Rd, vth=vth, lmfp=lmfp, OmegaK=OmegaK, 
                                    rhog=rhog, rhoint=rhoint)

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


def sfd_simple (msup, loc):
    """
    calculate the surface density at i based on the distance
    between i-1 and i+1
    """

    wdel = 0.5*np.log(loc[2:]/loc[:-2])

    #[22.08.19]added this hack to prevent the 
    #traffic jam at the inner boundary
    w0 = max(np.log(loc[1]/loc[0]), 0.5*np.log(loc[2]/loc[0]) )

    #extrapolate end points
    wdel = np.concatenate(([w0], wdel, [np.log(loc[-1]/loc[-2])]))

    sfd = msup /(2*np.pi *loc**2) /wdel
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






