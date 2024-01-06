import numpy as np
import cgs
import core 
import scipy.integrate as sciint
import scipy.optimize as sciop

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



def St_iterate (eta,v_K,v_th,lmfp,rho_g,
                Omega_K,R_d,rhoint=1.4,Sto=0.001, errorX=1e-4, nmax=100):
    """
    obtain Stokes number and drift velocity by fixed-point iteration
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
        v_r =-2*St/(St**2+1)*eta*v_K
        v_dg=np.abs(v_r)
        Rep=4*R_d*v_dg/v_th/lmfp
        #following Shibaike, eq....
        CD=24/Rep*(1+0.27*Rep)**0.43+0.47*(1-np.exp(-0.04*Rep**0.38))
        Stn=8/3/CD*rhoint*R_d/rho_g/v_dg*Omega_K

        #better to do relative error 
        error = abs(St/Stn-1)

        if error.max()<errorX:
            break
        else:
            St = Stn


    #print('# iterations:', niter)
    #if niter>50 and type(eta)==np.ndarray: import pdb; pdb.set_trace()

    return Stn, v_r


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


def epsilon_PA (PlanetsLoc,PlanetsMass,cross_p):
    """
    abnormal now
    
    different particles has different size and so different accretion rate, so maybe should change the total mass of particles???
    """
    # out = system.gas.get_key_disk_properties (PlanetsLoc, system.time)
    # disk = core.DISK (*out, PlanetsLoc, system.time)
    # disk.add_auxiliary ()
    # disk.user_difined ()

    # St=get_stokes_number(disk,system.time,cross_p[1],system.rhoint)
    # eta=disk.eta

    mus=PlanetsMass/cross_p.mcp

    # Hg= disk.Hg
    
    Hp=cross_p.Hg*(1+cross_p.St/ cross_p.alpha*(1+2*cross_p.St)/(1+cross_p.St))
    hp=Hp/PlanetsLoc

    delv_o_vK=0.52*(mus*cross_p.St)**(1/3)+cross_p.eta/(1+5.7*(mus/cross_p.eta**3*cross_p.St))
    
    P_eff=1/np.sqrt((0.32*np.sqrt(mus*delv_o_vK/ cross_p.St/ cross_p.eta**2))**(-2)+(0.39*mus/ cross_p.eta/hp)**(-2)) #Liu & Ormel 2018
    return P_eff




def M_critical (eta, St, mcp):

    M_critical=1/8*eta**3*St *mcp #Shibaike 2019
    return M_critical   
