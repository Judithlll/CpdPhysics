import numpy as np
import cgs
error=1e-8

def St_iterate(eta,v_K,v_th,lmfp,rho_g,
                Omega_K,R_d,rhoint=1.4,Sto=0.001):
        """
        obtain Stokes number and drift velocity by fixed-point iteration
        """
        #put in physics.py
    
        St=Sto
        nither=1
        while(nither<40):
            nither += 1
            v_r=-2*St/(St**2+1)*eta*v_K
            v_dg=np.abs(v_r)
            Rep=4*R_d*v_dg/v_th/lmfp
            #following Shibaike, eq....
            CD=24/Rep*(1+0.27*Rep)**0.43+0.47*(1-np.exp(-0.04*Rep**0.38))
            Stn=8/3/CD*rhoint*R_d/rho_g/v_dg*Omega_K
            delta_St=abs(St-Stn)
            
            if np.max(delta_St)<error:
                break
            else:
                St=Stn

        return Stn, v_r

def get_stokes_number(disk,t,ParticlePL,rhoint):
    r, mphy, mtot = ParticlePL
    Rd=(mphy/(rhoint*4/3*np.pi))**(1/3)

    eta=disk.eta(r,t)
    v_K=disk.vK(r,t)
    v_th=disk.vth(r,t)
    lmfp=disk.lmfp(r,t)
    rho_g=disk.rhog(r,t)
    Omega_K=disk.OmegaK(r,t)

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


def epsilon_PA (system,PlanetsLoc,PlanetsMass,cross_p):
    """
    abnormal now
    
    different particles has different size and so different accretion rate, so maybe should change the total mass of particles???
    """

    St=get_stokes_number(system.disk,system.time,cross_p,system.rhoint)
    eta=system.disk.eta(PlanetsLoc,system.time)

    mus=PlanetsMass/ system.disk.Mcp_t(system.time)

    Hg= system.disk.Hg(PlanetsLoc,system.time)
    
    Hp=Hg*(1+St/ system.disk.alpha*(1+2*St)/(1+St))
    hp=Hp/PlanetsLoc

    delv_o_vK=0.52*(mus*St)**(1/3)+eta/(1+5.7*(mus/eta**3*St))
    
    P_eff=1/np.sqrt((0.32*np.sqrt(mus*delv_o_vK/ St/ eta**2))**(-2)+(0.39*mus/ eta/hp)**(-2)) #Liu & Ormel 2018
    return P_eff

def M_critical(system,PlanetsLoc,cross_p):

    St=get_stokes_number(system.disk,system.time,cross_p,system.rhoint)
    eta=system.disk.eta(PlanetsLoc,system.time)

    M_critical=1/8*eta**3*St *system.disk.Mcp_t(system.time) #Shibaike 2019
    return M_critical   