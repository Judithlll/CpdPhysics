#all disk properties
import numpy as np
import cgs
import ode 
from scipy.integrate import quad

alpha = 1e-4
rinn = 5.89*cgs.RJ  #the location of Io
rout = 27*cgs.RJ
frac=0.2
ratio=0.01
# Mcp=0.4*cgs.MJ
sigmol=2e-15
rgg=1e-7
Mcp0=0.4*cgs.MJ

tdep= 3e6*cgs.yr2s  #depletion tiescale of gas (constant: 3e6 yr)
tgap= 1e6*cgs.yr2s  #time when the gap was built

def Mcp_t(t):
    """
    mass of central body with time (integration of dot_Mg)
    TBD: dot_Mg prescription has changed
    """

    Mcp=Mcp0+dotMg_gap()*tdep*(np.exp(tgap/tdep)-np.exp((tgap-t)/tdep))
    return Mcp

def dotMg_gap():
    """
    get the gas mass flow when the gap firstly biult

    fraction: the fraction of Jupiter mass which is accreted 
    """
    Mfg=frac*cgs.MJ/1e6/cgs.yr2s
    return Mfg

def dot_Mg(t):
    """
    accretion rate through the disk as fn of time (Shibaike+2019, eq.2)
    """
    if np.min(t)>tgap: 
        dotMg=dotMg_gap()*np.exp(-(t-tgap)/tdep)
    else:
        dotMg=dotMg_gap()  
    return dotMg   

def M_influx(t0,tEnd):
    """
    gas inflow into the CJD
    """
    Minflux,error=quad(dot_Mg,t0,tEnd)
    # Minflux=ratio*dotMg_gap()*tdep*(np.exp(tgap/tdep) -np.exp(-(t-tgap)/tdep))
    return Minflux

def dot_Md(dotMg):
    """
    get the solid mass flow

    ratio: the dust to gas ratio
    """
    Mfd=ratio*dotMg
    return Mfd

def v_K(r,t):
    v=np.sqrt(cgs.gC*Mcp_t(t)/r)
    return v

def Sigma_g(r,cs,OmegaK,dotMg):
    """
    get the surface density at loc and time 
    """
    Sg=dotMg/2/np.pi/rout*r**(3/2)/(alpha*cs**2/OmegaK)*(-2/9*r**(-1/2)+2/3*rout*r**(-3/2))
    return Sg

# def Sigma_g(r,t):
#     """
#     get the surface density at loc and time 
#     """
#     Sg=dot_Mg(t)/2/np.pi/rout*r**(3/2)/(alpha*c_s(r,t)**2/Omega_K(r,t))*(-2/9*r**(-1/2)+2/3*rout*r**(-3/2))
#     return Sg

# def T_d(r,t):
#     Td=(3*cgs.gC*Mcp_t(t)*dot_Mg(t)/8/np.pi/cgs.sigmaSB/r**3)**(1/4)
#     return Td

def key_disk_properties (r,t):
    
    #add a judgement if r is an array, because there is a comparison
    Mcp=Mcp_t(t)
    OmegaK=Omega_K(r,t,Mcp)
    dotMg=dot_Mg(t)
    if type(r)==np.ndarray:
        Ti=np.ones_like(r)
        n=1

        while n<20:
            if np.min(Ti)<0:
                print('Something wrong! T=',Ti,'<0')
                break

            #kapa=np.zeros_like(r)

            #kapa[Ti<160]=450*(Ti[Ti<160]/160)**2*rgg
            #kapa[Ti>=160]=450*rgg

            kapa = np.where(Ti<160, 450*(Ti/160)**2, 450) *rgg  
            cs = c_s(Ti)

            sigG = Sigma_g(r,cs,OmegaK,dotMg)
            tau = kapa*sigG

            g = (3/8*tau+1/4.8/tau)**(1/4)
            Td = (3*cgs.gC*Mcp_t(t)*dot_Mg(t)/8/np.pi/cgs.sigmaSB/r**3)**(1/4)*g #Shibaike 2019 equation (5)
            diff= abs(Ti/Td-1)
            Ti = Td
            # print(n,'diff_Td=',diff.max())
            if diff.max()<1e-4:
                break
            n+=1
    else:
        Ti=100
        n=1
        while n<20:
            kapa=np.zeros_like
            if 0<Ti<160:
                kapa=450*(Ti/160)**2*rgg
            elif Ti.any()<0:
                print('Something wrong! T=',Ti,'<0')
                break
            else:
                kapa=450*rgg

            cs=c_s(Ti)
            sigG = Sigma_g(r,cs,OmegaK,dotMg)
            tau = kapa*sigG

            g=(3/8*tau+1/4.8/tau)**(1/4)
            Td=(3*cgs.gC*Mcp_t(t)*dot_Mg(t)/8/np.pi/cgs.sigmaSB/r**3)**(1/4)*g #Shibaike 2019 equation (5)
            diff= abs(Ti-Td)
            Ti=Td
        #   print(n,'diff_Td=',diff)
            n+=1
    mu = 2.34*np.ones_like(sigG)
    return sigG,Td,mu

def T_d(sigmag,kapa,Mcp,dotMg,loc):
    '''
    Assume that the disk is viscously heated, calculate the midplane tempurature

    Parameters:
    Mcp: mass of the central planet
    dotMg: gas accration rate
    r: distance from the central planets
    rgg: the ratio of the surface dencity of grains that affect the temparature to gas surface dencity
    '''
    tau=kapa*sigmag
    g=(3/8*tau+1/4.8/tau)**(1/4)
    Td=(3*cgs.gC*Mcp*dotMg/8/np.pi/cgs.sigmaSB/loc**3)**(1/4)*g
    return Td

def c_s(T):
    return np.sqrt(cgs.kB*T/m_g())

def m_g():
    """
    get the mean molecular mass of disk
    """
    MeanMolecularMass=2.34*cgs.mp
    return MeanMolecularMass

# def c_s(r,t):
#     """
#     get sound speed at loc and time
#     """

#     cs=np.sqrt(cgs.kB*T_d(r,t)/m_g())

#     return cs

def Omega_K(r,t,Mcp):
    OK=np.sqrt(cgs.gC*Mcp/r**3)
    return OK

def H_g(cs,OmegaK):
    """
    get scale height
    """    
    GasScaleHeight=cs/OmegaK
    return GasScaleHeight

def viscosity(cs,Hg):
    """
    get viscosity
    """
    nu=alpha*cs*Hg
    return nu

def v_th(cs):
    """
    get the thermal velocity 
    """
    thermal_velocity=np.sqrt(8/np.pi)*cs
    return thermal_velocity

def l_mfp(rhog,mg):
    """
    get the mean free path

    """
    MeanFreePath=mg/sigmamol/rhog
    return MeanFreePath

def rho_g(Sigmag,Hg):
    """
    get the density of gas
    """
    rhogas=Sigmag/(2*np.pi)**0.5/Hg
    return rhogas


def eta(r,Mcp,dotMg,mg):
    """
    as in v=(1-eta)vK

    eta = 1/2 *n *h**2

    h = aspect ratio = Hgas/r

    n = d\logP /d\log r

    (This looks too complex)
    """

    e=0.477577*cgs.kB*(r**5.5-4.84615*r**4.5*rout)/cgs.gC/Mcp/mg/(r**4.5-3*r**3.5*rout)*(cgs.gC*Mcp*dotMg/r**3/cgs.sigmaSB)**(1/4)
    return e
