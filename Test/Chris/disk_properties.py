#all disk properties
import numpy as np
import cgs
import ode 
from scipy.integrate import quad

alpha = 1e-2
rinn = 6*cgs.RJ
rout = 27*cgs.RJ
frac=0.2
ratio=0.01
# Mcp=0.4*cgs.MJ
sigmamol=2e-15
rgg=1e-7
Mcp0=0.4*cgs.MJ

tdep= 3e6*cgs.yr2s  #depletion tiescale of gas (constant: 3e6 yr)
tgap= 1e6*cgs.yr2s  #time when the gap was built

def Mcp_t(t): #t should begin at tgap
    """
    mass of central body with time (integration of dot_Mg)
    TBD: dot_Mg prescription has changed
    """
    
    if type(t)==np.ndarray:
        if len(t.shape)==1:
            Mcp=np.array([])
            for time in t:
                Mcp=np.append(Mcp,Mcp0+quad(dot_Mg,0,time)[0])
        else:
            Mcp=np.zeros_like(t)
            for i in range(len(t)):
                for j in range(len(t[i])):
                    print([i,j])
                    Mcp[i,j]=Mcp0+quad(dot_Mg,0,t[i,j])[0]

    else:
        Mcp = Mcp0 + quad(dot_Mg,0,t)[0]
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


def dot_Mg(t):
    """
    accretion rate through the disk as fn of time (Shibaike+2019, eq.2)
    """
    if np.min(t)>tgap: 
        dotMg=dotMg_gap()*np.exp(-(t-tgap)/tdep)
    else:
        dotMg = dotMg_gap()
    return dotMg

def dot_Md(t):
    """
    get the solid mass flow

    ratio: the dust to gas ratio
    """
    Mfd=ratio*dot_Mg(t)
    return Mfd

def v_K(r,t):
    v=np.sqrt(cgs.gC*Mcp_t(t)/r)
    return v

def Sigma_g(r,t):
    """
    get the surface density at loc and time 
    """
    Sg=dot_Mg(t)/2/np.pi/rout*r**(3/2)/viscosity(r,t)*(-2/9*r**(-1/2)+2/3*rout*r**(-3/2))
    return Sg

def T_d(r,t):
    Td=(3*cgs.gC*Mcp_t(t)*dot_Mg(t)/8/np.pi/cgs.sigmaSB/r**3)**(1/4)
    return Td

# def T_d(r,t):
#     '''
#     Assume that the disk is viscously heated, calculate the midplane tempurature

#     Parameters:
#     Mcp: mass of the central planet
#     dotMg: gas accration rate
#     r: distance from the central planets
#     rgg: the ratio of the surface dencity of grains that affect the temparature to gas surface dencity
#     '''
#     Ti=100
#     n=1
#     while n<20:

#         if 0<Ti<160:
#             kapa=450*(Ti/160)**2*rgg
#         elif Ti<0:
#             print('Something wrong! T=',Ti,'<0')
#             break
#         else:
#             kapa=450*rgg
    
#         tau=kapa*Sigma_g(r,t)

#         g=(3/8*tau+1/4.8/tau)**(1/4)
#         Td=(3*cgs.gC*Mcp_t(t)*dot_Mg(t)/8/np.pi/cgs.sigmaSB/r**3)**(1/4)*g #Shibaike 2019 equation (5)
#         diff= abs(Ti-Td)
#         Ti=Td
#      #   print(n,'diff_Td=',diff)
#         n+=1

#     return Td

def m_g():
    """
    get the mean molecular mass of disk
    """
    MeanMolecularMass=2.34*cgs.mp
    return MeanMolecularMass

def c_s(r,t):
    """
    get sound speed at loc and time
    """

    cs=np.sqrt(cgs.kB*T_d(r,t)/m_g())

    return cs

def Omega_K(r,t):
    OK=np.sqrt(cgs.gC*Mcp_t(t)/r**3)
    return OK

def H_g(r,t):
    """
    get scale height
    """    
    GasScaleHeight=c_s(r,t)/Omega_K(r,t)
    return GasScaleHeight

def viscosity(r,t):
    """
    get viscosity
    """
    nu=alpha*c_s(r,t)*H_g(r,t)
    return nu

def v_th(r,t):
    """
    get the thermal velocity 
    """
    thermal_velocity=np.sqrt(8/np.pi)*c_s(r,t)
    return thermal_velocity

def l_mfp(r,t):
    """
    get the mean free path

    """
    MeanFreePath=m_g()/sigmamol/rho_g(r,t)
    return MeanFreePath

def rho_g(r,t):
    """
    get the density of gas
    """
    rhogas=Sigma_g(r,t)/(2*np.pi)**0.5/H_g(r,t)
    return rhogas


def eta(r,t):
    """
    v=(1-eta)vK
    """
    e=0.477577*cgs.kB*(r**5.5-4.84615*r**4.5*rout)/cgs.gC/Mcp_t(t)/m_g()/(r**4.5-3*r**3.5*rout)*(cgs.gC*Mcp_t(t)*dot_Mg(t)/r**3/cgs.sigmaSB)**(1/4)
    return e
