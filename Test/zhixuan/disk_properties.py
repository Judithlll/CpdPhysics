#all disk properties
import numpy as np
import cgs
import ode 
from scipy.integrate import quad
import physics
import parameters as pars

alpha = pars.dgasprop['alpha']
rinn = pars.dgasgrid['rinn']#the location of Io
rout = pars.dgasgrid['rout']
frac = pars.dgasprop['frac']
ratio = pars.dgasprop['dgratio']
#Mcp=0.4*cgs.MJ
sigmol= pars.dgasprop['sigmol']#
rgg= pars.dgasprop['rgg']
Mcp0= pars.dgasprop['Mcp0']
meanmol = pars.dgasprop['meanmol']

tdep= pars.dgasprop['tdep']  #depletion tiescale of gas (constant: 3e6 yr)
tgap= pars.dgasprop['tgap']  #time when the gap was built

def Mcp_t(t):
    """
    mass of central body with time (integration of dot_Mg)
    """

    Mcp=Mcp0+dotMg_gap()*tdep*(np.exp(tgap/tdep)-np.exp((tgap-t)/tdep))
    return Mcp

def user_add_fun ():
    """
    a list of functions to be added
    """
    return dot_Mg, beta_sigG, beta_P, dot_Md


def user_add_var ():
    """
    a list of attributes to be added
    """
    return {'alpha':alpha, 'rinn':rinn, 'rout':rout, 'ratio':ratio, 'beta_T':beta_T, 'beta_nu':beta_nu, 'beta_cs':beta_cs}


def user_add_eval ():
    return eta, mg

def dotMg_gap():
    """
    get the gas mass flow when the gap firstly biult

    fraction: the fraction of Jupiter mass which is accreted 
    """
    Mfg = frac*cgs.MJ/1e6/cgs.yr2s
    return Mfg

def mg(disk):
    return disk.mu* cgs.mp

def dot_Mg(t):
    """
    accretion rate through the disk as fn of time (Shibaike+2019, eq.2)
    """
    if t>tgap: 
        dotMg=dotMg_gap()*np.exp(-(t-tgap)/tdep)
    else:
        dotMg=dotMg_gap()  
    return dotMg   

def M_influx(t0,tEnd):
    """
    gas inflow into the CJD
    """
    Minflux=quad(dot_Mg,t0,tEnd)[0]
    # Minflux=ratio*dotMg_gap()*tdep*(np.exp(tgap/tdep) -np.exp(-(t-tgap)/tdep))
    return Minflux

def dot_Md(time):
    """
    get the solid mass flow

    ratio: the dust to gas ratio
    """
    Mfd=ratio*dot_Mg(time)
    return Mfd

def v_K(r,t):
    v=np.sqrt(cgs.gC*Mcp_t(t)/r)
    return v

def Sigma_g(r,cs,OmegaK,dotMg):
    """
    get the surface density at loc and time 
    """
    #Shibaike 2017 eq2
    Sg = dotMg/2/np.pi/rout*r**(3/2)/(alpha*cs**2/OmegaK)*(-2/9*r**(-1/2)+2/3*rout*r**(-3/2))
    #Sg = dotMg *OmegaK/3/np.pi/alpha/cs**2
    return Sg


def key_disk_properties (rad, t, dold=None):
    """
    This returns the key disk properties
        - surface density
        - temperature
        - mean molecular weight

    at location rad and time t
    """
    ## CWO: what would really help in the iteration is to provide estimated solutions

    #turn into array if necessary
    if type(rad) in [float, np.float64]:
        r = np.array([rad])
        returnfloat = True
    else:
        r = rad
        returnfloat = False
    
    #add a judgement if r is an array, because there is a comparison
    Mcp = Mcp_t(t)
    OmegaK = physics.Omega_K(r,Mcp)
    dotMg = dot_Mg(t)

    nmax = 10
    nn = 1

    #active indices
    ii = np.ones_like(r, dtype=bool)


    mu = meanmol *np.ones_like(r)

    #start from the guess solution
    if dold is not None and len(r)==len(dold['temp']):
        Td = dold['temp']
        sigG = dold['sigmaG']
    else:
        Td = 10*np.ones_like(r)
        sigG = np.ones_like(r)

    while nn<nmax:
        if np.min(Td)<0:
            print('Something wrong! T=', Ti, '<0')
            break

        #kapa=np.zeros_like(r)

        #kapa[Ti<160]=450*(Ti[Ti<160]/160)**2*rgg
        #kapa[Ti>=160]=450*rgg

        Ti = Td[ii]

        kapa = np.where(Ti<160, 450*(Ti/160)**2, 450) *rgg  

        cs = physics.c_s (Ti, mu[ii]) #
        #cs = c_s(Ti)

        sigG[ii] = Sigma_g(r[ii],cs,OmegaK[ii],dotMg)
        tau = kapa*sigG[ii]

        g = (3/8*tau+1/4.8/tau)**(1/4)
        Td[ii] = (3*cgs.gC*Mcp*dotMg/8/np.pi/cgs.sigmaSB/r[ii]**3)**(1/4)*g #Shibaike 2019 equation (5)

        diff = abs(Ti/Td[ii]-1)
        #Ti = Td
        # print(n,'diff_Td=',diff.max())
        #if diff.max()<1e-4: break

        #this inserts more False if condition is met
        ii[ii==True] = diff>1e-4

        if np.all(ii==False):
            break

        nn += 1

    #mu = MeanMolecularWeight()*np.ones_like(sigG)


    if returnfloat:
        return sigG[0],Td[0],mu[0]
    else:
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

#the dependence on location of some properties
beta_T = -3/4
beta_nu = beta_T+3/2
beta_cs = beta_T/2
def beta_sigG(loc):
    #from Shibaika.2017, just below the eq15
    beta_sigG = -beta_nu + loc/(loc-3*rout)
    return beta_sigG

def beta_P(loc):
    betasig = beta_sigG(loc)
    return betasig + beta_cs -3/2 

    
# def m_g():
#     """
#     get the mean molecular mass of disk
#     """
#     MeanMolecularMass=MeanMolecularWeight()*cgs.mp
#     return MeanMolecularMass

# def c_s(r,t):
#     """
#     get sound speed at loc and time
#     """

#     cs=np.sqrt(cgs.kB*T_d(r,t)/m_g())

#     return cs



# def H_g(cs,OmegaK):
#     """
#     get scale height
#     """    
#     GasScaleHeight=cs/OmegaK
#     return GasScaleHeight

# def viscosity(cs,Hg):
#     """
#     get viscosity
#     """
#     nu=alpha*cs*Hg
#     return nu

# def v_th(cs):
#     """
#     get the thermal velocity 
#     """
#     thermal_velocity=np.sqrt(8/np.pi)*cs
#     return thermal_velocity

# def l_mfp(rhog,mg):
#     """
#     get the mean free path

#     """
#     MeanFreePath=mg/sigmamol/rhog
#     return MeanFreePath

# def rho_g(Sigmag,Hg):
#     """
#     get the density of gas
#     """
#     rhogas=Sigmag/(2*np.pi)**0.5/Hg
#     return rhogas


def eta_old(r,Mcp,dotMg,mg):
    """
    v=(1-eta)vK
    """
    e=0.477577*cgs.kB*(r**5.5-4.84615*r**4.5*rout)/cgs.gC/Mcp/mg/(r**4.5-3*r**3.5*rout)*(cgs.gC*Mcp*dotMg/r**3/cgs.sigmaSB)**(1/4)
    return e

def eta (disk):
    """
    as in: vgas = (1-eta)vK

    This is a derived function, using parameters from the disk class
    """
    r = disk.loc
    #Mcp = disk.mcp
    #dotMg = dot_Mg(disk.time)
    Hg = disk.Hg
    # from Shibaika. 2019 eq13.
    # this involves the gradient of the surface density as defined above
    # which I calculated with mathematica
    e = -1/2 *(Hg/r)**2 *beta_P(r)
    #eo=0.477577*cgs.kB*(r**5.5-4.84615*r**4.5*rout)/cgs.gC/Mcp/mg/(r**4.5-3*r**3.5*rout)*(cgs.gC*Mcp*dotMg/r**3/cgs.sigmaSB)**(1/4)
    return e


