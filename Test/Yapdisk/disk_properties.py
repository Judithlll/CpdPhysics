import cgs 
from scipy.integrate import quad
from scipy.linalg import solve
from scipy.linalg import solve_banded
import matplotlib.pyplot as plt
import numpy as np 
import parameters as pars 
from scipy.optimize import root
from scipy.interpolate import interp1d

rinn = pars.dgasgrid['rinn']
#rout = pars.dgasgrid['rout']

alpha = 1e-3 
phi = 0.1
alphanu = alpha/(1+phi)
alphadw = alpha - alphanu
sigmol = 2e-19
gamma = 0.55 #from Yap.2024 or 0.69 for phi=10
#gamma = 0.69

#xi and gamma is defined artificially from Yap.2024
lamb = 3.5

xi = 1/4* (phi+1)*(np.sqrt(1+4*phi/(lamb-1)/(phi+1)**2)-1) #Yap.2024 eq.6 
fd = 0.01 #dg ratio 
kappa_d = 30 #opacity

def r_out(t,tacc0 = 0.3e6*cgs.yr):
    #eq 8 in Yap.2024 
    rout0 = pars.dgasgrid['rout']
    rout = rout0*(1+(t/(1+phi)/tacc0))**(1/(2-gamma))
    
    return rout

def Mcp_t(t):
    #seems they didn't consider the mass evolution of central star
    return cgs.Msun

def dot_Mg(t):
    return 0 

def M_influx(t0, tEnd):
    return 0.0

def user_add_var():
    #let's now see the midium viscous disk
    return {'alpha':alpha}

def user_add_fun():
    return [dot_Md]

def user_add_eval():
    return [eta, mu]

def temp_eq (x, p0, p1, p2, p3, p4):
    return p0 + p1*x + p2*x**2 + p3*x**3 + p4*x**4

def key_disk_properties(rad, t, dold=None, mode='numerical'):
    
    # globaly, if we do this numerically, we need to store the values in ghost cells also
    OmegaK = Omega_K(rad, Mcp_t(t))
    mu = 2.4*np.ones(len(rad))

    #get the total mass of disk and check the mass conservation 
    M_disk0 = 0.05*cgs.Msun 

    if t==0.:
        #if there's no old porperties, initiallize them.
        rout = r_out(0.)
        sigmaGc0 = M_disk0/2/np.pi/rout**2#ignore the beta for now
        sigmaG = sigmaGc0*(rad/rout)**(xi-gamma)*np.exp(-(rad/rout)**(2-gamma))

        #check the total mass of disk 
        dellogr = np.log10(rad[1]/rad[0])
        delr = (10**dellogr-1)*rad 

        mtot = np.sum(2*np.pi*rad*sigmaG*delr)

        global beta
        beta = 1/(mtot/0.05/cgs.Msun)

        #fix the initial sigmaG
        sigmaGc0 = beta*M_disk0/2/np.pi/rout**2
        sigmaG = sigmaGc0*(rad/rout)**(xi-gamma)*np.exp(-(rad/rout)**(2-gamma))

        # plt.figure()
        # plt.title('Surface density')
        # plt.xscale('log')
        # plt.yscale('log')
        # plt.plot(rad/cgs.au, sigmaG, 'x-')#convert the unit
        # plt.show()
        # import pdb;pdb.set_trace()
    elif mode == 'numerical': 
        ngrid = len(rad)
        dellogr = np.log10(rad[1]/rad[0]) 
        cs = np.sqrt(cgs.kB*dold['temp']/dold['muout']/cgs.mp)

        y = dold['sigmaG']*rad**(3/2)
        x = rad**(0.5)
        nu = np.mean(alphanu *cs**2 /OmegaK )
        mux = 3*nu/4/rad 


        #add the ghost cells for mu, y, and x 
        def add_ghost(arr):
            arr = np.append(0,arr)
            arr = np.append(arr, 0)
            return arr 

        x = add_ghost(x) 
        x[0] = np.sqrt(rad[0]*10**(-dellogr))
        x[-1] = np.sqrt(rad[-1]*10**(dellogr))
        y = add_ghost(y)
        mux = add_ghost(mux)

        delx = x*(10**(dellogr/2)-1)

        #set the boundary values
        y[0] = 0. 
        y[-1] = 0. 
        
        #extrapolate the ghost cells b/c the value of mux in ghost cells seems important
        f = interp1d (rad, mux[1:-1], kind='slinear', fill_value='extrapolate') 
        mux[0] = f(x[0]**2)
        mux[-1] = f(x[-1]**2)

        #the matrix for the term 'y_new' with ghost cells 
        E = np.eye(ngrid+2) 
    
        #the matrix for 2nd order derivative: 'd^2(mu y)/dx^2' with ghost cells
        Adiff = np.zeros((ngrid+2,ngrid+2))
        np.fill_diagonal(Adiff, -2)
        np.fill_diagonal(Adiff[1:], 1)
        np.fill_diagonal(Adiff[:,1:], 1)
        Adiff*=mux 
        Adiff*=1/delx**2 

        #try 0 influx first 
        afin = E - (t-dold['time'])*Adiff 

        #set the 0-flux boundary condition 
        afin[0,:3] = np.array([-2,0,2])*mux[:3]
        afin[-1,-3:] = np.array([2,0,-2])*mux[-3:] 

        #get the banded matrix 
        afin_band = np.zeros((5,ngrid+2))
        afin_band[0, 2:] = np.diag(afin,2)
        afin_band[1, 1:] = np.diag(afin,1)
        afin_band[2] = np.diag(afin)
        afin_band[3,:-1] = np.diag(afin,-1)
        afin_band[4,:-2] = np.diag(afin,-2)

        y_new = solve_banded((2,2), afin_band, y)

        sigmaG = y_new[1:-1]/x[1:-1]**3

        #if mtot != M_disk0/cgs.Msun:
        #import pdb;pdb.set_trace()


        try: 
            assert(np.all(sigmaG>=0))
        except:
            import pdb;pdb.set_trace()

    elif mode =='analytical':
        sigmaGc0 = beta*0.05*cgs.Msun/2/np.pi/r_out(0.)**2#ignore the beta for now
        cs = np.sqrt(cgs.kB*dold['temp']/dold['muout']/cgs.mp)
        tacc0 = r_out(0.)**2*OmegaK[-1]/3/(2-gamma)**2/cs[-1]**2/alpha 

        #import pdb;pdb.set_trace()

        sigmaGc = sigmaGc0 *(1+t/(1+phi)/tacc0)**((5+2*xi+phi)/2/(2-gamma))
        rout = r_out(t, tacc0)

        sigmaG = sigmaGc*(rad/rout)**(xi-gamma)*np.exp(-(rad/rout)**(2-gamma))

        print(tacc0/cgs.yr, dold['temp'][-1],sigmaGc, (5+2*xi+phi)/2/(2-gamma))
    else:
        raise ValueError('mode not recognized')
        

    #check the mass conservation 
    delr = (10**dellogr-1)*rad 
    mtot = np.sum(2*np.pi*rad*sigmaG*delr)/cgs.Msun

    print(mtot)
    #get the temperature 
    # p0 = -(7.3e27)*rad**(-12/7) 
    p1 = -27*fd*kappa_d*alphanu*cgs.kB*OmegaK*sigmaG**2/64/cgs.sigmaSB/mu/cgs.mp
    # p2 = 0. 
    # p3 = 0. 
    # p4 = 1. 
    Tvsic = (-p1)**(1/3)
    Tirr  = 150*(rad/cgs.au)**(-3/7)
    #temp = root(temp_eq, args=(p0,p1,p2,p3,p4), x0 = (-p1)**(1/3)).x
    temp = (Tirr**4 + Tvsic**4)**(1/4)

    return sigmaG, temp, mu


def Omega_K(loc, mcp):
    return np.sqrt(cgs.gC*mcp/loc**3)
def mu(disk):
    return 2.34*np.ones_like(disk.loc)


def dot_Md(time):
    return 0.01*dot_Mg(time)

def eta(disk):
    return 0.01*np.ones_like(disk.loc)

