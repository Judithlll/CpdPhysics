import cgs 
from scipy.integrate import quad
from scipy.linalg import solve
import matplotlib.pyplot as plt
import numpy as np 
import parameters as pars 

rinn = pars.dgasgrid['rinn']
#rout = pars.dgasgrid['rout']

alpha = 1e-3 
phi = 0.1
alphanu = alpha/(1+phi)
alphadw = alpha - alphanu
sigmol = 2e-19
#gamma = 0.55 #from Yap.2024 or 0.69 for phi=10
gamma = 0.69

def r_out(t,tacc0 = 0.5e6*cgs.yr):
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


def key_disk_properties(rad, t, dold=None, mode='numerical'):
    #xi and gamma is defined artificially from Yap.2024
    lamb = 3.5

    xi = 1/4* (phi+1)*(np.sqrt(1+4*phi/(lamb-1)/(phi+1)**2)-1) #Yap.2024 eq.6 
    fd = 0.01 #dg ratio 
    kappa_d = 30 #opacity
    
    OmegaK = Omega_K(rad, Mcp_t(t))
    mu = 2.34*np.ones(len(rad))

    if t==0.:
        #if there's no old porperties, initiallize them.
        rout = r_out(0.)
        sigmaGc0 = 0.05*cgs.Msun/2/np.pi/rout**2#ignore the beta for now
        sigmaG = sigmaGc0*(rad/rout)**(xi-gamma)*np.exp(-(rad/rout)**(2-gamma))

        #check the total mass of disk 
        dellogr = np.log10(rad[1]/rad[0])
        delr = (10**dellogr-1)*rad 

        mtot = np.sum(2*np.pi*rad*sigmaG*delr)

        global beta
        beta = 1/(mtot/0.05/cgs.Msun)

        #fix the initial sigmaG
        sigmaGc0 = beta*0.05*cgs.Msun/2/np.pi/rout**2#ignore the beta for now
        sigmaG = sigmaGc0*(rad/rout)**(xi-gamma)*np.exp(-(rad/rout)**(2-gamma))

        #plt.figure()
        #plt.title('Surface density')
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.plot(self.loc/cgs.au, self.sigmaG*10)#convert the unit
        #plt.show()
    elif mode == 'numerical': 
        #use the continuity equation to update the grid
        #eq. 4 in Yap.2024
        cs = np.sqrt(cgs.kB*dold['temp']/dold['muout']/cgs.mp)

        #equal in log space so get the del r like this
        dellogr = np.log10(rad[1]/rad[0])
        delr = (10**dellogr-1)*rad

        x = np.sqrt(rad)
        nu = alphanu*cs**2/OmegaK 
        u = 3*nu/4/rad 

        xout = x[-1]*np.sqrt(10**dellogr)
        delx = np.diff(np.append(x,xout))

        #build the matrix
        A = np.zeros((len(rad), len(rad)))
        np.fill_diagonal(A, -2*u/delx**2)
        np.fill_diagonal(A[1:], u[:-1]/delx[1:]**2)
        np.fill_diagonal(A[:,1:], u[1:]/delx[:-1]**2)

        #then modify the boundary
        A[0,0] = -u[0]/delx[0]**2 
        A[0,1] = u[1]/delx[0]**2 

        #the matrix of disk-wide term 
        z = u*phi/x  
        Adw = np.zeros((len(rad), len(rad)))
        np.fill_diagonal(Adw, -z/delx)
        np.fill_diagonal(Adw[:,1:], z[1:]/delx[:-1])
        
        #the mass loss term b/c of the disk wind.
        #eq. 3 in Yap.2024
        Aml = np.eye(len(rad))*u*phi/x**2

        E = np.eye(len(rad))

        afin = E - (t-dold['time'])*(A+Adw - Aml)

        sigmaG = solve(afin, dold['sigmaG'])



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
        
    temp = (27*fd*kappa_d*alphanu*cgs.kB*OmegaK*sigmaG**2/64/cgs.sigmaSB/mu/cgs.mp)**(1/3)

        
    
    return sigmaG, temp, mu


def Omega_K(loc, mcp):
    return np.sqrt(cgs.gC*mcp/loc**3)
def mu(disk):
    return 2.34*np.ones_like(disk.loc)


def dot_Md(time):
    return 0.01*dot_Mg(time)

def eta(disk):
    return 0.01*np.ones_like(disk.loc)

