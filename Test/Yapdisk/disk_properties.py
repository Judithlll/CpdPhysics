import cgs 
import matplotlib.pyplot as plt
import numpy as np 
import parameters as pars 

rinn = pars.dgasgrid['rinn']
#rout = pars.dgasgrid['rout']

alphanu = 1e-3/1.1
alphadw = 1e-3/11
phi=alphadw/alphanu
alpha = alphanu+alphadw
sigmol = 2e-19
gamma = 0.55 #from Yap.2024 or 0.69 for phi=10

def rout(t):
    #eq 8 in Yap.2024 
    rout0 = pars.dgasgrid['rout']
    tacc0 = 0.5e6*cgs.yr 
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


def key_disk_properties(rad, t, dold=None):
    #xi and gamma is defined artificially from Yap.2024
    lamb = 3.5

    xi = 1/4* (phi+1)*(np.sqrt(1+4*phi/(lamb-1)/(phi+1)**2)-1) #Yap.2024 eq.6 
    fd = 0.01 #dg ratio 
    kappa_d = 30 #opacity
    
    OmegaK = Omega_K(rad, Mcp_t(t))
    mu = 2.34*np.ones(len(rad))

    if t==0.:
        #if there's no old porperties, initiallize them.
        sigmaG0 = 0.05*cgs.Msun/2/np.pi/rout(t)**2 #ignore the beta for now
        sigmaG = sigmaG0*(rad/rout(t))**(xi-gamma)*np.exp(-(rad/rout(t))**(2-gamma))

        temp = (27*fd*kappa_d*alphanu*cgs.kB*OmegaK*sigmaG**2/64/cgs.sigmaSB/mu/cgs.mp)**(1/3)

        #plt.figure()
        #plt.title('Surface density')
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.plot(self.loc/cgs.au, self.sigmaG*10)#convert the unit
        #plt.show()
    else: 
        #use the continuity equation to update the grid
        #eq. 4 in Yap.2024
        cs = np.sqrt(cgs.kB*dold['temp']/dold['muout']/cgs.mp)
        nome1 = rad**2 * alphanu* dold['sigmaG']*cs**2 

        #so the boundary condition are:
        #1) at the out edge, the surface density is zero
        #2) at the inner edge, the dsigmadr is zero
        #nome1 = np.append(nome1, 0.0) 
        #rad = np.append(rad, rout(t))
        dnome1dr = (nome1[1:]-nome1[:-1])/(rad[1:]-rad[:-1])


        nome2 = dnome1dr/rad[:-1]/OmegaK[:-1]
        #nome2 = np.append(nome2[0],nome2)
        #rad = np.append(rinn-0.0024*cgs.au, rad)

        dnome2dr = (nome2[1:]-nome2[:-1])/(rad[2:]-rad[1:-1])
       
        sigmaG = dold['sigmaG'][1:-1] + (3/rad[1:-1]*dnome2dr)*(t-dold['time'])

        sigmaG_boundary = np.interp([rinn,rout(t)],rad[1:-1],sigmaG)
        sigmaG = np.pad(sigmaG, (1,1), 'constant', constant_values=sigmaG_boundary)
        temp = (27*fd*kappa_d*alphanu*cgs.kB*OmegaK*sigmaG**2/64/cgs.sigmaSB/mu/cgs.mp)**(1/3)

        try: 
            assert(np.all(sigmaG>=0))
        except:
            import pdb;pdb.set_trace()
        
    
    


    return sigmaG, temp, mu


def Omega_K(loc, mcp):
    return np.sqrt(cgs.gC*mcp/loc**3)
def mu(disk):
    return 2.34*np.ones_like(disk.loc)


def dot_Md(time):
    return 0.01*dot_Mg(time)

def eta(disk):
    return 0.01*np.ones_like(disk.loc)

