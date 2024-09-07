import cgs 
import numpy as np 
import parameters as pars 

rinn = pars.dgasgrid['rinn']
rout = pars.dgasgrid['rout']

sigmol = 2e-19

def Mcp_t(t):
    #seems they didn't consider the mass evolution of central star
    return cgs.Msun

def dot_Mg(t):
    return 0 

def M_influx(t0, tEnd):
    return 0.0

def user_add_var():
    #let's now see the midium viscous disk
    return {'alpha':1e-3}

def user_add_fun():
    return [dot_Md]

def user_add_eval():
    return [eta, mu]

def key_disk_properties(rad, t, dold=None):
    if dold is None:
        #if there's no old porperties, initiallize them.
        #xi and gamma is defined artificially from Yap.2024
        phi = 0.1 #the ratio alphadw/alphanu 
        gamma = 0.55 #from Yap.2024 or 0.69 for phi=10
        lamb = 3.5
        xi = 1/4* (phi+1)*(np.sqrt(1+4*phi/(lamb-1)/(phi+1)**2)-1) #Yap.2024 eq.6 
        sigmaG0 = 0.05*cgs.Msun/2/np.pi/self.rout**2 #ignore the beta for now
        sigmaG = sigmaG0*(rad/rout)**(xi-gamma)*np.exp(-(rad/rout)**(2-gamma))

        fd = 0.01 #dg ratio 
        kappa_d = 30 #opacity
        alphanu = 1e-3/1.1
        
        OmegaK = np.sqrt(cgs.gC*Mcp_t(0.0)/rad**3)
        mu = 2.34*np.ones(len(rad))
        temp = (27*fd*kappa_d*alphanu*cgs.kB*OmegaK*sigmaG**2/64/cgs.sigmaSB/mu/cgs.mp)**(1/3)

        #plt.figure()
        #plt.title('Surface density')
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.plot(self.loc/cgs.au, self.sigmaG*10)#convert the unit
        #plt.show()
    else: 
        #use the continuity equation to update the grid




        import pdb;pdb.set_trace()

    return sigmaG, temp, mu


def mu(disk):
    return 2.34*np.ones_like(disk.loc)

def eta(disk):
    return 0.01*np.ones_like(disk.loc)
