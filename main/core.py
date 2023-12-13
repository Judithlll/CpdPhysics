import numpy as np
# from scipy.integrate import odeint
import sys
import cgs
import ode
import matplotlib.pyplot as plt
import disk_properties as dp 
import planets_properties as pp
import functions as f
import copy
import parameters as pars
from gas import GAS

class System(object):

    ##CWO: are these parameters still used?
    ##they should anyway not be defined here

    rout=27*cgs.RJ
    #alpha=1e-4 # TBR
    fraction=0.02 
    dgratio=0.01  #dust to gas ratio
    #sigmamol=2e-15 #TBR
    rhoint = 1.4
    #deltaT=0.0
    mtot1 = 1e24 #total mass a single superparticle represents
    daction={}
    timestepn=10  #how many time points in every ODE solution process
    rhoPlanet=1.9

    #CWO: please work with default pars
    def __init__(self,Rdi=1.0,time=0.0,nini=10):

        #initialize parameter from txt file // disk.py


        self.Rdi=Rdi  #initial radius of particles
        self.nini=nini
        self.mini = 4*np.pi/3*self.rhoint *Rdi**3
        self.time=time  #initial time
        self.ntime = 0

        # define a disk class
        self.gas = self.init_gas ()

        #out = self.gas.get_key_disk_properties (np.array([4*cgs.rJup, 8*cgs.rJup]), time)
        #self.disk = Disk ()



        # define class of superparticles here
        #self.particles = Superparticles(nini,self.mini,self.disk.rinn,self.disk.rout,self.mtot1)

        self.particles = Superparticles(nini,self.mini,dp.rinn,dp.rout,self.mtot1)

        #the amount of solid mass that has crossed into the domain
        self.Minflux = 0
        #self.Mcp=self.disk.Mcp_t(self.time)


    def init_gas (self, gasL=None, dcomposL=None, dgrid={}):
        """
        init_gas is called at initialization. It adds an instance of the GAS class

        history
        [23.12.13]:this is copied from /NewLagrange code base...
                  :for the moment gasL and dcomposL are put to None

        """
        #also initializes a gas class (nested class)
        if pars.gasmodel=='grid':
            dgas = pars.dgasgrid
        else:
            dgas = {}

        dum = GAS (gasL, dcomposL, mode=pars.gasmodel, time=self.time, **dgas)

        return dum

    def update_particles (self,tEnd):
        """
        How to evolving the gas in the disk is still not sure

        Evolving system to the self.time

        """

        #time derivative of particles
        dydtP = self.particles.dY2d_dt(self.particles.Y2d,self.time,self.gas)
            
        #timescale
        tscaleArr = np.abs(self.particles.Y2d/dydtP)
        deltaT = np.nanmin(0.2*tscaleArr)
        # import pdb; pdb.set_trace()

        if self.time+deltaT>tEnd:
            deltaT = tEnd - self.time

        #update particle properties
        #make
        self.Y2dold = copy.deepcopy(self.particles.Y2d)
        Yt = self.particles.update(self.time,self.time+deltaT,self.gas,self.timestepn)

        self.deltaT = deltaT
        self.ntime += 1
        # stop here
        return Yt, deltaT
    
    def P_eff(self,Yt):
        """
        this function aims to get the pebble accretion rate of planets

        When the particles drift into the scale of the planets gravity, part of them will be accreted, but the number of the particles is not enough now, so we can only 
        1) add the initial number of particles and decrease the total mass of particles, but the code need very long time to run 
        2) change the total mass according to the accretion rate
        """
        self.Peff={}
        self.idx_Pars={}
        for i in range(len(Yt)):
            
            self.planets[:].time=self.time+self.deltaT/self.timestepn *i
            
            for j in range(len(self.planets)):
                self.idx_Pars[self.planets[j].time]=self.planets[j].get_effective_pebbles_index()
                self.Peff[self.planets[j].time]=self.planets[j].pebble_accretion_rate()

    
    def post_process (self):
        """
        returns indices of particles that cross boundaries 
        (which needs to be removed or added)
        """

        self.daction = {}

        loc = self.particles.Y2d[0]

        #particles that cross the inner disk edge
        idx, = (loc<dp.rinn).nonzero()
        if len(idx)>0:
            self.daction['remove'] = idx
        
        #particles that enter the domain
        Nadd = 0
        self.Minflux += dp.M_influx(self.time,self.time+self.deltaT)

        while self.Minflux>self.mtot1:
            Nadd += 1
            self.Minflux -= self.mtot1
        
        if Nadd>0:
            self.daction['add'] = Nadd

        
        #post_process particles
        if 'remove' in self.daction.keys():
            #remove the particles from Y2d!
            
            self.particles.remove_particles(self.daction['remove'])
            self.nini-=len(self.daction['remove'])
            # import pdb; pdb.set_trace()

        if 'add' in self.daction.keys():
            self.particles.add_particles(self.daction['add'])
            self.nini+=self.daction['add']

        #particles that are eaten by the planet
        # need to use pebble accretion rate 
        #....

        #particles that have become too big
        #....


def advance_planets (system):
    """
    [23.12.06]copied/edited from NewLagrange

    For now, planet migration and composition is not considered
    """
    for planet in system.planetL:

        #planet exists only after planet.time
        if planet.time<system.time:

            sploc = system.particles.Y2d[0]
            sploc_old = system.Y2dold[0]

            #particles that cross are those that
            #idx, = np.nonzero( (planet.loc<sp.loc) & (planet.loc>spN.loc) )

            idx, = np.nonzero( (planet.loc<sploc_old) & (planet.loc>sploc) )


            iterate = True
            niter = 0
            while iterate:

                crossL = []
                for ip in idx:
                    spi = system.Y2dold[:,ip] #makes a superparticle
                    crossL.append(spi)

                
                # if len(idx)!=0:
                #     import pdb; pdb.set_trace()
                #this user-defined function should detail how the
                #planet properties (location,mass,composition) change
                #with time and how the crossed particles are affected... 


                #this is about planet migration...
                #TBD later...
                if False:
                    loc_t, mass_t, fcomp_t = \
                    userfun.XY_planet (sim.time, planet.loc, planet.mass, planet.fcomp, 
                            crossL)
                else:
                    loc_t = -pp.planet_migration(system.gas,system.planetL[0].loc,system.planetL[0].mass,system.time,system.rhoPlanet)     #migration rate of planet
                    mass_t = 0.0    #gas accretion of planet
                    fcomp_t = 0.0   #how its composition changes

                #update planet properties from rates supplied by user
                planet_loc_nw = planet.loc + loc_t *system.deltaT

                planet_loc_nw = planet.loc


                #particles that cross are those that
                idxN, = np.nonzero( (planet.loc<sploc_old) & (sploc<planet_loc_nw) )


                if set(idxN)!=set(idx):
                    idx = idxN
                    niter += 1
                else:
                    iterate = False


            #update planet properties from rates supplied by user
            planet.loc += loc_t *system.deltaT
            planet.mass += mass_t *system.deltaT
            planet.fcomp += fcomp_t *system.deltaT


            #update s-particle properties from sp-crossings
            #assumes all particles are accreted (TBC!!)
            spN=system.particles

            for k, ip in enumerate(idxN):

                

                #mass that is being transferred (TBC!!)
                #need to calculate epsilon (PA efficiency)

                #calculate critical mass to verify if the pebble accretion can occur
                Mc=f.M_critical(system,planet.loc,crossL[k])
                if Mc<planet.mass:                    
                    epsilon = f.epsilon_PA(system,planet.loc,planet.mass,crossL[k])
                    delm = epsilon*crossL[k][2]#don't understand this line...
                    
                else:
                    "pebble accretion can not happen"
                    delm=0

                planet.mass += delm #increase mass (pebble accretion)
                # planet.fcomp += 0.  #TBD !!

                #spN -> system.particles.Y2d...
                spN.mtotL[ip] -= delm #decrease mass sp


#perhaps this class object is not necessary...
class DISK (object):
    """
    ther
    """

    def __init__ (self, sigmaG, temp, mu, loc, time):
        """
        initialize with key disk properties
        """
        self.loc = loc
        self.Sigmag = sigmaG
        self.temp = temp
        self.mu = mu
        self.time = time
        self.alpha = dp.alpha
        self.rout = dp.rout
        self.rinn = dp.rinn
        self.tgap=dp.tgap


    def add_auxiliary (self):
        self.Mcp = self.Mcp_t(self.time)  
        self.OmegaK=self.Omega_K(self.loc,self.time)       
        self.dotMg=self.dot_Mg(self.time)
        
        self.cs =  dp.c_s(self.temp) #dp.c_s(loc,time)
        self.vth = dp.v_th(self.cs)

        self.dotMd=self.dot_Md()
        self.vK=self.v_K(self.loc)
        self.vth=self.v_th()
        self.Hg=self.H_g()
        self.nu=self.viscosity()
        self.rhog=self.rho_g()
        self.mg=self.m_g()
        self.lmfp=self.l_mfp()
        self.eta=self.eta_cal(self.loc)


    def Mcp_t(self,time):
        return dp.Mcp_t(time)
    
    def dot_Mg(self,time):
        return dp.dot_Mg(time)
    
    def dot_Md(self):
        return dp.dot_Md(self.dotMg)
    
    def M_influx(self,t0,tEnd):
        return dp.M_influx(t0,tEnd)

    # def Sigmag (self,loc,time):
    #     c_s=self.cs(loc,time)
    #     Sg=dp.Sigma_g(loc,time,c_s)
    #     return Sg
    
    def Omega_K(self,loc,time):
        return dp.Omega_K(loc,time,self.Mcp)

    def v_K(self,loc):
        return self.OmegaK*loc
    
    def v_th(self):
        return dp.v_th(self.cs)
    
    def H_g(self):
        return dp.H_g(self.cs,self.OmegaK)
    
    def viscosity(self):
        return dp.viscosity(self.cs,self.Hg)
    
    def rho_g(self):
        return dp.rho_g(self.Sigmag,self.Hg)

    def m_g(self):
        return dp.m_g()

    def l_mfp(self):
        return dp.l_mfp(self.rhog,self.mg)
    
    def eta_cal(self,loc):
        return dp.eta(loc,self.Mcp,self.dotMg,self.mg)


class Superparticles(object):
    rhoint=1.4
    pi=np.pi
    Sto=0.0001  
    error=1e-8


    def __init__(self,nini,mini,rinn,rout,mtot1):
        """
        systems initial properties

        nini: initial number of the particles
        mini: initial mass for every particles
        rinn: inner edge of the disk
        rout: outer edge of the disk
        """
        self.nini=nini
        self.mini=mini #physical mass
        self.rinn=rinn
        self.rout=rout
        self.mtot1 = mtot1

        ndim = 3# location // mass // total mass

        self.locL=10**np.linspace(np.log10(rinn),np.log10(rout),nini)
        self.massL=[mini]*nini
        self.mtotL=[mtot1]*nini

        self.Y2d = np.empty((ndim,nini))
        self.Y2d[0] = self.locL
        self.Y2d[1] = self.massL
        self.Y2d[2] = self.mtotL


    def dY2d_dt (self,Y2d,t,gas):
        """
        input:
            Y2d -- state vector
            time -- time
            disk -- disk object
        """

        #unpack the state vector
        r, mphy, mtot = self.Y2d   #maybe the total mass needn't to be put in Y2d

        Rd=(mphy/(self.rhoint*4/3*np.pi))**(1/3)

        #get disk object instance from gas
        #maybe like this???

        out = gas.get_key_disk_properties (r, t)
        disk = DISK (*out, r, t)
        disk.add_auxiliary ()

        eta=disk.eta
        v_K=disk.vK
        v_th=disk.vth
        lmfp=disk.lmfp
        rho_g=disk.rhog
        Omega_K=disk.OmegaK
        H_g=disk.Hg
        dotMd=disk.dotMd


        #obtain Stokes number by iterating on drag law
        St,v_r = f.St_iterate(eta,v_K,v_th,lmfp,rho_g,Omega_K,Rd)

        #describe what's going on here
        v_dd = np.abs(v_r)/2
        H_d = H_g*(1+St/disk.alpha*(1+2*St)/(1+St))**(-0.5)

        #this becomes...
        H_d = disk.Hg *(1+St/disk.alpha*(1+2*St)/(1+St))**(-0.5)
            

        drdt = v_r
        #dR_ddt= v_dd*dot_M_d/4/pi**(3/2)/rho_int/H_d/r/v_r**2 *dr_dt

        sigD = dotMd /(-2*r*self.pi*v_r) #v_r<0
        
        dmdt=2*np.sqrt(np.pi)*Rd**2*v_dd/H_d*sigD  #eq. 5 of Shibaike et al. 2017

        Y2ddt = np.zeros_like(self.Y2d)
        Y2ddt[0] = drdt
        Y2ddt[1] = dmdt
        Y2ddt[2] = 0.0

        return Y2ddt 
    
    
    def update(self,t0,tFi,gas,nstep=10):
        """
        this integrate the particles until tFi
        -- d: disk object
        """

        tSpan=np.array([t0,tFi])
        tstep=(tFi-t0)/nstep #why 100?
    
        Y2copy = np.copy(self.Y2d)

        #integrates system to tFi
        Yt = ode.ode(self.dY2d_dt,Y2copy,tSpan,tstep,'RK5',gas)
        
        self.locL=Yt[-1,0,:]
        self.massL=Yt[-1,1,:]
        self.mtotL=Yt[-1,2,:]
        
        self.Y2d = Yt[-1,:,:]
        return Yt

    
    def remove_particles(self,remove_idx):
        self.Y2d = np.delete(self.Y2d, remove_idx, 1)     

    def add_particles(self,Nadd):

        ynew = np.array([[self.rout],[self.mini],[self.mtot1]])

        self.Y2d = np.append(self.Y2d, ynew,1)

        if Nadd!=1:
            print('[core]Error:can only add 1 particle // reduce timestep')
            sys.exit()

    
        #for moment put the implemented particles all at the outmost edge
        #lociL=np.linspace(self.rout,self.rout,add_number)
        #miL=np.linspace(self.mini,self.mini,add_number)
        #mtot=np.nanmean(self.Y2d[2])+add_number*self.mini
        #mtotL=np.linspace(mtot,mtot,add_number)
        #self.Y2d=np.append(self.Y2d,np.array([lociL,miL,mtotL]),1)
        #self.Y2d[-1]=mtot

        return


class PLANET ():
    """
    Copied over from /NewLagrange
    """
    def __init__(self, time, loc, mplanet, fcomp):
        self.loc = loc          #location
        self.time = time        #time when it appears
        self.mass = mplanet     #its mass
        self.fcomp = fcomp      #its composition
        self.spCrossTime = [0.0]   #list when particles cross
        self.spCrossMass = [0.0]   #corresponding mass
        self.spCrossTau = [-1]   #corresponding stokes number
        self.ncross = 0

    def record_cross_sp (self, time, sp, idxcross):
        for k in idxcross:
            self.spCrossTime.append(time)
            self.spCrossMass.append(sp.msup[k])
            self.spCrossTau.append(sp.tau[k])
            self.ncross += 1

    def calc_mdot (self, time, Nast=15):
        tlba = time -np.array(self.spCrossTime)[::-1]   #lookback time
        tast = tlba[:Nast].max()                            #characteristic timescale
        mdotarr = np.array(self.spCrossMass)[::-1] /tast
        wi = np.exp(-tlba/tast)     #weights
        mdot = np.sum(mdotarr *wi)  #mass flux through iceline
        return mdot


