import numpy as np
# from scipy.integrate import odeint
import sys
import cgs
import ode
import matplotlib.pyplot as plt
import disk_properties as dp 
import functions as ff
import copy
import parameters as pars
from gas import GAS
import scipy.optimize as sciop
import scipy.integrate as sciint
import time

class System(object):

    """
    SYSTEM: integrate the every part of the PPD (or CPD) and update them with time.
    """

    fraction=0.02 
    dgratio=0.01  #dust to gas ratio
    rhoint = 1.4
    mtot1 = 1e24 #total mass a single superparticle represents
    daction={}
    timestepn=3  #how many time points in every ODE solution process
    rhoPlanet=1.9
    tgap=dp.tgap

    def __init__(self,Rdi=0.01,time=0.0,nini=10,diskmass=0.01*cgs.MJ):
        
        #initialize parameter from txt file // disk.py
        self.Rdi=Rdi  #initial radius of particles
        self.nini=nini
        self.mini = 4*np.pi/3*self.rhoint *Rdi**3
        self.time=time  #initial time
        self.ntime = 0

        # define a disk class
        self.gas = self.init_gas ()

        #TBD: put this in dparticleprops
        self.diskmass=diskmass

        # define class of superparticles here
        #self.particles = Superparticles(nini,self.mini,self.disk.rinn,self.disk.rout,self.mtot1)
        #the amount of solid mass that has crossed into the domain
        self.Minflux = 0
        #self.Mcp=self.disk.Mcp_t(self.time)
    
    def init_particles(self, dparticleprops={}):
        """
        because we need to consider iceline, so separatly 
        initiallize the particles, for now just water iceline is considered  

        [23.12.30]CWO:instead of forwarding diskmass, I supply self.gas to the superparticle class
        """
        self.particles = Superparticles(dp.rinn,dp.rout,self.dcomposL,self.gas, **dparticleprops)


    def init_gas (self, gasL=None, dcomposL=None, dgrid={}):
        """
        init_gas is called at initialization. It adds an instance of the GAS class

        history
        [23.12.13]:this is copied from /NewLagrange code base...
                  :for the moment gasL and dcomposL are put to None

        """
        #also initializes a gas class (nested class)
        if pars.gasmodel=='gridstatic':
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

        #TBD: make an Y2d array...
        
        dydtP = self.particles.dY2d_dt(self.particles.Y2d,self.time,self.gas)

        #timescale
        tscaleArr = np.abs(self.particles.Y2d/dydtP)
        deltaT = np.nanmin(0.2*tscaleArr)
        # import pdb; pdb.set_trace()

        if self.time+deltaT>tEnd:
            deltaT = tEnd - self.time

        #update particle properties
        #make
        
        self.back_up_last_data()
        Yt = self.particles.update(self.time,self.time+deltaT,self.gas,self.timestepn)

        self.deltaT = deltaT
        self.ntime += 1
        # print(self.particles.locL, self.time)
        if self.time==np.nan or self.deltaT==np.nan:
            print('hello')
            import pdb;pdb.set_trace()
        return Yt, deltaT
    
    def back_up_last_data(self):
        self.locLold=copy.deepcopy(self.particles.locL)
        self.massLold=copy.deepcopy(self.particles.massL)
        self.mtotLold=copy.deepcopy(self.particles.mtotL)

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



def advance_iceline (system):
    """
    for now particles directly lose the mass of water without any other effect
    """

    sploc = system.particles.locL
    sploc_old = system.locLold
    for k,iceline in enumerate(system.icelineL):
        idx,=np.nonzero((iceline.loc<sploc_old) & (iceline.loc>sploc))

        ##[23.12.30]CWO: I have no idea what this does... can be removed, I think
        #frac=1
        #for i in range(k+1):
        #    frac-=system.icelineL[i].frac

        ic = pars.composL.index(iceline.species) #refers to species index
        if len(idx)!=0:     
            # for ix in idx:
            #     #system.particles.msup =


            fice = system.particles.fcomp[idx,ic]  #mass fraction in ice
            fremain = (1-fice)          #remain fraction
            fremain[fremain<1e-15] = 0  #loss of numbers (!!)
            system.particles.mtotL[idx] *= fremain    #reduce masses accordingly
            system.particles.massL[idx] *= fremain
            system.particles.fcomp[idx,ic] = 0.      #gone is the ice!

            #renormalize
            system.particles.fcomp[idx,:] = (system.particles.fcomp[idx,:].T /(system.particles.fcomp[idx,:].sum(1)+1e-100)).T



def advance_planets (system):
    import planets_properties as pp
    """
    [23.12.06]copied/edited from NewLagrange

    For now, planet migration and composition is not considered
    """
    for planet in system.planetL:

        #planet exists only after planet.time
        if planet.time<system.time:

            sploc = system.particles.locL
            sploc_old = system.locLold

            #particles that cross are those that
            idx, = np.nonzero( (planet.loc<sploc_old) & (planet.loc>sploc) )


            iterate = True
            niter = 0
            while iterate:

                crossL = []
                for ip in idx:
                    #makes a superparticle
                    spi = np.array([system.locLold[ip],system.massLold[ip],system.mtotLold[ip]])
                    crossL.append(spi)
                crossL=np.array(crossL)

                
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
                    mass_t = 0.0    #gas accretion of planet, TBD
                    fcomp_t = 0.0   #how its composition changes

                #update planet properties from rates supplied by user
                planet_loc_nw = planet.loc + loc_t *system.deltaT

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
                Mc = ff.M_critical(system,planet.loc,crossL[k])
                if Mc<planet.mass:                    
                    epsilon = ff.epsilon_PA(system,planet.loc,planet.mass,crossL[k])
                    delm = epsilon*crossL[k][2]#don't understand this line...
                    
                else:
                    "pebble accretion can not happen"
                    delm=0
                
                planet.mass += delm #in crease mass (pebble accretion)
                # planet.fcomp += 0.  #TBD !!
                
                # import pdb; pdb.set_trace()
                #spN -> system.particles.Y2d...
                spN.mtotL[ip] -= delm #decrease mass sp


#perhaps this class object is not necessary...
class DISK (object):
    """
    Disk object including every disk properperties use the methods definded in the disk_properties.py
    """

    def __init__ (self, sigmaG, temp, mu, loc, time):
        """
        initialize with key disk properties
        """
        self.loc = loc
        self.sigmaG = sigmaG
        self.temp = temp
        self.mu = mu
        self.time = time
        self.alpha = dp.alpha
        self.rout = dp.rout
        self.rinn = dp.rinn
        self.tgap=dp.tgap
        self.sigmol = dp.sigmol
        # self.tgap=dp.tgap


    def add_auxiliary (self):
        """
        this add auxiliary disk properties that directly follow
        from the key disk properties
        """
        self.Mcp = self.Mcp_t(self.time)  
        self.OmegaK = np.sqrt(cgs.gC *self.Mcp/self.loc**3)      
        
        self.cs =  np.sqrt(cgs.kB*self.temp/(self.mu*cgs.mp))
        self.vth = np.sqrt(8/np.pi)*self.cs 

        self.vK = self.loc *self.OmegaK
        self.Hg = self.cs/self.OmegaK 
        self.nu = self.alpha*self.cs*self.Hg
        self.rhog = self.sigmaG/(2*np.pi)**0.5/self.Hg
        self.lmfp = self.mu*cgs.mp/(self.sigmol*self.rhog)
        

    def user_difined(self):
        #move to user-defined
        self.dotMg = self.dot_Mg(self.time)
        self.dotMd = self.dot_Md()
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
    
    # def Omega_K(self,loc,time):
    #     return np.sqrt(cgs.gC *self.Mcp/loc**3)
        #return dp.Omega_K(loc,time,self.Mcp)

    def c_s (self):
        return np.sqrt(cgs.kB*self.temp/(self.mu*cgs.mp))

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
    
    def eta_cal (self,loc):
        return dp.eta(loc,self.Mcp,self.dotMg,self.mu*cgs.mp)


class Superparticles(object):

    #CWO: these parameters look ugly here..
    rhoint=1.4
    pi=np.pi
    Sto=0.0001  
    error=1e-8


    def __init__(self,rinn,rout,dcomposL,gas,nini=10,Rdi=0.1,
            initrule='equalmass'):
        """
        systems initial properties

        initrule:[equalmass,equallogspace]
                :how initial parameters are distributed

        gas     :gas object (needed to get gas surface density)   

        nini: initial number of the particles
        mini: initial mass for every particles
        rinn: inner edge of the disk
        rout: outer edge of the disk
        icelineLoc: the location of icelines from inner to outer
        ice_frac  : the coresponding ice fraction of every iceline
        """

        #import pdb; pdb.set_trace()

        self.nini=nini
        mini = Rdi**3 *4*np.pi/3 *1.0

        #[23.12.30]this was commented out; now uncommented
        self.mini=mini #physical mass

        self.rinn=rinn
        self.rout=rout

        #TBD this (=location, mphys, mtot, ... of initial particles) should become user defined...
        self.locL=10**np.linspace(np.log10(rinn),np.log10(rout),nini)
        self.massL=np.array([mini]*nini)
        # if icelineLoc==None:
        #     self.mtot1=diskmass/nini
        #     self.mtotL=self.mtot1*np.ones_like(self.locL)
        # if len(icelineLoc)==1:
        #     self.mtot1 = diskmass /(len(np.where(self.locL<icelineLoc))*ice_frac+len(np.where(self.locL>=icelineLoc)))
        #     self.mtotL = np.where(self.locL<icelineLoc[0], self.mtot1*0.5, self.mtot1)
        # if len(icelineLoc)==2:
        #     self.mtot1 = diskmass /()

        
        #[23.12.30]:copied from /NewLagrance
        def construct_farr (dcomposL):
            fnL = []; miL = []
            for dcompos in dcomposL: 
                fnL.append(dcompos['Z_init'])
                miL.append(dcompos['mask_icl'])

            def f_arr (rad):
                farr = [fn(rad) for fn in fnL]
                mirr = [mi(rad) for mi in miL]
                return np.array(farr) * np.maximum(0, mirr)
            
            return f_arr

        #function gives the initial abundance of species in disk
        f_arr = construct_farr (dcomposL)

        #[23.12.30]:copied from /NewLagrance
        def f_sample (r):
            #samples the initial solid mass distribution
            Zcom = f_arr(r)

            #get the initial surface density
            sigini, *dum = gas.get_key_disk_properties(r,0.0)
            return 2*np.pi*r*sigini *Zcom.sum()


        print('[core.Superparticles.init]:initialization superparticles under rule:', initrule)
        if initrule=='equallogspace':
            #put sp at equal distances in log space

            xarr = np.linspace(1/nini, 1, nini)
            radL = rinn *(rout/rinn)**xarr

            msup = np.zeros_like(radL)
            r0 = rinn
            for k,r1 in enumerate(radL):
                msup[k], err = sciint.quad(f_sample, r0, r1)
                r0 = r1

        elif initrule=='equalmass':
            #puts superparticles at location such that they have
            #equald mass

            Mtot, err = sciint.quad(f_sample, rinn, rout)
            t0 = time.time()
            print('[core.Superparticles.init]:calling rout.sample_equal... this may take a while')
            radL = ff.sample_equal (f_sample, rinn, rout, nini)
            radL = np.array(radL)
            radL = np.sqrt(radL[1:]*radL[:-1])
            t1 = time.time()
            print('[core.Superparticles.init]:sampling done ({:4.1f} sec)'.format(t1-t0))

            #the mass of the super-particles are equally spread through
            msup = np.ones_like(radL) *Mtot/nini

        #[23.12.30]:I don't know why/how mtot1 should be defined
        self.locL = np.array(radL)
        self.mtotL = np.array(msup)
        self.mtot1 = msup[-1] #??

        #[23.12.30]NEW:add composition data (fcomp)
        #[23.12.30]this looks a bit ugly...
        self.fcomp = np.empty((nini,len(pars.composL)))
        for k,rad in enumerate(radL):
            Zcomp = []
            for ic,scomp in enumerate(pars.composL):
                Zcomp.append(dcomposL[ic]['Z_init'](rad)*max(0,dcomposL[ic]['mask_icl'](rad)))
            Zcomp = np.array(Zcomp)

            #Zcomp = np.array([dcompos['Z_init'](rad)*max(0,dcompos['mask_icl'](rad)) for dcompos in dcomposL])
            self.fcomp[k,:] = Zcomp/sum(Zcomp)

        #[23.12.30]old stuff... this can be removed
        #TBR
        if False:
            if len(icelineL)!= 0:
                niceline=len(icelineL)
                fracs=[]
                nsum=0
                ns=[]
                for i in range(niceline):
                    n=len(self.locL[self.locL<icelineL[i].loc])-nsum
                    nsum+=n
                    ns.append(n)
                    frac=1
                    for j in range(i):
                        frac-=icelineL[-j].frac  #looking notes in Pad
                    fracs.append(frac)
                ns.append(len(self.locL)-nsum)
                f=1

                for k in range(niceline):
                    f-=icelineL[k].frac
                fracs.append(f)

                fracs.reverse()
                
                factors=np.array(ns)*np.array(fracs)
                self.mtot1=diskmass/np.sum(factors)
                self.mtotL=np.array([])
                
                for k in range(len(icelineL)+1):
                    self.mtotL=np.append(self.mtotL,np.ones(ns[k])*self.mtot1*fracs[k])
            else:
                self.mtot1=diskmass/nini
                self.mtotL=np.array([self.mtot1]*nini)


        self.generate_Y2d()   #get a Y2d used to integrate

        #TBR
        # self.Y2d = np.empty((ndim,nini))
        # self.Y2d[0] = self.locL
        # self.Y2d[1] = self.massL
        # self.Y2d[2] = self.mtotL

    def generate_Y2d(self):
        self.Y2d = np.array([self.locL, self.massL])

    def dY2d_dt (self,Y2d,t,gas):
        """
        input:
            Y2d -- state vector
            time -- time
            disk -- disk object
        """

        #unpack the state vector
        r, mphy = self.Y2d   #maybe the total mass needn't to be put in Y2d

        Rd=(mphy/(self.rhoint*4/3*np.pi))**(1/3)

        #get disk object instance from gas
        #maybe like this???

        #(1) provides key properties (T,mu,Sigma)
        #(2) make disk objects
        #(3) add additional properties, which follow from (1) to disk object
        #(4) add user-defined disk properties to disk class
        out = gas.get_key_disk_properties (r, t)
        disk = DISK (*out, r, t) #pro
        disk.add_auxiliary ()
        disk.user_difined ()

        eta=disk.eta
        v_K=disk.vK
        v_th=disk.vth
        lmfp=disk.lmfp
        rho_g=disk.rhog
        Omega_K=disk.OmegaK
        H_g=disk.Hg
        dotMd=disk.dotMd


        #obtain Stokes number by iterating on drag law
        St,v_r = ff.St_iterate(eta,v_K,v_th,lmfp,rho_g,Omega_K,Rd)

        #describe what's going on here
        v_dd = np.abs(v_r)/2
        H_d = H_g*(1+St/disk.alpha*(1+2*St)/(1+St))**(-0.5)

        #this becomes...
        H_d = disk.Hg *(1+St/disk.alpha*(1+2*St)/(1+St))**(-0.5)
            

        drdt = v_r
        #dR_ddt= v_dd*dot_M_d/4/pi**(3/2)/rho_int/H_d/r/v_r**2 *dr_dt

        sigD = dotMd /(-2*r*self.pi*v_r) #v_r<0
        
        dmdt = 2*np.sqrt(np.pi)*Rd**2*v_dd/H_d*sigD  #eq. 5 of Shibaike et al. 2017

        Y2ddt = np.zeros_like(self.Y2d)
        Y2ddt[0] = drdt
        Y2ddt[1] = dmdt
        # Y2ddt[2] = 0.0

        return Y2ddt 
    
    
    def update (self,t0,tFi,gas,nstep=10):
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
        # self.mtotL=Yt[-1,2,:] #no longer updated
        #TBR
        self.generate_Y2d()  #update Y2d for next step integration
        return Yt

    
    def remove_particles(self,remove_idx):
        #TBD: remove self.loc, self.mphys...
        self.mtotL =  np.delete(self.mtotL, remove_idx) 
        self.locL = np.delete(self.locL, remove_idx)
        self.massL = np.delete(self.massL, remove_idx)   
        
        self.generate_Y2d()  #get a new Y2d, update.

    def add_particles(self,Nadd):
        #TBD: no longer works with Y2d

        self.locL = np.append(self.locL, self.rout)
        self.massL = np.append(self.massL, self.mini)
        self.mtotL = np.append(self.mtotL, self.mtot1)

        self.generate_Y2d()  #update Y2d

        if Nadd!=1:
            print('[core]Error:can only add 1 particle // reduce timestep')
            sys.exit()


        # For the situation where not only one particles will be added, maybe useful in future
             
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

class ICELINE(object):
    def __init__(self,species,temp):
        self.species=species
        self.temp=temp
        #self.frac=frac
    
    def find_iceline (self,rad, time, gas):
        Tice = self.temp
        Tdisk = gas.get_key_disk_properties(rad,time)[1]
        return Tdisk -Tice

    def get_icelines_location(self,gas,time):
        """
        get location of iceline, whose temperature is assumped as 160K
        """
        locL = np.linspace(dp.rinn,dp.rout,1000)
        diffold = 1e5

        #change the bounds to make it general
        self.loc = sciop.brentq(self.find_iceline, cgs.RJ, 20*cgs.RJ, args=(time,gas))
        # import pdb; pdb.set_trace()


