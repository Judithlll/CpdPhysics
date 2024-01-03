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

class COPY (object):
    """
    this copies certain attributes over 
    """

    def __init__ (self, state, attributeL):

        for attr in attributeL:
            setattr(self, attr, copy.deepcopy(getattr(state,attr)))


class System(object):

    """
    SYSTEM: integrate the every part of the PPD (or CPD) and update them with time.
    """

    ## CWO: These parameters should be initialized with pars.dsystempars; not hard-coded like this..

    ## CWO: Many of these par also no longer used... please, clean up 
    fraction=0.02 #CWO -- what's this?
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

        #initiallize the old state
        self.oldstate=None
    
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


    def update_particles (self):
        """
        Integrate the particles forward by amount deltaT

        Evolving system to the self.time
        """

        #time derivative of particles

        #TBD: make an Y2d array...


        ## CWO: self.deltaT is now determined separately

        #update particle properties
        # self.back_up_last_data()
        Yt = self.particles.update(self.time,self.time+self.deltaT,self.gas,self.timestepn)

        #self.deltaT = deltaT
        # print(self.particles.locL, self.time)
        if self.time==np.nan or self.deltaT==np.nan:
            print('hello')
            import pdb;pdb.set_trace()

        ## CWO: why do we return stuff?
        return Yt

    def back_up_last_data(self):
        """
        copies present state to "old" 
        CWO: would it be good also to copy the time?
        """
        #LZX: a little confused about how to use this, maybe make the COPY a sub class?
        #CWO: OK, try it out
        self.oldstate = COPY (self, ['time', 'particles', 'planetL', 'icelineL'])



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

        ##CWO: TBD 
        #Make mtot1 smaller (greater) if total number of particles
        #falls below (exceeds) some limit
        
        #Ncrit = 100
        #N<90: mtot *= 0.9 && (wait some time?)


    def new_timestep (self, tEnd):
        """
        chooses a timestep
        """
        fp = 0.2 #The 0.2 should become a parameter...

        #organize the procedure a bit (for quasi-steady evolution... later!)
        mintimeL = []

        Y2d = self.particles.make_Y2d()
        Y2dp = self.particles.dY2d_dt(Y2d,self.time,self.gas)

        #timescale for the particles
        #I dont think there's need to use np.nanmin
        tpart = np.abs(Y2d/Y2dp)
        mintimeL.append({'name':'particles', 'tmin': fp*tpart.min(), 
                                'imin':np.unravel_index(tpart.argmin(),tpart.shape)})

        # if self.ntime>=1:
        #     import pdb; pdb.set_trace()
        
        #[24.01.02]CWO: I think it's possible to have an dm/dt & da/dt to planet objects?
        #               whay you do here is a bit ugly...
        if self.oldstate is not None:   
            #timescale for the planets (including migration and mass growth)
            PmassTscale=np.inf*np.ones_like(self.planetL)
            PlocaTscale=np.inf*np.ones_like(self.planetL) 
            for i in range(self.nplanet):
                if self.time>self.planetL[i].time:
                    # import pdb ;pdb.set_trace()
                        PmassTscale[i]=np.float64(self.planetL[i].mass)/abs(self.oldstate.planetL[i].mass-self.planetL[i].mass)*self.deltaT
                        PlocaTscale[i]=np.float64(self.planetL[i].loc)/abs(self.oldstate.planetL[i].loc-self.planetL[i].loc)*self.deltaT

            mintimeL.append({'name': 'planetsMigration', 'tmin': min(PlocaTscale)})
            mintimeL.append({'name': 'planetsGrowth', 'tmin': min(PmassTscale)})
            
            #timescale for the icelines
            IlocaTscale=np.inf*np.ones_like(self.icelineL)
            for i,iceline in enumerate(self.icelineL):
                tscale=np.float64(iceline.loc)/abs(self.oldstate.icelineL[i].loc-iceline.loc)*self.deltaT
                IlocaTscale[i]=tscale
            
            mintimeL.append({'name': 'icelineloca', 'tmin': min(IlocaTscale)})
        
        # put mintimeL into system object for now to check
        self.mintimeL=mintimeL

        deltaT = np.inf
        for ddum in mintimeL:
            if ddum['tmin'] < deltaT:
                deltaT = ddum['tmin']
            

        if self.time+deltaT>tEnd:
            deltaT = tEnd - self.time

        self.deltaT = deltaT



def advance_iceline (system):
    """
    for now particles directly lose the mass of water without any other effect
    """

    sploc = system.particles.locL
    sploc_old = system.oldstate.particles.locL
    for k,iceline in enumerate(system.icelineL):
        idx,=np.nonzero((iceline.loc<sploc_old) & (iceline.loc>sploc))

        ic = pars.composL.index(iceline.species) #refers to species index
        if len(idx)!=0:     
            
            for ix in idx:
                fice = system.particles.fcomp[ix,ic]  #mass fraction in ice
                fremain = (1-fice)          #remain fraction

                if fremain < 1e-15:
                    fremain=0 #loss of numbers (!!)
                system.particles.mtotL[ix] *= fremain    #reduce masses accordingly
                system.particles.massL[ix] *= fremain
                system.particles.fcomp[ix,ic] = 0.      #gone is the ice!

                #renormalize
                system.particles.fcomp[ix,:] = (system.particles.fcomp[ix,:].T /(system.particles.fcomp[ix,:].sum()+1e-100)).T


def advance_planets (system):
    import planets_properties as pp
    """
    [23.12.06]copied/edited from NewLagrange

    For now, planet migration and composition is not considered

    TBD:
        - add composition changes to planets
        - add migration rate
    """
    for planet in system.planetL:

        #planet exists only after planet.time
        if planet.time<system.time:

            sploc = system.particles.locL
            sploc_old = system.oldstate.particles.locL

            #particles that cross are those that
            idx, = np.nonzero( (planet.loc<sploc_old) & (planet.loc>sploc) )


            iterate = True
            niter = 0
            while iterate:


                crossL = []
                for ip in idx:
                    #makes a superparticle

                    ## CWO: perhaps nicer to make a particle object (instead of this array)
                    #       such that we can have spi.loc, spi.mphy, spi.msup, spi.fcomp
                    spi = np.array([system.oldstate.particles.locL[ip],system.oldstate.particles.massL[ip],system.oldstate.particles.mtotL[ip]])
                    #spi = Single (...)
                    #spi = particles.select(ip)
                    #spi = system.oldstate.particles.select(ip)
                    #spi = system.oldstate.particles.select_multi(idx)
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
            spN = system.particles
            
            for k, ip in enumerate(idxN):
                #mass that is being transferred (TBC!!)
                #need to calculate epsilon (PA efficiency)

                #calculate critical mass to verify if the pebble accretion can occur
                Mc = ff.M_critical(system,planet.loc,crossL[k])
                if Mc<planet.mass:                    
                    epsilon = ff.epsilon_PA(system,planet.loc,planet.mass,crossL[k])

                    #2nd index refers to total particle mass
                    #CWO: prefer to have smth like: delm = epsilon*crossL[k].mtot
                    delm = epsilon*crossL[k][2]
                    # import pdb; pdb.set_trace()

                else:
                    "pebble accretion can not happen"
                    delm=0
                
                planet.mass += delm #in crease mass (pebble accretion)
                # planet.fcomp += 0.  #TBD !!
                
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

        #user defined in disk_properties.py
        self.alpha = dp.alpha
        self.rout = dp.rout
        self.rinn = dp.rinn
        self.sigmol = dp.sigmol


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

        self.rhocompos=[]
        for compos in dcomposL:
            if compos['name']!= 'gas':
                self.rhocompos.append(compos['rhoint'])
        
        self.nini=nini

        #[23.12.30]this was commented out; now uncommented

        self.rinn=rinn
        self.rout=rout
        self.stokesOld = None
        
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

        #[24.01.01]this is a bit ugly... but necessary for adding particles
        self.fcompini = self.fcomp[-1]
        
        self.get_rhoint()
        
        self.massL = self.rhoint * 4/3*Rdi**3*np.pi
        self.mini = self.massL[-1]   #for adding particles

        self.generate_Y2d()   #get a Y2d used to integrate

    ## CWO: I'm not sure if we need to attach Y2d to self
    #       since Y2d is no longer fundamental...
    #       perhaps just "make_Y2d" is sufficient?
    def generate_Y2d(self):
        self.Y2d = np.array([self.locL, self.massL])

    def make_Y2d (self):
        return np.array([self.locL, self.massL])

    def get_rhoint(self):
        "get the true fcomp according to fcomp"
        Volume = np.zeros_like(self.locL)
        for i in range(len(self.rhocompos)):
            Volume[:] += self.fcomp[:,i]/self.rhocompos[i]
        
        self.rhoint=1/Volume

    def get_radius(self):
        "update the rhoint according to new fcomp"
        self.get_rhoint()
        return (self.massL/(self.rhoint*4/3*np.pi))**(1/3)
    
    def dY2d_dt (self,Y2d,time,gas):
        """
        input:
            Y2d -- state vector
            time -- time
            disk -- disk object
        """

        #unpack the state vector
        loc, mphy = Y2d

        Rd=self.get_radius()

        out = gas.get_key_disk_properties (loc, time)
        disk = DISK (*out, loc, time) #pro
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
        St,v_r = ff.St_iterate(eta,v_K,v_th,lmfp,rho_g,Omega_K,Rd,Sto=self.stokesOld)
        self.stokesOld = St

        #assume the relative velocity to be the half of radial velocity
        v_dd = np.abs(v_r)/2  

        #make the dust scale height to be user defined
        H_d = dp.H_d(H_g, St)     

        drdt = v_r
        #dR_ddt= v_dd*dot_M_d/4/pi**(3/2)/rho_int/H_d/r/v_r**2 *dr_dt

        ## CWO: surface density should follow from position of the Lagrangian particles...
        sigD = dotMd /(-2*loc*np.pi*v_r) #v_r<0
        
        ## CWO: this *could* become a userfun (b/c seems a bit specific)
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
        tstep=(tFi-t0)/nstep
    
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
        self.fcomp = np.delete(self.fcomp, remove_idx, axis=0) #[24.01.01] added

        
        #TBR: not necessary... only generate when you need it!
        self.generate_Y2d()  #get a new Y2d, update.

    def add_particles(self,Nadd):
        #TBD: no longer works with Y2d

        self.locL = np.append(self.locL, self.rout)
        self.massL = np.append(self.massL, self.mini)
        self.mtotL = np.append(self.mtotL, self.mtot1)  #mtot1 is needed here
        self.fcomp = np.append(self.fcomp, [self.fcompini], axis=0) #[24.01.01] added


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


