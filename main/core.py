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
import scipy.optimize as sciop
import physics
import userfun

class COPY (object):
    """
    This object used to back up the system state of last time point
    """

    def __init__ (self, state, attributeL):

        for attr in attributeL:
            setattr(self, attr, copy.deepcopy(getattr(state,attr)))

class Mintimes (object):
    """
    store the minimam times
    """
    def __init__(self, mintimedict):
        for i in range (len(mintimedict)):
            setattr(self, mintimedict[i]['name'], mintimedict[i]['tmin'])

        self.mintimelist(mintimedict)

    def mintimelist(self,mintimedict):
        self.tminarr = np.array([ddum['tmin'] for ddum in mintimedict])
        


class System(object):

    """
    SYSTEM: integrate the every part of the PPD (or CPD) and update them with time.
    """

    ## TBD-later: think about how to calculate rhoPlanet

    daction={}

    def __init__(self, timeini=0.0, rhoPlanet=1.9):
        """
        System initiallization parameters:
            
            timeini: [float]
                the initial time with default value 0.0(maybe not necessary)
            rhoPlanet: [float]
                internal density of planets with default value 1.9g/cm^3
        """
        self.rhoPlanet = rhoPlanet

        self.time=timeini
        self.ntime = 0  

        # define a disk class
        self.gas = self.init_gas ()

        #the amount of solid mass that has crossed into the domain
        self.dotMg = dp.dot_Mg(self.time)
        self.Minflux = 0
        self.Minflux_step = 0
        #self.Mcp=self.disk.Mcp_t(self.time)

        #initiallize the old state
        self.oldstate=None

        #self.planetMassData=[]
        self.njump = 0
        self.njumptime = 0

        self.milestones = {}
        #d = self.milestones
        #d[0] = 'initialization'
        #d[time1] = 'insert planet'
        #d[time2] = 'insert planet'
        #d[tmax] = 'final time'
        #
        #
        #timeM = d.keys()
        #eventsM = d.values()
        #for k,v in d.items(): if v=='event'

    
    def init_particles(self, dparticleprops={}):
        """
        because we need to consider iceline, so separatly 
        initiallize the particles, for now just water iceline is considered  

        history:
        [23.12.30]CWO: instead of forwarding diskmass, I supply self.gas to the superparticle class
        [24.01.08]LZX: mtot1 is generated from Superparticles, but is used much in post_process, so get this from superparticles for now
        """
        self.particles = Superparticles(dp.rinn,dp.rout,self.dcomposL,self.gas, **dparticleprops)

        
        self.mtot1 = self.particles.mtot1


    def init_gas (self, gasL=None, dcomposL=None, dgrid={}):
        """
        init_gas is called at initialization. It adds an instance of the GAS class

        history:
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


    def update_particles (self, timestepn=3, **kwargs):
        """
        Integrate the particles forward by amount deltaT

        Evolving system to the self.time
        """
        Yt = self.particles.update(self.time,self.time+self.deltaT,self.gas,timestepn)


        if self.time==np.nan or self.deltaT==np.nan:
            print('hello')
            import pdb;pdb.set_trace()

        return Yt

    def back_up_last_data(self):
        """
        copies present state to "old" 
        """

        self.oldstate = COPY (self, ['time', 'particles', 'planetL', 'icelineL', 'Minflux_step', 'dotMg'])



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


        delmgasIn, error = sciint.quad(dp.dot_Mg,self.time,self.time+self.deltaT)
        self.dotMg = dp.dot_Mg(self.time)
            
        #self.Minflux_step = dp.M_influx(self.time,self.time+self.deltaT)
        self.Minflux_step = delmgasIn
        self.Minflux += self.Minflux_step

        while self.Minflux>self.mtot1:
            Nadd += 1
            self.Minflux -= self.mtot1
        
        if Nadd>0:
            self.daction['add'] = Nadd
        
        #post_process particles
        if 'remove' in self.daction.keys():
            #remove the particles from Y2d!
            
            self.particles.remove_particles(self.daction['remove'])


        if 'add' in self.daction.keys():
            self.particles.add_particles(self.daction['add'])


        #[24.01.04]
        #it is really difficult to stabalize particle numbers, b/c 
        #of the huge lag... I think the below algorith accomplishes smth
        #but is a bit ugly

        #to be brief..
        part = self.particles

        Np = len(part.massL)

        nch = 4
        if Np==part.ninit: part.Nplevel = part.ninit #reset

        #particle level is decreasing...
        if Np<part.Nplevel -nch and Np<part.ninit:
            part.Nplevel -= nch
            eps = abs(part.ninit -part.Nplevel) /part.ninit
            self.mtot1 *= 1 -eps

        #particle level is decreasing, but above ninit: modest decrease mtot1
        elif Np<part.Nplevel -nch and Np>part.ninit:
            part.Nplevel -= nch
            eps = nch /part.ninit
            self.mtot1 *= 1 -eps

        #particle level is increasing, but below ninit: modest increase mtot1
        elif Np>part.Nplevel +nch and Np<part.ninit:
            part.Nplevel += nch
            eps = nch /part.ninit
            self.mtot1 *= 1 +eps

        #particle level is increasing...
        elif Np>part.Nplevel +nch and Np>part.ninit:
            part.Nplevel += nch
            eps = abs(part.ninit -part.Nplevel) /part.ninit
            self.mtot1 *= 1 +eps

        #Ncrit = 100
        #N<90: mtot *= 0.9 && (wait some time?)
        self.particles.generate_Y2d()


    def new_timestep (self, tEnd, deltaTfraction = 0.2, afterjump = False, **kwargs):
        """
        chooses a timestep
        """
        #organize the procedure a bit (for quasi-steady evolution... later!)
        mintimeL = []
        
        
        Y2d = self.particles.make_Y2d()
        Y2dp = self.particles.dY2d_dt(Y2d,self.time,self.gas)

        #timescale for the particles
        #I dont think there's need to use np.nanmin
        tpart = np.abs(Y2d/Y2dp)
        mintimeL.append({'name':'particles', 'tmin': deltaTfraction*tpart.min(), 
                                'imin':np.unravel_index(tpart.argmin(),tpart.shape)})
        
        #after the jump, we only have particles drift timescale
        if afterjump:
            self.minTimes = Mintimes(mintimeL)
            self.mintimeL = mintimeL
            self.deltaT = deltaTfraction*tpart.min() #min([ob['tmin'] for ob in mintimeL])

            return

        #[23.01.19]LZX: maybe we don't need this, because we can always get the accurate value from the user-defined function
        #central mass growth timescale
        #Mcpnew=dp.Mcp_t(self.time)  
        #Mcpold=dp.Mcp0
        #McTscale = np.inf
        #if self.time > 0:
        #    McTscale = Mcpnew/ abs(Mcpold - Mcpnew) *self.deltaT
        #    Mcpold = Mcpnew
        #    mintimeL.append({'name': 'CentralMassGrowth', 'tmin': McTscale})
            # import pdb; pdb.set_trace()
        
       
        #Mass influx timescale
        if self.time > 0:  # better to get rid of the first step one, but the effect is tiny
            #[24.01.07]CWO: in cases there is no influx, this should evaluate the infinity...
            #InfluxTscale = (1e-100 + np.float64(self.Minflux_step)) / abs(self.oldstate.Minflux_step - self.Minflux_step) *self.deltaT
            mdotgTscale = (1e-100 + np.float64(self.dotMg)) / abs(self.oldstate.dotMg - self.dotMg) *self.deltaT
            mintimeL.append({'name': 'Mass_Influx', 'tmin': mdotgTscale})

        
        if self.oldstate is not None:   
            #timescale for the planets 
            # (including migration and mass growth)
            
            if pars.doPlanets:
                PmassTscale = np.inf*np.ones(self.nplanet)
                PlocaTscale = np.inf*np.ones(self.nplanet) 
                PcompTscale = np.inf*np.ones(self.nplanet)
                thre_jump_max = 1e-3  #threshold when getting the max jumpT
                for i in range(self.nplanet):
                    planet = self.planetL[i]
                    if self.time>planet.starttime:

                            #store mass data first
                        if self.oldstate.planetL[i].mass != planet.mass:
                            # self.masstime=
                            planet.planetMassData.append([self.time , planet.mass])


                        #then try to fit the mass to a curve
                        #planet.loc_t
                        PlocaTscale[i]=np.float64(self.planetL[i].loc)/abs(self.oldstate.planetL[i].loc-self.planetL[i].loc)*self.deltaT
                        
                        #consider the data is not large enough to make the fit

                        Npts = len(planet.planetMassData)
                        Nfit = 10
                        if Npts >= Nfit:
                            #better way to do
                            timedots, massdots = np.log(np.array(planet.planetMassData).T)
                            #timedots = np.log10([planet.planetMassData[j][0] for j in range(len(planet.planetMassData))])
                            def mass_fit(logt,p,b):
                                logm = p*logt+b
                                return logm 
                            def jac_mass_fit (logt,p,b):
                                jac = np.array([logt, np.ones_like(logt)]).T
                                return jac

                            #massdots = np.log10([planet.planetMassData[j][1] for j in range(len(planet.planetMassData))])

                            relp = np.inf #relative error/uncertainty in pwl index
                            popt = np.array([1,1]) #initial guess (p0)
                            while Nfit<=Npts:
                                ##[24.01.20]CWO: this may help a bit, but can still be much faster b/c of linear square root
                                #popt, pcov = sciop.curve_fit(mass_fit, timedots[-Nfit:], massdots[-Nfit:])
                                popt, pcov = sciop.curve_fit(mass_fit, timedots[-Nfit:], massdots[-Nfit:], p0=popt, jac=jac_mass_fit)

                                #line = '{:7.4f} {:10.3e} {:10.3e} {:6d}'.format(popt[0], np.sqrt(pcov[0,0]), np.sqrt(pcov[0,0])/popt[0], Nfit)
                                #print(line)
                                if np.sqrt(pcov[0,0])/np.abs(popt[0])<relp:
                                    pidx = popt[0]
                                    psig = np.sqrt(pcov[0,0])
                                    relp = psig/np.abs(pidx)
                                else:
                                    break
                                Nfit = int(Nfit*1.2) +1
                            planet.relp_mass = relp
                            #if Npts>=10 and self.ntime>1000:
                            #    import pdb; pdb.set_trace()
                            #plt.scatter(timedots, massdots)
                            #t_list=np.linspace(timedots[0], timedots[-1], 30)
                            #plt.plot(t_list, mass_fit(t_list, *popt))
                            #plt.savefig('/home/lzx/CpdPhysics/Test/Zhixuan/test.jpg')
                            
                            planet.dmdt = pidx *planet.mass/(self.time) 
                            planet.dmdt_err = abs(psig *planet.mass/self.time) ##LZX: please check expressions
                            PmassTscale[i] = 1/abs(pidx)*(np.exp(timedots[-1]) - planet.starttime)
                            #planet.tmass_err = abs(psig/popt[0] *PmassTscale[i])
                            
                            #jump time is limited by uncertainty in the fit
                            denom = (planet.dmdt_err - planet.dmdt*thre_jump_max)
                            if denom<0:
                                planet.max_jumpT = np.inf
                            else:
                                planet.max_jumpT = thre_jump_max*planet.mass /denom

                        else:
                            planet.dmdt = 0
                            planet.relp_mass = np.inf 
                            planet.max_jumpT = 0.0 #np.nan
                        
                        #get the planet composition change time scale
                        #TBD: do the composition 'crude' for the moment, assuming the comp of added mass during the jump is the same as the nearest particles
                        if False:
                            if self.oldstate.planetL[i].fcomp[0] != planet.fcomp[0]:
                                planet.compData.append([self.time, planet.fcomp[0]])
                            
                            if len(planet.compData) >10:
                                timedots, compdots = np.array(planet.compData).T
                                timedots = np.log(timedots)
                                def comp_fit(t,a,b):
                                    c = a*t+b
                                    return c
                                
                                relp = np.inf
                                while Nfit<=Npts:
                                    popt, pcov = curve_fit(mass_fit, timedots[-Nfit:], compdots[-Nfit:])
                                    #line = '{:7.4f} {:10.3e} {:10.3e} {:6d}'.format(popt[0], np.sqrt(pcov[0,0]), np.sqrt(pcov[0,0])/popt[0], Nfit)
                                    #print(line)
                                    if np.sqrt(pcov[0,0])/np.abs(popt[0])<relp:
                                        pidx = popt[0]
                                        psig = np.sqrt(pcov[0,0])
                                        relp = psig/np.abs(pidx)
                                    else:
                                        break
                                    Nfit = int(Nfit*1.2) +1
                                planet.relp_comp = relp
                                
                                PcompTscale[i] = (np.exp(timedots[-1])-planet.starttime)* compdots[-1]/abs(popt[0]) 
                                planet.comp_tc = PcompTscale[i]

                            else:
                                planet.relp_comp = np.nan
                        
                mintimeL.append({'name': 'planetsMigration', 'tmin': min(PlocaTscale)})
                mintimeL.append({'name': 'planetsGrowth', 'tmin': min(PmassTscale)})
                mintimeL.append({'name': 'planetsComp', 'tmin': min(PcompTscale)})

            #timescale for the icelines
            if pars.doIcelines:
                IlocaTscale=np.inf*np.ones_like(self.icelineL)
                for i,iceline in enumerate(self.icelineL):
                    if not np.isnan(iceline.loc):
                        tscale = np.float64(iceline.loc)/abs(self.oldstate.icelineL[i].loc-iceline.loc)*self.deltaT
                        iceline.loc_tc = tscale
                        IlocaTscale[i] = tscale
            
                mintimeL.append({'name': 'icelineloca', 'tmin': min(IlocaTscale)})
        
        # put mintimeL into system object for now to check
        self.mintimeL=mintimeL
        #make a class to use this mintimeL, change the name 
        self.minTimes = Mintimes(mintimeL)
        deltaT = np.inf
        for ddum in mintimeL:
            if ddum['tmin'] < deltaT:
                deltaT = ddum['tmin']
            
        if self.time+deltaT>tEnd:
            deltaT = tEnd - self.time

        self.deltaT = deltaT


    def query_system_jump(self, jumpfrac):
        """
        investigates w/r we can jump and by how much
        ...

        returns True/False, {jump properties}
        """
        Tscale_ratio = []
        for t in self.minTimes.tminarr[1:]:
            Tscale_ratio.append( t/self.minTimes.particles)
        #print(self.minTimes.tminarr) 
        if len(self.minTimes.tminarr[1:])==0:
            return False, {}

        #jump time is given by the min evol. timescales, excluding those of the particles
        max_tevol = jumpfrac*min(self.minTimes.tminarr[1:])
        
        #(cannot jump over "important events" -> Milestones) 
        #adjust jumpT according to milestones
        timepoints = np.sort( np.array(list(self.milestones.keys())) )

        ii = timepoints >self.time
        if sum(ii)==0:
            max_tms = np.inf
        else:
            max_tms = min(timepoints[ii]-self.time)
        
         
        #reached_ms = np.where( (timepoints - self.time - jumpT<0) & (timepoints - self.time>0) )[0] 
        #we cannot jump over milestones
        #if len(reached_ms) != 0:
        #   jumpT = timepoints[reached_ms[0]] - self.time
        
        #jumpT cannot exceed max_jumpT
        #NOTE: planet may not exist... hence the for/if construct
        max_tpl = np.inf
        for planet in self.planetL:
            if self.time>planet.starttime:
                max_tpl = min(max_tpl, planet.max_jumpT)


        #the jump time is the minimum
        #[24.01.21]I wrote this a bit more clearly...
        tjumpkeys = ['system-evol', 'milestones', 'planet-growth-fit']
        tjumparr = np.array([max_tevol, max_tms, max_tpl])
        jumpT = min(tjumparr)

        #evol.timescales >> drift timescales && ensure significant jump
        con0 = self.ntime >self.njumptime +100 
        con1 = (min(Tscale_ratio) > 1e2) and (jumpT/self.deltaT>100) #should have more judgements

        if False:
            # the fit goodness of planets grow >0.95, if the r2 is np.inf, the planets isn't added, if the planet is already added but r2 is still np.nan, then the mass data points are not enough 
            relp_massL = np.array([p.relp_mass for p in self.planetL if p.starttime < self.time])
            if len(relp_massL) != 0:
                con2 = relp_massL.min() < 2e-3
            else:
                con2 = False
            #if self.time > 50*cgs.yr:
            #    import pdb; pdb.set_trace()
            # consider the fit goodness of the planet's composition change, similar with the mass
            relp_compL = np.array([abs(p.relp_comp) for p in self.planetL if p.starttime < self.time])
            if len(relp_compL) != 0:
                con3 = relp_compL.min() < 2e-3
            else:
                con3 = False

            jumptf = con1 & con2 &con3
        else:
            jumptf = con0 & con1

        djump = {'jumpT':jumpT, 'tjumpkeys':tjumpkeys, 'tjumparr':tjumparr}
        #if jumptf: import pdb; pdb.set_trace()

        return jumptf, djump


    def reset_after_jump(self):
        #reset the planet mass data
        for planet in self.planetL:
            planet.planetMassData = []
            planet.relp_mass = np.inf #Chris: never assign "nan" please
            planet.max_jumpT = 0.0
            #planet.compData = []
            #planet.relp_comp = np.nan
        
    def system_jump(self, djump):
        """
        execute the system jump
        """
        self.jumpT = djump['jumpT']
        
        # parameters needs to be upda`ted:
        # planets: location and mass and composition(this maybe very complex)
        # icelines :location
        
        if pars.doPlanets:
            for planet in self.planetL:
                if self.time > planet.starttime:

                    planet.loc += planet.dlocdt *self.jumpT
                    #planet.loc -= planet.loc/ self.minTimes.planetsMigration *jumpT
                    #planet.mass += planet.mass/ self.minTimes.planetsGrowth *jumpT
                    planet.mass += planet.dmdt *self.jumpT

                    #TBD: generalize this. Perhaps best way is to make planet.dmdt a vector
                    #       planet.dmdt = [dmdt comp 1, dmdt comp 2, ...]
                    sil_comp = planet.fcomp[0]-planet.fcomp[0]/self.minTimes.planetsComp *self.jumpT
                    planet.fcomp -=[sil_comp, 1-sil_comp]

        if pars.doIcelines:
            for iceline in self.icelineL:
                iceline.loc -= iceline.loc/self.minTimes.icelineloca *self.jumpT
        
        self.njump +=1
        self.njumptime = self.ntime

        im = djump['tjumparr'].argmin()
        #maybe interesting to store and plot which factor limits the jumpT
        self.jump_limitation = djump["tjumpkeys"][im]

        nameL = [d['name'] for d in self.mintimeL]
        tminarr = np.array([d['tmin'] for d in self.mintimeL])
        imin = 1 +tminarr[1:].argmin() #minimum evolution

        print(f'[system.jump]:at {self.time/cgs.yr:5.2f} yr jumped by {self.jumpT/cgs.yr:5.2f} yr')
        print(f'[core.system_jump]:jump time limited by: {self.jump_limitation}')
        print(f'[core.system_jump]:min. evolution time ({nameL[imin]}) {tminarr[imin]/cgs.yr:8.2e} yr')
        
        if self.planetL[0].loc <dp.rinn:
            print('[core.system.jump]: the first planet migrates across the inner edge')
            import pdb;pdb.set_trace()
        #"erase" previous planet.crossL OR record the jump time to planet.
        #such that new fit for dm/dt starts w/ N=0 particles




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
        
        #a little time comsuming, don't consider this for now
        if system.time > 1e6 *cgs.yr:
            loc_pv = system.oldstate.icelineL[k].loc
            iceline.get_icelines_location(system.gas,system.time,guess=loc_pv)



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
        if planet.starttime<system.time:

            sploc = system.particles.locL
            sploc_old = system.oldstate.particles.locL

            #particles that cross are those that
            idx, = np.nonzero( (planet.loc<sploc_old) & (planet.loc>sploc) )


            iterate = True
            niter = 0
            while iterate:


                crossL = []
                for ip in idx:

                    spi = system.oldstate.particles.select_single(ip)

                    crossL.append(spi)

                #crossL=np.array(crossL)

                
                # if len(idx)!=0:
                #                     #this user-defined function should detail how the
                #planet properties (location,mass,composition) change
                #with time and how the crossed particles are affected... 


                #this is about planet migration...
                #TBD later...
                if False:
                    ## CWO: switch to the user-defined approach
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
            planet.dlocdt = loc_t

            planet.time = system.time  ## add this time to planet

            #planet.mass += mass_t *system.deltaT
            #planet.fcomp += fcomp_t *system.deltaT

            #update s-particle properties from sp-crossings
            #assumes all particles are accreted (TBC!!)
            spN = system.particles

            #mass*composition that is acccreted
            delmcomp = np.zeros((spN.fcomp.shape[1]))
            
            for k, ip in enumerate(idxN):

                #calculate critical mass to verify if the pebble accretion can occur
                 
                #[24.01.05]
                ## CWO: I dont think it's a good idea to forward the entire system
                #       class to such simple functions..
                #
                #       (for the moment I ignore it... we can discuss)


                spk = crossL[k]
                Mc = userfun.M_critical(spk.eta, spk.St, spk.mcp)

                delmass = np.zeros_like(crossL[0].fcomp)

                #TBD-later: I dont like the need for an M_critical userfun..
                #   instead, we can incorporate this into userfun.epsilon_PA
                #
                #On the other hand, particles may get stuck in pressure bump   
                #and fail to accrete and drift. This behavior would be nice to capture
                #but this requires more thinking..
                #
                #TBD-later: in reality particles may be "stuck" in pressure bump
                #           incorporate in planet object?
                if Mc<planet.mass:                    


                    epsilon = userfun.epsilon_PA(planet.loc,planet.mass,spk)

                    #accreted mass by composition
                    delmcomp += epsilon*crossL[k].fcomp *crossL[k].mtotL

                    #delm = epsilon*crossL[k].mtotL
                    #for i in range(len(crossL[k].fcomp)):
                    #    delmass[i] = epsilon* crossL[k].fcomp[i]*crossL[k].mtotL
                    
                    #planet.mass*planet.fcomp
                    
                    #planet.fcomp = [ (delmass[i]+planet.mass*planet.fcomp[i]) / (planet.mass+delmass.sum()) for i in range(len(delmass))]
                    #planet.mass += delmass.sum() #increase mass (pebble accretion)
                else:
                    print('[core]: pebble accretion can not happen')
                    import pdb; pdb.set_trace()


                # planet.fcomp += 0.  #TBD !!
                
                #spN -> system.particles.Y2d...
                spN.mtotL[ip] -= delmass.sum() #decrease mass sp
            masscomp = planet.fcomp*planet.mass +delmcomp
            planet.mass += delmcomp.sum()
            planet.fcomp = masscomp /planet.mass


class SingleSP(object):
    """
    Used to get single Superparticle
    """
    def __init__ (self,**kwargs):
        for key,val in kwargs.items():
            setattr(self,key,val)


class Superparticles(object):

    def __init__(self, rinn, rout, dcomposL, gas, nini=10, Rdi=0.1, 
            initrule='equalmass'):
        """
        systems initial properties

        initrule:[equalmass,equallogspace]
                :how initial parameters are distributed

        gas     :gas object (needed to get gas surface density)   

        nini    : initial number of the particles
        rinn    : inner edge of the disk
        rout    : outer edge of the disk
        dcomposL: the composition list

        """

        self.rhocompos=[]
        for compos in dcomposL:
            if compos['name']!= 'gas':
                self.rhocompos.append(compos['rhoint'])
        
        self.num = nini
        self.ninit =nini
        self.Nplevel = nini

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
                #after change the mask_icl getting location, there will be a strange warning, by set the limit=100 can remove this warning
                msup[k], err = sciint.quad(f_sample, r0, r1, limit =100)
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

        self.locL = np.array(radL)
        self.mtotL = np.array(msup)
        self.mtot1 = msup[-1] #for adding new particles

        #[23.12.30]NEW:add composition data (fcomp)
        #[23.12.30]this looks a bit ugly...
        #[24.01.09]LZX:the initial fcomp also need to consider the iceline location, TBD
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


    def generate_Y2d(self):
        self.Y2d = np.array([self.locL, self.massL])


    def select_single(self, ix):

        kwargs = {}
        # select the properties that are list or numpy.ndarray
        propL = [attr for attr in dir(self) if not attr.startswith('__') and isinstance(getattr(self, attr), list) or isinstance(getattr(self, attr), np.ndarray)]   
        # propL = ['locL','massL','mtotL','fcomp','St','eta'] maybe just select properties artificially is better
        # select the properties that are float
        propSol = [attr for attr in dir(self) if not attr.startswith('__') and isinstance(getattr(self, attr), float)]

        for prop in propL:
            if len(getattr(self,prop)) > len(self.rhocompos):
                kwargs[prop] = getattr(self,prop)[ix]
        
        for prop in propSol:
            kwargs[prop] = getattr(self,prop)
        

        spi = SingleSP (**kwargs)
        return spi


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

    
    def dY2d_dt (self,Y2d,time,gas, returnMore=False):
        """
        input:
            Y2d -- state vector
            time -- time
            disk -- disk object
        """

        #unpack the state vector
        loc, mphy = Y2d

        #get radii of the particles
        Rd = self.get_radius()


        #we need these 2 things to initalize the class object
        mcp = dp.Mcp_t(time)  
        out = gas.get_key_disk_properties (loc, time)

        disk = physics.DISK (*out, loc, time, mcp) #pro
        disk.add_auxiliary ()
        
        userparL = disk.add_uservar (dp.user_add_var())    #variables
        disk.add_userfun (dp.user_add_fun())    #functions only

        userevalL = disk.add_user_eval (dp.user_add_eval()) #evaluations

        #obtain Stokes number by iterating on drag law
        St, v_r = ff.St_iterate (disk.eta,
                                 disk.vK,
                                 disk.vth,
                                 disk.lmfp,
                                 disk.rhog,
                                 disk.OmegaK,
                                 Rd,
                                 Sto=self.stokesOld)
        self.stokesOld = St


        #assume the relative velocity to be the half of radial velocity
        #v_dd = np.abs(v_r)/2    

        #make the dust scale height to be user defined
        Hd = physics.H_d(disk.Hg, St)     

        drdt = v_r
        #dR_ddt= v_dd*dot_M_d/4/pi**(3/2)/rho_int/H_d/r/v_r**2 *dr_dt

        ## CWO: surface density should follow from position of the Lagrangian particles...
        sigD = disk.dotMd /(-2*loc*np.pi*v_r) #v_r<0
        rhoD = sigD /np.sqrt(np.pi*2) /Hd

        #relative velocity may depend on: alpha, cs, St, rhod/rhog, ..
        delv = userfun.del_v (St, rhoD, disk)
        
        ## CWO: this *should* become a userfun (b/c seems a bit specific)
        dmdt = userfun.dm_dt (Rd, delv, Hd, sigD)

        Y2ddt = np.zeros_like(self.Y2d)
        Y2ddt[0] = drdt
        Y2ddt[1] = dmdt

        #[24.01.05]:also return additional particle properties
        if returnMore:
            #[24.01.07]CWO: alpha cannot be returned here, b/c some disk don't have it!
            #
            dMore = {'Rd':Rd, 'St':St, 'v_r':v_r, 'mcp':mcp, 'Hg':disk.Hg} 
            for key in userparL+userevalL:
                dMore[key] = getattr(disk,key)

            return Y2ddt, dMore

        else:
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

        self.locL = Yt[-1,0,:]
        self.massL =Yt[-1,1,:]

        # self.mtotL=Yt[-1,2,:] #no longer updated

        #[24.01.05]CWO
        #after the integration, extract the particle properties
        #for future use
        dum, daux = self.dY2d_dt (Yt[-1], tSpan[-1], gas, returnMore=True)

        for key, val in daux.items():
            setattr(self, key, val)

        return Yt

    
    def remove_particles(self,remove_idx):
        self.mtotL =  np.delete(self.mtotL, remove_idx) 
        self.locL = np.delete(self.locL, remove_idx)
        self.massL = np.delete(self.massL, remove_idx)   
        self.fcomp = np.delete(self.fcomp, remove_idx, axis=0) #[24.01.01] added
        self.num -= len(remove_idx)


    def add_particles (self,Nadd):

        self.locL = np.append(self.locL, self.rout)
        self.massL = np.append(self.massL, self.mini)
        self.mtotL = np.append(self.mtotL, self.mtot1)  #mtot1 is needed here
        self.fcomp = np.append(self.fcomp, [self.fcompini], axis=0) #[24.01.01] added
        self.num += 1 

        if Nadd!=1:
            #[24.01.07]CWO: I dont understand why we need to add 2 particles sometimes?
            print('[Super.add_particles]WARNING:can only add 1 particle // reduce timestep')


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
        self.starttime = time        #time when it appears
        self.mass = mplanet     #its mass
        self.fcomp = fcomp      #its composition
        self.spCrossTime = [0.0]   #list when particles cross
        self.spCrossMass = [0.0]   #corresponding mass
        self.spCrossTau = [-1]   #corresponding stokes number
        self.ncross = 0
        self.planetMassData = [[time, mplanet]]
        self.relp_mass = np.inf
        #self.compData = [[time,fcomp[0]]]
        self.relp_comp = np.nan
        self.max_jumpT = 0.0

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

    def get_icelines_location (self,gas,time,bounds=None,guess=None):
        """
        get location of iceline, whose temperature is assumped as 160K
        """

        if guess is not None:
            dsol = sciop.root_scalar(self.find_iceline, x0=guess, args=(time,gas), 
                            method='secant', rtol=1e-6)
            self.loc = dsol.root
            return

        if bounds==None:
            bounds = (dp.rinn, dp.rout)
        #change the bounds to make it general
        #[23.01.08]LZX: if we change alpha larger,there exist possibility that can't find the iceline location
        #               so for now make a try-exception here, if can't find, then set it to np.nan
        try:
            self.loc = sciop.brentq(self.find_iceline, *bounds, args=(time,gas))
        except:
            self.loc = np.nan
