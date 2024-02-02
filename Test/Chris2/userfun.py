import numpy as np
import matplotlib.pyplot as pl
import cgs
import datetime
import csv
import copy
import physics

def del_v (St, rhoD, disk):
    vr = 2*disk.eta *disk.vK *St/(1+St**2)

    #take half of the velocity...
    return np.abs(vr)/2

def dm_dt(Rd, delv, Hd, sigD):
    """
    the time derivetive of particles's mass, determine how particles grow
    """
    return 2*np.sqrt(np.pi)*Rd**2*delv/Hd*sigD   #eq. 5 of Shibaike et al. 2017


def H_d (St, disk):
    return physics.H_d(disk.Hg, St, disk.alpha) 


def epsilon_PA (PlanetsLoc,PlanetsMass,sp):
    """
    Get the pebble accretion rate
    """

    mus=PlanetsMass/sp.mcp

    Hp = physics.H_d(sp.Hg, sp.St, sp.alpha)
    # Hp=sp.Hg*(1+sp.St/ sp.alpha*(1+2*sp.St)/(1+sp.St))
    hp=Hp/PlanetsLoc

    delv_o_vK=0.52*(mus*sp.St)**(1/3)+sp.eta/(1+5.7*(mus/sp.eta**3*sp.St))
    
    P_eff=1/np.sqrt((0.32*np.sqrt(mus*delv_o_vK/ sp.St/ sp.eta**2))**(-2)+(0.39*mus/ sp.eta/hp)**(-2)) #Liu & Ormel 2018
    return P_eff


def M_critical (eta, St, mcp):

    M_critical=1/8*eta**3*St *mcp #Shibaike 2019
    return M_critical


def init_planets ():
    """
    this initializes the planets 

    return planet properties: 
        - time it appears
        - location
        - mass
        - composition

    [23.12.05]:it's probably better to format this as dictionaries (TBD!)
    """
    #composition should follow what you have defined in parameters
    #and be normalized

    #[23.12.05]:let's worry about composition issues later...
    #fcomp = np.ones_like(pars.composL, dtype=float)
    #fcomp = fcomp /sum(fcomp)

    #return lists for the N-planets we have 
    return [1*cgs.yr, 1e3*cgs.yr], [7*cgs.rJup, 10*cgs.rJup], [3e23, 3e23], [1.0, 1.0]

def init_icelines():
    """
    these quantities are not sure
    """
    return ['H20', 'CO'], [160, 150], [0.5, 0.1]

def init_compos (material):
    dcompos = {}
    if material=='silicates':
        dcompos['Zinit'] = 0.01
    elif material=='H2O':
        dcompos['name'] = 'H2O'
        dcompos['Zinit'] = 0.01
        dcompos['iceline'] = True
        dcompos['rhoint'] = 1.0
        dcompos['iceline_temp'] = 160
        dcompos['iceline_init'] = 2*cgs.RJ

    return dcompos

def do_stuff (system, init=False, final=False):
    #data class is available...
    # import pdb; pdb.set_trace()
    if init:
        #data.data_process(system.particles.locL,system.particles.massL,system.particles.mtotL,system.time,system.daction,system.planetL)
        #initialize your data class
        #data = 
        pass
    else:
        #data.data_process(system.particles.locL,system.particles.massL,system.particles.mtotL,system.time,system.daction,system.planetL)
        #data object should be available...

        data.data_add(system)
        
        tminarr = np.array([ddum['tmin'] for ddum in system.mintimeL])
        isort = tminarr.argsort()

        if len(tminarr)>1:
            tevol = tminarr[isort[1]]
        else:
            tevol = -1

        sfmt = '{:7d} {:5d} {:10.2e} {:10.2e} {:10.2e}'
        line = sfmt.format(system.ntime, len(system.particles.massL), system.deltaT, tevol/cgs.yr, system.time/cgs.yr)
        #print(line)

        if system.time>1.1*data.time_lastplot and system.ntime>100+data.ntime_lastplot or final:
            data.data_plot(copy.copy(system.time), copy.copy(system.ntime))


class Planet(object):

    def __init__(self, time, mass, loc):
        self.massL = [mass]
        self.locL = [loc]
        self.timeL = [time]

    def add_tpoint (self, time, mass, loc):
        self.massL.append(mass)
        self.timeL.append(time)
        self.locL.append(loc)


class Data(object):
    """
    To store data and process data (convert\plot)
    """

    def __init__(self):

        self.timeL=[]
        self.planetL=[]
        self.time_lastplot = 0
        self.ntime_lastplot = 0
        
    def update_cumulative_change(self,daction):
        if 'remove' in daction.keys():
            self.cumulative_change['remove']+=list(daction['remove'])
        if 'add' in daction.keys():
            self.cumulative_change['add']+=daction['add']
    
    def data_add (self,system):
        """
        time: the system.time 
        Y2d: the particles properties' list
        """

        self.timeL.append(system.time)
        for k, planet in enumerate(system.planetL):
            if k<len(self.planetL):
                self.planetL[k].add_tpoint(system.time, planet.mass, planet.loc)
            else:
                self.planetL.append(Planet(system.time, planet.mass, planet.loc))


    def data_plot(self, time, ntime):

        fg, ax = pl.subplots(1,1)
        for k,planet in enumerate(self.planetL):
            rad = np.array(planet.locL)
            mass = np.array(planet.massL)

            ax.loglog(rad, mass, '-')


        
        pl.savefig('test.jpg')
        pl.close(fg)
        self.time_lastplot = time
        self.ntime_lastplot = ntime


    def plot_disk(self,time,disk):
        r_Span=np.linspace(disk.rinn,disk.rout)
        Sigmag=disk.Sigma_g(r_Span,time)
        Td=disk.T_d(r_Span,time)

    def plot_planets_accretion(self,planet,system):
            
        keysL=list(self.radD.keys())

        for i in range(len(self.radD)):
            if planet.time/cgs.RJ<keysL[i]:
                plt.figure()
                plt.ylim((6,28))
                plt.figure('pebble accretion')

                particle_index=np.argwhere(np.isnan(self.radD[keysL[i]])==False)
                pn= particle_index.size
                particles=np.linspace(0,pn,pn)
                sizeL=self.mtotD[keysL[i]][particle_index]/system.mtot1 *5
                plt.scatter(particles, self.radD[keysL[i]][particle_index],s=sizeL,label='totally'+str(pn)+'are in disk')
                
                plt.plot(particles,np.linspace(planet.loc,planet.loc,pn)/cgs.RJ,linewidth=1,label='Planets Location')

                plt.legend(fontsize=5)
                print (pn)
                # import pdb ; pdb.set_trace()
                plt.savefig('/home/lzx/CpdPhysics_Chris/Test/Zhixuan/planets&pebbles/'+str(self.timeL[i])+'.jpg') 
                plt.close()

data = Data() #sets it up
