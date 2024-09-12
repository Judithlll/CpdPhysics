import numpy as np
import subprocess as sp
import shutil
import cv2
import random
import fileio
import matplotlib.pyplot as plt
import cgs
import datetime
import csv
import copy
import parameters as pars     
#import imageio.v2 as imageio
import os
import pandas as pd
import parameters as pars
import disk_properties as dp
import physics
from scipy.optimize import curve_fit
from matplotlib.colors import LinearSegmentedColormap

def PIM():
    return 1.48e40
def Stokes_number(delv, Rd, vth, lmfp, OmegaK, rhog, rhoint = 1.4):
    v_dg=np.abs(delv)
    Rep=4*Rd*v_dg/vth/lmfp
    #following Shibaike, eq....
    CD=24/Rep*(1+0.27*Rep)**0.43+0.47*(1-np.exp(-0.04*Rep**0.38))
    St=8/3/CD*rhoint*Rd/rhog/v_dg*OmegaK
    return St

def magneto_radius (B_mag, rho, dotMg, gas, r, m, t=0.0):
    """
    Calculate the magneto radius of central body. 

    Input parameters:
    B_mag: The magnetic field strength of central body. 
    rho: density of central body. 
    dotMg: the gas mass inflow of disk. 
    gas: The gas object used for getting the temperature. 
    r: the radius of central body. 
    m: the mass of central body.
    t: the time of system evolution.

    Output:
    If the temperature is larger than the ionization temperature, the 
    there exists the cavity, otherwise there's not.
    """
    T_crit = 1e5 #Temporarily use the ionization Temperature of H 

    #standard Lamb ... cavity opening (assumes ionized)
    r_cav = (B_mag**4 *r**(12)/4/cgs.gC/m/dotMg**2)**(1/7)

    Tinn = gas.get_key_disk_properties(r_cav, t)[1]

    if Tinn <= T_crit:
        r_cav = r 
        print('''The inner disk is not hot enough, the r_cav is now set to 
              radius of central body ''')

    return r_cav 

#LZX [24.08.04]: not used now, but put it here for preperation
def St_fragment(alpha, cs, eta, loc, Omega_K, v_crit, icelineloc):
    #following Shibaike 2019 eq14
    sil_part = np.argwhere(np.ndarray.flatten(loc<icelineloc))
    icy_part = np.argwhere(np.ndarray.flatten(loc>=icelineloc))

    sq_part = np.append(np.sqrt(9*alpha**2*cs[sil_part]**4 +4*eta[sil_part]**2*loc[sil_part]**2*Omega_K[sil_part]**2*v_crit['silicates']**2),
                        np.sqrt(9*alpha**2*cs[icy_part]**4 +4*eta[icy_part]**2*loc[icy_part]**2*Omega_K[icy_part]**2*v_crit['icy']**2))
    
    St = (-3*alpha*cs**2+sq_part)/(2*eta**2*loc**2*Omega_K**2)
    return St

def add_planet(system, mode):
    """
    The add planet here should have some options, but both of them can
    apply the add_planet under System, which add planet from a candidate 
    list, so here we should modify the candidate list.

    The options are:
    1. given: the add from the given information
    2. random: the add from the planetesimal generating rate
    """
    #First modify the candidate 
    ##LZX[24.08.27]:maybe check the planetesimal formation models
    
    #Then add the planet to the planet list
    if mode == 'capture':
        #maybe also some newly captured seeds
        #pnew = core.PLANET(time, mcp0, rho)
        #system.planet_candidate.append(pnew) 
        
        if len(system.planet_candidate)>0:
            for planet in system.planet_candidate:
                if planet.starttime <= system.time:
                    system.add_planet(planet)
                    system.planet_candidate.remove(planet)
    elif mode == 'insitu':
        pass

    #Finally we should sort the planets according to the location.

    return system

def planet_migration (gas,planetLoc,planetMass,time,rhopl):
    #aerodynamic drag migration
    out = gas.get_key_disk_properties (planetLoc, time)

    mcp = dp.Mcp_t(time)

    disk = physics.DISK (*out, planetLoc, time, mcp) #pro
    disk.add_auxiliary ()
    disk.add_uservar (dp.user_add_var())    #variables
    disk.add_userfun (dp.user_add_fun())
    disk.add_user_eval (dp.user_add_eval())

    try:
        rinn = disk.rinn
    except:
        rinn = 0.0
    
    if planetLoc < rinn:
        vt1 = 0.0
    else:
        
        #kley-nelson 2011 eq14
        torq0 = (planetMass/mcp)**2 *(disk.Hg/planetLoc)**(-2)*disk.sigmaG *planetLoc**4*disk.OmegaK**2
        p_phi = planetMass*np.sqrt(cgs.gC*mcp*planetLoc)
        
        torq_tot = -(1.36-0.62*disk.beta_sigG(planetLoc)- 0.43*disk.beta_T) *torq0
        vt1 = 2*planetLoc*torq_tot/p_phi
        if vt1 >0:
            import pdb;pdb.set_trace()
    # v_mig=vt1+vd

    #if planetMass > 5e26: import pdb; pdb.set_trace()
    #if vt1< -0.04:
    #    import pdb;pdb.set_trace()
    return vt1 # v_mig

def del_v (St, disk):

    #take half of the velocity...
    return np.abs(disk.v_r)/2


def H_d (St, disk):
    return physics.H_d(disk.Hg, St, disk.alpha) 
    

def dm_dt(particles):
    """
    the time derivetive of particles's mass, determine how particles grow
    if icy fraction == 0, then come to fragmentation (some assumption here)
    """
    #TBD: make a smoother function to do this 
    
    Rd = particles.Rd 
    sigD = particles.sfd
    St = particles.St 
    fcomp = particles.fcomp
    delv = del_v(St, particles)
    Hd = H_d(St, particles)

    vc = pars.vc['silicates']*fcomp[:,0]+ pars.vc['icy']*fcomp[:,1]
    Fcomp = np.ones_like(Rd)
    Fcomp = np.where(delv/vc>1, -1, 1) 

    return Fcomp *2*np.sqrt(np.pi)*Rd**2*delv/Hd*sigD   #eq. 5 of Shibaike et al. 2017

def epsilon_PA (planetLoc,planetMass,cross_p):
    """
    Get the pebble accretion rate
    #"""

    St = cross_p.St
    eta = cross_p.eta

    #prefer to call this qp
    mus = planetMass/cross_p.mcp
    Hp = physics.H_d(cross_p.Hg, St, cross_p.alpha)
    # Hp=cross_p.Hg*(1+cross_p.St/ cross_p.alpha*(1+2*cross_p.St)/(1+cross_p.St))
    hp = Hp/planetLoc

    #expression from eq.18 S+19
    delv_o_vK = 0.52*(mus*St)**(1/3)+eta/(1+5.7*(mus/eta**3*cross_p.St))
    
    P_eff = 1/np.sqrt( (0.32*np.sqrt(mus*delv_o_vK/ St/ eta**2))**(-2)
                      +(0.39*mus/ eta/hp)**(-2)) #Liu & Ormel 2018

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

    [23.12.05]:it's probably better to format this as dictionaries (TBD:!)
    """
    #composition should follow what you have defined in parameters
    #and be normalized

    #[23.12.05]:let's worry about composition issues later...
    #fcomp = np.ones_like(pars.composL, dtype=float)
    #fcomp = fcomp /sum(fcomp)

    #return lists for the N-planets we have 
    timeL = [0.4e6*cgs.yr, 1.5e6*cgs.yr, 1.7e6*cgs.yr, 1.9e6*cgs.yr] 
    #some things wrong with the initial location is set to the out edge
    #about particles number
    locationL = [50*cgs.rJup, 50*cgs.rJup, 50*cgs.rJup, 50*cgs.rJup] 
    massL = [3e23, 3e23, 3e23, 3e23] 
    compoL = np.array([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])

    return timeL, locationL, massL, compoL

def init_compos (material):
    dcompos = {}
    if material=='silicates':
        dcompos['Zinit'] = 0.0008  #cwo made smaller by factor 10
    elif material=='H2O':
        dcompos['name'] = 'H2O'
        dcompos['Zinit'] = 0.0008
        dcompos['iceline'] = True
        dcompos['rhoint'] = 1.0
        dcompos['iceline_temp'] = 160
        dcompos['iceline_init'] = 15*cgs.RJ

    return dcompos

def do_stuff (system, init=False, final= False):
    #data class is available...
    # import pdb; pdb.set_trace()
    if init:
        data.data_process(system)
        data.gas = system.gas
        #initialize your data class
    elif final:
        data.data_process(system)
        data.get_plot_list(doParticles = False)
        #store system components as pickles
        fileio.store_class(system, 'system')
        fileio.store_class(data, 'data')##CWO: this and all stuff below does not seem to be general. Move to do_stuff perhaps

    else:
        data.data_process(system)
        sp.run('tail -n1 log/system_evol.log', shell=True)
        #data object should be available...
        #tminarr = system.minTimes.tminarr 
        #sfmt = '{:6d} {:5d} {:10.2e} {:3d} {:7.2f} {:22s} {:7.2e}'
        #line = sfmt.format(system.ntime, system.particles.num, system.deltaT, tminarr.argmin(), system.time/cgs.yr, nameM, tmin)
        #print(line)  

def plot_massfit(planetMassData):
    def mass_fit(t,a,b):
        m = a*t+b
        return m 

        
    plt.figure()
    if len(planetMassData)>2:
        timedots=np.log([planetMassData[j][0] for j in range(len(planetMassData))])
        massdots=np.log([planetMassData[j][1] for j in range(len(planetMassData))])
        popt, pcov = curve_fit(mass_fit, timedots, massdots)
        plt.scatter(timedots, massdots)
        t_list=np.linspace(timedots[0], timedots[-1], 30)
        plt.plot(t_list, mass_fit(t_list, *popt))
    plt.savefig('/home/lzx/CpdPhysics/Test/Zhixuan/plot/test.jpg')
    plt.close()

   

class Data(object):
    """
    To store data and process data (convert\plot)
    """

    def __init__(self):

        self.timeL=[]
        self.radD={}
        self.mD={}
        self.mtotD={}
        self.fcompD={}
        self.v_rD = {}
        self.StD = {}
        self.sfd = {}
        self.cumulative_change={'remove':[],'add':0}
        self.planetsmass = {}
        self.planetsloc = {}
        self.icelinesloc = {}
        self.jumpstuff = []
        self.planetsfcomp = {}

        self.planetsdmdt = {}

        self.dactionD ={}
        
    def update_cumulative_change(self,daction):
        for daction in self.dactionD.values():
            if 'remove' in daction.keys():
                self.cumulative_change['remove']+=list(daction['remove'])
            if 'add' in daction.keys():
                self.cumulative_change['add']+=daction['add']
    
    def data_process(self, system):
        """
        time: the system.time 
        Y2d: the particles properties' list
        """
        
        time = system.time
        
        # store time series 
        self.timeL.append(time)
        
        # store particles data
        self.radD.setdefault(time,system.particles.locL)
        self.mD.setdefault(time, system.particles.massL)
        self.mtotD.setdefault(time, system.particles.mtotL)
        self.fcompD.setdefault(time, system.particles.fcomp)
        self.sfd.setdefault(time, system.particles.sfd)
        self.dactionD.setdefault(time, system.daction)
        try: 
            self.v_rD[time]=system.particles.v_r
            self.StD[time] = system.particles.St 
        except:
            self.v_rD[time] =np.array([])
            self.StD[time] =np.array([])

        #store palnets' data
        if pars.doPlanets:
            # for now we assume the planets are in order 
            try:
                lengt = max(system.planetD.keys()) +1
            except:
                lengt = 0
            pmassL = [np.nan] * lengt
            plocL  = [np.nan] * lengt
            pcompL = [[np.nan]*2] *lengt #TBD:not general now 
            pdmdtL = [np.nan] * lengt
            for k,v in system.planetD.items():
                pmassL[k] = v.mass
                plocL [k] = v.loc
                pcompL[k] = v.fcomp
                pdmdtL[k] = v.dmdt

            self.planetsmass.setdefault(system.time,pmassL)
            self.planetsloc.setdefault(system.time,plocL )
            self.planetsfcomp.setdefault(system.time, pcompL)
            self.planetsdmdt.setdefault(system.time, pdmdtL)
        #store icelines' data
        if pars.doIcelines:
            self.icelinesloc.setdefault(system.time, [iceline.loc for iceline in system.icelineL])
        
        #if jump is done, then shore something about jump
        if system.doJump:
            stuff = {'njump': system.njump, 'njumptime': system.njumptime, 'jumptime':system.time-system.jumpT, 'jumpT': system.jumpT, 'jump_limitation':system.jump_limitation}
            self.jumpstuff.append(stuff)

        if 'remove' in system.daction.keys():
            self.cumulative_change['remove'].append(system.daction['remove'])
        if 'add' in system.daction.keys():
            self.cumulative_change['add']+=system.daction['add']


    def get_plot_list(self, doParticles = False):
        """
        process data and make it in order

        doParticles: decide to process particles' data or not
        """
        # if particles ware removed, then put a np.nan into this location. 
        # for now, I can just use the number of the removed particles to add the np.nan to the head of every list. I think more reasonable is to put in the np.nan according to the removed index, but it seems too complex.
        if doParticles:
            self.num = [len(v) for v in self.radD.values()]
            complen = len(self.fcompD[self.timeL[0]][0])
            for tt, daction in self.dactionD.items():
                idx = self.timeL.index(tt)
                if 'remove' in daction.keys():
                    added_L = np.array([np.nan]*len(daction['remove']))
                    fcomp_adL = np.array([[np.nan]*complen]*len(daction['remove']))
                    for re_time in self.timeL[idx:]:
                        self.StD[re_time] = np.append(added_L, self.StD[re_time])
                        self.v_rD[re_time] = np.append(added_L, self.v_rD[re_time]) 
                        self.radD[re_time] = np.append(added_L, self.radD[re_time])
                        self.mD[re_time] = np.append(added_L, self.mD[re_time])
                        self.mtotD[re_time] = np.append(added_L, self.mtotD[re_time])
                        self.fcompD[re_time] = np.append(fcomp_adL, self.fcompD[re_time],0)
               
                if 'add' in daction.keys():
                    added_L = np.array([np.nan]*daction['add'])
                    fcomp_adL = np.array([[np.nan]*complen]*daction['add'])
                    for ad_time in self.timeL[:idx]:
                        self.StD[ad_time] = np.append(self.StD[ad_time], added_L)
                        self.v_rD[ad_time] = np.append(self.v_rD[ad_time], added_L)
                        self.radD[ad_time] = np.append(self.radD[ad_time], added_L)
                        self.mD[ad_time] = np.append(self.mD[ad_time], added_L)
                        self.mtotD[ad_time] = np.append(self.mtotD[ad_time], added_L)
                        self.fcompD[ad_time] = np.append(self.fcompD[ad_time], fcomp_adL, 0)
            
            self.radL=np.array(list(self.radD.values()))
            self.mL=  np.array(list(self.mD.values()))
            self.mtotL=np.array(list(self.mtotD.values()))
            self.fcompL = np.array(list(self.fcompD.values()))
        
                

        if False:
            if len(self.cumulative_change['remove'])>0:
                rL=np.insert(locL,0,np.full(len(self.cumulative_change['remove']),np.nan))
                mL=np.insert(massL,0,np.full(len(self.cumulative_change['remove']),np.nan))
                mtL=np.insert(mtotL,0,np.full(len(self.cumulative_change['remove']),np.nan))
               
                insertv = np.full((len(self.cumulative_change['remove']),len(fcompL[0])), np.nan)
                fL=np.append(insertv, fcompL, axis =0)
            else:

                rL=locL
                mL=massL
                mtL=mtotL
                fL = fcompL

            max_len=max(len(v) for v in self.radD.values())

            # want to make the data dict in the same length
            for k,v in self.radD.items():
                
                #if particles were added, then add np.nan to the former list to make the length the same.
                if len(v)< max_len:
                    self.radD[k]=np.pad(v, (0, max_len - len(v)), constant_values=np.nan)
                    self.mD[k]=np.pad(self.mD[k], (0, max_len - len(v)), constant_values=np.nan)
                    self.mtotD[k]=np.pad(self.mtotD[k], (0, max_len - len(v)), constant_values=np.nan)
                    
                    apdv = np.full((max_len-len(v), len(fcompL[0])), np.nan)
                    self.fcompD[k] = np.append(self.fcompD[k], apdv, axis = 0)


        pmL = copy.deepcopy(list(self.planetsmass.values()))
        plL = copy.deepcopy(list(self.planetsloc.values()))
        pdmdtL = copy.deepcopy(list(self.planetsdmdt.values()))
        pfL = copy.deepcopy(list(self.planetsfcomp.values()))
        max_len = max([len(l) for l in pmL])
        
        # maybe can use the number of planet
        
        for i in range (len(pmL)):
            if len(pmL[i]) != max_len:
                length = max_len - len(pmL[i])
                pmL[i].extend([np.nan]*length)
                plL[i].extend([np.nan]*length)
                pdmdtL[i].extend([np.nan]*length)
                pfL[i].extend([[np.nan]*2]*length)
        self.planetsmassL = np.array(pmL)
        self.planetslocL = np.array(plL)
        self.planetsfcompL = np.array(pfL)
        self.planetsdmdtL = np.array(pdmdtL)

        self.icelineslocL = np.array(list(self.icelinesloc.values()))


    def plot_stuff(self):
        
        #[time,loc]=np.meshgrid(self.timeL,np.linspace(pars.dgasgrid['rinn']/cgs.RJ,pars.dgasgrid['rout']/cgs.RJ,len(self.timeL)))
        #sigmag=gas.get_key_disk_properties(loc,time)[0]

        # import pdb; pdb.set_trace()
        plt.figure(figsize=(12,18))
        plt.subplot(211)
        plt.title('Particles Location')
        plt.xlabel('time [year]')
        plt.ylabel('Distance from central planet')
        
        

        plt.subplot(212)
        plt.title('Particles Mass')
        plt.xlabel('time [year]')
        plt.ylabel('Mass [g]')   


        for i in range(len(self.radL[0])):
            plt.subplot(211)
            plt.plot(np.array(self.timeL)/cgs.yr ,self.radL[:,i]/cgs.RJ)
            plt.subplot(212)
            plt.plot(np.array(self.timeL)/cgs.yr ,self.mL[:,i])
        
        plt.subplot(211)
        #plt.contourf(time,loc,sigmag,alpha=0.3)
        #plt.colorbar()

        plt.savefig('./plot/test.jpg')

    def plot_disk(self,time):
        r_Span=np.linspace(1*cgs.RJ,pars.dgasgrid['rout'])
        plt.figure(figsize=(16,6))
        if type(time) == float or type(time) == np.float64:
            Sigmag,Td=self.gas.get_key_disk_properties(r_Span,time)[0:2]
            plt.subplot(121)
            plt.ylabel('Surface density $[g/cm^2]$', )
            plt.xlabel('Location [$R_J$]')
            plt.xscale('log')
            plt.yscale('log')
            plt.plot(r_Span/cgs.RJ,Sigmag,label=str(time/cgs.yr))
            plt.subplot(122)
            plt.ylabel('Midplane Temperature $[K]$')
            plt.xlabel('Location [$R_J$]')
            plt.xscale('log')
            plt.yscale('log')
            plt.plot(r_Span/cgs.RJ,Td,label=str(time/cgs.yr))
            plt.legend()
            plt.savefig('./plot/diskproperties.jpg')
        else:
            plt.subplot(121)
            plt.ylabel('Surface density $[g/cm^2]$')
            plt.xlabel('Location [$R_J$]')
            plt.xscale('log')
            plt.yscale('log')
            
            plt.subplot(122)
            plt.ylabel('Midplane Temperature $[K]$')
            plt.xlabel('Location [$R_J$]')
            plt.xscale('log')
            plt.yscale('log')
            
            labels=[r'$0yr$', r'$3\times 10^6yr$', r'$1\times 10^7yr$', r'$3\times 10^7yr$']
            for i in range(len(time)):
                Sigmag, Td = self.gas.get_key_disk_properties(r_Span, time[i])[0:2]
                plt.subplot(121)
                plt.plot(r_Span/cgs.RJ,Sigmag,label=labels[i])
                plt.subplot(122)
                plt.plot(r_Span/cgs.RJ,Td,label=labels[i])

            plt.subplot(122)
            plt.axhline(160, label='iceline', color = 'skyblue', linestyle = 'dashed')
            plt.legend()
            plt.savefig('./plot/diskproperties.jpg')
            plt.close()
    
    def plot_particles_number(self):
        plt.figure()
        plt.title('Particles number')
        plt.xlabel('time [yr]')
        plt.ylabel('number')
        self.num = [len(v) for v in self.radD.values()]
        plt.plot(np.array(self.timeL)/cgs.yr, self.num, color = 'red')
        plt.savefig('./plot/particleNum.jpg')
        plt.close()

    def plot_peff_log(self,logpath='./log/'):
        l=0
        planet_num = 4
        timeL = []
        data = []
        massL = [[],[],[],[]]
        peffL = [[],[],[],[]]
        with open(logpath+'planets.log', 'r') as file:

            columns = file.readline().strip().split()
            # Read the rest of the lines
            while True:
                line = file.readline()
                if not line:
                    break  # End of file
                lm = line.strip().split()
                
                timeL.append(float(lm[0]))

                num = (len(lm) -1.)/5.

                for j in range(int(num)):
                    massL[j].append(float(lm[4+5*j]))
                    peffL[j].append(float(lm[5+5*j]))
                for k in range(int(num),planet_num):
                    massL[k].append(np.nan)
                    peffL[k].append(np.nan)

                data.append(line.strip().split())    

        plt.figure()
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('Satellite mass [g]')
        plt.ylabel('Accretion effieiency [%]')
#yticks_values = [0.0001, 0.01, 0.02, 0.03, 0.04, 0.05]
#yticks_labels = ['0.0001', '', '', '', '', '0.05']
#plt.yticks(yticks_values, yticks_labels)

        xticks_values = [1e23, 1e24, 1e25, 1e26]
        xticks_labels = ['$1$', '1e24', '1e25', '1e26']
#plt.xticks(xticks_values, xticks_labels)
        for i in range(planet_num):
            plt.plot(massL[i], np.array(peffL[i])*100)
        plt.savefig('./plot/peff.jpg')
        plt.close()

    def plot_growth_timescale(self):
        plt.figure()
        plt.xlabel('Satellite mass [g]')
        plt.ylabel('Growth timescale [yr]')
        plt.yscale('log')
        plt.xscale('log')
        idx=[]
        for j in self.jumpstuff:
            idx.append(j['njumptime'])

        for i, dmdtL in enumerate(self.planetsdmdtL.T):
            plt.plot(self.planetsmassL.T[i][idx], self.planetsmassL.T[i][idx]/dmdtL[idx]/cgs.yr, label = 'Satellite'+str(i+1))

        plt.legend()
        plt.savefig('./plot/tgrowth.jpg')
        plt.close()

    def plot_sfd(self,locL,sfd,time,imin,deltaT):
        plt.figure(figsize=(12,8))
        plt.subplot(211)
        plt.xlim(0,100)
        plt.ylim(0.001, 50)
        plt.yscale('log')
        plt.title('Surface density profile at {:.2f}yr'.format(time/cgs.yr))
        plt.plot(locL/cgs.RJ, sfd, 'x-', label=str(imin)+'\n'+'{:.2f}'.format(deltaT))
        plt.scatter(locL[imin[1]]/cgs.RJ, sfd[imin[1]], c= 'red')
        plt.axvline(5.89, linestyle='dashed', color='black', linewidth = 1)
        plt.legend(loc='lower right')

        plt.subplot(212)
        plt.xlim(time/cgs.yr-500,time/cgs.yr+500)
        plt.xticks([time/cgs.yr], ['{:.2f}'.format(time/cgs.yr)])
        plt.ylim(1,1e8)
        plt.yscale('log')
        plt.plot(np.array(self.timeL)/cgs.yr, np.append(np.diff(self.timeL), deltaT) )
        plt.scatter(time/cgs.yr, deltaT, c='red')

        plt.savefig('./sfdevol/{:.2f}.png'.format(time))
        plt.close()

    def plot_planets_accretion(self,planet,system):
            
        keysL=list(self.radD.keys())

        for i in range(len(self.radD)):
            if planet.time/cgs.RJ<keysL[i]:
                plt.figure(figsize=(6,15))
                plt.ylim((6,28))
                # plt.yticks([5,10,15,20,25,30])
                plt.title('pebble accretion')

                particle_index=np.argwhere(np.isnan(self.radD[keysL[i]])==False)
                pn= particle_index.size
                particles=np.linspace(pn,pn,pn)
                lines=np.linspace(0,2*pn,pn)
                sizeL=self.mtotD[keysL[i]][particle_index]/system.mtot1 *0.3
                plt.scatter(particles, self.radD[keysL[i]][particle_index],s=6,c=sizeL,cmap='rainbow',label='totally'+str(pn)+'are in disk')
                for iceline in system.icelineL:
                    plt.plot(lines, iceline.loc*np.ones(pn)/cgs.RJ,linewidth=3,label='water iceline')
                    plt.plot(lines, iceline.loc*np.ones(pn)/cgs.RJ,linewidth=3,label='CO iceline')
                
                for planet in system.planetL:
                    plt.plot(lines, planet.loc*np.ones(pn)/cgs.RJ,linewidth=3,label='Planets Location')

                plt.legend(fontsize=12)
                print (pn)
                # import pdb ; pdb.set_trace()
                plt.savefig('./plot/planets&pebbles/'+str(self.timeL[i])+'.jpg') 
                plt.close()
    
    def plot_deltaT(self):
        plt.figure()
        plt.title('Time step variation')
        plt.xlabel('System Time [yr]')
        plt.ylabel('Time Step [yr]')
        plt.plot(self.timeL[:-1],np.diff(self.timeL))
        plt.savefig('./plot/Delta_t.jpg')
        plt.close()

    def plot_satepart(self):
        folder_path ='./plot/satepart/'

        #check if the folder exists 
        if os.path.exists(folder_path):
            # if exists, clear this folder
            for filename in os.listdir(folder_path):
                file_path = os.path.join(folder_path, filename)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print('failed to delete %s. reason: %s' % (file_path, e))
        else:
            #if not exists, create
            os.makedirs(folder_path)

        mass_0 = 1e-1
        Pmass_0 = 1e23

        tst1 = init_planets()[0][0]
        ik = np.argmin(np.abs(np.array(self.timeL)- tst1))
        time = self.timeL[ik]
        while ik<=len(self.timeL):
            time = self.timeL[ik]

            locL = np.array([5.89, 9.38, 15.0, 26.3])*cgs.RJ
            massL = np.array([0.893, 0.480, 1.48, 1.08])*1e26 
            #fcompL = np.array([0., 0.08, 0.47, 0.52])
            nameL = ['I', 'E', 'G', 'C']


            fig, ax1 = plt.subplots(figsize=(10,9))
            fig.subplots_adjust(bottom=0,top=0.95)
            #fig.subplots_adjust(right=0.95)
            ax1.set_title('Time: {:.2f}Myr'.format(time/cgs.yr/1e6),loc='left')
            ax1.set_xlabel('Location $[R_J]$')
            ax1.set_ylabel(r'Satellites Mass $[M_{\oplus}]$')
            ax1.set_ylim(1.e23/cgs.Mea, 1.e27/cgs.Mea)
            ax1.set_yscale('log')
            ax1.set_xscale('log')
            

            #define the colorbar 
            cmap = LinearSegmentedColormap.from_list(
                'custom_cmap', [(0,'gray'), (0.35,'gray'), (1,'blue')], 
                N=256
            )
            #plot the Gealian Satellites 
            dotsize = 16*(massL/Pmass_0)**(1/3)
            for i,loc in enumerate(locL):
                ax1.axvline(loc/cgs.RJ, color = 'black', linestyle='dashed', alpha = 0.3)
                ax1.text(loc/cgs.RJ, massL[i]/cgs.Mea, nameL[i], fontsize=20, ha='left', va='center')
            #sca=ax1.scatter(locL/cgs.RJ, massL/cgs.Mea, s = dotsize, c =fcompL, cmap =cmap ,vmin=0.0,vmax=0.5, alpha =1)
            #ax1.scatter(locL/cgs.RJ, massL/cgs.Mea, s = 3*dotsize, facecolor ='none', edgecolor='black' )

            #radL = self.radD[time]
            #redgeL = np.array([dp.rinn]+ [np.sqrt(radL[i]*radL[i+1]) for i in range(len(radL)-1)]+[dp.rout])
            #get the number of dots, and normalized by the minimum number
            #numL =np.array((self.mtotD[time]/self.mD[time]))
            #numL_nom = ((numL/numL.min())**(1/2)).astype(int) *100
            #for i,rad in enumerate(radL):
                #the particles number is too large...
            #    random_loc_parti = np.random.uniform(redgeL[i], redgeL[i+1],numL_nom[i])
            #    vert_parti = 10**(np.linspace(23,27,numL_nom[i]))
            #    plt.scatter(random_loc_parti/cgs.RJ,vert_parti, s=(self.mD[time][i]/mass_0)**(1/3),color='gray',alpha=0.3, edgecolors='none')
                #ax.axvline(rad/cgs.RJ, linewidth= (self.mD[time][i]/mass_0)**(1/3), color = 'gray', linestyle='dotted', alpha =0.3)
            #import pdb;pdb.set_trace()

            #2ed y axis
            ax2 = ax1.twinx()
            ax2.set_ylabel('Stokes number')
            ax2.set_ylim(1e-6,1)
            ax2.set_xscale('log')
            ax2.set_xlim(4.,100.)
            ax2.set_yscale('log')
            ax2.plot(self.radD[time]/cgs.rJup, self.StD[time], color='black',alpha=0.3, linewidth=1)
            ax2.scatter(self.radD[time]/cgs.RJ, self.StD[time], s=np.where(self.radD[time]>self.icelinesloc[time][0], 4,2), c=np.where(self.radD[time]>self.icelinesloc[time][0], 'blue', 'gray'))
            ax2.set_xticks([dp.rinn/cgs.RJ,self.icelinesloc[time][0]/cgs.RJ, 10.,20.,50.],
                           ['{:.2f}'.format(dp.rinn/cgs.RJ),
                            '{:.2f}'.format(self.icelinesloc[time][0]/cgs.RJ), '','20','50'])

            ax1.tick_params(axis='both', which='major', length=7, width=1) 
            #plot the special locations in the disk:[inner egde, iceline]
            ax1.axvline(dp.rinn/cgs.RJ, color = 'black', linewidth = 1, linestyle= 'dotted', label = 'inner edge')
            ax1.axvline(self.icelinesloc[time][0]/cgs.RJ, color = 'blue', linestyle = 'dashed',label = 'Iceline')
            
            #to make a legend and colorbar with a totoally transparent dot 
            sca=ax1.scatter(0,0, c=0.1, cmap=cmap,vmin=0.0,vmax=0.5)
            if not np.isnan(self.planetsloc[time]).all():
                dotsize = 16*(np.array(self.planetsmass[time])/Pmass_0)**(1/3) # need to be modified
                sca=ax1.scatter(np.array(self.planetsloc[time])/cgs.RJ, np.array(self.planetsmass[time])/cgs.Mea, s = dotsize, c =np.array(self.planetsfcomp[time])[:,1], cmap =cmap ,vmin=0.0,vmax=0.5, alpha =1)


            fig.colorbar(sca, label=r'Water Fraction [%]',orientation='horizontal',pad =0.1)
            ax1.legend(loc='upper right')
            fig.savefig('./plot/satepart/{}.jpg'.format(time/cgs.yr))
            plt.close()

            ik+=100 
            print([ik,time/cgs.yr])

    def plot_St_t(self):
        folder_path ='./plot/St_t/'

        #check if the folder exists 
        if os.path.exists(folder_path):
            # if exists, clear this folder
            for filename in os.listdir(folder_path):
                file_path = os.path.join(folder_path, filename)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print('failed to delete %s. reason: %s' % (file_path, e))
        else:
            #if not exists, create
            os.makedirs(folder_path)

        ik = np.argmin(np.abs(np.array(self.timeL)-1e6*cgs.yr))
        time = self.timeL[ik]
        while ik <= len(self.timeL):
            time = self.timeL[ik]
        #for time in self.timeL:
            plt.figure(figsize=(8,6))
            plt.title('Time: {:.2f}'.format(time/cgs.yr),loc='left')
            plt.xlabel('Location $[R_J]$')
            plt.ylabel('Stokes number')
            plt.ylim(1e-6,1)
            plt.xlim(4.,100.0)
            plt.yscale('log')
            plt.xscale('log')

            plt.plot(self.radD[time]/cgs.rJup, self.StD[time], '.-',label = "{:.2f}".format(time/cgs.yr))
            
            plt.axvline(self.icelinesloc[time][0]/cgs.RJ, linestyle = 'dashed',color='gray')

            plt.savefig('./plot/St_t/{:2f}.png'.format(time/cgs.yr))
            plt.close()

            ik+=100
            print('Plotting St:',[ik,time/cgs.yr])
        
    def plot_planet_migration(self):
        plt.figure()
        plt.title('Planet migration')
        plt.xlabel('Planets location [$R_{Jup}$]' )
        plt.ylabel('System time [yr]')
        plt.yscale('log')
        plt.xscale('log')
        loclist = self.planetslocL.T
        time = np.array(list(self.planetsloc.keys()))

        for jump in self.jumpstuff:
            plt.axhspan(jump['jumptime']/cgs.yr, (jump['jumptime']+jump['jumpT'])/cgs.yr, alpha = 0.3)
        for i,loc in enumerate(loclist):
            plt.plot(loc/cgs.RJ, time/cgs.yr, label = 'Satellite'+str(i+1))


            plt.axhline((jump['jumptime']+jump['jumpT'])/cgs.yr, color = 'green', linewidth = 0.2)
        plt.axvline(dp.rinn/cgs.RJ, color = 'gray', linewidth = 0.5, label = 'inner edge')
        plt.legend()
        plt.savefig('./plot/planet_migration.jpg',dpi=600)
        plt.close()   

    def plot_planet_evolution(self):
        plt.figure()
        plt.title('Satellites evolution')
        plt.xlabel('Distance from the planet [$R_{J}$]' )
        plt.ylabel('Evolution Time [yr]')
        plt.yscale('log')
        plt.xscale('log')
        loclist = self.planetslocL.T
        masslist = self.planetsmassL.T
        time = np.array(list(self.planetsloc.keys()))
        
        planetst = 0.8e6*cgs.yr
        plt.ylim(planetst/cgs.yr,20e6)
        stidx = np.argwhere(time>planetst)[0][0]
        
        dotssize = masslist/np.nanmin(masslist)*0.1
        #cmap = LinearSegmentedColormap.from_list("my_colormap", ["y", "royalblue"])

        for jump in self.jumpstuff:
            plt.axhspan((jump['jumptime'])/cgs.yr, 
                        (jump['jumptime']+jump['jumpT'])/cgs.yr, alpha = 0.3)
        for i,loc in enumerate(loclist):
            plt.scatter(loc[stidx:]/cgs.RJ, (time[stidx:])/cgs.yr, s = dotssize[i][stidx:], c =self.planetsfcompL[stidx:,i][:,1], cmap ='Spectral', alpha =1 )

        for i,loc in enumerate(loclist):
            plt.plot(loc/cgs.RJ, (time)/cgs.yr, color='black', alpha=0.5,linewidth=0.5)

            #plt.axhline((jump['jumptime']+jump['jumpT'])/cgs.yr, color = 'green', linewidth = 0.2)
        plt.colorbar(label = "Water Fraction")
        plt.axvline(dp.rinn/cgs.RJ, color = 'gray', linewidth = 1, label = 'inner edge')
        plt.plot(self.icelineslocL[stidx:,0]/cgs.RJ, (time[stidx:])/cgs.yr, color = 'blue', linestyle = 'dashed',label = 'Iceline')
        plt.legend()
        plt.xticks([5.89,10,20,50],['5.89','10','20','50'])
        plt.savefig('./plot/planet_evolution.jpg',dpi=600)
        plt.close()
   

    def plot_iceline(self):
        plt.figure()
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel('iceline location')
        plt.xlabel('time')
        plt.xlim(1e3,20e6)
        plt.yticks([14.8,15],['14.8','15'])
        plt.axhline(15.)
        plt.plot(np.array(self.timeL)/cgs.yr, self.icelineslocL.T[0]/cgs.RJ, 'x-')
        plt.savefig('./plot/iceline.jpg')
        plt.close()

    def plot_planet_massloc(self):
        plt.figure()
        plt.title('PLanet mass-location')
        plt.xlabel('Planets location [$R_J$]' )
        plt.ylabel('Planets mass [$R_J$]')
        loclist = self.planetslocL.T
        masslist = self.planetsmassL.T
        
        for i,loc in enumerate(loclist):
            plt.plot(loc/cgs.RJ, masslist[i]/cgs.MJ, label = 'planet'+str(i+1), linewidth = 3)
        import pdb;pdb.set_trace()
        #for jump in self.jumpstuff:
        #    plt.axhspan(jump['jumptime']/cgs.yr, (jump['jumptime']+jump['jumpT'])/cgs.yr, alpha = 0.5)
        #    plt.axhline((jump['jumptime']+jump['jumpT'])/cgs.yr, color = 'green', linewidth = 0.2)

        plt.axvline(dp.rinn/cgs.RJ, color = 'gray', linewidth = 0.5, label = 'inner edge')
        plt.legend()
        plt.savefig('./plot/planet_massloc.jpg',dpi=600)
        plt.close()
    
    def plot_disk_profile(self):
        """
        a little annoying, let's leave it here for now TBD 
        """

        locL =np.linspace(dp.rinn, dp.rout, 200)
        tL = np.array([50, 1e6, 3e6]) *cgs.yr
        colors = ['r', 'b','g']
        pld = {tL[i]:colors[i] for i in range(3)}
        fig, ax1 = plt.subplots()
        plt.title ('Disk Properties')

        ax1.set_ylabel('Surface Density $[g/cm^2]$')
        ax1.set_xlabel('Location $[R_J]$')
        ax2 = ax1.twinx()
        ax2.set_ylabel('Temperature [K]')
        #plot the surface density
        for t,c in pld.items():

            sigmaGL,temp = self.gas.get_key_disk_properties(locL,t)[0:2]
            ax1.plot(locL/cgs.RJ, sigmaGL, label = 'Surface density', color = c, linestyle = 'dotted')

            ax2.plot(locL/cgs.RJ, temp, label = 'Temperature', color = c, linestyle = 'dashed')
        
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        l1 = [lines[0].set_color('black'), lines[1].set_color('black')]
        la1 = ['Surface density', 'Temperature']
        la2 = ['50', '1e6', '3e6']
        import pdb;pdb.set_trace() 

        ax2.legend(l1+l2, la1+la2, loc='upper right')
        plt.savefig('./plot/disk_properties.jpg')
        plt.close()
    
    def plot_pebble_Sigma(self, tL, mode = 'grid'):
        timeL =tL
        
        plt.figure()
        plt.title('pebble surface density')
        
        if mode == 'grid':

            grids = np.logspace(np.log10(dp.rinn), np.log10(dp.rout), 10000)
            width = np.diff(grids)

            plt.yscale('log')
            plt.xscale('log')
            plt.ylim(1e-3,10)
            for time in timeL:
                tidx = np.argmin(np.abs(self.timeL-time))
                ti = self.timeL[tidx]
                mtot = self.mtotD[ti]
                loc = self.radD[ti] #ordered
                sfd = self.sfd[ti]

                boundaries = np.sqrt(loc[1:]*loc[:-1])
                boundaries = np.append(dp.rinn,boundaries)
                boundaries = np.append(boundaries,dp.rout)
                warr = np.diff(boundaries)
                sigma = mtot /(2*np.pi*loc*warr)
                    
                #plt.plot(loc/cgs.rJup, sigma, 'x-', label = str('{:7.2f}'.format(time/cgs.yr)))
                #dotMd = dp.dot_Mg(ti)*dp.ratio
                #v_r = self.v_rD[ti]


                #sigmaP = dotMd/(2*np.pi*loc*(-v_r))/len(fcomp[0])*np.count_nonzero(fcomp,axis=1)
                plt.plot(loc/cgs.rJup, sfd, label = 'particles'+"{:.2f}".format(time/cgs.yr) )
                #sigmaP = np.array([])
                #for i in range(len(loc)-1):
                #    invo_idx = np.argwhere((grids>loc[i]) &(grids<loc[i+1])) 
                #    
                #    invo_grids = grids[invo_idx]
                #    invo_width = width[invo_idx]
                #    
                #    m_split= mtot[i]/len(invo_grids)
                #    
                #    sigmaP =np.append(sigmaP, m_split/(2*np.pi*invo_grids*invo_width))
                #plt.plot(grids[-len(sigmaP):]/cgs.rJup, sigmaP, label = str(ti/cgs.yr))
                #for r in loc:
                #    plt.axvline(r/cgs.rJup)
            plt.legend()
            plt.savefig('./plot/sigmaP.jpg')

        elif mode == 'particle':
            plt.yscale('log')
            for time in timeL:
                tidx = np.argmin(np.abs(self.timeL-time))
                ti = self.timeL[tidx]

                dotMd = dp.dot_Mg(ti)*dp.ratio
                loc = self.radD[ti]
                v_r = self.v_rD[ti]


                sigmaP = dotMd/(2*np.pi*loc*(-v_r))
                plt.plot(loc/cgs.rJup, sigmaP, label = "{:.2f}".format(ti/cgs.yr) )
                plt.legend()
                plt.savefig('./plot/sigmaP.jpg')
    
    def plot_St(self, tL, savepath='./plot/'):
        timeL =tL
        plt.figure()
        plt.title('Stokes number')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e-6, 1)
        plt.xlabel('Distance from the planet [$R_J$]')
        plt.ylabel('Stokes number')
        color = ['blue', 'm', 'green']
        labels=[r'$4\times 10^6  yrs$',r'$5\times 10^6  yrs$',r'$6\times 10^6  yrs$']
        for i,time in enumerate(timeL):
            tidx = np.argmin(np.abs(self.timeL-time))
            ti = self.timeL[tidx]
            
            St = self.StD[ti]
            loc = self.radD[ti]
            
            #q = -dp.beta_T
            #loc_new = np.linspace(1*cgs.RJ,51*cgs.RJ, 1000)
            #p = -dp.beta_sigG(loc_new)
            #sigG, T = self.gas.get_key_disk_properties(loc_new, time)[0:2]
            #St_new = 0.23*(2/(3+2*p+q))**(4/5)*(10/(18-39*q))**(2/5)*(dp.ratio/0.003)**(2/5)*(dp.alpha/1e-4)**(1/5)*(T/160)**(-2/5)*(dp.Mcp_t(time)/cgs.MJ)**(2/5)*(loc_new/10/cgs.RJ)**(-2/5)
            
            plt.plot(loc/cgs.rJup, St, '.-',label = "{:.2f}".format(time/cgs.yr))
            #plt.plot(loc_new/cgs.rJup, St_new, '--',label = "{:.2f}_powerlaw".format(ti/cgs.yr))
            #for i,loc in enumerate(self.planetsloc[ti]): 
            #    plt.axvline(loc/cgs.rJup, linestyle='dashed', color='gray')
            
            plt.axvline(self.icelinesloc[ti][0]/cgs.RJ, linestyle = 'dotted')
        plt.legend()
        plt.savefig(savepath+'St{:.2f}.jpg'.format(time/cgs.yr))
        plt.close()

    def plot_vr(self, tL):
        timeL =tL
        plt.figure()
        plt.title('Raidal velocity')
        plt.xscale('log')
        for time in timeL:
            tidx = np.argmin(np.abs(self.timeL-time))
            ti = self.timeL[tidx]
            
            vr = np.abs(self.v_rD[ti])
            loc = self.radD[ti]
            
            plt.plot(loc/cgs.rJup, vr, label = "{:.2f}".format(time/cgs.yr))
        #plt.axhline(pars.vc['icy'], label = 'Fragmentation velocity [icy]')
        plt.axhline(pars.vc['silicates'], label = 'Fragmentation velocity [sil]')
            
        plt.legend()
        plt.savefig('./plot/v_r.jpg')
        plt.close()


    def plot_jumpT(self):
        plt.figure()
        plt.xlabel('System time [yr]')
        plt.ylabel('Jump time [yr]')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(7e5,1e7)
        for jump in self.jumpstuff:
            #plt.axvspan((jump['jumptime'])/cgs.yr, 
            #           (jump['jumptime']+jump['jumpT'])/cgs.yr, alpha = 0.3,color = 'gray')
            plt.axvline((jump['jumptime']+jump['jumpT'])/cgs.yr, color = 'black', alpha =0.5, linewidth = 0.3, linestyle='dotted')
        for i, st in enumerate([0.8e6, 1.1e6, 1.4e6, 1.7e6]):
            plt.axvline(st, color='black', linestyle='dashed', label = 'Seed'+str(i+1), linewidth=2)
        jumpTlist = [f['jumpT']/cgs.yr for f in self.jumpstuff]
        timelist = [f['jumptime']/cgs.yr for f in self.jumpstuff]
        plt.plot(timelist, jumpTlist, '.-',label='jump time', color = 'blue')
        #plt.scatter(timelist, jumpTlist)
        plt.legend()
        plt.savefig('./plot/jumpT.jpg')
        plt.close()
        

def make_animation(mp4name, path='./plot/satepart_splitmerge'):
    pic_list = []
    pics=os.listdir(path)
    pics_sorted=sorted(pics, key=lambda x: float(x[:-4]))
    frame = cv2.imread(os.path.join(path,pics_sorted[0]))
    height, width, layers = frame.shape
    video_name =  mp4name
    fps=10
    video_codec = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(video_name, video_codec, fps, (width, height))
    for pic in pics_sorted:
        video.write(cv2.imread(os.path.join(path,pic)))

    cv2.destroyAllWindows()
    video.release()

    #cv2.VideoWriter(video_name, video_codec, fps, (width,height))
    # import pdb;pdb.set_trace()
    #for pic in pics_sorted:
    #    im = imageio.imread(path+"/"+pic)
    #    pic_list.append(im)
    #imageio.mimsave(save_name_gif, pic_list, 'GIF', loop=0)


data = Data() #sets it up
