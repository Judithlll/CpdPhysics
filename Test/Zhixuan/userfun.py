import numpy as np
import matplotlib.pyplot as plt
import cgs
import datetime
import csv
import copy
import parameters as pars     
import imageio.v2 as imageio
import os
import pandas as pd
import parameters as pars
import disk_properties as dp
import physics
from scipy.optimize import curve_fit
from matplotlib.colors import LinearSegmentedColormap

def PIM():
    return 1.48e26
def Stokes_number(v_r, R_d, v_th, lmfp, Omega_K, rho_g, rhoint = 1.4):
    v_dg=np.abs(v_r)
    Rep=4*R_d*v_dg/v_th/lmfp
    #following Shibaike, eq....
    CD=24/Rep*(1+0.27*Rep)**0.43+0.47*(1-np.exp(-0.04*Rep**0.38))
    St=8/3/CD*rhoint*R_d/rho_g/v_dg*Omega_K
    return St

def St_fragment(alpha, cs, eta, loc, Omega_K, v_crit, icelineloc):
    #following Shibaike 2019 eq14
    sil_part = np.argwhere(np.ndarray.flatten(loc<icelineloc))
    icy_part = np.argwhere(np.ndarray.flatten(loc>=icelineloc))

    sq_part = np.append(np.sqrt(9*alpha**2*cs[sil_part]**4 +4*eta[sil_part]**2*loc[sil_part]**2*Omega_K[sil_part]**2*v_crit['silicates']**2),
                        np.sqrt(9*alpha**2*cs[icy_part]**4 +4*eta[icy_part]**2*loc[icy_part]**2*Omega_K[icy_part]**2*v_crit['icy']**2))
    
    St = (-3*alpha*cs**2+sq_part)/(2*eta**2*loc**2*Omega_K**2)
    return St

# seems very complex...
def St2Rd(St,v_r,v_th, lmfp, Omega_K, rho_g, rhoint=1.4):
    v_dg = np.abs(v_r)
    Rep=4*R_d*v_dg/v_th/lmfp
    CD=24/Rep*(1+0.27*Rep)**0.43+0.47*(1-np.exp(-0.04*Rep**0.38))
    return St*rho_g*v_dg*Omega_K*CD*3/8/rhoint

def planet_migration (gas,planetLoc,planetMass,time,rhopl):
    #aerodynamic drag migration
    out = gas.get_key_disk_properties (planetLoc, time)

    mcp = dp.Mcp_t(time)

    disk = physics.DISK (*out, planetLoc, time, mcp) #pro
    disk.add_auxiliary ()
    disk.add_uservar (dp.user_add_var())    #variables
    disk.add_userfun (dp.user_add_fun())
    disk.add_user_eval (dp.user_add_eval())
    #disk = core.DISK (*out, planetLoc, time)
    #disk.add_auxiliary ()
    #disk.user_difined ()

    # eta=disk.eta
    # v_K=disk.vK
    # v_th=disk.vth
    # lmfp=disk.lmfp
    # rho_g=disk.rhog
    # Omega_K=disk.OmegaK
    # dotMg=disk.dotMg
    # Mcp=disk.Mcp
    # mg=disk.mu*cgs.mp
    # Sigmag=disk.sigmaG
    # cs=disk.cs
    
    # Rpl=(planetMass/(4/3*np.pi*rhopl))**(1/3)

    # Stpl,vd=f.St_iterate(eta,v_K,v_th,lmfp,rho_g,Omega_K,Rpl)

    try:
        rinn = disk.rinn
    except:
        rinn = 0.0
    
    if planetLoc < rinn:
        vt1 = 0.0
    else:
    #Type I migration
        
        #qr=-0.14493*disk.dot_Mg(disk.time)*cgs.gC*disk.mcp*(-0.206349*planetLoc**(5.5)+planetLoc**4.5*disk.rout)/disk.alpha/planetLoc**(8.5)/disk.rout/np.sqrt(cgs.kB*(disk.dot_Mg(time)*cgs.gC*disk.mcp/planetLoc**3/cgs.sigmaSB)**(1/4)/disk.mg) #pressure gradient
        
        #kley-nelson 2011 eq14
        ## TBD: make more check: this migration rate is much larger than old one, and seems will be positive sometimes
        torq0 = (planetMass/mcp)**2 *(disk.Hg/planetLoc)**(-2)*disk.sigmaG *planetLoc**4*disk.OmegaK**2
        p_phi = planetMass*np.sqrt(cgs.gC*mcp*planetLoc)
        
        torq_tot = -(1.36-0.62*disk.beta_sigG(planetLoc)- 0.43*disk.beta_T) *torq0
        vt1 = 2*planetLoc*torq_tot/p_phi
        #CI=0.1
        #bt1=CI*(2.7+1.1*qr)   #a constant Ogihara 2014
        #vt1o=bt1*(planetMass/disk.mcp)*(disk.sigmaG*planetLoc**2/disk.mcp)*(disk.vK/disk.cs)**2*disk.vK
        if vt1 >0:
            import pdb;pdb.set_trace()
    # v_mig=vt1+vd

    #if planetMass > 5e26: import pdb; pdb.set_trace()

    return vt1 # v_mig

def del_v (St, disk):
    vr = 2*disk.eta *disk.vK *St/(1+St**2)

    #take half of the velocity...
    return np.abs(vr)/2


def H_d (St, disk):
    return physics.H_d(disk.Hg, St, disk.alpha) 
    

def dm_dt(Rd, delv, Hd, sigD, fcomp):
    """
    the time derivetive of particles's mass, determine how particles grow
    if icy fraction == 0, then come to fragmentation (some assumption here)
    """
    #TBD calculate/define fragm. velocity based on composition
    vc = pars.vc['silicates']*fcomp[:,0]+ pars.vc['icy']*fcomp[:,1]
    Fcomp = np.ones_like(Rd)
    Fcomp = np.where(delv/vc>1, -1, 1) 

    return Fcomp *2*np.sqrt(np.pi)*Rd**2*delv/Hd*sigD   #eq. 5 of Shibaike et al. 2017

def epsilon_PA (planetLoc,planetMass,cross_p):
    """
    Get the pebble accretion rate
    """

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
    timeL = [0.8e6*cgs.yr, 1.3e6*cgs.yr, 1.8e6*cgs.yr, 2.3e6*cgs.yr] 
    #some things wrong with the initial location is set to the out edge
    #about particles number
    locationL = [50*cgs.rJup, 50*cgs.rJup, 50*cgs.rJup, 50*cgs.rJup] 
    massL = [3e23, 3e23, 3e23,3e23] 
    compoL = np.array([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])

    return timeL, locationL, massL, compoL

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
        dcompos['iceline_init'] = 15*cgs.RJ

    return dcompos

def do_stuff (system, init=False, final= False):
    #data class is available...
    # import pdb; pdb.set_trace()
    if init:
        data.data_process(system)
        data.gas = system.gas
        #initialize your data class
    else:
        data.data_process(system)
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
        self.cumulative_change={'remove':[],'add':0}
        self.planetsmass = {}
        self.planetsloc = {}
        self.icelinesloc = {}
        self.jumpstuff = []
        self.planetsfcomp = {}

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
            for k,v in system.planetD.items():
                pmassL[k] = v.mass
                plocL [k] = v.loc
                pcompL[k] = v.fcomp

            self.planetsmass.setdefault(system.time,pmassL)
            self.planetsloc.setdefault(system.time,plocL )
            self.planetsfcomp.setdefault(system.time, pcompL)
        #store icelines' data
        if pars.doIcelines:
            self.icelinesloc.setdefault(system.time, [iceline.loc for iceline in system.icelineL])
        
        #if jump is done, then shore something about jump
        if system.doJump:
            stuff = {'njump': system.njump, 'njumptime': system.njumptime, 'jumptime':system.time-system.jumpT, 'jumpT': system.jumpT, 'jump_limitation':system.jump_limitation}
            self.jumpstuff.append(stuff)


    def get_plot_list(self, doParticles = False):
        """
        process data and make it in order

        doParticles: decide to process particles' data or not
        """
        #TBD: here process the particles data with dactionD
        
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
        pfL = copy.deepcopy(list(self.planetsfcomp.values()))
        max_len = max([len(l) for l in pmL])
        
        # maybe can use the number of planet
        
        for i in range (len(pmL)):
            if len(pmL[i]) != max_len:
                length = max_len - len(pmL[i])
                pmL[i].extend([np.nan]*length)
                plL[i].extend([np.nan]*length)
                pfL[i].extend([[np.nan]*2]*length)
        self.planetsmassL = np.array(pmL)
        self.planetslocL = np.array(plL)
        self.planetsfcompL = np.array(pfL)

        self.icelineslocL = np.array(list(self.icelinesloc.values()))


    def data_store (self,path=os.getcwd()):
        
        self.get_plot_list()

        #store particles data
        df_rad = pd.DataFrame(self.radL)
        df_mass = pd.DataFrame(self.mL)
        df_mtot = pd.DataFrame(self.mtotL)
        df_fcomp = pd.DataFrame(self.fcompL)
        
        df_rad.index = self.timeL
        df_mass.index = self.timeL
        df_mtot.index = self.timeL
        df_fcomp.index = self.timeL

        ## write CSV
        import pdb; pdb.set_trace()
        writer = pd.ExcelWriter('particles_data.xlsx', engine='xlsxwriter')
        df_rad.to_excel (writer, sheet_name= 'location data')
        df_mass.to_excel (writer, sheet_name= 'mass data')
        df_mtot.to_excel (writer, sheet_name= 'total mass data')
        df_fcomp.to_excel (writer, sheet_name= 'composition data')
        writer.close()

        #store planets data
        df_plmass = pd.DataFrame(self.planetsmassL)
        df_plloc = pd.DataFrame(self.planetslocL)
        df_plmass.index = self.timeL
        df_plloc.index = self.timeL

        writer = pd.ExcelWriter('planets_data.xlsx', engine='xlsxwriter')
        df_plmass.to_excel (writer, sheet_name= 'mass data')
        df_plloc.to_excel (writer, sheet_name= 'location data')     
        
        writer.close()

        #store iceline data
        df_illoc = pd.DataFrame(self.icelineslocL)
        df_illoc.index = self.timeL

        writer = pd.ExcelWriter('icelines_data.xlsx', engine='xlsxwriter')
        df_illoc.to_excel (writer, sheet_name= 'location data')     
        
        writer.close()



        #with open(path+str(datetime.datetime.now())+'data_particles.csv', 'w', newline='') as csvfile:
        #    writer=csv.DictWriter(csvfile,fieldnames=self.timeL)
        #    writer.writeheader()
        #    writer.writerows([self.radD,self.mD,self.mtotD])

        # with open(str(datetime.datetime.now())+'data_mass.csv', 'w', newline='') as csvfile:
        #     writer=csv.DictWriter(csvfile,fieldnames=self.timeL)
        #     writer.writeheader()
        #     writer.writerows(self.mD)       

    
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
            plt.title('Surface dencity $[g/cm^2]$', )
            plt.xlabel('Location [$R_J$]')
            plt.xscale('log')
            plt.yscale('log')
            plt.plot(r_Span/cgs.RJ,Sigmag,label=str(time/cgs.yr))
            plt.subplot(122)
            plt.title('Midplane Temperature $[K]$')
            plt.xlabel('Location [$R_J$]')
            plt.xscale('log')
            plt.yscale('log')
            plt.plot(r_Span/cgs.RJ,Td,label=str(time/cgs.yr))
            plt.legend()
            plt.savefig('./plot/diskproperties.jpg')
        else:
            plt.subplot(121)
            plt.title('Surface dencity $[g/cm^2]$')
            plt.xlabel('Location [$R_J$]')
            plt.xscale('log')
            plt.yscale('log')
            
            plt.subplot(122)
            plt.title('Midplane Temperature $[K]$')
            plt.xlabel('Location [$R_J$]')
            plt.xscale('log')
            plt.yscale('log')

            for i in range(len(time)):
                Sigmag, Td = self.gas.get_key_disk_properties(r_Span, time[i])[0:2]
                plt.subplot(121)
                plt.plot(r_Span/cgs.RJ,Sigmag,label=str(time[i]/cgs.yr))
                plt.subplot(122)
                plt.plot(r_Span/cgs.RJ,Td,label=str(time[i]/cgs.yr))

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
        plt.plot(self.timeL, self.num, color = 'red')
        plt.savefig('./plot/particleNum.jpg')
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
        plt.title('Planet evolution')
        plt.xlabel('Planets location [$R_{Jup}$]' )
        plt.ylabel('System time [yr]')
        plt.yscale('log')
        plt.ylim(1e4,30e6)
        plt.xscale('log')
        loclist = self.planetslocL.T
        masslist = self.planetsmassL.T
        time = np.array(list(self.planetsloc.keys()))
        
        planetst = 0.8e6*cgs.yr
        stidx = np.argwhere(time>planetst)[0][0]
        
        dotssize = masslist/np.nanmin(masslist)*0.1
        #cmap = LinearSegmentedColormap.from_list("my_colormap", ["y", "royalblue"])

        for jump in self.jumpstuff:
            if jump['jumptime'] > planetst:
                plt.axhspan((jump['jumptime']-planetst)/cgs.yr, 
                            (jump['jumptime']+jump['jumpT']-planetst)/cgs.yr, alpha = 0.3)
        for i,loc in enumerate(loclist):
            plt.plot(loc[stidx:]/cgs.RJ, (time[stidx:]-planetst)/cgs.yr)
            plt.scatter(loc[stidx:]/cgs.RJ, (time[stidx:]-planetst)/cgs.yr, s = dotssize[i][stidx:], c =self.planetsfcompL[stidx:,i][:,1], cmap ='Spectral', alpha =1 )


            #plt.axhline((jump['jumptime']+jump['jumpT'])/cgs.yr, color = 'green', linewidth = 0.2)
        plt.colorbar(label = "Water Fraction [%]")
        plt.axvline(dp.rinn/cgs.RJ, color = 'gray', linewidth = 0.5, label = 'inner edge')
        plt.plot(self.icelineslocL[stidx:,0]/cgs.RJ, (time[stidx:]-planetst)/cgs.yr, color = 'blue', linestyle = 'dashed',label = 'Iceline')
        plt.legend()
        plt.xticks([5.89,10,14.8,15,20,50],['5.89','10','14.8','','20','50'])
        plt.savefig('./plot/planet_evolution.jpg',dpi=600)
        plt.close()
   

    def plot_iceline(self):
        plt.figure()
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel('iceline location')
        plt.xlabel('time')
        plt.yticks([14.8,15],['14.8','15'])
        plt.axhline(15.)
        plt.plot(np.array(self.timeL)/cgs.yr, self.icelineslocL.T[0]/cgs.RJ, 'x-')
        plt.savefig('./plot/iceline.jpg')
        plt.close()
        import pdb;pdb.set_trace()
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
    
    def plot_pebble_Sigma(self, tL, mode = 'particle'):
        timeL =tL + dp.tgap
        
        plt.figure()
        plt.title('pebble surface density')
        
        if mode == 'grid':

            grids = np.logspace(np.log10(dp.rinn), np.log10(dp.rout), 10000)
            width = np.diff(grids)

            plt.yscale('log')
            plt.xscale('log')
            plt.ylim(1e-5,10)
            for time in timeL:
                tidx = np.argmin(np.abs(self.timeL-time))
                ti = self.timeL[tidx]
                mtot = self.mtotD[ti]
                loc = self.radD[ti] #ordered

                boundaries = np.sqrt(loc[1:]*loc[:-1])
                boundaries = np.append(dp.rinn,boundaries)
                boundaries = np.append(boundaries,dp.rout)
                warr = np.diff(boundaries)
                sigma = mtot /(2*np.pi*loc*warr)
                    
                plt.plot(loc/cgs.rJup, sigma, 'x-', label = str('{:7.2f}'.format(time/cgs.yr)))
                #dotMd = dp.dot_Mg(ti)*dp.ratio
                #loc = self.radD[ti]
                #v_r = self.v_rD[ti]


                #sigmaP = dotMd/(2*np.pi*loc*(-v_r))
                #plt.plot(loc/cgs.rJup, sigmaP, label = 'particles'+"{:.2f}".format(ti/cgs.yr) )
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
    
    def plot_St(self, tL):
        timeL =tL
        plt.figure()
        plt.title('Stokes number')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e-6, 1)
        plt.xlabel('Location [$R_J$]')
        plt.ylabel('Stokes number')
        for time in timeL:
            tidx = np.argmin(np.abs(self.timeL-time))
            ti = self.timeL[tidx]
            
            St = self.StD[ti]
            import pdb;pdb.set_trace()
            loc = self.radD[ti]
            
            plt.plot(loc/cgs.rJup, St, 'x-',label = "{:.2f}".format(ti/cgs.yr))
            for i,loc in enumerate(self.planetsloc[ti]): 
                plt.axvline(loc/cgs.rJup, linestyle='dashed', color='gray', label='planet'+str(i))
            
        plt.axvline(self.icelineslocL[-1]/cgs.RJ, linestyle = 'dashed', color = 'gray', label = 'water iceline')
        plt.legend()
        plt.savefig('./plot/St.jpg')
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
            plt.axhline(pars.vc['icy'], label = 'Fragmentation velocity [icy]')
            plt.axhline(pars.vc['silicates'], label = 'Fragmentation velocity [sil]')
            
            plt.legend()
            plt.savefig('./plot/v_r.jpg')


    def plot_jumpT(self):
        plt.figure()
        plt.title('Jump time')
        plt.xlabel('System time [yr]')
        plt.ylabel('Jump time [yr]')
        jumpTlist = [f['jumpT']/cgs.yr for f in self.jumpstuff]
        timelist = [f['jumptime']/cgs.yr for f in self.jumpstuff]
        plt.plot(timelist, jumpTlist)
        plt.scatter(timelist, jumpTlist)
        plt.savefig('./plot/jumpT variation')
        plt.close()

def make_animation(path='pebbles&planets'):
    save_name_gif =  "Cpd.gif"
    pic_list = []
    pics=os.listdir(path)
    pics_sorted=sorted(pics, key=lambda x: float(x[:-4]))
    # import pdb;pdb.set_trace()
    for pic in pics_sorted:
        im = imageio.imread(path+"/"+pic)
        pic_list.append(im)
    imageio.mimsave(save_name_gif, pic_list, 'GIF', loop=0)

def load_data(path=os.getcwd()):
    #TBD: load data from excel
    
    #load particles data
    filename = '/particles_data.xlsx'
    df_rad = pd.read_excel (path + filename, sheet_name = 'location data', engine = 'openpyxl')
    df_mass = pd.read_excel (path + filename, sheet_name = 'mass data', engine = 'openpyxl')
    df_mtot = pd.read_excel (path + filename, sheet_name = 'total mass data', engine = 'openpyxl')
    df_fcomp = pd.read_excel (path + filename, sheet_name = 'composition data', engine = 'openpyxl')

    #load planets data
    filename = '/planets_data.xlsx'
    df_plloc = pd.read_excel (path + filename, sheet_name = 'location data', engine = 'openpyxl')
    df_plmass = pd.read_excel (path + filename, sheet_name = 'mass data', engine = 'openpyxl')

    #load icelines data
    filename = '/icelines_data.xlsx'
    df_illoc = pd.read_excel (path + filename, sheet_name = 'location data', engine = 'openpyxl')
    
    loaddata = Data()
    loaddata.radL = df_rad.iloc[:,1:]
    loaddata.mL=df_mass.iloc[:,1:]
    loaddata.mtotL = df_mtot.iloc[:,1:]
    loaddata.fcomL = df_fcomp.iloc[:,1:]
    loaddata.planetsmass = df_plmass.iloc[:,1:]
    loaddata.planetsloc = df_plloc.iloc[:,1:]
    loaddata.icelinesloc = df_illoc.iloc[:,1:]
    loaddata.timeL = df_illoc.iloc[:,0]

    return loaddata


data = Data() #sets it up
