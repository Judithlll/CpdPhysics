import numpy as np
import matplotlib.pyplot as plt
import cgs
import datetime
import csv
import copy
import parameters as pars     
import imageio.v2 as imageio
import os
import glob
import pandas as pd
import parameters as pars
import disk_properties as dp
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

def epsilon_PA (PlanetsLoc,PlanetsMass,cross_p):
    """
    Get the pebble accretion rate
    """

    mus=PlanetsMass/cross_p.mcp

    Hp= physics.H_d(cross_p.Hg, cross_p.St)
    # Hp=cross_p.Hg*(1+cross_p.St/ cross_p.alpha*(1+2*cross_p.St)/(1+cross_p.St))
    hp=Hp/PlanetsLoc

    delv_o_vK=0.52*(mus*cross_p.St)**(1/3)+cross_p.eta/(1+5.7*(mus/cross_p.eta**3*cross_p.St))
    
    P_eff=1/np.sqrt((0.32*np.sqrt(mus*delv_o_vK/ cross_p.St/ cross_p.eta**2))**(-2)+(0.39*mus/ cross_p.eta/hp)**(-2)) #Liu & Ormel 2018
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
    timeL = [1*cgs.yr, 1e3*cgs.yr] 
    locationL = [7*cgs.rJup, 10*cgs.rJup] 
    massL = [3e23, 3e23] 
    compoL = np.array([[1.0, 0.0], [1.0, 0.0]])

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
        dcompos['iceline_init'] = 2*cgs.RJ

    return dcompos

def do_stuff (sys, init=False):
    #data class is available...
    system=copy.deepcopy(sys)
    # import pdb; pdb.set_trace()
    if init:
        data.data_process(system)
        #initialize your data class
    else:
        data.data_process(system)        
        #data object should be available...
        tminarr = np.array([ddum['tmin'] for ddum in system.mintimeL])

        sfmt = '{:5d} {:5d} {:10.2e} {:3d} {:7.2f}'
        line = sfmt.format(system.ntime, system.particles.num, system.deltaT, tminarr.argmin(), system.time/cgs.yr)
        print(line) #TBD: print more things 


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
        self.cumulative_change={'remove':[],'add':0}
        self.planetsmass = {}
        self.planetsloc = {}
        self.icelinesloc = {}
        
    def update_cumulative_change(self,daction):
        if 'remove' in daction.keys():
            self.cumulative_change['remove']+=list(daction['remove'])
        if 'add' in daction.keys():
            self.cumulative_change['add']+=daction['add']
    
    def data_process(self, system):
        """
        time: the system.time 
        Y2d: the particles properties' list
        """
        
        locL = system.particles.locL
        massL = system.particles.massL
        mtotL = system.particles.mtotL
        fcompL = system.particles.fcomp
        daction = system.daction
        time = system.time
        self.update_cumulative_change(daction)
        
        # store time series 
        self.timeL.append(time/cgs.yr2s)
        
        # store particles data
        # if particles ware removed, then put a np.nan into this location. 
        # for now, I can just use the number of the removed particles to add the np.nan to the head of every list. I think more reasonable is to put in the np.nan according to the removed index, but it seems too complex.
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

        self.radD.setdefault(time/cgs.yr2s,rL/cgs.RJ)
        self.mD.setdefault(time/cgs.yr2s,mL)
        self.mtotD.setdefault(time/cgs.yr2s,mtL)
        self.fcompD.setdefault(time/cgs.yr2s, fL)

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


        #store palnets' data
        if pars.doPlanets:
            self.planetsmass.setdefault(system.time, [planet.mass for planet in system.planetL])
            self.planetsloc.setdefault(system.time, [planet.loc for planet in system.planetL])

        #store icelines' data
        if pars.doIcelines:
            self.icelinesloc.setdefault(system.time, [iceline.loc for iceline in system.icelineL])


    def get_plot_list(self):
        self.radL=np.array(list(self.radD.values()))
        self.mL=  np.array(list(self.mD.values()))
        self.mtotL=np.array(list(self.mtotD.values()))
        self.fcompL = np.array(list(self.fcompD.values()))
        
        self.planetsmassL = np.array(list(self.planetsmass.values()))
        self.planetslocL = np.array(list(self.planetsloc.values()))

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

    
    def plot_stuff(self,gas):
        
        [time,loc]=np.meshgrid(self.timeL,np.linspace(pars.dgasgrid['rinn']/cgs.RJ,pars.dgasgrid['rout']/cgs.RJ,len(self.timeL)))
        sigmag=gas.get_key_disk_properties(loc,time)[0]

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


        for i in range(len(self.radL)):
            plt.subplot(211)
            plt.plot(self.timeL,self.radL[i])
            plt.subplot(212)
            plt.plot(self.timeL,self.mL[i])
        
        plt.subplot(211)
        plt.contourf(time,loc,sigmag,alpha=0.3)
        plt.colorbar()

        plt.savefig('test.jpg')

    def plot_disk(self,time,gas):
        r_Span=np.linspace(pars.dgasgrid['rinn'],pars.dgasgrid['rout'])
        plt.figure(figsize=(24,9))
        if type(time) == float or type(time) == np.float64:
            Sigmag=gas.get_key_disk_properties(r_Span,time)[0]
            Td=gas.get_key_disk_properties(r_Span,time)[1]
            plt.subplot(121)
            plt.title('Surface dencity $[g/cm^2]$')
            plt.plot(r_Span/cgs.RJ,Sigmag,label=str(time/cgs.yr))
            plt.subplot(122)
            plt.title('Midplane Temperature $[K]$')
            plt.plot(r_Span/cgs.RJ,Td,label=str(time/cgs.yr))
            plt.legend()
            plt.savefig('diskproperties.jpg')
        else:
            for i in range(len(time)):
                Sigmag, Td = gas.get_key_disk_properties(r_Span, time[i])[0:2]
                plt.subplot(121)
                plt.title('Surface dencity $[g/cm^2]$')
                plt.plot(r_Span/cgs.RJ,Sigmag,label=str(time[i]/cgs.yr))
                plt.subplot(122)
                plt.title('Midplane Temperature $[K]$')
                plt.plot(r_Span/cgs.RJ,Td,label=str(time[i]/cgs.yr))
            plt.legend()
            plt.savefig('diskproperties.jpg')


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
                plt.savefig('planets&pebbles/'+str(self.timeL[i])+'.jpg') 
                plt.close()
    
    def plot_deltaT(self):
        plt.figure()
        plt.title('Time step variation')
        plt.xlabel('System Time [yr]')
        plt.ylabel('Time Step [yr]')
        plt.plot(self.timeL[:-1],np.diff(self.timeL))
        plt.savefig('Delta_t.jpg')
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
