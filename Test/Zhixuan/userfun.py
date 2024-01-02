import numpy as np
import matplotlib.pyplot as plt
import cgs
import datetime
import csv
import copy
import parameters as pars
import imageio.v2 as imageio
import os

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
        data.data_process(system.particles.locL,system.particles.massL,system.particles.mtotL,system.time,system.daction,system.planetL)
        #initialize your data class
    else:
        data.data_process(system.particles.locL,system.particles.massL,system.particles.mtotL,system.time,system.daction,system.planetL)
        #data object should be available...



class Data(object):
    """
    To store data and process data (convert\plot)
    """

    def __init__(self):

        self.timeL=[]
        self.radL=np.array([])
        self.mL=np.array([])
        self.mtotL=np.array([])
        self.radD={}
        self.mD={}
        self.mtotD={}
        self.cumulative_change={'remove':[],'add':0}
        self.planetL=[]
        
    def update_cumulative_change(self,daction):
        if 'remove' in daction.keys():
            self.cumulative_change['remove']+=list(daction['remove'])
        if 'add' in daction.keys():
            self.cumulative_change['add']+=daction['add']
    
    def data_process(self,locL,massL,mtotL,time,daction,planetL):
        """
        time: the system.time 
        Y2d: the particles properties' list
        """
        
        self.update_cumulative_change(daction)

        self.timeL.append(time/cgs.yr2s)

        # if particles ware removed, then put a np.nan into this location. 
        # for now, I can just use the number of the removed particles to add the np.nan to the head of every list. I think more reasonable is to put in the np.nan according to the removed index, but it seems too complex.
        if len(self.cumulative_change['remove'])>0:
            rL=np.insert(locL,0,np.full(len(self.cumulative_change['remove']),np.nan))
            mL=np.insert(massL,0,np.full(len(self.cumulative_change['remove']),np.nan))
            mtL=np.insert(mtotL,0,np.full(len(self.cumulative_change['remove']),np.nan))

        else:

            rL=locL
            mL=massL
            mtL=mtotL

        self.radD.setdefault(time/cgs.yr2s,rL/cgs.RJ)
        self.mD.setdefault(time/cgs.yr2s,mL)
        self.mtotD.setdefault(time/cgs.yr2s,mtL)

        max_len=max(len(v) for v in self.radD.values())

        # want to make the data dict in the same length
        for k,v in self.radD.items():
            
            #if particles were added, then add np.nan to the former list to make the length the same.
            if len(v)< max_len:
                self.radD[k]=np.pad(v, (0, max_len - len(v)), constant_values=np.nan)
                self.mD[k]=np.pad(self.mD[k], (0, max_len - len(v)), constant_values=np.nan)
                self.mtotD[k]=np.pad(self.mtotD[k], (0, max_len - len(v)), constant_values=np.nan)
        
        self.get_particles_plot_list()
        self.planetL=planetL


    def get_particles_plot_list(self):
        self.radL=np.array(list(self.radD.values())).T
        self.mL=  np.array(list(self.mD.values())).T
        self.mtotL=np.array(list(self.mtotD.values())).T

    def data_store(self,path):
        with open(path+str(datetime.datetime.now())+'data_particles.csv', 'w', newline='') as csvfile:
            writer=csv.DictWriter(csvfile,fieldnames=self.timeL)
            writer.writeheader()
            writer.writerows([self.radD,self.mD,self.mtotD])

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
        Sigmag=gas.get_key_disk_properties(r_Span,time)[0]
        Td=gas.get_key_disk_properties(r_Span,time)[1]
        plt.figure(figsize=(24,9))
        plt.subplot(121)
        plt.title('Surface dencity $[g/cm^2]$')
        plt.plot(r_Span/cgs.RJ,Sigmag,label=str(time/cgs.yr))
        plt.subplot(122)
        plt.title('Midplane Temperature $[K]$')
        plt.plot(r_Span/cgs.RJ,Td,label=str(time/cgs.yr))
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

data = Data() #sets it up