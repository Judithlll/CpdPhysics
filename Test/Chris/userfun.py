import numpy as np
import matplotlib.pyplot as plt
import cgs
import datetime
import csv

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
    return [10*cgs.yr, 1e3*cgs.yr], [7*cgs.rJup, 10*cgs.rJup], [1e23, 1e23], [1.0, 1.0]


class data_process(object):
    """
    To store data and process data (convert\plot)
    """

    def __init__(self):

        self.timeL=[]
        self.radL=np.array([])
        self.mL=np.array([])
        self.radD={}
        self.mD={}
        self.cumulative_change={'remove':[],'add':0}
        
    def update_cumulative_change(self,daction):
        if 'remove' in daction.keys():
            self.cumulative_change['remove']+=list(daction['remove'])
        if 'add' in daction.keys():
            self.cumulative_change['add']+=daction['add']
    
    def data_process(self,Y2d,time,daction):
        """
        time: the system.time 
        Y2d: the particles properties' list
        """
        
        self.update_cumulative_change(daction)

        self.timeL.append(time/cgs.yr2s)

        # if particles ware removed, then put a np.nan into this location. 
        # for now, I can just use the number of the removed particles to add the np.nan to the head of every list. I think more reasonable is to put in the np.nan according to the removed index, but it seems too complex.
        if len(self.cumulative_change['remove'])>0:
            rL=np.insert(Y2d[0],0,np.full(len(self.cumulative_change['remove']),np.nan))
            mL=np.insert(Y2d[1],0,np.full(len(self.cumulative_change['remove']),np.nan))

        else:

            rL=Y2d[0]
            mL=Y2d[1]

        self.radD.setdefault(time/cgs.yr2s,rL/cgs.RJ)
        self.mD.setdefault(time/cgs.yr2s,mL)
        max_len=max(len(v) for v in self.radD.values())

        # want to make the data dict in the same length
        for k,v in self.radD.items():
            
            #if particles were added, then add np.nan to the former list to make the length the same.
            if len(v)< max_len:
                self.radD[k]=np.pad(v, (0, max_len - len(v)), constant_values=np.nan)
                self.mD[k]=np.pad(self.mD[k], (0, max_len - len(v)), constant_values=np.nan)
        self.get_particles_plot_list()
        

    def get_particles_plot_list(self):
        self.radL=np.array(list(self.radD.values())).T
        self.mL=  np.array(list(self.mD.values())).T

    def data_store(self):
        with open(str(datetime.datetime.now())+'data_location.csv', 'w', newline='') as csvfile:
            writer=csv.DictWriter(csvfile,fieldnames=self.timeL)
            writer.writeheader()
            writer.writerows([self.radD,self.mD])

        # with open(str(datetime.datetime.now())+'data_mass.csv', 'w', newline='') as csvfile:
        #     writer=csv.DictWriter(csvfile,fieldnames=self.timeL)
        #     writer.writeheader()
        #     writer.writerows(self.mD)       

    
    def plot_stuff(self,disk):
        
        [time,loc]=np.meshgrid(self.timeL,np.linspace(disk.rinn/cgs.RJ,disk.rout/cgs.RJ,len(self.timeL)))
        sigmag=disk.Sigmag(loc,time)


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

    def plot_disk(self,time,disk):
        r_Span=np.linspace(disk.rinn,disk.rout)
        Sigmag=disk.Sigma_g(r_Span,time)
        Td=disk.T_d(r_Span,time)

        

