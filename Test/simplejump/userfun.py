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
    #make this a smaller value, to stop the accretion and allow the jump
    return 1e10
def Stokes_number(delv, Rd, vth, lmfp, OmegaK, rhog, rhoint = 1.4):
    v_dg=np.abs(delv)
    Rep=4*Rd*v_dg/vth/lmfp
    #following Shibaike, eq....
    CD=24/Rep*(1+0.27*Rep)**0.43+0.47*(1-np.exp(-0.04*Rep**0.38))
    St=8/3/CD*rhoint*Rd/rhog/v_dg*OmegaK
    return St
def M_critical (*args):

    return 0.0

def planet_migration (*args):
    return -0.02

def del_v (St, disk):

    return 0.0


def H_d (St, disk):
    return physics.H_d(disk.Hg, St, disk.alpha) 
    

def dm_dt(Rd, delv, Hd, sigD, fcomp):
    #don't consider the growth of particles
    return 0.0 


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
    timeL = [0.1e6*cgs.yr, 0.2e6*cgs.yr] 
    #some things wrong with the initial location is set to the out edge
    #about particles number
    locationL = [50*cgs.rJup, 50*cgs.rJup] 
    massL = [3e23, 3e23] 
    compoL = np.array([[1.0, 0.0], [1.0, 0.0]])

    return timeL, locationL, massL, compoL


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
        fileio.store_class(data, 'data')
    else:
        data.data_process(system)
        sp.run('tail -n1 log/system_evol.log', shell=True)



   

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




data = Data() #sets it up
