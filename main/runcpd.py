#!/usr/bin/env python
import sys, os
getwcd = os.getcwd()    #these lines are needed 
sys.path.append(getwcd) #on my desktop (I don't know why)
import parameters as pars
sys.path.append('../../main/')
import core
import numpy as np
import init
import cgs
import userfun
import time 
import functions as f

start=time.time()
argL = sys.argv #maybe use later
calldir = init.init_default_pars (argL[0]) #directory from which this is called (maybe we need later)

#this initializes the system...
system, gasL = init.sim_init (calldir, pars.dsystempars)
system.init_particles(pars.dparticleprops)

#initialize userfun's data class
userfun.do_stuff(system, init=True)

while system.time<pars.tmax:

    #[24.01.01]:determines the timestep for this step
    #[24.01.02],LZX: don't understand why this should be here
    system.new_timestep (pars.tmax, pars.evolvingparas['deltaTfraction'])  #get deltaT through comparing some time scales
    # if system.time > system.planetL[0].time:
    #     import pdb; pdb.set_trace()
    system.back_up_last_data()       #back up the data of last step
                        
    #integrate the super particles
    Yt = system.update_particles (pars.evolvingparas['timestepn'])
    
    #change planet and super particle properties
    #due to crossings and intrinsic evolution
    if system.nplanet>0:
        core.advance_planets (system)

    if system.niceline>0:
        ## CWO: removed "idx"
        core.advance_iceline(system)



    system.post_process()
    
    system.time += system.deltaT
    system.ntime += 1

    userfun.do_stuff(system)

    end = time.time()
    runTime = end-start

print('[runcpd]:finished')

