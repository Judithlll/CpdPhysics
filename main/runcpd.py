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

argL = sys.argv #maybe use later
calldir = init.init_default_pars (argL[0]) #directory from which this is called (maybe we need later)

#this initializes the system...
system = init.sim_init (pars.dsystempars)


#initialize userfun's data class
userfun.do_stuff(system, init=True)


while system.time<pars.tmax:

    #integrate the super particles
    Yt = system.update_particles (pars.tmax)

    #import pdb; pdb.set_trace()

    #change planet and super particle properties
    #due to crossings and intrinsic evolution
    if system.nplanet>0:
        core.advance_planets (system)

    system.post_process()
    #TBD: postprocess particles (add/remove)
    
    #TBD: change system.time
    system.time += system.deltaT

    userfun.do_stuff(system)
    system.iceline_loc()
    idx=np.array([])
    idx=core.advance_iceline(system)
    

    #print('hello', system.time/cgs.yr)

