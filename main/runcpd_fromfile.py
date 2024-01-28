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
import fileio


#read data from files
pickle_path = './pickles/'
pickles = f.find_pickles(pickle_path)
#need to judge whether the pickles is complete

completeness = f.check_pickles(pickles)
if not completeness:
    print('[runcpd_fromfile]: the pickles file is not complete' )
    sys.exit()
classes = {}
for pcs in pickles:
    classes[pcs.rsplit('.',1)[0]] = fileio.load_class(pickle_path, pcs)

#put the classes into sim_init
#maybe: 
#   sim_init(.., from_file=True)
#need to be changed: sim_init, init_particles, 
start=time.time()
argL = sys.argv #maybe use later
calldir = init.init_default_pars (argL[0]) #directory from which this is called (maybe we need later)

#this initializes the system...
system = init.sim_init (calldir, pars.dsystempars, fromfile = True , classesdict=classes)
system.init_particles(fromfile= True, particles_fromfile=classes['particles'])
system.milestones[pars.tmax] = 'Ending_time'

userfun.do_stuff(system, init=True)
#import pdb; pdb.set_trace()
#get the initial deltaT
system.new_timestep (pars.tmax, **pars.dtimesteppars)  #get deltaT through comparing some time scales
#backup the initial data
system.back_up_last_data()       #back up the data of last step

while system.time<pars.tmax:

    #[24.01.01]:determines the timestep for this step
    #[24.01.02],LZX: don't understand why this should be here
        # if system.time > system.planetL[0].time:
    #     import pdb; pdb.set_trace()
                        
    #integrate the super particles
    Yt = system.update_particles (**pars.dtimesteppars)
    
    #change planet and super particle properties
    #due to crossings and intrinsic evolution
    if system.nplanet>0:
        core.advance_planets (system)

    if system.niceline>0:
        core.advance_iceline(system)


    doJump, djump = system.query_system_jump(jumpfrac=0.2)
    
    if doJump:
        #1)do jump,
        #2)update the system time with jumpT
        #3)then the deltaT should be the jumpT
        #4)finally, generate the deltaT for next step
        system.system_jump(djump) #changes system.time
        system.time += djump['jumpT']
        #system.deltaT = system.jumpT
        
        #import pdb;pdb.set_trace()
        system.new_timestep(pars.tmax, afterjump = True, **pars.dtimesteppars)
        system.reset_after_jump()
    else:
        #1)update the system time 
        #2)get the new deltaT 
        system.time +=system.deltaT
        system.new_timestep(pars.tmax, **pars.dtimesteppars)
    #do this again? LZX: maybe not, post_process mainly for particles
    #   system.post_process()
    system.post_process()
    
    system.back_up_last_data()       #back up the data of last step
    system.ntime += 1
    
    userfun.do_stuff(system, doJump = doJump)


    end = time.time()
    runTime = end-start

#store  finally state of system, not sure if this is userful
fileio.store_class(system.particles, 'particles')
fileio.store_class(system.gas, 'gas')
for i in range(system.nplanet):
    fileio.store_class(system.planetL[i], 'planet'+str(i+1))

for i in range(system.niceline):
    fileio.store_class(system.icelineL[i], 'iceline'+str(i+1))

#put some necessary properties into store_systemclass and store it
nece_pd = {'time':system.time, 'jumpT': system.jumpT, 'ntime':system.ntime, 'rhoPlanet':system.rhoPlanet, 'nplanet': system.nplanet, 'niceline':system.niceline, 'milestones':system.milestones}

system_store = core.store_system(nece_pd)
fileio.store_class(system_store, 'system_store')

