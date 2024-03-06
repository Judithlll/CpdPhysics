#!/usr/bin/env python
import sys, os
getwcd = os.getcwd()    #these lines are needed 
sys.path.append(getwcd) #on my desktop (I don't know why)
import parameters as pars
sys.path.append('../../main/')
import core
import init
import userfun
import time 
import fileio
import time
import matplotlib.pyplot as plt
import numpy as np
import cgs

start=time.time()
argL = sys.argv #maybe use later
if 'fromfile' in argL:
    system = fileio.load_class('./pickles/', 'system.pickle')
    if pars.tmax <= system.time:
        print('\033[31m WARNING\033[0m : the ending time is smaller than the current system time [runcpd]')
        sys.exit()
    print('\033[32m [runcpd]: run from pickle file at {system.time/cgs.yr} \033[0m')
    time.sleep(2)  
else:
    calldir = init.init_default_pars (argL[0]) #directory from which this is called (maybe we need later)

#this initializes the system...
    system, gasL = init.sim_init (calldir, pars.dsystempars)
    system.init_particles(pars.dparticleprops)
    system.milestones[pars.tmax] = 'Ending_time'
#initialize userfun's data class
    #userfun.do_stuff(system, init=True)
#import pdb; pdb.set_trace()
#get the initial deltaT
    system.new_timestep (pars.tmax, **pars.dtimesteppars)  #get deltaT through comparing some time scales
#backup the initial data
    system.back_up_last_data()       #back up the data of last step
    print('\033[32m [runcpd]: run from the beginning \033[0m')
    time.sleep(2)  

while True:

    #[24.01.01]:determines the timestep for this step
    #[24.01.02],LZX: don't understand why this should be here
        # if system.time > system.planetL[0].time:
    #     import pdb; pdb.set_trace()
                        
    #integrate the super particles
    if len(system.planet_candidate)>0:
        for planet in system.planet_candidate:
            if planet.starttime <= system.time:
                system.add_planet(planet)
                system.planet_candidate.remove(planet)
    
    Yt = system.update_particles (**pars.dtimesteppars)
    
    #change planet and super particle properties
    #due to crossings and intrinsic evolution
    if system.nplanet>0:
        core.advance_planets (system)

    if system.niceline>0:
        core.advance_iceline(system)


    ## cwo: resonance crossing criteria
    if False:
        #check if any planets crossed resonances
        for iplanet in range(1,system.nplanet):
            prat = (system.yvec[iplanet,0]/system.yvec[iplanet-1,0])**1.5
            inxt = (dres['prat']<prat).argmax()

            pinfoi = system.planetinfo[iplanet]

            #it locks into a resonance
            #if inxt>pinfoi['inxt'] and pinfoi['resS']!='R':



    djump = system.query_system_jump(jumpfrac=0.2)
    #system.doJump = False
    if system.doJump:
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
    
    if system.deltaT <= 0:
        print('[runcpd]: something wrong with deltaT')
        import pdb; pdb.set_trace()

    system.back_up_last_data()       #back up the data of last step
    system.ntime += 1
    
    #[24.01.30]cwo: dont append other arguments to "do_stuff"
    #so, find a way to append "djump" to system

    #[24.02.01]cwo: added stopping condition // we could add more
    final = system.time>=pars.tmax

    userfun.do_stuff(system, final=final)
    
    print ([p.dlocdt for p in system.planetL], [p.loc/cgs.RJ for p in system.planetL], system.time/cgs.yr)
    if final: 
        end = time.time()
        runTime = end-start
        break
#store system components as pickles
fileio.store_class(system, 'system')
fileio.store_class(userfun.data, 'data')
userfun.data.plot_planet_migration()
userfun.data.plot_jumpT()
userfun.data.plot_stuff(system.gas)
import pdb;pdb.set_trace()
plt.figure()
plt.plot(userfun.data.timeL, userfun.data.planetslocL.T[1]/userfun.data.planetslocL.T[0])
plt.plot(userfun.data.timeL, 2*np.ones_like(userfun.data.timeL))
plt.savefig("pratio.jpg")
print('[runcpd]:finished')


