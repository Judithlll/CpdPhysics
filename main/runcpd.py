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
        #data = fileio.load_class('./pickles/', 'data.pickle')
        doEvo = False
    else:
        print('\033[32m [runcpd]: run from pickle file at {system.time/cgs.yr} \033[0m')
        doEvo = True
        time.sleep(1)  
else:
    doEvo = True
    calldir = init.init_default_pars (argL[0]) #directory from which this is called (maybe we need later)

    #this initializes the system...
    system, gasL = init.sim_init (calldir, pars.dsystempars)
    system.init_particles(pars.dparticleprops)
    system.milestones[pars.tmax] = 'Ending_time'
    #initialize userfun's data class
    userfun.do_stuff(system, init=True)
    #import pdb; pdb.set_trace()
    #get the initial deltaT
    system.new_timestep (pars.tmax, jumpfracD=pars.jumpfracD, **pars.dtimesteppars)  #get deltaT through comparing some time scales
    #backup the initial data
    system.back_up_last_data()       #back up the data of last step

    # write temp data to log files
    system.update_log(init = True)
    print('\033[32m[runcpd]: run from the beginning \033[0m')
    time.sleep(1)  

#[24/11/28]LZX: tbr when these are not needed 
outflux =[]
removenum = 0
addnum = 0 
plotnum= 0

while doEvo:

    #[24.01.01]:determines the timestep for this step
    #[24.01.02],LZX: don't understand why this should be here
        # if system.time > system.planetL[0].time:
    #     import pdb; pdb.set_trace()
                        
    #[24.08.06]LZX: put the add_planet() into userfun
    if pars.doPlanets:
        system = userfun.add_planet(system, pars.planetaddmode)
    
    #integrate the super particles
    Yt = system.update_particles (**pars.dtimesteppars)

    #change planet and super particle properties
    #due to crossings and intrinsic evolution
    #NOTE:order may matter a lot!!
    if system.niceline>0:
        core.advance_iceline(system)

    if system.nplanet>0:
        core.advance_planets (system)



    #do this again? LZX: maybe not, post_process mainly for particles
    #   system.post_process()
    system.post_process()
    #print([system.daction, system.particles.mtot1, system.Minflux])
    

    ## cwo: resonance crossing criteria
    if False:
        #check if any planets crossed resonances
        for iplanet in range(1,system.nplanet):
            prat = (system.yvec[iplanet,0]/system.yvec[iplanet-1,0])**1.5
            inxt = (dres['prat']<prat).argmax()

            pinfoi = system.planetinfo[iplanet]

            #it locks into a resonance
            #if inxt>pinfoi['inxt'] and pinfoi['resS']!='R':

    djump = system.query_system_jump()

    if system.doJump:
        #1)do jump,
        #2)update the system time with jumpT
        #3)then the deltaT should be the jumpT
        #4)finally, generate the deltaT for next step
        system.system_jump(djump) #changes system.time
        system.time += djump['jumpT']
        #system.deltaT = system.jumpT
        system.get_auxiliary(system.time)
        
        #import pdb;pdb.set_trace()
        system.new_timestep(pars.tmax, afterjump = True, jumpfracD = pars.jumpfracD, **pars.dtimesteppars)
        system.reset_after_jump()
    else:
        ## All steps have been taken, update ALL disk and particle properties to the new time
        #1)update the system time 
        #2)get the new deltaT 
        system.time += system.deltaT

        system.get_auxiliary(system.time)
        system.new_timestep(pars.tmax, jumpfracD=pars.jumpfracD, **pars.dtimesteppars)
        system.re_sample()
        #system.query_splitmerge ()
        #system.split_merge (..) 
    
    system.timeL.append(system.time)


    system.back_up_last_data()       #back up the data of last step
    system.ntime += 1

    #LZX[24.09.06]:tbr when these are not needed
    outflux.append(system.Moutflux)
    if 'remove' in system.daction.keys():
        removenum += len(system.daction['remove'])
    if 'add' in system.daction.keys():
        addnum += system.daction['add']



    #[24.02.01]cwo: added stopping condition // we could add more
    final = system.time>=pars.tmax

    # wrrite key information into logs and print things
    #[24.04.21]cwo:no absolute paths, please!
    system.update_log(djump = djump)
    userfun.do_stuff(system, final=final)

    #tbr
    # plot the surface density profile
    if system.time/cgs.yr > plotnum: #plot every 1 yr
        userfun.plot_sfd(system.particles.locL, system.particles.sfd, system.time, system.minTimes.dpart['imin'], system.deltaT, system.timeL, system.resam_time)
        plotnum += 1
        
    # print ([p.dlocdt for p in system.planetL], [p.loc/cgs.RJ for p in system.planetL], system.time/cgs.yr)
    if final: 
        end = time.time()
        runTime = end-start
        print ('[runcpd]: Congrats!! You finished a sucessful run which consume {:.2f} seconds'.format(runTime))
        break

    else:
        if system.deltaT <= 0:
            print('[runcpd]: something wrong with deltaT')
            import pdb;pdb.set_trace()



#tbr[24.09.06]
plt.figure()
plt.xscale('log')
plt.plot(system.timeL, outflux)
plt.savefig('outflux_'+pars.sfdmode+'.png')
plt.close()


print('remove number:', removenum)
print('add number:', addnum)
print('[runcpd]:finished')


