import sys
sys.path.append('/home/lzx/CpdPhysics/main')

import init
import core
import matplotlib.pyplot as plt
import numpy as np
import cgs
import userfun
import planets_properties as pp
import functions as f


#initialize the simulation
config_path='/home/lzx/CpdPhysics/'
pars=init.init_default_pars(config_path)

#initial size pf particles
Rdi=pars.dsystempars['Rdi']
#initial number of particles
nini=pars.dsystempars['nini']
#initial time
tini=pars.dsystempars['tini']
#initial planets properties
Ploca=pars.PlanetsLoca
Pmass=pars.PlanetsMass
Ptime=pars.PlanetsTime


system=init.sim_init(Rdi,nini,tini)

system.time=tini

tmax = pars.tmax
import pdb; pdb.set_trace()
data=userfun.data_process()

data.data_process(system.particles.Y2d,system.time,system.daction)

St0=[]
v_r0=[]
deltaT=[]
ntime=1
while system.time<tmax:

    Yt = system.update(tmax)

    # import pdb; pdb.set_trace()
    data.data_process(system.particles.Y2d,system.time,system.daction)
    v_r0.append(system.particles.get_stokes_number(system.disk,system.time)[1])
    deltaT.append(system.deltaT)

    # if 'add' in system.daction:
    #     import pdb; pdb.set_trace()
    
    # if system.ntime%10==0:
    #     data.plot_stuff(system.disk)
        # import pdb; pdb.set_trace()
    ntime+=1

TimeScale=f.dotMgTscale(data.radL,deltaT)

planets=pp.Planets(system.disk,system.particles,7*cgs.RJ,3e23)

# data.plot_stuff(system.disk)

