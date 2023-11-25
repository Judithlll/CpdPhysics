import core
import matplotlib.pyplot as plt
import numpy as np
import cgs
import userfun
import planets_properties as pp

#initial size pf particles
Rdi=0.01  
#initial number of particles
nini=10

system = core.System(Rdi,nini)


time0=1e6
tEnd=7e9

system.time=0.

tmax = 6e9
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


def dotMgTscale(radL,deltaTL):
    v=np.diff(radL*cgs.RJ,axis=1)/deltaTL
    vTimesr=v*radL[:,0:-1]*cgs.RJ
    first=np.diff(vTimesr,axis=1)/deltaTL[:-1]
    second=v[:,:-1]*np.diff(vTimesr,axis=1)/np.diff(radL,axis=1)[:,:-1]/cgs.RJ
    Tscale=1/((first-second)/v[:,:-1]/radL[:,:-2]/cgs.RJ)
    return Tscale

TimeScale=dotMgTscale(data.radL,deltaT)

planets=pp.Planets(system.disk,system.particles,7*cgs.RJ,3e23,system.time)
data.plot_stuff(system.disk)

