#!/usr/bin/env python
import sys
sys.path.append('../../main/')
import fileio
import numpy as np
import cgs
import matplotlib.pyplot as plt
import parameters as pars
import core

argL = sys.argv
data=[]
if len(argL) > 1:
    for arg in argL[1:]:
        data .append( fileio.load_class('./pickles/',arg))
else:
    data .append( fileio.load_class('./pickles/', 'data.pickle'))
print('[plot]: data.pickle has been loaded')
system = fileio.load_class('./pickles/','system.pickle')

data[0].plot_iceline()
data[0].plot_planet_evolution()
data[0].plot_St(tL =np.array(np.array([1e6,5e6,6e6,10e6,15e6])*cgs.yr))
data[0].plot_disk(np.array([0.,3.e6, 10.e6,30.e6])*cgs.yr)
import pdb;pdb.set_trace()
#plot St in different dust-to-gas ratio 
#plt.figure(figsize=(10,7))
#plt.ylabel('Stokes number')
#plt.xlabel('Distance from central planet [$R_J$]')
#plt.yscale('log')
#plt.xscale('log')
#labelL=[r'$\dot{M}_d/\dot{M}_g$ = 1',r'$\dot{M}_d/\dot{M}_g$ = 0.1',r'$\dot{M}_d/\dot{M}_g$ = 0.01',r'$\dot{M}_d/\dot{M}_g$ = 0.001']
#for data in dataL:
#    plt.plot(data.radL[99]/cgs.RJ, data.StD

#plot St w and w/o frag:
plt.figure(figsize=(10,7))
#plt.title('Stokes number')
plt.ylabel('Stokes number')
plt.xlabel('Distance from central planet [$R_J$]')
plt.yscale('log')
plt.xscale('log')
labelL=['without fragmentation', 'with fragmentation']
for i,s in enumerate(data):
    plt.plot(s.particles.locL/cgs.RJ, s.particles.St, 'x-',label = labelL[i])
plt.axvline(data[0].icelineL[0].loc/cgs.RJ, color='gray', linestyle='dashed', label='water iceline')
plt.axvline(5.89, color='gray', linestyle='dotted', label='inner edge')
plt.xticks([5.89,6,10,15,20,27],['5.89','','10','15','20','27'])
plt.legend()
plt.savefig('./plot/St_compare.jpg')
#data.plot_pebble_Sigma(tL =np.array([0]), mode='particle')
#data.plot_pebble_Sigma(tL =np.array([1e6*cgs.yr, 2e6*cgs.yr,3e6*cgs.yr]),mode='particle')
#data.plot_particles_number()
#data.plot_disk(np.array([50,1e6,3e6])*cgs.yr)
#data.plot_planet_evolution()

#data.plot_planet_migration()
#data.plot_jumpT()
#data.plot_disk_profile()
#plt.figure()
#plt.plot(data.timeL, data.planetslocL.T[1]/data.planetslocL.T[0])
#plt.plot(data.timeL, 2*np.ones_like(userfun.data.timeL))
#plt.savefig("pratio.jpg")

