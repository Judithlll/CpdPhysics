#!/usr/bin/env python
import sys
sys.path.append('../../main/')
import fileio
import numpy as np
import cgs
import matplotlib.pyplot as plt
import parameters as pars
import core
import userfun


userfun.make_animation('sfdevol.mp4','./sfdevol/')
argL = sys.argv
data=[]
if len(argL) > 1:
    for arg in argL[1:]:
        data .append( fileio.load_class('./pickles/',arg))
else:
    data .append( fileio.load_class('./pickles/', 'data.pickle'))
print('[plot]: data.pickle has been loaded')
#system = fileio.load_class('./pickles/','system.pickle')


import pdb;pdb.set_trace()
data[0].plot_deltaT()
data[0].plot_planet_evolution()
data[0].plot_peff_log()
#data[0].plot_iceline()
data[0].plot_growth_timescale()
#data[0].plot_satepart()
#data[0].plot_disk(np.array([0.0,1e6,10e6])*cgs.yr)
data[0].plot_pebble_Sigma(tL =np.array(np.array([1e6,5e6,9e6])*cgs.yr))


data[0].plot_St(tL =np.array(np.array([1e6,5e6,9e6])*cgs.yr))
data[0].make_animation(mp4name= 'satepart.mp4', path='./plot/satepart')
#data[0].plot_St_t()
#data[0].make_animation(mp4name= 'St_t.mp4',path='./plot/St_t')
data[0].plot_jumpT()
data[0].plot_vr(tL =np.array(np.array([5e6,6e6,10e6])*cgs.yr))
data[0].plot_disk(np.array([0,3e6,10e6,30e6])*cgs.yr)

plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Iceline location [R_J]')
plt.xlabel('Time [yr]')
plt.xlim(1e3,20e6)
labelL=[r'$r_{gg} = 1.7\times 10^{-6}$',r'$r_{gg} = 1.7\times 10^{-7}$',r'$r_{gg} = 1.7\times 10^{-8}$',]
for i,d in enumerate(data):
    plt.plot(np.array(d.timeL)/cgs.yr, d.icelineslocL.T[0]/cgs.RJ, '.-', label = labelL[i])

for loc in [5.89, 9.38, 15.0, 26.3]:
    plt.axhline(loc, linestyle='dashed', color ='gray')

plt.yticks([5.89, 9.38, 15.0, 26.3],['5.89', '9.38', '15.0', '26.3'])

plt.legend()
plt.savefig('./plot/iceline.jpg')
plt.close()
#data[0].plot_disk(np.array([0.,3.e6, 10.e6,30.e6])*cgs.yr)
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
#data.plot_planet_evolution()

#data.plot_planet_migration()
#data.plot_jumpT()
#data.plot_disk_profile()
#plt.figure()
#plt.plot(data.timeL, data.planetslocL.T[1]/data.planetslocL.T[0])
#plt.plot(data.timeL, 2*np.ones_like(userfun.data.timeL))
#plt.savefig("pratio.jpg")

