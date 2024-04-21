#!/usr/bin/env python
import sys
sys.path.append('../../main/')
import fileio
import numpy as np
import cgs
import matplotlib.pyplot as plt

argL = sys.argv

if len(argL) > 1:
    data = fileio.load_class('./pickles/',argL[1])
else:
    data = fileio.load_class('./pickles/', 'data.pickle')
print('[plot]: data.pickle has been loaded')
system = fileio.load_class('./pickles/','system.pickle')
import pdb;pdb.set_trace()
#data.plot_disk([0,3e6*cgs.yr,10e6*cgs.yr])
#data.plot_pebble_Sigma(tL =np.array([0]), mode='particle')
data.plot_St(tL =np.array([1e6*cgs.yr, 2e6*cgs.yr,3e6*cgs.yr]))
data.plot_pebble_Sigma(tL =np.array([1e6*cgs.yr, 2e6*cgs.yr,3e6*cgs.yr]),mode='particle')
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

