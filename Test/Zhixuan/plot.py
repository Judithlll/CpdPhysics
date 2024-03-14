#!/usr/bin/env python
import sys
sys.path.append('../../main/')
import fileio

argL = sys.argv

if len(argL) > 1:
    data = fileio.load_class('./pickles/',argL[1])
else:
    data = fileio.load_class('./pickles/', 'data.pickle')


#data.plot_particles_number()
#data.plot_disk(np.array([50,1e6,3e6])*cgs.yr)
data.plot_planet_evolution()

#data.plot_planet_migration()
#data.plot_jumpT()
#data.plot_disk_profile()
import pdb;pdb.set_trace()
plt.figure()
plt.plot(data.timeL, data.planetslocL.T[1]/data.planetslocL.T[0])
plt.plot(data.timeL, 2*np.ones_like(userfun.data.timeL))
plt.savefig("pratio.jpg")

