import os
import sys 
sys.path.insert(0, '../../main/')
import subprocess
import runcpd as rc 

execfile('../../main/runcpd.py')

doEvo = rc.doEvo
data = rc.userfun.data
system = rc.userfun.system
if doEvo and False:
    data = userfun.data
    data.get_plot_list(doParticles = False)
    #store system components as pickles
    fileio.store_class(system, 'system')
    fileio.store_class(data, 'data')##CWO: this and all stuff below does not seem to be general. Move to do_stuff perhaps

