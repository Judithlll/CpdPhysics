import numpy as np
import sys
import parameters as pars
sys.path.append('../../main/')
import core
import init
import cgs

argL = sys.argv #maybe use later
calldir = init.init_default_pars (argL[0]) #directory from which this is called (maybe we need later)

#this initializes the system...
system = init.sim_init (calldir,pars.dsystempars)


while system.time<pars.tmax:

    Yt = system.update (pars.tmax)

    print('hello', system.time/cgs.yr)

