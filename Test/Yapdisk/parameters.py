import numpy as np 
import cgs

composL = ['silicates', 'H2O']
intrhoL = [3.0, 1.0]

#parameters fed into System
dsystempars = {}
dparticleprops = {'Rdi':0.01,  #initial size pf particles
                 'nini':100,  #initial number of particles
                 'initrule':'equallogspace'} #how particles are distributed

dtimesteppars = {'deltaTfraction':0.001, 'itgmethod':'RK4'}
dgasgrid = {'rinn':0.1*cgs.au,'rout':15*cgs.au}
gasmodel = 'gridstatic'

tmax = 3e6 *cgs.yr

doPlanets = False 
doIcelines = False 
doJump = False 
doResonance = False 

resampleMode = 'splitmerge'
dresample = {'fdelS':0.04, 'fdelM':0.01} 
