import cgs

#put all parameters in this file
composL = ['silicates']
intrhoL = [3.0]

#parameters fed into System
dsystempars = {}
dparticleprops = {'Rdi':0.0001,  #initial size of particles in cm
                'nini':1000,  #initial number of particles
                'initrule':'equallogspace'} #how particles are distributed

dtimesteppars = {'deltaTfraction':0.8, 'itgmethod':'Euler'}

dgasgrid = {'rinn':0.1*cgs.au,'rout':1e3*cgs.au}
gasmodel = 'prescribed'
dragmodel = 'Epstein'

tmax = 1.01e7 *cgs.yr

doJump = False
doPlanets = False
doIcelines = False

resampleMode = 'global_resample'
dresample = {'fdelS':0.02, 'fdelM':0.002, 'fdelDM':1e-7}
