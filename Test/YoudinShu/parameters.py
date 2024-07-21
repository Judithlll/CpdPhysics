import cgs

#put all parameters in this file
composL = ['silicates']
intrhoL = [3.0]

#parameters fed into System
dsystempars = {}
dparticleprops = {'Rdi':0.1,  #initial size of particles in cm
                'nini':1000,  #initial number of particles
                'initrule':'equallogspace'} #how particles are distributed

dtimesteppars = {'deltaTfraction':0.01}

dgasgrid = {'rinn':0.1*cgs.au,'rout':800*cgs.au}
gasmodel = 'prescribed'
dragmodel = 'Epstein'

tmax = 2.1e5 *cgs.yr

doJump = False
doPlanets = False
doIcelines = False

resampleMode = None
dresample = {'fdelM':1e-3, 'fdelDM':1e-5}
