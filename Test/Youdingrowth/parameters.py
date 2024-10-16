import cgs

#put all parameters in this file
composL = ['silicates']
intrhoL = [3.0]

#parameters fed into System
dsystempars = {}
dparticleprops = {'Rdi':0.01,  #initial size of particles in cm
                'nini':2000,  #initial number of particles
                'initrule':'equallogspace'} #how particles are distributed
dtimesteppars = {'deltaTfraction':0.5, 'itgmethod':'RK4'}

dgasgrid = {'rinn':0.1*cgs.au,'rout':800*cgs.au}
gasmodel = 'prescribed'
dragmodel = 'Epstein'

tmax = 1.01e6 *cgs.yr

doJump = False
doPlanets = False
doIcelines = False

resampleMode = 'splitmerge'
dresample = {'fdelS':1.04, 'fdelM':1e-7, 'fdelDM':3e-3}
