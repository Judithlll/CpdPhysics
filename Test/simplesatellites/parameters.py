import cgs

#put all parameters in this file
composL = ['silicates']
intrhoL = [3.0]

#parameters fed into System
dsystempars = {}
dparticleprops = {'Rdi':0.01,  #initial size pf particles
                'nini':100,  #initial number of particles
                'initrule':'equallogspace'} #how particles are distributed

dgasgrid = {'rinn':0.1*cgs.au,'rout':100*cgs.au}
gasmodel = 'prescribed'

tmax = 15e6 *cgs.yr

doJump = False
doPlanets = True
doResonance = True
doIcelines = False

resampleMode = 'splitmerge'
dresample = {'fdelS':0.04, 'fdelM':0.01}
