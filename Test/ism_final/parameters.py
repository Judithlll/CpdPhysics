import cgs

#put all parameters in this file
composL = ['silicates']
intrhoL = [3.0]

#parameters fed into System
dsystempars = {}
dparticleprops = {'Rdi': 0.01,  #initial size of particles in cm
                'nini':5000,  #initial number of particles
                'initrule':'equallogspace'} #how particles are distributed

dtimesteppars = {'deltaTfraction':0.2, 'itgmethod':'RK4'}

dgasgrid = {'rinn':10*cgs.au,'rout':2000*cgs.au}
gasmodel = 'prescribed'

dgasprop = {'alpha' : 10e-5, 'Mcp0':cgs.Msun}
dragmodel = 'Epstein'

tmax = 1e7 *cgs.yr

fixed_St = 10**(-3)

planetaddmode = 'capture'

doJump = False
doPlanets = True 
doIcelines = False
doResonance = False

resampleMode = 'global_resample'
dresample = {'fdelS':1e99, 'fdelM':1e-3, 'fdelDM':1e-5}

vc = {'silicates': 100}
