import cgs

#put all parameters in this file
composL = ['silicates', 'H2O']
intrhoL = [3.0, 1.0]


#parameters fed into System
dsystempars = {}
dparticleprops = {'Rdi':0.01,  #initial size pf particles
                'nini':100,  #initial number of particles
                'initrule':'equallogspace'} #how particles are distributed

dtimesteppars = {'deltaTfraction':0.2, 'itgmethod':'Heun'}

dgasgrid = {'rinn':0.1*cgs.au,'rout':100*cgs.au}
gasmodel = 'prescribed'


tmax = 1e6 *cgs.yr

doJump = False
doPlanets = False
doIcelines = True

sfdmode = 'special'
resampleMode = 'fixed_resample' #'global_resample4'# None#

dresample = {'fdelS':0.04, 'fdelM':0.01}
