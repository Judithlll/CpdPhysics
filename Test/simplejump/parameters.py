import cgs

#put all parameters in this file
composL = ['silicates']
intrhoL = [3.0]

#parameters fed into System
dsystempars = {}
#make the initial size larger to let them drift faster
dparticleprops = {'Rdi': 0.05,  #initial size pf particles
                'nini':100,  #initial number of particles
                'initrule':'equallogspace'} #how particles are distributed

dtimesteppars = {'itgmethod': 'Heun',
                'deltaTfraction': 0.5}
#use the size of typical CPD's, to make the con3 for jump can be met
dgasgrid = {'rinn': 6*cgs.RJ,'rout':100*cgs.RJ}
gasmodel = 'prescribed'

tmax = 1e6 *cgs.yr

doJump = True
doPlanets = True
doResonance = True
doIcelines = False

resampleMode = 'splitmerge'
dresample = {'fdelS':0.04, 'fdelM':0.01}
