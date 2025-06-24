import cgs

#put all parameters in this file
composL = ['silicates', 'H2O']
intrhoL = [3.0, 1.0]


#parameters fed into System
dsystempars = {}
dparticleprops = {'Rdi':0.01,  #initial size pf particles
                'nini':1000,  #initial number of particles
                'initrule':'equallogspace'} #how particles are distributed

dtimesteppars = {'deltaTfraction':0.2, 'itgmethod':'Heun'}

dgasgrid = {'rinn':0.1*cgs.au,'rout':100*cgs.au}
gasmodel = 'prescribed'


vc = {'icy':1.e2,
      'silicates':5.e2}

tmax = 1e6 *cgs.yr

doJump = False
doPlanets = False
doIcelines = True

sfdmode = 'sfd_kernel'
#resampleMode = 'face_splitmerge' #'global_resample4'# None#
resampleMode = 'new_splitmerge_zxl'
dresample = {'Xspecial':10}
