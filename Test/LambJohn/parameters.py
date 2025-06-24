import cgs

#put all parameters in this file
composL = ['silicates']
intrhoL = [3.0]

#parameters fed into System
dsystempars = {}
dparticleprops = {'Rdi':0.000001,  #initial size of particles in cm
                'nini':250,  #initial number of particles
                'initrule':'equallogspace'} #how particles are distributed

dtimesteppars = {'deltaTfraction':0.8, 'itgmethod':'Heun'}

dgasgrid = {'rinn':0.1*cgs.au,'rout':2e3*cgs.au}
gasmodel = 'prescribed'
dragmodel = 'Epstein'

tmax = 1.01e7 *cgs.yr

doJump = False
doPlanets = False
doIcelines = False

sfdmode = 'sfd_kernel' #'sfd_face'#'sfd_zxl'#'special' #'fixed_bin'# 'sfd_spline'# 'simple'#

resampleMode =  'face_splitmerge'#'new_splitmerge_zxl'# 'local_splitmerge'#  'global_resample'# 'global_resample4'#  'splitmerge'#
dresample = {'fdelS':0.05, 'fdelM':0.01, 'fdelDM':1e-7}
