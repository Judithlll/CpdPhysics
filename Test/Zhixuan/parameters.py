import cgs

#put all parameters in this file
composL = ['silicates', 'H2O']
intrhoL = [3.0, 1.0]

#parameters fed into System
dsystempars = {'timeini': 0.0,
               'rhoPlanet': 1.9,
               'timestepn': 3}

dparticleprops = {'Rdi':0.01,  #initial size pf particles
                'nini':100,
                'initrule':'equallogspace'}  #initial number of particles

dgasgrid = {'rinn':6*cgs.RJ,'rout':27*cgs.RJ}
gasmodel = 'prescribed'

tmax = 1e2 *cgs.yr

doPlanets = True
doIcelines = True
