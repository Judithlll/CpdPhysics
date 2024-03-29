import cgs

#put all parameters in this file
composL = ['silicates', 'H2O']
intrhoL = [3.0, 1.0]

#parameters fed into System
dsystempars = {'timeini': 0.0,
               'rhoPlanet': 1.9}

dparticleprops = {'Rdi':0.01,  #initial size pf particles
                'nini':100,
                'initrule':'equallogspace'}  #initial number of particles

#[24.01.08]LZX: don't know whether this is a good idea
dtimesteppars = {'deltaTfraction': 0.2,
                 'timestepn': 1} #belong to integrator

#TBD: implement integrator properties
dintegratorpars = {'name':'RK5', 'timestepn':3}

dgasgrid = {'rinn':6*cgs.RJ,'rout':27*cgs.RJ}
gasmodel = 'prescribed'

tmax = 3.1e6 *cgs.yr

doPlanets = True
doIcelines = True
doResonance = True
