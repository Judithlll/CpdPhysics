import cgs

#put all parameters in this file
composL = ['silicates', 'H2O']
intrhoL = [3.0, 1.0]

#parameters fed into System
dsystempars = {'timeini': 0.0,
               'rhoPlanet': 1.9}

dparticleprops = {'Rdi':0.01,  #initial size pf particles
                'nini':700,
                'initrule':'equallogspace'}  #initial number of particles

#[24.01.08]LZX: don't know whether this is a good idea
dtimesteppars = {'deltaTfraction': 0.2,
                 'timestepn': 1} #belong to integrator

#TBD: implement integrator properties
dintegratorpars = {'name':'RK5', 'timestepn':3}

#give the jumpfraction here, 
#which means how accrute the jump you want to have, the smaller, the more accrute
#and you can also define jumpfraction for specific timescale
jumpfracD={'general': 0.2,
          'PlanetsRes': 0.5,}

dgasgrid = {'rinn':5.89*cgs.RJ,'rout':100*cgs.RJ}
gasmodel = 'prescribed'

dgasprop = {'alpha': 1e-4, 
            'frac': 0.2/1.5, # the accretion mass with unit: M_J/Myr
            'dgratio': 0.002,
            'rgg': 1.7e-7,
            'Mcp0':0.4*cgs.MJ,
            'meanmol': 2.34,
            'sigmol':2e-15,
            'tdep': 2.9e6*cgs.yr,
            'tgap': 1.e6*cgs.yr
            }

tmax = 12018782.555918233*cgs.yr

#critical velocity for fragmentation
vc = {'icy': 5.e3,
      'silicates': 5.e2}

resampleMode = 'splitmerge'
doJump = True
doPlanets = True
doIcelines = True
doResonance = True
