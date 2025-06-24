import cgs

#put all parameters in this file
composL = ['silicates', 'H2O']
intrhoL = [3.0, 1.0]

#parameters fed into System
dsystempars = {'timeini': 0.0,
               'rhoPlanet': 1.9}

dparticleprops = {'Rdi':0.01,  #initial size pf particles
                'nini':200,
                'initrule':'equallogspace'}  #initial number of particles

#[24.01.08]LZX: don't know whether this is a good idea
dtimesteppars = {'itgmethod': 'Heun',
                 'deltaTfraction': 0.2,
                 'coltimefrac':4.0,
                 'timestepn': 1, 
                 } #belong to integrator


#give the jumpfraction here, 
#which means how accrute the jump you want to have, the smaller, the more accrute
#and you can also define jumpfraction for specific timescale
jumpfracD={'general': 0.2,
          'PlanetsRes': 0.5,}

dgasgrid = {'rinn':5.89*cgs.RJ,'rout':100*cgs.RJ}
gasmodel = 'prescribed'

dragmodel = 'Epstein'

dgasprop = {'alpha':5e-05, 
            'frac': 0.2/1.5, # the accretion mass with unit: M_J/Myr
            'dgratio': 0.0016,
            'rgg': 0.88e-7,
            'Mcp0':cgs.MJ,
            'Rcp0': 2.*cgs.RJ,
            'meanmol': 2.34,
            'sigmol':2e-15,
            'tdep': 3.e6*cgs.yr,
            'tgap': 1.e6*cgs.yr
            }

fraginit = False#True
#tmax = 347140156224193.44 
tmax = 3e4*cgs.yr
#critical velocity for fragmentation
vc={'icy': 1.e2,
    'silicates': 1.e2, 
    }

sfdmode = 'simple'
resampleMode = 'new_splitmerge_zxl' #'Nplevel'
dresample = {'fchange':0.9,'fdelDM':1e-3, 'Xspecial':10}
#dresample = {'fdelS':1.,'fdelM':0.001,'fdelDM':1e-7, 'fspec':0.002}
doJump = False
doPlanets = True
planetaddmode = 'capture'
doIcelines = True
doResonance = True
