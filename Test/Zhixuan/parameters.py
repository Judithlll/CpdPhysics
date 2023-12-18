import cgs

#put all parameters in this file

#parameters fed into System
dsystempars = {'Rdi':0.01,  #initial size pf particles
                'nini':100,
                'ice_frac':0.5}  #initial number of particles

dgasgrid = {'rinn':6*cgs.RJ,'rout':27*cgs.RJ}
gasmodel = 'gridstatic' #'prescribed'

tmax = 1e2 *cgs.yr

doPlanets = True
