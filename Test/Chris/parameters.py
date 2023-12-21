import cgs

#put all parameters in this file

#parameters fed into System
dsystempars = {'Rdi':0.01,  #initial size pf particles
                'nini':100}  #initial number of particles

dgasgrid = {'rinn':6*cgs.RJ,'rout':27*cgs.RJ}
gasmodel = 'prescribed'

tmax = 1e2 *cgs.yr

doPlanets = True
doIcelines = False

