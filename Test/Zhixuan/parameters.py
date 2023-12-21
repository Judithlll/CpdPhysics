import cgs

#put all parameters in this file

#parameters fed into System
dsystempars = {'Rdi':0.01,  #initial size pf particles
                'nini':12,
                'ice_frac':0.5,
                'diskmass':0.01*cgs.MJ}  #initial number of particles

dgasgrid = {'rinn':5.89*cgs.RJ,'rout':27*cgs.RJ}
gasmodel = 'gridstatic' #'prescribed'

tmax = 1.2e2 *cgs.yr

doPlanets = True
