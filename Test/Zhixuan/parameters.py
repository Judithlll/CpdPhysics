import cgs

#put all parameters in this file

#parameters fed into System
dsystempars = {'Rdi':0.01,  #initial size pf particles
                'nini':20,
                'time':0.0
                }
#               'PlanetsLoca':[10*cgs.RJ], 'PlanetsMass':[3e23], 'PlanetsTime':[3e6*cgs.yr2s]}  #initial number of particles

tmax = 1e2 *cgs.yr
