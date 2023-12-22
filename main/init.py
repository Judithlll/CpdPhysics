import sys
#sys.path.append('/home/lzx/CpdPhysics/Test/Zhixuan/') ## This is not a good idea...
import parameters as pars 
import functions
import core
import userfun

"""
We need the parameters in defaults.txt, which is a dict or list.

if it's a dict, should use ';' to split different elements, use ':' to split keys and values.

if it's a list,
"""
def init_default_pars (fcall):
    """
    iniitialize parameters to None if they have not been defined 

    [23.12.03] -- This file is copied from Chris' /NewLagrance
    """
    ix = fcall.rfind('/')
    calldir = fcall[:ix+1]
    fname = calldir+'../config/defaults.txt'
    with open(fname,'r') as f:
        line = f.readline()
        while line:
            if line[0]!='#' and len(line)>=5:
                key, sep, val, *dumL = line.split()  

                #dictionaries
                if val[0]=='{':
                    dout = {}
                    sss = val[1:-1]
                    subL = sss.split(';')
                    for sub in subL:
                        key1, val1 = sub.split(':')
                        dout[key1] = functions.determine_type(val1)

                    #update it
                    if hasattr(pars, key):
                        dout.update(getattr(pars,key))

                    setattr(pars, key, dout)

                #if it's already there... dont use defaults
                elif hasattr(pars, key) is False:

                    out = functions.determine_type (val)
                    setattr(pars, key, out)


            line = f.readline()

    #in pars there must be a 'addgasL'
    if hasattr(pars, 'addgasL') is False:
        pars.addgasL = []

    return pars  #why return calldir? why not pars


def sim_init (dsystempars={},*args):
    """
    The idea is that the system object is defined here... TBD
    """

    #let's System be initialized this way with keyword parameters... 
    #please don't change it again
    system = core.System(**dsystempars)

    #add the planets
    if pars.doPlanets is None:
        pars.doPLanets = False
    elif pars.doPlanets:
        #TBD: perhaps better to work with dictionaries
        tarr, radarr, mplanarr, fcomparr = userfun.init_planets ()

        #TBD: check w/r fcompL is normalized

        nplan = len(radarr)

        #planet list
        planL = []
        for k in range(nplan):
            planL.append(core.PLANET (tarr[k], radarr[k], mplanarr[k], fcomparr[k]))

        system.planetL = planL
        system.nplanet = nplan
    else:
        # the planet list
        system.planetL = [] 
        system.nplanet = 0

    if pars.doIcelines is None:
        pars.doIcelines = False
    elif pars.doIcelines:
        speciesarr,temparr = userfun.init_icelines()
        icelineL = []
        icelineLocL=[]
        niceline = len(speciesarr)
        for i in range (niceline):
            iceline=core.ICELINE (speciesarr[i], temparr[i])
            iceline.get_icelines_location(system.gas,system.time)
            icelineL.append(iceline)

        system.icelineL = icelineL
        system.niceline = niceline
        
    else:
        system.icelineL = []
        system.niceline = 0
    # core.ICELINE.get_icelines_location(system)

    #TBD: add icelines to system...
    #
    #returns: icelines species, condensation temperature...
    #
    #diceline = userfun.init_icelines ()


    #Then assign particle properties

    return system

