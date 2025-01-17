import sys
#sys.path.append('/home/lzx/CpdPhysics/Test/Zhixuan/') ## This is not a good idea...
import parameters as pars 
import functions
import core
import userfun
import os
import numpy as np

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
                        if ':' in sub:#cwo: sub can be empty
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

    return calldir


def Z_fixed (val):
    def Z_compos (rad):
        return val *np.ones_like(rad)
    return Z_compos


def construct_mask_icl (rice):
    def mask_icl (rad):
        return np.where(rad>rice, 1.0, -1.0)
    return mask_icl


def load_composdata (calldir, composL):
    """
    loads the composition data as dictionary

        dcompos

    The default is to load parameters from /config
    BUT this can be overwritten by:

        dcompos = userfuncs.init_compos()

    In particular, iceline [locations] must be explicitly specified by the user

    Finally, we construct a function Z_init(rad) which provides the initial composition. 
    The default is a constant, based on dcompos['Zinit']
    But this can also be specified by userfuncs.init_compos()
    """

    dcomposL = []
    for compos in composL:
        #look/load materials file
        fname = calldir+'../config/'+compos+'.txt'
        if os.path.isfile(fname):
            dcompos = functions.load_dict_from_file (fname)
        else:
            dcompos = {}

        if 'iceline' not in dcompos:
            dcompos['iceline'] = False #default

        #now look for user initializations...
        #that will overwrite the default parameters
        if hasattr(userfun, 'init_compos'):
            dadd = userfun.init_compos(compos)
            dcompos.update(dadd)

        
        #Z_init is some function that determines the initial abundance
        if 'Z_init' not in dcompos:
            try:
                dcompos['Z_init'] = Z_fixed(dcompos['Zinit'])
            except:
                print('[init.py]:initialization of >>', compos.upper(), '<< abundances unclear. Specify in userfuncs?' )
                print('[init.py]:aborting')
                sys.exit()

        dcomposL.append(dcompos)

    return dcomposL




def sim_init (calldir, dsystempars={},*args):
    """
    The idea is that the system object is defined here

    ...
    it also defines a composition list of dictionaries
        >> dcomposL <<

    which consists of:
        [species in composL] + [noncondensible gas] + [species in addgasL]

    """

    #let's System be initialized this way with keyword parameters... 
    #please don't change it again
    
    system = core.System(**dsystempars)


    #load composition data
    #this is being copied from /NewLagrange
    dcomposL = load_composdata (calldir, pars.composL)

    #next, my /NewLagrange initialization has some routines
    #that initialize the surface density, also accounting for icelines 
    #(actually this should be merged with your superparticles initialization)

    #now add gas and additional vapor components to the composition object
    dcomposL += load_composdata (calldir, ['gas'] +pars.addgasL)
    #import pdb; pdb.set_trace()
    #construct the gas species we need to follow
    gasL = [dcomp['name'] for dcomp in dcomposL if dcomp['iceline']==True]
    gasL.append('gas')      #always
    gasL += pars.addgasL    #specified by user


    system.gasL = gasL

    #add the planets
    if pars.doPlanets is None:
        pars.doPLanets = False
    elif pars.doPlanets: #from file?
        #TBD: perhaps better to work with dictionaries
        tarr, radarr, mplanarr, fcomparr = userfun.init_planets ()

        nplan = len(radarr)

        #planet list
        planL = []
        for k in range(nplan):
            planL.append(core.PLANET (tarr[k], radarr[k], mplanarr[k], fcomparr[k]))
            system.milestones[tarr[k]] = 'insert_planet'
        system.planet_candidate = planL
    else:
        # the planet list
        system.planet_candidate = [] 
    system.planetL = []
    system.nplanet = 0
    system.planetD = {}
    #[24.02.02] cwo
    #for simplicity, let's always enter this block ...
    #such that we always have a "construct_mask_icl" function

    icelineL = []
    for dcompos in dcomposL:
        if dcompos['iceline']==True:
            iceline = core.ICELINE (dcompos['name'], dcompos['iceline_temp'])  
            iceline.get_icelines_location(system.gas,system.time, bounds=[system.rinn, system.rout])      
            icelineL.append(iceline)
           #(initial) position of iceline
            rice = iceline.loc
            
        elif dcompos['iceline'] == 'None': #all vapor
            rice = np.inf
        elif dcompos['iceline'] == False: #all refractory
            rice = 0.0

        dcompos['mask_icl'] = construct_mask_icl (rice)
        dcompos['rice_init'] = rice #[25.01.15] added this info

    system.icelineL = icelineL
    system.niceline = len(icelineL)
    system.dcomposL = dcomposL

    return system,gasL

