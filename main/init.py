import parameters as pars
import functions
import core

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

    if hasattr(pars, 'addgasL') is False:
        pars.addgasL = []

    return calldir


def sim_init (calldir, dsystempars={}):
    """
    The idea is that the system object is defined here... TBD
    """

    system = core.System(**dsystempars)

    return system

