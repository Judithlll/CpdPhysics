import pickle
import os

def store_class(class_name, filename):
    if not os.path.exists('./pickles/'):
        os.makedirs('./pickles')
    with open('./pickles/'+filename+'.pickle','wb') as f:
        pickle.dump(class_name, f)
    print('[fileio]: '+filename+' class has been stored into system.pickle')


def load_class(path,filename):
    with open(path+filename,'rb') as f:
        class_object = pickle.load(f)
    return class_object

class store_system(object):
    def __init__(self, SystemPd):
        for k, v in SystemPd.items() :
            setattr(self, k, v)

def store_state(system):
    """
    this function store the system components respectively
    """
#store  finally state of system, not sure if this is userful
    store_class(system.particles, 'particles')
    store_class(system.gas, 'gas')
    for i in range(system.nplanet):
        store_class(system.planetL[i], 'planet'+str(i+1))

    for i in range(system.niceline):
        store_class(system.icelineL[i], 'iceline'+str(i+1))

#put some necessary properties into store_systemclass and store it
    nece_pd = {'time':system.time, 'jumpT': system.jumpT, 'ntime':system.ntime, 'rhoPlanet':system.rhoPlanet, 'nplanet': system.nplanet, 'niceline':system.niceline, 'milestones':system.milestones}

    system_store = store_system(nece_pd)
    store_class(system_store, 'system_store')

## SAVE the current state (-> pickle) --> run1234.pickle

## LOAD a state.


## write formated tables (??)

## simple plots (??)
