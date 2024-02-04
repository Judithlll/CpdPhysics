import pickle
import os
import sys
import parameters as pars

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
    Not used now!
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

def find_pickles(path):
    files = os.listdir(path)
    pickles = []
    for file in files:
        if '.pickle' in file:
            pickles.append(file)
    return pickles

def check_pickles(pickles):
    pk1 = 'particles.pickle' in pickles
    pk2 = 'gas.pickle' in pickles
    pk3 = 'planet1.pickle' in pickles
    pk4 = ('iceline1.pickle' in pickles) or (pars.doIcelines == False)
    pk5 = ('system_store.pickle' in pickles) or (pars.doPlanets == False)
    if pk1 *pk2*pk3*pk4*pk5:
        tf = True
    else:
        tf = False

    return tf

def load_state(pickle_path):
    """
    not used now!
    """
    if not os.path.exists(pickle_path):
        print("WARNING: you don't have pickle directory [runcpd_fromfile]")
        sys.exit()

    pickles = find_pickles(pickle_path)
    #need to judge whether the pickles is complete

    completeness = check_pickles(pickles)
    if not completeness:
        print('[runcpd_fromfile]: the pickles file is not complete' )
        sys.exit()

    classes = {}
    for pcs in pickles:
        classes[pcs.rsplit('.',1)[0]] = load_class(pickle_path, pcs)

    return classes
## SAVE the current state (-> pickle) --> run1234.pickle

## LOAD a state.


## write formated tables (??)

## simple plots (??)
