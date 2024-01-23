import pickle

def store_system(system):
    with open('system.pickle','wb') as f:
        pickle.dump(system, f)
    print('system class has been stored into system.pickle')


def load_system(path):
    with open(path+'system.pickle','rb') as f:
        system = pickle.load(f)
    return system

## SAVE the current state (-> pickle) --> run1234.pickle

## LOAD a state.


## write formated tables (??)

## simple plots (??)
