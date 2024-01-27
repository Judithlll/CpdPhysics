import pickle

def store_class(class_name, filename):
    with open('./pickles/'+filename+'.pickle','wb') as f:
        pickle.dump(class_name, f)
    print('[fileio]: '+filename+' class has been stored into system.pickle')


def load_class(path,filename):
    with open(path+filename,'rb') as f:
        class_object = pickle.load(f)
    return class_object

## SAVE the current state (-> pickle) --> run1234.pickle

## LOAD a state.


## write formated tables (??)

## simple plots (??)
