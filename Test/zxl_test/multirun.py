import sys 
sys.path.insert(0,'../../main/')
import subprocess as sp 
import cgs
import os 
import parameters as pars 
import numpy as np 

def change_par(parname, parval, filename='parameters.py'):
    # check the type of the parameter we want to change
    if ':' in parname:
        idx = parname.find(':')
        valname = parname[idx+1:]
    else:
        valname = parname 
        idx = None 

    parvalo = getattr(pars, parname[0:idx])
    

    with open(filename, 'r') as f:
        lines = f.readlines()

    for i,line in enumerate(lines):
        
        if valname in line:
            if isinstance(parvalo, list) or isinstance(parvalo, float) or isinstance(parvalo, np.ndarray):
                lines[i] = valname + ' = ' + str(parval)+' \n'
             
            if isinstance(parvalo, dict):
                idx_colon = line.find(':')
                lines[i] = line[0:idx_colon] + ':' + str(parval) +', '+'\n'
    
    with open(filename, 'w') as f:
        f.writelines(lines)

    print('[multirun]: the parameter {} has been modified to {}'.format(parname, parval))
            
    

parname = 'tmax'
parL = [11.e6*cgs.yr]

for par in parL: 
    #1. change the parameters and write in parameters.py 
    change_par(parname , par)

    #2. run ../../main/runcpd.py 
    print('[multirun]: begin to run the system')
    import pdb;pdb.set_trace()
    sp.run(['python', '../../main/runcpd.py'])

    #3. change the name of pickles according to the parameters used 
    dirname = './pickles/{}+{}'.format(parname, par)
    sp.run(['mkdir', dirname])
    sp.run(['mv', './pickles/data.pickle', dirname+'/data.pickle'])
    sp.run(['mv', './pickles/system.pickle', dirname+'/system.pickle'])
    sp.run(['cp', './parameters.py', dirname+'/parameters.py'])
    print('[multirun]: the pickles and parameters have been put in the specifical directory' )

    #plot the multi parameter figure
    #sp.run(['python', 'plot.py'])
#TBD: change some parameters at one time 
#TBD: plot multi parameters figures. 

