#!/usr/bin/env python
import sys
sys.path.append('../../main/')
import fileio
import numpy as np
import cgs
import matplotlib.pyplot as plt
import pandas as pd
argL = sys.argv

#pl_data = pd.read_csv('../../main/log/planets.log', delim_whitespace=True, header=0, engine = 'python-fwf')
l=0
planet_num = 4
timeL = []
data = []
massL = [[],[],[],[]]
peffL = [[],[],[],[]]
with open('./log/planets.log', 'r') as file:

    columns = file.readline().strip().split()
    # Read the rest of the lines
    while True:
        line = file.readline()
        if not line:
            break  # End of file
        lm = line.strip().split()
        
        timeL.append(float(lm[0]))

        num = (len(lm) -1.)/5.

        for j in range(int(num)):
            massL[j].append(float(lm[4+5*j]))
            peffL[j].append(float(lm[5+5*j]))
        for k in range(int(num),planet_num):
            massL[k].append(np.nan)
            peffL[k].append(np.nan)

        data.append(line.strip().split())    

import pdb;pdb.set_trace()
plt.figure()
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Satellite mass [g]')
plt.ylabel('Accretion effieiency [%]')
#yticks_values = [0.0001, 0.01, 0.02, 0.03, 0.04, 0.05]
#yticks_labels = ['0.0001', '', '', '', '', '0.05']
#plt.yticks(yticks_values, yticks_labels)

xticks_values = [1e23, 1e24, 1e25, 1e26]
xticks_labels = ['$1$', '1e24', '1e25', '1e26']
#plt.xticks(xticks_values, xticks_labels)
for i in range(planet_num):
    plt.plot(massL[i], np.array(peffL[i])*100, label=str(i))
plt.legend()
plt.savefig('./plot/peff.jpg')

