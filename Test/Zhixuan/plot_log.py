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
timeL = []
data = []
massL = []
peffL = []
with open('./log/planets.log', 'r') as file:

    columns = file.readline().strip().split()
    # Read the rest of the lines
    while True:
        line = file.readline()
        if not line:
            break  # End of file
        lm = line.strip().split()
        if len(lm)>=6:
            timeL.append(lm[0])
            massL.append(lm[4])
            peffL.append(lm[5])
        data.append(line.strip().split())    

plt.figure()
plt.yscale('log')
plt.xscale('log')
# 设置 y 轴刻度和标签
yticks_values = [0.0001, 0.01, 0.02, 0.03, 0.04, 0.05]
yticks_labels = ['0.0001', '0.01', '0.02', '0.03', '0.04', '0.05']
plt.yticks(yticks_values, yticks_labels)

# 设置 x 轴刻度和标签
xticks_values = [1e23, 1e24, 1e25, 1e26]
xticks_labels = ['1e23', '1e24', '1e25', '1e26']
plt.xticks(xticks_values, xticks_labels)
plt.scatter(massL, peffL)
import pdb;pdb.set_trace()

