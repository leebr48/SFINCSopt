#!/usr/bin/env python

import numpy as np

iters = []
DOFvals = []
Lvals = []
Jvals = []
with open('log.txt', 'r') as f:
    lines = f.readlines()[1:]
    for line in lines:
        splitline = line.strip().split(',')
        
        iters.append(int(splitline[0].strip()))
        
        DOFvals.append([float(item) for item in splitline[1].strip().split()])
        
        Lvals.append([float(item) for item in splitline[2].strip().split()])

        Jvals.append(float(splitline[3].strip()))

minJind = np.argmin(Jvals)
minJ = Jvals[minJind]
minIter = iters[minJind]

print('The minimum objective value of {} is achieved on iteration {}.'.format(minJ, minIter))
