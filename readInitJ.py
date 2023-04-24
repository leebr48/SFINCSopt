#!/usr/bin/env python

from os import environ
environ['SLURM_JOB_ID'] = '-666' # dummy value for other functions to work properly

import numpy as np
from main import LTargets
from main import targetWeights
from main import extract_value_subdirs

initDirname = 'prep'

flatValues = extract_value_subdirs(initDirname).flatten()
flatTargets = np.array(LTargets).flatten()
flatWeights = np.array(targetWeights).flatten()
    
out = np.sum(flatWeights * (flatValues - flatTargets)**2)

print('The initial value of the objective function is {}.'.format(out))
