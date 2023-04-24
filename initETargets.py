#!/usr/bin/env python

from os import environ
environ['SLURM_JOB_ID'] = '-666' # dummy value for other functions to work properly

import numpy as np
from main import extract_value_subdirs

initDir = 'prep'

vals = extract_value_subdirs(initDir)
L11e = vals[0::2][:,0]

np.set_printoptions(precision=100)
print('The initial values of L11e, in flux surface order, are {}.'.format(L11e))
