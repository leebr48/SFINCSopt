#!/usr/bin/env python

import h5py
from os.path import join

def read_transportMatrix(dirname):
    with h5py.File(join(dirname, 'sfincsOutput.h5'),'r') as f:
        ret1 = float(f['transportMatrix'][()][0][0])
        ret2 = float(f['transportMatrix'][()][2][0])
        ret = [ret1, ret2]
    return ret

if __name__=="__main__":
    import sys
    if len(sys.argv) > 1:
        dirname = sys.argv[1]
    else:
        dirname = '.'
    ret = read_transportMatrix(dirname)
    print(ret, flush=True)
