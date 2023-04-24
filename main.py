#!/usr/bin/env python

####################################################################################################
# User inputs
####################################################################################################

nconcur = 6 # Number of SFINCS calculations to run in parallel
nNodes = 3 # Number of nodes per SFINCS job; (total number of nodes / nconcur) if you wish to allocate all the nodes to SFINCS
nprocesses = 120 # Processes per SFINCS job - you need to request (nconcur * nprocesses) processors in the slurm file
memPerCpu = '2G' # Memory available to each CPU; if this is set incorrectly, the jobs may not run in parallel
surfs = [0.25, 0.50, 0.75] # = STELLOPT rho = SFINCS rN
speciesParams = [[-1, 1], [5.446170214863400e-04, 2.496634561495773]] # [list of Zs, list of mHats]; each sublist should be in species order and use appropriate SFINCS units
nHats = [[1.925759290562555, 1.925759290562555], [1.8000000003021746, 1.8000000003021746], [1.5952669602019312, 1.5952669602019312]] # Species densities in SFINCS units; first index is for the surface, second index is for the species
dnHatdrHats = [[-0.21746074335787643, -0.21746074335787643], [-0.26725895961109203, -0.26725895961109203], [-1.1496241130928588, -1.1496241130928588]] # Derivatives of species densities wrt rHat; first index is for the surface, second index is for the species
THats = [[14.925059121976009, 14.925059121976009], [10.000000001797774, 10.000000001797774], [4.999998821932152, 4.999998821932152]] # Species temperatures in SFINCS units; first index is for the surface, second index is for the species
dTHatdrHats = [[-10.625892808407443, -10.625892808407443], [-9.954507300393985, -9.954507300393985], [-9.929753618960675, -9.929753618960675]] # Derivatives of species temperatures wrt rHat; first index is for the surface, second index is for the species
Ers = [15, 15, 15] # Radial electric field in SFINCS units; index indicates flux surface
LTargets = [[[-0.00016516411764600346, 0], [0, 0]], [[-0.0005159076261538864, 0], [0, 0]], [[-0.0015055649421042677, 0], [0, 0]]] # Target thermal transport coefficient values; first index is for the surface, second index is for the species, third index is for the coefficient (L11 and L31, in that order)
targetWeights = [[[1, 1], [1, 1]], [[1, 1], [1, 1]], [[1, 1], [1, 1]]] # Weights for the terms in the objective function; the meanings of each dimension are the same as in <LTargets>
unfixMajorRadius = False # If True, the VMEC parameter RBC(0,0) will be included in the optimization space
unfixPHIEDGE = False # If True, the VMEC parameter PHIEDGE will be included in the optimization space
autobound = 0.1 # Domain of valid values for variables in the optimization space will be [(1 - autobound)*x, (1 + autobound)*x], where x denotes an arbitrary variable 
DEkwargs = {'args':(),
            'strategy':'best1bin',
            'maxiter':1000,
            'popsize':1,
            'tol':1e-2,
            'mutation':(0.6, 1),
            'recombination':0.7,
            'seed':None,
            'disp':False,
            'callback':None,
            'polish':True,
            'init':'halton',
            'atol':0,
            'updating':'deferred',
            'workers':1,
            'constraints':(),
            'x0':None} # NOTE: These were set using a mix of SciPy defaults, SciPy recommendations, and experience with STELLOPT - modifications may be necessary
BIG_NUM = 1e3 # Large number that is returned to the optimizer in place of L11 if SFINCS fails to run
fullVMECInputFileName = "input.vmec" # This VMEC input file specifies the initial equilibrium for the optimization
strippedVMECInputFileName = "stripped_input.vmec" # The should be the same as <fullVMECInputFileName>, but with all the boundary parameters deleted
baseSFINCSInputFileName = "input.namelist" # This SFINCS input file specifies the resolution parameters, calculation 'mode', etc.
logfileName = "log.txt" # Name for the optimization log file - the default is recommended

####################################################################################################
# Code begins below
####################################################################################################

# Import external libraries
import numpy as np
import os
from glob import glob
from shutil import copy
from scipy.optimize import differential_evolution
from simsopt.geo.surfacerzfourier import SurfaceRZFourier

# Import internal functions
from run_compute import run_compute
from modify_sfincs_input import changevar
from read_transportMatrix import read_transportMatrix

# Define some useful functions
def extract_value(dirname):
    try:
        ret = read_transportMatrix(dirname)
    except FileNotFoundError:
        # NOTE: could increase SFINCS resolution here
        print("Failed to extract data from: " + dirname, flush=True)
        ret = BIG_NUM
    
    return ret

def readPhiEdge(dirname, fileName='input.vmec'):
    inputFile = os.path.join(dirname, fileName)
    with open(inputFile, 'r') as f:
        
        for line in f.readlines():
            if 'phiedge' in line.lower():
                phiedge = float(line.split('!', 1)[0].split('=', 1)[-1].strip())
                return phiedge
        
        else: # If we fail
            raise IOError('No value could be found for PHIEDGE.')

def dirname_from_indices(i,n):
    dirname = str(i).zfill(5) + "_" + str(n).zfill(3)
    
    return dirname

def prepare_next_simulation(x, n = 0):
    global my_iter
    global strippedVMECInputFileName

    dirname = dirname_from_indices(my_iter,n)
    
    if unfixPHIEDGE:
        surf.x = x[:-1]
    else:
        surf.x = x

    boundary_string = surf.get_nml()
    boundary_string = boundary_string.rsplit('\n',2)[0]
    boundary_string = boundary_string.split('\n',3)[-1] + "\n"
    
    if not os.path.exists(dirname):
        os.mkdir(dirname)
        copy("input.namelist",dirname)
        new_vmec_input = dirname + "/input.vmec" 
        copy(strippedVMECInputFileName,new_vmec_input)
         
        # modify vmec file
        with open(new_vmec_input,'r') as f:
            lines = f.readlines()
        Nlines= len(lines)
        if lines[Nlines-1].strip() == "/":
            lines.insert(Nlines-1,boundary_string)
        else:
            lines.insert(Nlines-2,boundary_string)
        
        with open(new_vmec_input,'w') as f:
            f.write("".join(lines))
        if unfixPHIEDGE:
            changevar(dirname, 'INDATA', 'PHIEDGE', x[-1], inputname = "input.vmec")
    
    else:
        print("dir already exists: " + dirname, flush=True)
    
    return dirname

def extract_value_subdirs(dirname):
    subsubdirs = sorted(glob(dirname + "/*/*/"))
    values = []
    for d in subsubdirs:
        values.append(extract_value(d))
    
    return np.array(values)

def opt(x):
    global my_iter
    global logfileName
    
    dirname = prepare_next_simulation(x)
    dirnames = [dirname]
    run_compute(dirnames, surfs, specArray, nHats, dnHatdrHats, THats, dTHatdrHats, Ers, nconcur, nprocesses=nprocesses, memPerCpu=memPerCpu, nNodes=nNodes, sfincs_input=baseSFINCSInputFileName, vmec_output=woutFileName)
   
    flatValues = extract_value_subdirs(dirname).flatten()
    flatTargets = np.array(LTargets).flatten()
    flatWeights = np.array(targetWeights).flatten()
    
    if len(flatValues) != len(flatTargets): # Some run or read error occurred
        flatValues = np.repeat(BIG_NUM, len(flatTargets))

    out = np.sum(flatWeights * (flatValues - flatTargets)**2)
    
    print('Iteration number: ', my_iter, flush=True)
    print('L values: ', flatValues, flush=True)
    print('Objective function value: ', out, flush=True)
    print('\n', flush=True)

    with open(logfileName,'a') as f:
        f.write(", ".join(( str(my_iter), np.array_str(x, max_line_width=np.inf, precision=None).replace('\n', ',').replace('[', '').replace(']', ''), np.array_str(flatValues, max_line_width=np.inf, precision=None).replace('\n', ',').replace('[', '').replace(']', ''), str(out) + "\n")))
    my_iter = my_iter + 1
    
    return out

# Perform the optimization
if __name__ == "__main__": # This structure supposedly sometimes helps Python MPI applications run correctly
    
    # Initialize some optimization-relevant variables
    my_iter = 0
    specArray = np.array(speciesParams)
    Nspecies = specArray.shape[1]
    woutFileName = 'wout_'+fullVMECInputFileName.split('.')[-1].strip()+'.nc'
    surf = SurfaceRZFourier.from_vmec_input(fullVMECInputFileName)

    # Organize variables in the optimization space
    if not unfixMajorRadius:
        surf.fix("rc(0,0)")
    
    if unfixPHIEDGE:
        x = np.zeros((len(surf.x) + 1))
        x[:len(surf.x)] = surf.x
        x[-1] = readPhiEdge('.', fileName=fullVMECInputFileName)
        DOF_NAMES = surf.dof_names + ['phiedge']
    else:
        x = surf.x
        DOF_NAMES = surf.dof_names

    # Save DOF names straight away in case the optimization does not complete
    np.save("dof_names.npy", DOF_NAMES)

    # Create the log file
    header = 'iter, DOFvals, Lvals, J'
    with open(logfileName,'w') as f:
        f.write(header + '\n')

    # Set up the bounds for the global optimizer
    boundsList = [(float(xval * (1 - autobound)), float(xval * (1 + autobound))) for xval in x]

    # Carry out the optimization
    args = [opt, boundsList]
    res = differential_evolution(*args, **DEkwargs)

    # Save some relevant outputs
    np.save("x.npy", res.x)
    np.save("fun.npy", res.fun)
    np.save("my_iter.npy", my_iter)

    # Print some relevant outputs in the log file
    print("success = " + str(res.success), flush=True)
    print("x =  " + str(res.x), flush=True)
    print("status = " + str(res.message), flush=True)
    print("f = " + str(res.fun), flush=True)
    print("OPTIMIZATION COMPLETE.", flush=True)
