#!/usr/bin/env python

import subprocess
import os
from shutil import copy

from slurm_utils import ErrOut, getLatestJobIDInDir

from modify_sfincs_input import changevar

jobID = os.environ['SLURM_JOB_ID']
outname = "out." + jobID
errname = "err." + jobID

def create_sfincs_dirs(vmec_dirnames, surfs, speciesArray, nHats, dnHatdrHats, THats, dTHatdrHats, Er, sfincs_input = "input.namelist", vmec_output = "wout_vmec.nc"):
    nsurfs = len(surfs)
    nSpecies = speciesArray.shape[1]

    # speciesArray should be a numpy array with two rows (Zs and mHats), and one column for each species
    
    # create dirnames as a flat list
    sfincs_dirnames = []
    for d in vmec_dirnames:
        # we assume the latest jobID is correct
        status = ErrOut(d).status
        print('For directory {}, the VMEC status is: {}'.format(d, status), flush=True)
        if status == "STEPLIMIT":
            # instantly abort job
            os.write(2, b"srun: error: Step limit reached for this job\n")
            exit(1)

        elif status == "VMECDONE":
            for i in range(nsurfs):
                dirname0 = d + "/" + str(i).zfill(3)
                if not os.path.exists(dirname0):
                    os.mkdir(dirname0)
                for j in range(nSpecies):
                    dirname = os.path.join(dirname0, str(j).zfill(2))
                    if not os.path.exists(dirname):
                        os.mkdir(dirname)
                        sfincs_dirnames.append(dirname)
                        copy(d + "/" + sfincs_input, dirname)
                        copy(d + "/" + vmec_output, dirname)
                        changevar(dirname, "geometryParameters","rN_wish", surfs[i])
                        changevar(dirname, "geometryParameters", "equilibriumFile", os.path.abspath(os.path.join(dirname, 'wout_vmec.nc')))
                        changevar(dirname, "speciesParameters", "Zs", speciesArray[0,j])
                        changevar(dirname, "speciesParameters", "mHats", speciesArray[1,j])
                        changevar(dirname, "speciesParameters", "nHats", nHats[i][j])
                        changevar(dirname, "speciesParameters", "dnHatdrHats", dnHatdrHats[i][j])
                        changevar(dirname, "speciesParameters", "THats", THats[i][j])
                        changevar(dirname, "speciesParameters", "dTHatdrHats", dTHatdrHats[i][j])
                        changevar(dirname, "physicsParameters", "Er", Er[i])
                    
    return sfincs_dirnames
    
    
def check_VMEC_to_run(dirnames):
    to_run_dirnames = []
    for d in dirnames:
        # assume latest jobID is correct
        status = ErrOut(d).status
        if (status == "NOOUT") or (status == "TIME") or (status == "CANCELLED") or (status == "STEPLIMIT"):
            to_run_dirnames.append(d)
        else:
            print("Skipping running VMEC in directory '" + d + "'.", flush=True)
    return to_run_dirnames


def run_task_parallel(task, nconcur, dirnames, nprocesses = 40, memPerCpu='2G', nNodes=None, timeout=None):
    if nNodes is None:
        cmd = ['srun',"--exclusive","-n", str(nprocesses), "--mem-per-cpu={}".format(memPerCpu), task]
    else:
        cmd = ['srun',"--exclusive","-n", str(nprocesses), "-N", str(nNodes), "--mem-per-cpu={}".format(memPerCpu), task]
    nremaining = len(dirnames)
    i = 0
    while nremaining > 0:
        processes =[]
        if nconcur > nremaining:
            n = nremaining
        else:
            n = nconcur
        for j in range(n):
            dirname = dirnames[i]
            #dir should already be created and ready
            if not os.path.exists(dirname):
                raise ValueError("Directory '" + dirname + "' does not exist!")
            
            with open(dirname + "/" + outname, 'a') as stdout, open(dirname + "/" + errname, 'a') as stderr:
            # Popen() runs in parallel, unlike run()
                # need to add memory limit to actually run in parallel , "--mem-per-cpu=200MB" 
                process = subprocess.Popen(cmd, cwd = dirname, text=True, stdout = stdout, stderr = stderr)
                processes.append(process)
                # it may be a problem that stderr, stdout are close()'d here
            i = i + 1
            
        try:
            exit_codes = [p.wait(timeout=timeout) for p in processes]
        except subprocess.TimeoutExpired:
            pass
        
        nremaining = nremaining - n

def run_compute(dirnames, surfs, speciesArray, nHats, dnHatdrHats, THats, dTHatdrHats, Er, nconcur, nprocesses = 40, memPerCpu='2G', nNodes=None, sfincs_input="input.namelist", vmec_output="wout_vmec.nc", VMECtimeout=None):
    run_VMEC_dirs = check_VMEC_to_run(dirnames)

    run_task_parallel("../taskVMEC.sh", nconcur, run_VMEC_dirs, nprocesses = nprocesses, memPerCpu=memPerCpu, nNodes=nNodes, timeout=VMECtimeout)

    # this also checks if sfincs needs to be run
    sfincs_dirnames = create_sfincs_dirs(dirnames, surfs, speciesArray, nHats, dnHatdrHats, THats, dTHatdrHats, Er, sfincs_input=sfincs_input, vmec_output=vmec_output)
    run_task_parallel("../../../taskSFINCS.sh", nconcur, sfincs_dirnames, nprocesses = nprocesses, memPerCpu=memPerCpu, nNodes=nNodes)
    return sfincs_dirnames
