#!/usr/bin/env python
import subprocess
import numpy as np 

def changevar(dirname, group, var, value, inputname = "input.namelist"):

    filename = dirname + "/" + inputname
    # Warning: this command will fail silently if the pattern is not found. Sorry about that.
    # Warning: case insensitive
    if type(value) == bool:
        if value == True:
            value = ".true."
        else:
            value = ".false."
    elif type(value) == str:
        #strings must be enclosed in "" in namelists
        #may be wise to see if the string contains citation marks...
        if (value.find("'") != -1) or (value.find('"') != -1):
            print("Warning! String to changevar contains a ' or \" character.", flush=True)
        value = '"' + value + '"'
        # escape slashes
        value = value.replace("/","\/")
    elif (type(value) == list) or (type(value) == np.ndarray):
        # arrays are space seperated
        delimiter=' '
        value_temp = '' 
        for val in value:
            value_temp =  value_temp + str(val) + delimiter
        value = value_temp.rsplit(delimiter,1)[0]
    else:
        value = str(value)
    subprocess.call(r"sed -i -e '/\&"+group+r"/I,/\&/{ s/^\(\s*"+var+r"\s*=\s*\).*/\1"+value+r"/I } ' " + filename, shell=True)



if __name__ == "__main__":
    import sys
    import os
    if len(sys.argv) > 1:
        dirname = argv[1]
    else:
        dirname = '.'
    changevar(dirname,"speciesParameters", "Zs", [1.0, 666.0])
    changevar(dirname, "geometryParameters", "equilibriumFile", os.path.abspath(os.path.join(dirname, 'wout_vmec.nc')))
