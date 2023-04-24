# SFINCSopt
This repo contains simple Python scripts for performing neoclassical optimizations with SFINCS. The code was originally written for a slightly different purpose by [Stefan Buller](https://github.com/daringli). Minor modifications were made by [Brandon Lee](https://github.com/leebr48/). In its current form, the code allows users to optimize the boundary shape of a VMEC equilibrium to achieve desired values of the thermal transport coefficients (typically written as $L_{ab}^{s}$, with $a,b=1,2,3$ and $s$ a species). This can be done for an arbitrary number of magnetic surfaces and species.

## How to Use
The code must be given input files for [SFINCS](https://github.com/landreman/sfincs) and [VMEC](https://github.com/PrincetonUniversity/STELLOPT), as well as a "stripped" VMEC input file with the boundary coefficients removed. These files primarily serve to set resolution parameters and so forth for the SFINCS and VMEC runs. They should live in the same directory as the Python files. Examples have been included in this repo. In `main.py`, the user can change a variety of options related to the numerical settings of each calculation, plasma profiles, optimization settings, and so forth. It is typical to prepare for an optimization by creating a `prep/` subdirectory with the same structure as a run subdirectory (`prep/fluxSurfaceNumber/speciesNumber/runFiles`) and running SFINCS for each combination of flux surface and species. This information can be used to set the $L$ target values appropriately given the initial condition. The `initETargets.py` script, or a slightly modified version of it, may help with this. The `readInitJ.py` script can also print the initial value of the objective function. Note that the resource allocation parameters in `main.py` should be compatible with those in `job`. Keep in mind that the computational requirements (especially related to memory) for tasks inside the optimization loop are typically greater than those for the same tasks run outside the optimization loop. To start an optimization, one can run `sbatch job`. Optimization subdirectories (for each iteration) will be created in the same repository where the python scripts and supporting files live.

After an optimization is complete, the user can identify the iteration with the lowest value of the objective function with `readLog.py`. If the user wishes to delete the files create by an optimization, they can call `cleandir.sh`.
