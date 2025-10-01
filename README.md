# Boolean

This repository contains tools for analyzing and optimizing Boolean networks using various methods, including CellNOpt [1], MEIGO [2], and caspo [3].

[1] Terfve, C., Cokelaer, T., Henriques, D. et al. CellNOptR: a flexible toolkit to train protein signaling networks to data using multiple logic formalisms. BMC Syst Biol 6, 133 (2012). https://doi.org/10.1186/1752-0509-6-13

[2] Egea, J.A., Henriques, D., Cokelaer, T. et al. MEIGO: an open-source software suite based on metaheuristics for global optimization in systems biology and bioinformatics. BMC Bioinformatics 15, 136 (2014). https://doi.org/10.1186/1471-2105-15-136

[3] Santiago Videla, Julio Saez-Rodriguez, Carito Guziolowski, Anne Siegel, caspo: a toolbox for automated reasoning on the response of logical signaling networks families, Bioinformatics, Volume 33, Issue 6, March 2017, Pages 947–950, https://doi-org.proxy.bnl.lu/10.1093/bioinformatics/btw738

## Load HPC Modules
To run the code on a High-Performance Computing (HPC) cluster, you need to load the necessary modules. Here are the commands to do so:

```bash
salloc -N 1 -n 1 --exclusive --time=4:00:00
salloc -N 1 -n 1 --exclusive -C skylake
module load lang/R/4.4.1-gfbf-2023b
module load lang/Python/3.11.5-GCCcore-13.2.0
module load lang/SciPy-bundle/2023.11-gfbf-2023b
source ${HOME}/.virtualenvs/boolean/bin/activate
module load lang/Python/3.11.5-GCCcore-13.2.0
module load math/MATLAB/2024a-r6
```
Since the architectures of the clusters are different, and the partitions in each cluster have different architectures, you may need to adjust the module versions accordingly.

```bash
module load lang/SciPy-bundle/2023.11-gfbf-2023b
# python -m pip install --user virtualenv # For local installation
virtualenv --system-site-packages ${HOME}/.virtualenvs/boolean
```
The --system-site-packages option inherits all packages already installed in your system site-packages directory. To use the environment, just load the module and the environment:
```bash
module load lang/SciPy-bundle/2023.11-gfbf-2023b
source ${HOME}/.virtualenvs/boolean/bin/activate
deactivate
```

## Installation
Please follow the instructions in the [installation.md](doc/installation.md) file for installation.
```bash
make install
```
Note: the version of python should be 3.11 or higher and should install `clasp` and `gringo`.

## Usage
```bash
Rscript Step_01_Synthesis_data.R # For TCell network generating experiment data

# Run each dataset
make tcell 
make toy
make dream 
```
