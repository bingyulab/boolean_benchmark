
## Installation

CaSQ is available on PyPI and can be installed using pip:

```bash
pip install -r requirements.txt
```

### R Packages Installation
```bash
Rscript Install/install.R
```

Note, in CellNOpt, we need download CPLEX solver from IBM. You can download it from [IBM CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) and install it in the `~/CPLEX_Studio2211` directory. Make sure to set the environment variable `CPLEX_PATH` to point to the CPLEX installation directory.

## Matlab Installation
```bash
https://github.com/sysbiolux/optPBN.git
https://github.com/gingproc-IIM-CSIC/MEIGO64.git
```

```bash
unzip optPBN_stand_alone_v2.2.3
cd optPBN/optPBN_stand_alone_v2.2.3/
matlab -nodisplay -nosplash -r "run('$HOME/boolean/Install/install_optpbn.m'); exit"
```

```bash
cd ~/MEIGO64/MEIGO
matlab -nodisplay -nosplash -r "\
  addpath('$HOME/MEIGO64/MEIGO');               % add MEIGO64 folder
  addpath(genpath('$HOME/MEIGO64/MEIGO'));
  install_MEIGO;                   % run the installer script
  savepath('~/matlab_startup/pathdef.m');                        % save the updated pathdef
  exit" 
```

```bash
mkdir -p ~/matlab_startup

cd ~/matlab_startup && matlab -nodisplay -nosplash
run('~/matlab_startup/startup.m')
```

### Installation for Genetic Algorithms

#### Linux
```bash
# Install dependencies
apt-get install -y libtinfo5
apt-get install graphviz libgraphviz-dev pkg-config

# Install Python packages
pip install git+https://github.com/DEAP/deap@master
pip install git+https://github.com/hklarner/pyboolnet@master
pip install pydot igraph cairosvg pygraphviz
```

#### macOS
```bash
# Install dependencies
brew install libtinfo5
brew install graphviz libgraphviz-dev pkg-config

# Install Python packages
pip install git+https://github.com/DEAP/deap@master
pip install git+https://github.com/hklarner/pyboolnet@master
pip install pydot igraph cairosvg pygraphviz
```


#### Remove Large Files from Git History

1. Download ![BFG Repo-Cleaner](https://rtyley.github.io/bfg-repo-cleaner/)
2. From your repo root, run:
```bash
java -jar bfg.jar --delete-files cplex_studio2211.linux_x86_64.bin
git reflog expire --expire=now --all && git gc --prune=now --aggressive
git push origin main --force
```