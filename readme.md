# HARP - Hierarchical Atomic Resolution Perception
HARP is an advanced, model-based, physics-informed machine learning algorithm capable of perceiving the local resolution of a biomolecule in its imaged density. Using the [Bayesian inference-based shape calculation framework](https://bayes-shape-calc.github.io/) it can be used to validate the trustworthiness of (*e.g.*, cryoEM-derived) structural models.

## Quick start
1. Install HARP

2. Activate the HARP environment (`conda activate harp` or `source harpenv/bin/activate`; see above)

3. use the CLI (using PDB ID 6J6J as an example; see tutorial on documentation website):
``` bash
harpcalc -id 6j6j
```

``` 
Loading 6j6j
N_adfs = 10, N_blobs = 20
Using default weights: [0.05 1.   1.   1.   2.   2.  ]
Chain: A, 119 residues, <P_res> = 0.7948, t = 1.40 sec, <t/res> = 11.73 msec
Chain: B, 119 residues, <P_res> = 0.8067, t = 1.34 sec, <t/res> = 11.28 msec
Chain: C, 119 residues, <P_res> = 0.8062, t = 1.39 sec, <t/res> = 11.65 msec
Chain: D, 119 residues, <P_res> = 0.8045, t = 1.44 sec, <t/res> = 12.08 msec
<P_res>: 0.8030
Finite: 476 / 476
```

## Install
HARP comes as a Python module called `harp` that has a command line interface (CLI) program called `harpcalc`. Installing the `harp` module into a separate Python environment is recommended practice, because it allows you to maintain required library dependencies. You can use `conda` or `venv` to do this. 

### conda
`conda` is a Python environment manager. You'll have to [download and install conda](https://conda.io/docs/user-guide/install/) or the smaller [miniconda](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html). With that installed, you need to get the HARP code, make an environment for HARP, and then install HARP. Here is an example of how to do that from a terminal window:

``` bash
git clone https://github.com/bayes-shape-calc/HARP.git
cd HARP
conda env create
conda activate harp
pip install ./
```

### venv + pip
If you don't want to use `conda`, you can use `venv` to make an environment and `pip` to handle installing all the libraries. Here is an exmaple of how to do that from a terminal window:

``` bash
python -m venv harpenv
source harpenv/bin/activate
git clone https://github.com/bayes-shape-calc/HARP.git
cd HARP
pip install ./
```

### Dependencies
These are not the absolute minimum required versions. They should be automatically installed by the `pip` command above (step 4).

| Library    | Version| Required? | Use                                    |
| ---------  | ------ |-----------|--------------------------------------- |
| python     | 3.7    |  Yes      | Programming                            |
| mrcfile    | 1.3    |  Yes      | For loading density maps               |
| numba      | 0.55   |  Yes      | For fast model building                |
| numpy      | 1.22   |  Yes      | Math                                   |
| gemmi      | 0.5.5  |  No       | (optional) For X-ray SF loading        |


## [Optional] Optimization: Faster model building with a C library
The functions to create Gaussian blob density models at certain points in space are in `HARP/harp/models.py`. The default is a python version that is JIT compiled with Numba for fast calculations. There is also a version in C that must first be compiled into a shared library. It is than faster the Numba version.

### Requirements
You'll need `gcc` and `make` for this to work. If you are running Ubuntu, you can try the following to get `gcc` and `make`:
``` bash
sudo apt update
sudo apt install build-essential
```

If you are on a mac, you can try first installing [homebrew](https://brew.sh), and then running:
``` bash
brew install gcc
brew install make
```

### Compiling
A makefile exists that you can use to easily compile this library if you are on a linux/mac computer:
``` bash
cd HARP
make
```

### Usage
In the harp module, the C functions are called through `ctypes`. You can switch between the Numba and C versions using:
``` python
harp.models.use_c()
harp.models.use_python()
```

Alternatively, if you are using `harpcalc`, you can use the argument `--use_c` or `--use_python`
