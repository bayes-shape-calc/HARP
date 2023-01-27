# HARP - Hierarchical Atomic Resolution Perception
HARP is an advanced, model-based, physics-informed machine learning algorithm capable of perceiving the local resolution of a biomolecule in its imaged density. Using the [Bayesian inference-based shape calculation framework](https://bayes-shape-calc.github.io/) it can be used to validate the trustworthiness of (*e.g.*, cryoEM-derived) structural models.


## Quick start
1. Download this repository

2. Navigate to the HARP folder
```bash
cd <path to the location of the HARP folder>
cd HARP
```

3. [*Optional*] To compile the (fast) C library, run the command:
```bash
make
```

4. To install the harp module, run the command:
``` bash
pip install ./
```

5. Option 1: Write python scripts:
``` python
import harp
```

6. Option 2: use the installed command line interface (using PDB ID 6J6J as an example):
``` bash
harpcalc 6j6j
```


### Isolated install
Follow these to use instructions to use and install in an isolated virtual environment

``` bash
git clone https://github.com/bayes-shape-calc/HARP.git
python -m venv harpenv
source harpenv/bin/activate
cd HARP
make
pip install ./
cd ..
harpcalc 6j6j --use_c
deactivate
```

and you should see output that looks like:
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

## Library versions
Note: these are not the minimum required versions. They should be automatically installed by the `pip` command above (step 4).

| Library    | Version| Required? | Use                                    |
| ---------  | ------ |-----------|--------------------------------------- |
| gemmi      | 0.5.5  |  No       | (optional) For X-ray SF loading        |
| mrcfile    | 1.3    |  Yes      | For loading density maps               |
| numba      | 0.55   |  Yes      | For fast model building                |
| numpy      | 1.22   |  Yes      | Math                                   |
| python     | 3.7    |  Yes      | Programming                            |



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
