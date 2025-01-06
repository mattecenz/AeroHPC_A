# AeroHPC_A

## Overview

\< Available soon \>

This is part of a project-work for High Performance Scientific Computing In Aerospace course [@Polimi](https://www.polimi.it/)

### Authors
Project developed by:
- \< Available soon \>

### Problem description, Implementation decisions and Result analysis

Can be found into report - \< Available soon \>

### File Structure

* \< Available soon \>

### Build dependencies

In order to build all targets the following software libraries must be installed:

* `OpenMPI` - [Reference site](https://www.open-mpi.org/) - tested version: 5.0.x
* `FFTW` - [Reference site](https://fftw.org/) - tested version: 3.3.x

The program makes also use of a modified version of the following code (already provided in the repo):
* `C2Decomp` - [Reference site](https://github.com/emathew1/2Decomp_C)

Since installing methods can variate between linux distribution we do not provide any installation tutorial for the above libraries.

### How to build

In order to build the executable, from the root folder run the following commands:

```bash
$ mkdir build
$ cd build
$ cmake .. _FLAGS_ _TEST_FLAGS_
$ make -j %n
```

`_FLAGS_` are optional, they can be:
* `-D FFTW_LIBRARY_PATH='/path/to/libfftw*.so'` - specify fftw library path
* `-D FFTW_INCLUDE_PATH='/path/to/fftw_headers` - specify fftw headers path
* `-D USE_FLOAT=(0 or 1)` - specify if you want to use single or double precision (default double) \[make sure you have the necessary fftw library object\]
* `-D MAKE_TEST=(0 or 1)` - specify whether to build the test main file or the release one

`_TEST_FLAGS_` are optional, they can be:
* `-D DISABLE_PRESSURE=1'` - disable pressure terms
* `-D DEBUG_PRINT_BUFFERS=1` - enable print of numerical method buffers (keep attention, may generate lot of files)

output executables will be:
* `AeroHPC_A` - Main executable

### How to run solver

```bash
$ mpirun -n N_PROCESSORS ./AeroHPC_A /path/to/testCaseInfo.txt
```
where:
* `N_PROCESSORS` - is the number of processors you want to use (make sure they are compatible with TestCaseInfo)
* `path/to/testCaseInfo.txt` - path to test case input data file

### Supported format of test case input data file

```txt
nTimeSteps
deltaT
Nx Ny Nz
Npy Npz
TestCaseID
```
where:
* `nTimeSteps` - Is an integer, defines the number of solver timesteps
* `deltaT` - Is a float, defines the size of time variation between timesteps 
* `Nx Ny Nz` - Are integers, define number of discretization points in X Y and Z directions
* `Npy Npz` - Are integers, define number of domain partition on Y and Z direction (one partition will be assigned to each processor)
* `TestCaseID` - Is the identification number of the test case, based on this number the solver will apply all remaining problem parameters (as Boundary condition, domain size etc.)
