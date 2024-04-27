## Introduction

The code is a finite element algorithm for linear elasticity with surface stress and surface stiffness theory, or non-linear saint Venant-Kirchhoff solid in large deformations. The doc folder contains a markdown documentation of the input parameters. This code is not yet standalone and requires the user to provide a structure with the mesh connectivity and mechanical properties.

## Installation

First, make sure that GCC and the external library Eigen is correctly installed in your system. Then move to the source folder and compile the project with CMake. An example of compiler command is given in compile.sh.
```sh
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
```

The code alone will not run, you must add to solver.cpp the required code that builds the data structure used as input for the solver, this latter contains the information about the mesh and boundary conditions as described in the documentation. The executable is generated in the main folder and can be executed by
```sh
./NanoFEM
```
