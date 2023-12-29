## Introduction

Code developed for obtaining the **advanced Master's degree in Nanotechnology**. The code is a finite element algorithm for linear elasticity with surface stress and surface stiffness theory, or non-linear saint Venant-Kirchhoff elasticity for large deformations. The doc folder contains a compilable documentation of the input parameters. This code is not standalone and requires an input in solver. [Thesis report here](https://hdl.handle.net/2268/293576).

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

## Author

* Martin Lacroix
