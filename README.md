## Introduction

Code developed for obtaining the **advanced Master's degree in Nanotechnology**. The code is a finite element algorithm for linear elasticity with surface stress and surface stiffness theory, or non-linear saint Venant-Kirchhoff elasticity for large deformations. The doc folder contains a compilable documentation of the input parameters. This code is not standalone and requires an input in solver.cpp. [Thesis report here](https://hdl.handle.net/2268/293576).

## Installation

First, make sure that GCC and the external library Alglib is correctly installed in your system. Then move to the main folder and compile the project by providing the path of Alglib to the compiler. An example of command-line compilation is provided in compile.sh.
```sh
g++ -I path-to-Alglib/src -o NanoFEM -O3 -static -static-libgcc
-static-libstdc++ path-to-Alglib/src/*.cpp source/*.cpp
```

The code alone will not run, you must add to solver.cpp the required code that builds the data structure used as input for the solver, this latter contains the information about the mesh and boundary conditions as described in the documentation. The executable is generated in the main folder and can be executed by
```sh
./Mechanics
```

## Author

* Martin Lacroix
