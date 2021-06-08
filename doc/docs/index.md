# <img src="Cpp.svg" width="60"/> Data Structure
--------------------------------------------------

This structure is located in `Mesh` and contains all the information required for building the mesh and computing the finite element solution. This is the most important input of the algorithm.

<br />

```cpp
    int step;                   // Number of load steps
    int order;                  // Order of the quadrature and the elements
    double tol;                 // Tolerance for Newton-Raphson
```

The `step` must be a strictly positive integer, this gives the nulmber of intermediate equilibrium configurations computed in the SVK model. The `order` must be a strictly positive integer, but is generally fixed at 1 or 2. The `tol` is the value of the relative correction of the displacement under which the convergence of the Newton-Raphson algorithm ends, generally 1e-6.

<br />

```cpp
    std::vector<array3d> nXYZ;          // Coordinates of the nodes
```

The vector `nXYZ` has the dimension **(n,3)** where **n** is the number of nodes in the mesh. Each `array3d` stores the 3 coordinates **(x,y,z)** of the node assigned with the corresponding index.

<br />

```cpp
    std::vector<ivector> eNode;         // Nodes of each element
```

The vector `eNode` has the dimension **(e,k)** where **e** is the number of element in the mesh and **k** is the number of nodes in each of these element, an element of order **n** contains **(n+1)^3** nodes. Each `ivector` contains the list of the nodes of the element assigned with the corresponding index. The nodes are stored in the order **z -> y -> x** increasing coordinate in the local space.

<br />

```cpp
    std::vector<ivector> eSurf;         // Free surfaces of each element
```

The vector `eSurf` has the dimension **(e,k)** where **e** is the number of element in the mesh and **k** is the number of free surface in each of these element, an element can have up to 6 free surfaces. Each `ivector` contains a list of the indices of the free surfaces of the element. The indices of the surfaces of an element are noted **(-z,+z,-y,+y,-x,+x) = (0,1,2,3,4,5)** with **-z** denoting the face perpendicular to the **z** axis and located at **z=-1** in the local space. This vector is empty is the element has no free surfaces.

<br />

```cpp
    std::vector<ivector> neuFace;       // Faces for Neumann boundary conditions
```

The vector `neuFace` has the dimension **(f,k)** where **f** is the number of faces where Neumann boundary conditions are applied and **k** is the number of nodes in each of these faces, a face of order **n** contains **(n+1)^2** nodes. Each `ivector` contains the list of the nodes of the face assigned with the corresponding index. The nodes are stored in the order y->x increasing coordinate in the local space.

<br />

```cpp
    std::vector<darray> neuVal;         // Forces for Neumann boundary conditions
```

The vector `neuVal` has the dimension **(f,3)** where **f** is the number of faces where Neumann boundary conditions are applied, the indices correspond to the face described in the `neuFace` vector. Each `darray` contains the 3 components **(x,y,z)** of the applied force on the corresponding face.

<br />

```cpp
    ivector dirNode[3];         // Nodes for Dirichlet boundary conditions
```

The vector `dirNode` has the dimension **(3,n)** where **n** is the number of nodes where Dirichlet boundary conditions are applied, the first dimension of length 3 corresponds to the direction in space **(x,y,z)** of the imposed displacement. Each `ivector` contains the list of node indices where the Dirichlet BC is imposed.

<br />

```cpp
    dvector dirVal[3];          // Displacment for Dirichlet boundary conditions
```

The vector `dirVal` has the dimension **(3,n)** where **n** is the number of nodes where Dirichlet boundary conditions are applied, the first dimension of length 3 corresponds to the direction in space **(0,1,2) = (x,y,z)** of the imposed displacement. Each `ivector` contains the list of values of the Dirichlet BC in the corresponding dimension.

<br />

```cpp
    std::vector<ivector> coupNode[3];       // List of coupled pack of nodes
```

The vector `coupNode` has the dimension **(3,k,n)** where **k** is the number of packs of coupled nodes, **n** is the number of nodes in a pack. The first dimension of length 3 corresponds to the direction in space **(0,1,2) = (x,y,z)** of the coupled displacement. Each `ivector` contains the list of coupled node indices. All coupled nodes of a same pack will have the same displacement in the corresponding dimension.

<br />

```cpp
    std::vector<std::pair<int,int>> deltaNode[3];       // List of delta pair of nodes
```

The vector `deltaNode` has the dimension **(3,k,2)** where **k** is the number of change of variables. The first dimension of length 3 corresponds to the direction in space **(0,1,2) = (x,y,z)** of the change of variable. Each `pair` contains the two nodes indices for the change of variable **(u1,u2) -> (u1-u2,u2)** such that the first variable (nodal value) gives now the relative displacement between the two nodes in the pair in the corresponding dimension, while the second variable remains unchanged.

<br />

```cpp
    std::vector<array3d> LmR;           // Bulk element parameters
    std::vector<array3d> LmS;           // Surface element parameters
```

Both `LmR` and `LmS` have the dimension **(n,3)** where **n** is the number of elements in the mesh. Each `array3d` contains the first and second Lamé parameters for the bulk or the surface material. The last entry if the density for `LmR` and the surface tension for the `LmS` vector.

<br />

# <img src="Cpp.svg" width="60"/> Shape Structure
-----------------------------------------------------

This structure is computed automatically by the algorithm according to the patrameters located in the `Data` structure. It contains the informations about the Legendre shape functions and Legendre quadrature rule for numerical integration. This structure is the same for all elements and is strus stored only once.

<br />

```cpp
    int gLen;               // Number of integration points
```

The value of `gLen` is positive and is equal to **(n+1)^d** where **n** is the order of the quadrature rule and **d** is the dimension. This integer must be strictly positive.

<br />

```cpp
    matrix N;               // Nodal shape function
```

The matrix `N` has the size **(n,g)** where **n** is the number of shape functions (so the number of nodes in an element) and **g** is the number of integration points. This matrix contains the value of each shape function evaluated at each integration points.

<br />

```cpp
    dvector weight;         // Nodal shape function
```

The vector `weight` has the size **g** which is the number of integration points. This vector contains the weights corresponding to the integration points located in the vector `gRST` with corresponding indices.

<br />

```cpp
    std::vector<matrix> dN;         // Nodal shape function
```

The vector `dN` has the size **(d,n,g)** where **d** is the dimension (1, 2 or 3), **n** is the number of shape functions of the elements and **g** is the number of integration points. It contains **d** matrices, each `matrix` contains the value of derivative of the shape functions wity respect to the dimension **d** evaluated at the corresponding integration point in local coordinates.

<br />

```cpp
    std::vector<dvector> gRST;          // Nodal shape function
```

The vector `gRST` has the size **(g,d)** where **g** is the number of integration points and **d** is the dimension. Each `dvector` contains the coordinates of the integration points **(r,s)** or **(r,s,t)** in the local space.

<br />

# <img src="Cpp.svg" width="60"/> Quad Structure
-----------------------------------------------------

This structure is computed automatically by the algorithm according to the patrameters located in the `Data` structure. It is a temporary structure used to construct `Shape` and contains only ine informations about the quadrature rule.

<br />

```cpp
    int gLen;               // Number of integration points
```

The value of `gLen` is positive and is equal to **(n+1)^d** where **n** is the order of the quadrature rule and **d** is the dimension. This integer must be strictly positive.

<br />

```cpp
    dvector weight;         // Nodal shape function
```

The vector `weight` has the size **g** which is the number of integration points. This vector contains the weights corresponding to the integration points located in the vector `gRST` with corresponding indices.

<br />

```cpp
    std::vector<dvector> gRST;          // Nodal shape function
```

The vector `gRST` has the size **(g,d)** where **g** is the number of integration points and **d** is the dimension. Each `dvector` contains the coordinates of the integration points **(r,s)** or **(r,s,t)** in the local space.

<br />

# <img src="Cpp.svg" width="60"/> Mesh Class
-----------------------------------------------------

This class if the main class of the algoirithm and contains all the methods, objects and structures required for the computation of the global matrices and global vector used in the linear system of equation representing the problem.

<br />

```cpp
    dataStruct data;                // main information about the mesh
    shapeStruct shape3D;            // Shape functions for 3D hexahedrons
    shapeStruct shape2D;            // Shape functions for 2D quadrangles
    shapeStruct shapeS[6];          // Shape functions 3D -> 2D for element faces
```

The structure `data` is the one presented in the first section. The structures `shape3D` and `shape2D` are the ones presented in the second section, they are computed respectively for the 3D hexahedrons used for elemental matruces, and the 2D quadrangles ujused in Neumann boundary conditions. The array `shape3D` contains the 2D shape functions of the 6 element faces defined in the 3D local space of the element. The indices of the faces of an element are noted **(-z,+z,-y,+y,-x,+x) = (0,1,2,3,4,5)** with **-z** denoting the face perpendicular to the **z** axis and located at **z=-1** in the local space. 

<br />

```cpp
    int nLen;           // Number of nodes
    int eLen;           // Number of elements
    int fLen;           // Number of Neumann BC
```

The integer `nLen` is strictly positive and denotes the number of nodes in the mesh. The integer `eLen` is strictly positive and denotes the number of Lagrange hexahedrons in the mesh. The integer `fLen` is positive and denotes the number of faces on which Neumann boundary conditions are applied.

<br />

# <img src="Cpp.svg" width="60"/> Element Class
-----------------------------------------------------

This class represent a single finite element composing the mesh. They are stored in a vector into the `Mesh` class, but are regularily cleaned and rebuilt in order to save memory.

<br />