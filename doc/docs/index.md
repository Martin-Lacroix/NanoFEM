# <img src="Cpp.svg" width="60"/> Data Structure
--------------------------------------------------

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

The vector `nXYZ` has the dimension **(n,3)** where **n** is the number of nodes in the mesh. Each `array3d` stores the 3 coordinates (x,y,z) of the node assigned with the corresponding index.

<br />

```cpp
    std::vector<ivector> eNode;         // Nodes of each element
```

The vector `eNode` has the dimension **(e,k)** where **e** is the number of element in the mesh and **k** is the number of nodes in each of these element, an element of order **n** contains **(n+1)^3** nodes. Each `ivector` contains the list of the nodes of the element assigned with the corresponding index. The nodes are stored in the order z->y->x increasing coordinate in the local space.

<br />

```cpp
    std::vector<ivector> eSurf;         // Free surfaces of each element
```

The vector `eSurf` has the dimension **(e,k)** where **e** is the number of element in the mesh and **k** is the number of free surface in each of these element, an element can have up to 6 free surfaces. Each `ivector` contains a list of the indices of the free surfaces of the element. The indices of the surfaces of an element are noted [-z,+z,-y,+y,-x,+x] = [0,1,2,3,4,5], where **-z** denotes the face perpendicular to the **z** axis and located at **z=-1** in the local space. This vector is empty is the element has no free surfaces.

<br />

```cpp
    std::vector<ivector> neuFace;       // Faces for Neumann boundary conditions
```

The vector `neuFace` has the dimension **(f,k)** where **f** is the number of faces where Neumann boundary conditions are applied and **k** is the number of nodes in each of these faces, a face of order **n** contains **(n+1)^2** nodes. Each `ivector` contains the list of the nodes of the face assigned with the corresponding index. The nodes are stored in the order y->x increasing coordinate in the local space.

<br />

```cpp
    std::vector<darray> neuVal;         // Forces for Neumann boundary conditions
```

The vector `neuVal` has the dimension **(f,3)** where **f** is the number of faces where Neumann boundary conditions are applied, the indices correspond to the face described in the `neuFace` vector. Each `darray` contains the 3 components (x,y,z) of the applied force on the corresponding face.

<br />

```cpp
    ivector dirNode[3];         // Nodes for Dirichlet boundary conditions
```

The vector `dirNode` has the dimension **(3,n)** where **n** is the number of nodes where Dirichlet boundary conditions are applied, the first dimension of length 3 corresponds to the direction in space (x,y,z) of the imposed displacement. Each `ivector` contain the list of node indices where the Dirichlet BC is imposed.

<br />

```cpp
    dvector dirVal[3];          // Displacment for Dirichlet boundary conditions
```

The vector `dirVal` has the dimension **(3,n)** where **n** is the number of nodes where Dirichlet boundary conditions are applied, the first dimension of length 3 corresponds to the direction in space (x,y,z) of the imposed displacement. Each `ivector` contain the list of values of the Dirichlet BC in the corresponding dimension.

<br />

# <img src="Cpp.svg" width="60"/> Polynomial Basis
-----------------------------------------------------

A polynomial basis can be constructed by with an exponent table and a coefficient matrix. The parameter `d` is the dimension, `m` the number of monomials and `p` is the number of polynomials.

```python
    poly = Polynomial(expo,coef,csr=0)              # Class of olynomial basis
```

| Input             | Type                  | Description                       |
|-------------------|-----------------------|-----------------------------------|
| *expo*            | *(d,m) array*         | *exponent table*                  |
| *coef*            | *(p,m) array*         | *coefficients matrix*             |
| *csr*             | *bool*                | *for sparse csr format*           |

<br />

The class of polynomial basis contains some in-built methods such as the evaluation of its polynomials at an array of `n` point:

```python
    V = poly.eval(point)                                # Computes the Vandermonde matrix
    poly.clean(index)                                   # Selects the relevant polynomials
    poly.trunc(order)                                   # Truncates the polynomial basis
```

| Input             | Type                  | Description                                   |
|-------------------|-----------------------|-----------------------------------------------|
| *point*           | *(n,d) array*         | *points to evaluate polynomials*              |
| *index*           | *(-) array*           | *indices of the polynomials to keep*          |
| *order*           | *int*                 | *order of truncation*                         |

<br />

An orthonormal polynomial basis with respect to a sample can be generated by Gram Schmidt process. The parameter `weight` must be provided if the points are quadrature nodes, otherwise a Monte Carlo integration is assumed.

```python
    poly = gschmidt(order,point,weight=0,trunc=1)               # Gram-Schmidt process
```

| Input             | Type                  | Description                                   |
|-------------------|-----------------------|-----------------------------------------------|
| *order*           | *int*                 | *maximum order of the polynomials*            |
| *point*           | *(n,d) array*         | *points for the quadrature*                   |
| *weight*          | *(n) array*           | *weights associated to the points*            |
| *trunc*           | *float*               | *hyperbolic truncation q-norm*                |

<br />

Similarly, orthogonal polynomials can be constructed using a three terms recurrence coefficients relation related to a well-knows distribution or joint distribution.

```python
    poly = polyrecur(order,dist,trunc=1)                    # Three terms recurrence relation
```

| Input             | Type                          | Description                                           |
|-------------------|-------------------------------|-------------------------------------------------------|
| *order*           | *int*                         | *maximum order of the polynomials*                    |
| *dist*            | *oject or (-) array*          | *joint distribution or list of distributions*         |
| *trunc*           | *float*                       | *hyperbolic truncation norm*                          |

<br />

# <img src="Cpp.svg" width="60"/> Quadrature Rules
-----------------------------------------------------

A quasi-Monte Carlo quadrature rule can be generated with a low-discrepancy sequence. The parameter `pdf` must have the same behaviour as the `pdf` method of a distribution class.

```python
    point,weight = qmcquad(nbrPts,dom,pdf=0,seq='halton')       # Quasi-Monte Carlo quadrature
```

| Input             | Type                  | Description                                   |
|-------------------|-----------------------|-----------------------------------------------|
| *nbrPts*          | *int*                 | *number pf points to generate*                |
| *dist*            | *(d,2) array*         | *multidimentional integration domain*         |
| *pdf*             | *callable*            | *callable probability density function*       |
| *seq*             | *string*              | *name of the low-discrepancy sequence*        |

<br />

A tensor product quadrature rule with respect to the probability density function of a well-known distribution or joint disctrubution can be generated by

```python
    point,weight = tensquad(order,dist)           # Tensor product quadrature
```

| Input             | Type                          | Description                                           |
|-------------------|-------------------------------|-------------------------------------------------------|
| *order*           | *int*                         | *order of the quadrature rule*                        |
| *dist*            | *oject or (-) array*          | *joint distribution or list of distributions*         |

<br />

Different sparse quadrature rules for a polynomial basis can be generated from a Monte Carlo integration set. The original number of points must be larger than the number of polynomials in the basis.

```python
    index,weight = simquad(point,poly)                 # Revised simplex algorithm
    index,weight = fekquad(point,poly)                 # Approximate Fekete points
    index,weight = nulquad(point,poly,weight)          # Positive quadrature with null space
    index,weight = newquad(point,poly,weight)          # Positive quadrature with Newton
```

| Input             | Type                  | Description                           |
|-------------------|-----------------------|---------------------------------------|
| *point*           | *(n,d) array*         | *original set of points*              |
| *poly*            | *object*              | *polynomial basis object*             |
| *weight*          | *(n) array*           | *previous weights of the points*      |

<br />

# <img src="Cpp.svg" width="60"/> Expansion Coefficients
-----------------------------------------------------------

The polynomial chaos coefficients can be computed by spectral projection of least squares regression. The parameter `n` is the number of points and `d` the dimension. The parameter `weight` must be provided if the points are quadrature nodes, otherwise a Monte Carlo integration is assumed.

```python
    coef = spectral(resp,poly,point,weight=0)              # Spectral projection
    coef = colloc(resp,poly,point,weight=0)                # Point collocation
```

| Input             | Type                  | Description                                   |
|-------------------|-----------------------|-----------------------------------------------|
| *resp*            | *(n,-) array*         | *response of thr model at the points*         |
| *poly*            | *object*              | *polynomial basis object*                     |
| *point*           | *(n,d) array*         | *quadrature points*                           |
| *weight*          | *(n) array*           | *weights of the points*                       |

<br />

The least angle regression algorithm selects the relevant polynomials in addition to compute their coefficients. The parameter `index` is the indices of the selected polynomials, if `it` is not provided, all the available polynomials are selected.

```python
    coef,index = lars(resp,poly,point,weight=0,it=np.inf)           # Least angle regression
    coef,index = lasso(resp,poly,point,weight=0,it=np.inf)          # Least shrinkage operator
```

| Input             | Type                  | Description                                   |
|-------------------|-----------------------|-----------------------------------------------|
| *resp*            | *(n,-) array*         | *response of thr model at the points*         |
| *poly*            | *object*              | *polynomial basis object*                     |
| *point*           | *(n,d) array*         | *original set of points*                      |
| *weight*          | *(n) array*           | *weights of the points*                       |
| *it*              | *int*                 | *maximum number of iterations*                |

<br />

# <img src="Cpp.svg" width="60"/> Surrogate Model
----------------------------------------------------

The polynomial chaos model can be constructed by invoking the constructor of the expansion class. The parameter `p` is the number of polynomials in the basis. The class provides different built-in methods:

```python
    model = Expansion(coef,poly)                # Class containing the surrogate model
    resp = model.eval(point)                    # Evaluates the surrogate model
    mean = model.mean                           # Returns the mean of the output
    var = model.var                             # Return the variance of the output
```

| Input             | Type                  | Description                                   |
|-------------------|-----------------------|-----------------------------------------------|
| *coef*            | *(p,-) array*         | *polynomial chaos coefficients*               |
| *poly*            | *object*              | *polynomial basis object*                     |
| *point*           | *(n,d) array*         | *points at which evaluate the model*          |

<br />

In addition, a simple polynomial mapping between two one-dimensional random variables can be computed by

```python
    mapping = transfo(invcdf,order,dist)            # Transforms a distribution
    y = mapping(x)                                  # Mapping from the x space to the y space
```

| Input             | Type                  | Description                                       |
|-------------------|-----------------------|---------------------------------------------------|
| *invcdf*          | *callable*            | *inverse cumulative distribution function*        |
| *order*           | *int*                 | *order of the expansion*                          |
| *dist*            | *object*              | *distribution object*                             |

<br />

# <img src="Cpp.svg" width="60"/> Other Functions
----------------------------------------------------

A class of principal components analysis whitening can be created for a linearly correlated sample. The parameter `n` is the number of points and `d` is the dimension.

```python
    mapping = Pca(sample)                   # Class of PCA whitening
```

| Input             | Type                  | Description                           |
|-------------------|-----------------------|---------------------------------------|
| *sample*          | *(n,d) array*         | *reference sample of points*          |

<br />

The class will act as a mapping function between the whitened and the original random vectors.

```python
    whitened = mapping.white(sample)            # Whitens a sample of points
    sample = mapping.corr(whitened)             # Recovers the original sample
```

| Input             | Type                  | Description                                               |
|-------------------|-----------------------|-----------------------------------------------------------|
| *sample*          | *(n,d) array*         | *sample from the same distributio as the reference*       |
| *whitened*        | *(n,d) array*         | *white noise sample of points*                            |

<br />

The Sobol sensitivity indices of a model can be directly obtained from the polynomial chaos coefficients by

```python
    sobol = anova(coef,poly)           # Computes the Sobol sensitivity indices
```

| Input             | Type                  | Description                               |
|-------------------|-----------------------|-------------------------------------------|
| *coef*            | *(p,-) array*         | *polynomial chaos coefficients*           |
| *poly*            | *object*              | *polynomial basis object*                 |

<br />

For dependent random variables, the analysis of covariance indices can be obtained from the polynomial chaos model with

```python
    index,ancova = ancova(model,point,weight=0)         # Computes the ancova indices
```

| Input             | Type                  | Description                           |
|-------------------|-----------------------|---------------------------------------|
| *model*           | *object*              | *expansion object*                    |
| *point*           | *(n,d) array*         | *quadrature points*                   |
| *weight*          | *(n) array*           | *weights of the points*               |

<br />

Finally, any objects generated by Chaoslib can be saved in a pickle file:

```python
    save(item,name)                # Saves an object in a file.pickle
```

| Input             | Type                  | Description                               |
|-------------------|-----------------------|-------------------------------------------|
| *item*            | *object*              | *any object to be saved*                  |
| *name*            | *string*              | *desired name for the object*             |
