#include "math.h"
#ifndef ELEM_H
#define ELEM_H

// -------------------------------------------------------------------|
// Structure storing the shape functions in the element local space   |
// -------------------------------------------------------------------|

struct shapeStruct{

    // Local derivatives of the shape functions at Gauss nodes

    matrix N;
    matrix drN;
    matrix dsN;
    matrix dtN;
};

// --------------------------------------------------|
// Class of 8-node hexahedron linear element in 3D   |
// --------------------------------------------------|

class Elem{

    public:

    // Functions available in the elem.cpp file

    Elem(std::vector<dvector> nXYZ,shapeStruct shape);
    darray stress(quadStruct quad,matrix D,darray u);
    matrix selfS(quadStruct quad,dvector xyz);
    matrix selfK(quadStruct quad,matrix D);

    // Parameters specific to each element

    std::vector<dvector> gXYZ;
    darray detJ;
    matrix dxN;
    matrix dyN;
    matrix dzN;
    int nLen;
};

// --------------------------------------------------|
// Class of 4-node quadrangle linear element in 2D   |
// --------------------------------------------------|

class Face{

    public:

    // Functions available in the elem.cpp file

    Face(std::vector<dvector> nXYZ,shapeStruct shape);
    matrix selfM(shapeStruct shape,quadStruct quad);

    // Parameters specific to each face

    darray dJ2D;
    int nLen;
};

#endif