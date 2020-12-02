#include "math.h"
#ifndef ELEM_H
#define ELEM_H

struct shapeStruct{

    int nLen;
    int gLen;

    // Shape functions aty Gauss points
    // Local derivatives of the shape functions at Gauss nodes

    matrix N;
    matrix drN;
    matrix dsN;
    matrix dtN;
};

class Elem{

    public:

    Elem(std::vector<dvector> nXYZ,shapeStruct shape);
    matrix selfS(quadStruct quad,dvector xyz);
    matrix selfK(quadStruct quad,matrix D);

    // Parameters specific to each element

    std::vector<dvector> gXYZ;
    darray detJ;
    matrix dxN;
    matrix dyN;
    matrix dzN;
    matrix N;
    int nLen;
};

class Face{

    public:

    Face(std::vector<dvector> nXYZ,shapeStruct shape);
    matrix selfN(shapeStruct shape,quadStruct quad);

    // Parameters specific to each face

    darray detJ;
    matrix dxN;
    matrix dyN;
    int nLen;
};

#endif