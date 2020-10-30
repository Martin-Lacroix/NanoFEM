#include "func.h"
#ifndef ELEM_H
#define ELEM_H

struct shapeStruct{

    int nLen;
    int gLen;

    matrix N;
    matrix drN;
    matrix dsN;
    matrix dtN;
};

class Elem{

    public:

    Elem(std::vector<darray> nXYZ,shapeStruct shape);
    matrix selfS(quadStruct quad,darray xyz);
    matrix selfK(quadStruct quad,matrix D);

    std::vector<darray> gXYZ;
    darray detJ;
    matrix dxN;
    matrix dyN;
    matrix dzN;
    matrix N;
    int nLen;
};

class Face{

    public:

    Face(std::vector<darray> nXYZ,shapeStruct shape);
    matrix selfN(shapeStruct shape,quadStruct quad);

    darray detJ;
    matrix dxN;
    matrix dyN;
    int nLen;
};

#endif