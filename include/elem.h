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
    matrix localK(quadStruct quad,matrix D);

    int nLen;
    matrix N;
    darray detJ;
    matrix dxN;
    matrix dyN;
    matrix dzN;
    matrix gXYZ;
};

class Face{

    public:

    Face(std::vector<darray> nXYZ,shapeStruct shape);
    matrix localN(shapeStruct shape,quadStruct quad);

    int nLen;
    darray detJ;
    matrix dxN;
    matrix dyN;
};

#endif