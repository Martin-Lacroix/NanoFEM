#include "math.h"
#ifndef ELEM_H
#define ELEM_H

// -------------------------------------------------------------------|
// Structure storing the shape functions in the element local space   |
// -------------------------------------------------------------------|

struct shapeStruct{

    int gLen;
    matrix N;
    dvector weight;

    // Local derivatives of the shape functions at Gauss nodes

    std::vector<matrix> dN;
    std::vector<dvector> gRST;
};

// --------------------------------------------------|
// Class of 8-node hexahedron linear element in 3D   |
// --------------------------------------------------|

class Elem{

    public:

    // Constructor and functions available in the elem.cpp file

    Elem(std::vector<array3d> nXYZ,ivector surface);
    darray stress(shapeStruct &shape,array3d EvR,darray u);
    dvector surfaceJ(shapeStruct &shape,ivector &node,int index);
    void updateJ(shapeStruct &shape);

    // Functions of elemental stiffness and mass matrices

    matrix selfM(shapeStruct &shape,double rho);
    matrix selfK(shapeStruct &shape,array3d EvR);
    std::pair<matrix,darray> selfKB(shapeStruct &shape,shapeStruct (&shapeS)[6],array3d EvS);

    // Parameters specific to each element

    std::vector<array3d> nXYZ;
    std::vector<matrix> J;
    ivector surface;
    dvector detJ;
    matrix dN[3];
    int sLen;
    int nLen;
};

// --------------------------------------------------|
// Class of 4-node quadrangle linear element in 2D   |
// --------------------------------------------------|

class Face{

    public:

    // Constructor and functions available in the elem.cpp file

    Face(std::vector<array3d> nXYZ,shapeStruct &shape);
    darray selfB(shapeStruct &shape,darray F);
    dvector detJ2D;
    int nLen;
};

#endif