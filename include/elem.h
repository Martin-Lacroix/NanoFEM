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
    std::vector<matrix> jacobian(shapeStruct shape);
    darray stress(shapeStruct shape,matrix D,darray u);
    void derivative(shapeStruct shape,std::vector<matrix> J);

    // Functions of elemental stiffness and mass matrices

    matrix selfG(shapeStruct shape);
    matrix selfS(shapeStruct shape);
    matrix selfK(shapeStruct shape,matrix D);
    matrix selfM(shapeStruct shape,double rho);
    matrix selfKs(shapeStruct shape[6],matrix D);

    // Parameters specific to each element

    std::vector<array3d> nXYZ;
    ivector surface;
    dvector detJ;
    matrix dN[3];
    int nLen;
};

// --------------------------------------------------|
// Class of 4-node quadrangle linear element in 2D   |
// --------------------------------------------------|

class Face{

    public:

    // Constructor and functions available in the elem.cpp file

    Face(std::vector<array3d> nXYZ,shapeStruct shape);
    darray selfB(shapeStruct shape,darray F);
    darray dJ2D;
    int nLen;
};

#endif