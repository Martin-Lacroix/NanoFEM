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

    void clean();
    void updateJ(shapeStruct &shape);
    void updateS(shapeStruct (&shape)[6]);
    void updateF(shapeStruct &shape,darray u);

    // Stress extraction functions

    double stress(shapeStruct &shape,array3d LmR);
    double stress(shapeStruct &shape,array3d LmR,darray u);

    // Constructor and functions available in the file

    Elem(std::vector<array3d> nXYZ,ivector surface);
    dvector surfaceJ(shapeStruct &shape,ivector &node,int index);

    // Functions of elemental stiffness and mass matrices

    matrix selfM(shapeStruct &shape,double rho);
    matrix selfK(shapeStruct &shape,array3d LmR);
    darray selfFS(shapeStruct (&shape)[6],array3d LmS);
    matrix selfKS(shapeStruct (&shape)[6],array3d LmS);

    // Functions for large deformation elasticity

    matrix selfKN(shapeStruct &shape,array3d LmR);
    matrix selfKL(shapeStruct &shape,array3d LmR);
    darray selfFX(shapeStruct &shape,array3d LmR);

    matrix selfKNS(shapeStruct (&shape)[6],array3d LmR);
    matrix selfKLS(shapeStruct (&shape)[6],array3d LmR);

    // Parameters specific to each element

    std::vector<array3d> nXYZ;
    ivector surface;
    int nLen;

    // Parameters at quadrature nodes

    std::vector<darray3d> norm[6];
    std::vector<matrix3d> F;
    std::vector<darray> E;
    dvector detJ2D[6];
    matrix dNs[6][3];
    dvector detJ;
    matrix dN[3];
};

// --------------------------------------------------|
// Class of 4-node quadrangle linear element in 2D   |
// --------------------------------------------------|

class Face{

    public:

    // Constructor and functions available in the elem.cpp file

    Face(std::vector<array3d> nXYZ,shapeStruct &shape);
    darray selfFT(shapeStruct &shape,darray F);
    dvector detJ2D;
    int nLen;
};

#endif