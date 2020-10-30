#include "elem.h"
#ifndef MESH_H
#define MESH_H

struct meshStruct{

    // Nodes indices of each element [m,8]
    // Node coordinates int [n,3]
    // Stiffness tensor [6,6]
    // Quadrature order

    std::vector<iarray> eNode;
    std::vector<darray> nXYZ;
    int order;
    matrix D;
};

struct bcStruct{

    // Values of Neumann BC [n,3]
    // Face-node indices int [n,4]

    std::vector<darray> fVal;
    std::vector<iarray> fNode;

    // Values of Dirichlet BC [3,m]
    // Node indices int [3,m]

    std::vector<darray> nVal;
    std::vector<iarray> nIdx;
};

class Mesh{

    public:

    sparse localK();
    darray neumann();
    sparse nonLocalK();
    matrix elemK(int idx);
    matrix totalS(darray xyz);
    void dirichlet(darray &B);
    void dirichlet(sparse &K);
    shapeStruct shape(int dim);
    Mesh(meshStruct mesh,bcStruct bc);

    bcStruct bcParam;
    quadStruct eQuad;
    quadStruct fQuad;
    shapeStruct eShape;
    shapeStruct fShape;
    meshStruct meshParam;
    std::vector<Elem> eList;
    std::vector<Face> fList;
};

#endif