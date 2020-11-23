#include "elem.h"
#ifndef MESH_H
#define MESH_H

struct meshStruct{

    // Quadrature order
    // Nodes indices of each element [m,8]
    // Nodes indices of each face [m,4]
    // Node coordinates int [n,3]
    // Stiffness tensor [6,6]

    int order;
    std::vector<darray> nXYZ;
    std::vector<iarray> eNode;
    std::vector<iarray> fNode;
    matrix D;
};

struct bcStruct{

    // Values of the applied stress on the faces

    std::vector<darray> neumann;

    // Dimension first and and node index second

    std::vector<std::vector<int>> dirichlet;
    
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