#include "elem.h"
#ifndef MESH_H
#define MESH_H

struct meshStruct{

    // Nodes indices of each element [m,8]
    // Node coordinates [n,3]
    // Stiffness tensor [6,6]

    std::vector<iarray> eNode;
    std::vector<darray> nXYZ;
    matrix D;
};

struct bcStruct{

    // Values of Neumann BC [n,3]
    // Face-node indices [n,4]

    std::vector<darray> fVal;
    std::vector<iarray> fNode;

    // Values of Dirichlet BC [3,m]
    // Node indices [3,m]

    std::vector<darray> nVal;
    std::vector<iarray> nIdx;
};

class Mesh{

    public:

    sparse localK();
    darray neumann();
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