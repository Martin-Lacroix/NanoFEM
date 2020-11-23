#include "elem.h"
#ifndef MESH_H
#define MESH_H

struct meshStruct{

    // Nodes indices of each element
    // Nodes indices of each face
    // Global node coordinates

    std::vector<darray> nXYZ;
    std::vector<iarray> eNode;
    std::vector<iarray> fNode;
};

struct otherStruct{

    // Stiffness tensor and quadrature order
    // Neumann = values of the applied stress on the faces
    // Dirichlet = dimension first and and node index second

    matrix D;
    int order;
    std::vector<darray> neumann;
    std::vector<std::vector<int>> dirichlet;
    
};

class Mesh{

    public:

    // Functions for computing the system matrix

    sparse localK();
    darray neumann();
    sparse nonLocalK();
    matrix elemK(int idx);
    matrix totalS(darray xyz);
    void dirichlet(darray &B);
    void dirichlet(sparse &K);
    shapeStruct shape(int dim);
    Mesh(meshStruct mesh,otherStruct bc);

    // Internal variables of the FEM system

    quadStruct eQuad;
    quadStruct fQuad;
    shapeStruct eShape;
    shapeStruct fShape;
    meshStruct meshData;
    otherStruct otherData;
    std::vector<Elem> eList;
    std::vector<Face> fList;
};

#endif