#include "elem.h"
#ifndef MESH_H
#define MESH_H

struct meshStruct{

    // D = [element] [stiffness tensor]
    // fNode = [face index] [face nodes]
    // order = order of the quadrature rule
    // nXYZ = [node index] [node coordinates]
    // eNode = [element index] [element nodes]
    // fElem = [element index of the free face]

    int order;
    ivector fElem;
    std::vector<matrix> D;
    std::vector<dvector> nXYZ;
    std::vector<ivector> eNode;
    std::vector<ivector> fNode;

    // neuFace = [face index]
    // dirNode = [dimension] [node index]
    // dirValue = [dimension] [displacement]
    // neuValue = [neuFace index] [traction vector]
    // perNode = [antiperiodic axis] [node index pair]
    
    ivector neuFace;
    std::vector<darray> neuValue;
    std::vector<ivector> dirNode;
    std::vector<dvector> dirValue;
    std::vector<pivector> perNode;
};

class Mesh{

    public:

    // Functions for computing the system matrix

    sparse localK();
    darray neumann();
    sparse nonLocalK();
    matrix elemK(int idx);
    shapeStruct shape(int dim);
    matrix totalS(dvector xyz);
    void periodic(sparse&K,darray &B);
    void dirichlet(sparse&K,darray &B);
    void complete(darray &u);
    Mesh(meshStruct &mesh);

    // Internal variables of the system

    meshStruct mesh;
    quadStruct quad3D;
    quadStruct quad2D;
    shapeStruct shape3D;
    shapeStruct shape2D;

    // List of elements and faces

    std::vector<Elem> eList;
    std::vector<Face> fList;
};

#endif