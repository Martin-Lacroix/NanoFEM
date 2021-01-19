#include "elem.h"
#ifndef MESH_H
#define MESH_H

struct meshStruct{

    // nXYZ = [node index] [node coordinates]
    // eNode = [element index] [element nodes]
    // fNode = [face index] [face nodes]

    std::vector<dvector> nXYZ;
    std::vector<ivector> eNode;
    std::vector<ivector> fNode;

    // Order = order of the quadrature rule
    // fElem = [element index of the free face]
    // D = [element] [stiffness tensor]
    
    int order;
    ivector fElem;
    std::vector<matrix> D;

    // dirNode = [dimension] [node index]
    // dirValue = [dimension] [displacement]
    
    std::vector<dvector> dirValue;
    std::vector<darray> neuValue;

    // neuFace = [face index]
    // neuValue = [neuFace index] [traction vector]

    ivector neuFace;
    std::vector<ivector> dirNode;

    // perNode = [antiperiodic axis] [node index pair]

    std::vector<pvector> perNode;
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