#include "elem.h"
#ifndef MESH_H
#define MESH_H

struct meshStruct{

    // nXYZ = node index then coordinates
    // eNode = element index then nodes
    // fNode = face index then nodes

    std::vector<dvector> nXYZ;
    std::vector<ivector> eNode;
    std::vector<ivector> fNode;

    // Order = order of the quadrature rule
    // fElem = element indices corresponding to the faces
    // D = stiffness tensor of each element
    
    int order;
    ivector fElem;
    std::vector<matrix> D;

    // dirNode = dimension then node index
    // dirValue = dimension then displacement
    
    std::vector<dvector> dirValue;
    std::vector<darray> neuValue;

    // neuFace = face index
    // neuValue = face index then traction vector

    ivector neuFace;
    std::vector<ivector> dirNode;
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
    void dirichlet(sparse&K,darray &B);
    Mesh(meshStruct &mesh1);

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