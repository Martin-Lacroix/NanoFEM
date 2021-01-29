#include "elem.h"
#ifndef MESH_H
#define MESH_H

struct meshStruct{

    // Main informations about the mesh

    int order;
    std::vector<matrix> D;
    std::vector<dvector> nXYZ;
    std::vector<ivector> eNode;
    std::vector<ivector> eSurf;

    // Neumann boundary conditions
    
    std::vector<ivector> neuNode;
    std::vector<darray> neuVal;

    // Dirichlet boundary conditions

    std::vector<ivector> dirNode;
    std::vector<dvector> dirVal;

    // Periodic boundary conditions

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