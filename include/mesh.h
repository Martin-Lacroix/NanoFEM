#include "elem.h"
#ifndef MESH_H
#define MESH_H

// ----------------------------------------|
// Structure storing the mesh parameters   |
// ----------------------------------------|

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

    // Periodic and coupled boundary conditions

    std::vector<std::vector<ivector>> coupNode;
};

// -------------------------------|
// Class of finite element mesh   |
// -------------------------------|

class Mesh{

    public:

    // Functions available in the mesh.cpp file

    sparse localK();
    darray neumann();
    sparse nonLocalK();
    matrix elemK(int idx);
    shapeStruct shape(int dim);
    matrix totalS(dvector xyz);
    void coupling(sparse&K,darray &B);
    void dirichlet(sparse&K,darray &B);
    void complete(darray &u);
    void update(darray &u);
    Mesh(meshStruct &mesh);

    // Internal variables of the system

    meshStruct mesh;
    quadStruct quad3D;
    quadStruct quad2D;
    shapeStruct shape3D;
    shapeStruct shape2D;
};

#endif