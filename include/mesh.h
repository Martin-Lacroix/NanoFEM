#include "elem.h"
#ifndef MESH_H
#define MESH_H

// ----------------------------------------|
// Structure storing the mesh parameters   |
// ----------------------------------------|

struct meshStruct{

    int ordElem;
    int ordQuad;

    // Main informations about the mesh

    std::vector<dvector> Ev;
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
    dvector stress(darray &u);
    matrix totalS(dvector xyz);
    shapeStruct shape(int dim,int order);

    // Operations on the total stiffness matrix
    
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