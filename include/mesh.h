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
    
    std::vector<ivector> neuFace;
    std::vector<darray> neuVal;

    // Dirichlet boundary conditions

    ivector dirNode[3];
    dvector dirVal[3];

    // Nodes with the same displacement or Δu

    std::vector<ivector> coupNode[3];
    std::vector<std::pair<int,int>> deltaNode[3];

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
    matrix totalS(dvector xyz);
    shapeStruct shape(int dim,int order);
    std::vector<darray> stress(darray &u);

    // Operations on the total stiffness matrix
    
    Mesh(meshStruct &mesh);
    void update(darray &u);
    void complete(darray &u);
    void delta(sparse &K,darray &B);
    void coupling(sparse&K,darray &B);
    void dirichlet(sparse&K,darray &B);

    // Internal variables of the system

    meshStruct mesh;
    quadStruct quad3D;
    quadStruct quad2D;
    shapeStruct shape3D;
    shapeStruct shape2D;

    // Number of nodes, elements and faces

    int nLen;
    int eLen;
    int fLen;
};

#endif