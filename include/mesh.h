#include "elem.h"
#ifndef MESH_H
#define MESH_H

// ----------------------------------------|
// Structure storing the mesh parameters   |
// ----------------------------------------|

struct meshStruct{

    // Main informations about the mesh

    int order;
    std::vector<array3d> EvR;
    std::vector<std::array<double,2>> EvS;
    std::vector<array3d> nXYZ;
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

    sparse totalK();
    sparse totalM();
    darray neumann();
    shapeStruct shape(quadStruct quad);
    std::vector<darray> stress(darray &u);

    template<class type>type shape2(quadStruct quad);

    // Constructor and operations on the stiffness matrix
    
    Mesh(meshStruct &mesh);
    void update(darray &u);
    void complete(darray &u);
    void delta(sparse &K,darray &B);
    void coupling(sparse&K,darray &B);
    void dirichlet(sparse&K,darray &B);

    // Internal variables of the system

    meshStruct mesh;
    shapeStruct shape3D;
    shapeStruct shape2D;
    shapeStruct shapeS[6];

    // Number of nodes, elements and faces
    
    int nLen;
    int eLen;
    int fLen;
};

#endif