#include "elem.h"
#ifndef MESH_H
#define MESH_H

// ----------------------------------------|
// Structure storing the mesh parameters   |
// ----------------------------------------|

struct dataStruct{

    int step;
    int order;
    double tol;

    // Main informations about the mesh

    std::vector<array3d> nXYZ;
    std::vector<ivector> eNode;
    std::vector<ivector> eSurf;

    // Neumann boundary conditions
    
    std::vector<ivector> neuFace;
    std::vector<darray3d> neuVal;

    // Dirichlet boundary conditions

    ivector dirNode[3];
    dvector dirVal[3];

    // Nodes with same displacement or (u1,u2) => (u1-u2,u2)

    std::vector<ivector> coupNode[3];
    std::vector<std::pair<int,int>> deltaNode[3];

    // Volume [λ μ ρ] and Surface [λs μs τ]

    std::vector<array3d> LmR;
    std::vector<array3d> LmS;
};

// -------------------------------|
// Class of finite element mesh   |
// -------------------------------|

class Mesh{

    public:
    Mesh(dataStruct &&data);

    // Functions for finite element matrices
    
    void totalB(darray &B);
    void totalM(sparse &M);
    void totalKB(sparse &K,darray &B);
    void totalKT(sparse &K,darray &B,darray &u);

    // COmputes the shape functions and stress
    
    shapeStruct shape(quadStruct &quad);
    dvector stress(darray &u,int svk);

    // Boundary conditions on the total system
    
    void update(darray &u);
    void neumann(darray &B);
    void complete(darray &u);
    void delta(sparse &K,darray &B);
    void coupling(sparse &K,darray &B);
    void dirichlet(sparse &K,darray &B);

    // Internal variables of the system

    dataStruct data;
    shapeStruct shape3D;
    shapeStruct shape2D;
    shapeStruct shapeS[6];
    std::vector<Elem> elem;

    // Number of nodes, elements and faces
    
    int nLen;
    int eLen;
    int fLen;
};

#endif