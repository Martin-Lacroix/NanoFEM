#include "..\include\mesh.h"
using namespace std;

// -------------------------------------------------------------|
// Builds the elements list, shape functions and quadratures    |
// -------------------------------------------------------------|

Mesh::Mesh(meshStruct &input){
    
    mesh = input;
    int fLen = mesh.neuNode.size();
    int eLen = mesh.eNode.size();

    // Stores the quadrature rules and shape functions

    quad3D = math::legendre(3,mesh.order);
    quad2D = math::legendre(2,mesh.order);
    shape3D = shape(3);
    shape2D = shape(2);

    // Stores the mesh elements into a vector

    for(int i=0; i<eLen; i++){

        vector<dvector> eXYZ;
        int nLen = mesh.eNode[i].size();
        for(int j=0; j<nLen; j++){eXYZ.push_back(mesh.nXYZ[mesh.eNode[i][j]]);}
        Elem elem(eXYZ,shape3D);
        eList.push_back(elem);
    }

    // Stores the mesh faces into a vector

    for(int i=0; i<fLen; i++){

        vector<dvector> fXYZ;
        int nLen = mesh.neuNode[i].size();
        for(int j=0; j<nLen; j++){fXYZ.push_back(mesh.nXYZ[mesh.neuNode[i][j]]);}
        Face face(fXYZ,shape2D);
        fList.push_back(face);
    }
}

// -------------------------------------------------------------|
// Returns the structure of 2D and 3D linear shape functions    |
// -------------------------------------------------------------|

shapeStruct Mesh::shape(int dim){

    int nLen,gLen;
    shapeStruct shape;
    if(dim==3){gLen = quad3D.gRST.size(); nLen = 8;}
    if(dim==2){gLen = quad2D.gRST.size(); nLen = 4;}
    vector<ivector> n(nLen,ivector(dim));
    shape.gLen = gLen;
    shape.nLen = nLen;

    // Performs the memory allocation

    shape.drN.setlength(nLen,gLen);
    shape.dsN.setlength(nLen,gLen);
    shape.dtN.setlength(nLen,gLen);
    shape.N.setlength(nLen,gLen);

    if(dim==3){

        // Sets the coordinates of the 3D nodes

        n[0]={-1,-1,-1}; n[1]={1,-1,-1}; n[2]={1,1,-1}; n[3]={-1,1,-1};
        n[4]={-1,-1,1}; n[5]={1,-1,1}; n[6]={1,1,1}; n[7]={-1,1,1};

        for(int i=0; i<gLen; i++){

            // Extracts the 3D coordinates of the Gauss points

            double r = quad3D.gRST[i][0];
            double s = quad3D.gRST[i][1];
            double t = quad3D.gRST[i][2];

            for(int j=0; j<nLen; j++){

                // Stores the shape functions evaluated at the Gauss points

                shape.N(j,i) = (1+n[j][0]*r)*(1+n[j][1]*s)*(1+n[j][2]*t)/8;
                shape.drN(j,i) = n[j][0]*(1+n[j][1]*s)*(1+n[j][2]*t)/8;
                shape.dsN(j,i) = n[j][1]*(1+n[j][0]*r)*(1+n[j][2]*t)/8;
                shape.dtN(j,i) = n[j][2]*(1+n[j][0]*r)*(1+n[j][1]*s)/8;
            }
        }
    }

    if(dim==2){

        // Sets the coordinates of the 2D nodes

        n[0]={-1,-1}; n[1]={1,-1};
        n[2]={1,1}; n[3]={-1,1};

        for(int i=0; i<gLen; i++){

            // Extracts the 2D coordinates of the Gauss points

            double r = quad2D.gRST[i][0];
            double s = quad2D.gRST[i][1];

            for(int j=0; j<nLen; j++){

                // Stores the shape functions evaluated at the Gauss pojnts

                shape.N(j,i) = (1+n[j][0]*r)*(1+n[j][1]*s)/4;
                shape.drN(j,i) = n[j][0]*(1+n[j][1]*s)/4;
                shape.dsN(j,i) = n[j][1]*(1+n[j][0]*r)/4;
            }
        }
    }
    return shape;
}

// -----------------------------------------------------|
// Computes the total stiffness matrix for local FEM    |
// -----------------------------------------------------|

sparse Mesh::localK(){

    sparse K;
    double eLen = eList.size();
    double nLen = mesh.nXYZ.size();
    alglib::sparsecreate(3*nLen,3*nLen,K);

    // Computes the elemental K matrices

    for(int i=0; i<eLen; i++){

        int nbr = eList[i].nLen;
        ivector eNode = mesh.eNode[i];
        matrix K1 = eList[i].selfK(quad3D,mesh.D[i]);

        // Inserts the elemental matrices into the global K matrix

        for(int j=0; j<nbr; j++){
            for(int k=0; k<nbr; k++){
                for(int m=0; m<3; m++){
                    for(int n=0; n<3; n++){

                        double val = K1[j+m*nbr][k+n*nbr];
                        alglib::sparseadd(K,eNode[j]+m*nLen,eNode[k]+n*nLen,val);
                    }
                }
            }
        }
    }
    return K;
}

// -----------------------------------------------------------|
// Computes the total stiffness matrix K for non-local FEM    |
// -----------------------------------------------------------|

sparse Mesh::nonLocalK(){

    sparse K;
    double eLen = eList.size();
    double nLen = mesh.nXYZ.size();
    alglib::sparsecreate(3*nLen,3*nLen,K);

    // Computes the elemental K matrices

    for(int i=0; i<eLen; i++){

        matrix K1 = elemK(i);
        int nbr = eList[i].nLen;
        ivector eNode = mesh.eNode[i];

        // Inserts the elemental matrices into the global K matrix

        for(int k=0; k<3*nLen; k++){
            for(int j=0; j<nbr; j++){
                for(int m=0; m<3; m++){

                    double val = K1[j+m*nbr][k];
                    alglib::sparseadd(K,eNode[j]+m*nLen,k,val);
                }
            }
        }
    }
    return K;
}

// ---------------------------------------------------------------------|
// Evaluates the total S matrix at a point (x,y,z) for non-local FEM    |
// ---------------------------------------------------------------------|

matrix Mesh::totalS(dvector xyz){

    matrix S;
    double eLen = eList.size();
    double nLen = mesh.nXYZ.size();
    S.setlength(6,3*nLen);
    math::zero(S);

    // Computes the elemental S matrices

    for(int i=0; i<eLen; i++){

        int nbr = eList[i].nLen;
        matrix S1 = eList[i].selfS(quad3D,xyz);

        // Inserts the elemental matrices into the global S matrix

        for(int j=0; j<6; j++){
            for(int k=0; k<nbr; k++){
                for(int n=0; n<3; n++){
                    S(j,mesh.eNode[i][k]+n*nLen) += S1(j,k+n*nbr);
                }
            }
        }
    }
    return S;
}

// ---------------------------------------------------------------|
// Computes the elemental stiffness matrix K for non-local FEM    |
// ---------------------------------------------------------------|

matrix Mesh::elemK(int idx){

    matrix B,K;
    Elem& elem = eList[idx];
    double nLen = elem.nLen;
    double gLen = quad3D.gRST.size();
    K.setlength(3*nLen,3*mesh.nXYZ.size());
    B.setlength(3*nLen,6);
    math::zero(B);
    math::zero(K);

    // Performs the numerical integration

    for(int i=0; i<gLen; i++){
        for(int j=0; j<nLen; j++){

            // Computes the shape functions derivative matrix B
            
            B(j,0) = B(j+nLen,3) = B(j+2*nLen,4) = elem.dxN(j,i);
            B(j,3) = B(j+nLen,1) = B(j+2*nLen,5) = elem.dyN(j,i);
            B(j,4) = B(j+nLen,5) = B(j+2*nLen,2) = elem.dzN(j,i);
        }

        // Computes K by Gauss-Legendre quadrature

        matrix S = totalS(elem.gXYZ[i]);
        matrix K1 = math::prod(quad3D.weight[i],B,mesh.D[idx]);
        matrix K2 = math::prod(elem.detJ[i],K1,S);
        math::add(1,1,K2,K);
    }
    return K;
}

// -------------------------------------------------------------|
// Computes the vector of applied stresses B from Neumann BC    |
// -------------------------------------------------------------|

darray Mesh::neumann(){

    darray B;
    double nLen = mesh.nXYZ.size();
    double fLen = fList.size();
    B.setlength(3*nLen);
    math::zero(B);

    // Computes the applied surface tractions on the 2D elements

    for(int i=0; i<fLen; i++){

        int nbr = fList[i].nLen;
        matrix N = fList[i].selfN(shape2D,quad2D);
        darray B1 = math::prod(1,N,mesh.neuVal[i]);

        // Inserts the subvectors into the global vector

        for(int j=0; j<nbr; j++){
            for(int k=0; k<3; k++){
                B(mesh.neuNode[i][j]+k*nLen) += B1(j+k*nbr);
            }
        }
    }
    return B;
}

// ----------------------------------------------------|
// Edits the matrix K and vector B for Dirichlet BC    |
// ----------------------------------------------------|

void Mesh::dirichlet(sparse &K,darray &B){

    int nLen = mesh.nXYZ.size();
    matStruct mat = math::mapclean(K);

    // Gets the number of nodes in the dirichlet list

    for(int n=0; n<3; n++){
        int dLen = mesh.dirNode[n].size();

        // Gets the index of the node in the solution vector

        for(int i=0; i<dLen; i++){
            int idx = mesh.dirNode[n][i]+n*nLen;

            // Substract the node coefficients in K to B and cancels this row

            for(int j:mat.col[idx]){

                B(j) -= alglib::sparseget(K,j,idx)*mesh.dirVal[n][i];
                alglib::sparseset(K,j,idx,0);
            }

            // Also cancels the corresponding column in k

            for(int j:mat.row[idx]){
                alglib::sparseset(K,idx,j,0);
            }
        }
    }

    // Places the dirichlet value the B vector and 1 in the diagonal of K

    for(int n=0; n<3; n++){
        for(int m=0; m<mesh.dirNode[n].size(); m++){
            
            int idx = mesh.dirNode[n][m]+n*nLen;
            alglib::sparseset(K,idx,idx,1);
            B(idx) = mesh.dirVal[n][m];
        }
    }
}

// ------------------------------------------------------------|
// Edits the matrix K and vector B for periodic and flat BC    |
// ------------------------------------------------------------|

void Mesh::periodic(sparse &K,darray &B){

    double val;
    int nLen = mesh.nXYZ.size();
    matStruct mat = math::mapclean(K);

    // Gets the number of nodes in the periodic list

    for(int n=0; n<3; n++){
        for(int i=0; i<mesh.perNode[n].size(); i++){
            int dLen = mesh.perNode[n][i].size();

            // Gets the index of the curent and final coupled nodes

            for(int j=0; j<dLen-1; j++){

                int idx1 = mesh.perNode[n][i][j]+n*nLen;
                int idx2 = mesh.perNode[n][i].back()+n*nLen;

                // Adds the current row to the final row of K and cancels this row

                for(int k:mat.row[idx1]){
                    val = alglib::sparseget(K,idx1,k);

                    if(val!=0){

                        if(alglib::sparseget(K,idx2,k)==0){mat.col[k].push_back(idx2);}
                        alglib::sparseadd(K,idx2,k,val);
                        alglib::sparseset(K,idx1,k,0);
                    }
                }

                // Adds the current column to the final column of K and cancels this column

                for(int k:mat.col[idx1]){
                    val = alglib::sparseget(K,k,idx1);
                    
                    if(val!=0){

                        if(alglib::sparseget(K,k,idx2)==0){mat.row[k].push_back(idx2);}
                        alglib::sparseadd(K,k,idx2,val);
                        alglib::sparseset(K,k,idx1,0);
                    }
                }

                // Removes the curent value in B and adds it to the final one

                alglib::sparseset(K,idx1,idx1,1);
                B(idx2) += B(idx1);
                B(idx1) = 0;
            }
        }
    }
}

// ----------------------------------------------------|
// Completes the solution with periodic and flat BC    |
// ----------------------------------------------------|

void Mesh::complete(darray &u){

    int nLen = mesh.nXYZ.size();

    for(int n=0; n<3; n++){
        for(int i=0; i<mesh.perNode[n].size(); i++){
            for(int j=0; j<mesh.perNode[n][i].size(); j++){

                // Copy the value of the last node into the coupled ones

                int idx1 = mesh.perNode[n][i][j]+n*nLen;
                int idx2 = mesh.perNode[n][i].back()+n*nLen;
                u[idx1] = u[idx2];
            }
        }
    }
}