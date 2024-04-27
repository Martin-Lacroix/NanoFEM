#include "mesh.h"
using namespace std;

// -------------------------------------------------------------|
// Builds the elements list, shape functions and quadratures    |
// -------------------------------------------------------------|

Mesh::Mesh(dataStruct &&input) : data{move(input)}{
    
    quadStruct quad;
    nLen = data.nXYZ.size();
    eLen = data.eNode.size();
    fLen = data.neuFace.size();
    int sLen = 0.1+pow(data.order+1,3);

    // Stores the quadrature rules and shape functions

    quad = math::legendre(3,data.order);
    shape3D = shape(quad);

    quad = math::legendre(2,data.order);
    vector<quadStruct> quadS(6,quad);
    shape2D = shape(quad);

    // Stores the quardature rules for surface FEM

    for(int i = 0; i < quad.gLen; i++){
        for(int j = 1; j < 3; j++){

            quadS[j-1].gRST[i] = {quad.gRST[i][0],quad.gRST[i][1],pow(-1,j)};
            quadS[j+1].gRST[i] = {quad.gRST[i][0],pow(-1,j),quad.gRST[i][1]};
            quadS[j+3].gRST[i] = {pow(-1,j),quad.gRST[i][0],quad.gRST[i][1]};
        }
    }

    // Stores the shape functions for surface FEM

    for(int i = 0; i < 6; i++){
        shapeS[i] = shape(quadS[i]);
    }

    // Creates the list of elements with node coordinates

    for(int i = 0; i < eLen; i++){
        
        vector<array3d> eXYZ(sLen);
        for(int j = 0; j < sLen; j++) eXYZ[j] = data.nXYZ[data.eNode[i][j]];
        Elem hexa(eXYZ,data.eSurf[i]);
        elem.push_back(hexa);
    }
}

// -------------------------------------------------------------|
// Returns the structure of 2D and 3D linear shape functions    |
// -------------------------------------------------------------|

shapeStruct Mesh::shape(quadStruct &quad){

    shapeStruct shape;
    dvector node(data.order+1);
    int dim = quad.gRST[0].size();
    int nLen = 0.1+pow(data.order+1,dim);

    // Copies the quadrature in the shape structure

    shape.gLen = quad.gLen;
    shape.gRST = quad.gRST;
    shape.weight = quad.weight;

    // Memory allocation and 1D node list

    shape.dN.resize(dim);
    shape.N.resize(nLen,quad.gLen);
    for(int i = 0; i <= data.order; i++) node[i] = i*2.0/data.order-1;
    for(int i = 0; i < dim; i++) shape.dN[i].resize(nLen,quad.gLen);

    // Sets the coordinates of the Gauss points

    for(int i = 0; i < shape.gLen; i++){

        vector<dvector> dN(dim);
        dvector N = math::lagrange(-1,node,quad.gRST[i]);
        for(int j = 0; j < dim; j++) dN[j] = math::lagrange(j,node,quad.gRST[i]);

        // Stores the shape functions evaluated at the Gauss points

        for(int j = 0; j < nLen; j++){
            shape.N(j,i) = N[j];
            
            for(int k = 0; k < dim; k++){
                shape.dN[k](j,i) = dN[k][j];
            }
        }
    }
    return shape;
}

// -----------------------------------------------------|
// Computes the total stiffness matrix K for SST/SET    |
// -----------------------------------------------------|

void Mesh::totalKB(sparse &K,darray &B){

    int sLen = shape3D.N.rows();
    for(int i = 0; i < eLen; i++){

        // Computes the elemental matrices and vectors

        elem[i].updateS(shapeS);
        elem[i].updateJ(shape3D);
        matrix Ke = elem[i].selfK(shape3D,data.LmR[i]);
        matrix Ks = elem[i].selfKS(shapeS,data.LmS[i]);
        darray Fs = elem[i].selfFS(shapeS,data.LmS[i]);
        elem[i].clean();

        // Inserts the elemental vector into the global B vector

        for(int k = 0; k < 3; k++){
            for(int j = 0; j < sLen; j++){

                int idx = data.eNode[i][j]+k*nLen;
                B(idx) += Fs(j+k*sLen);
            }
        }
        
        // Inserts the elemental matrix into the global K matrix

        for(int m = 0; m < 3; m++){
            for(int n = 0; n < 3; n++){
                for(int j = 0; j < sLen; j++){
                    for(int k = 0; k < sLen; k++){

                        int row = data.eNode[i][j]+m*nLen;
                        int col = data.eNode[i][k]+n*nLen;

                        // Adds the element only in the upper triangle

                        if(row <= col){

                            double val = Ke(j+m*sLen,k+n*sLen)+Ks(j+m*sLen,k+n*sLen);
                            K.coeffRef(row,col) += val;
                        }
                    }
                }
            }
        }
    }
}

// -------------------------------------------------|
// Computes the total stiffness matrix K for SVK    |
// -------------------------------------------------|

void Mesh::totalKT(sparse &K,darray &B,darray &u){

    int sLen = shape3D.N.rows();
    darray ue(3*sLen);

    for(int i = 0; i < eLen; i++){

        // Stores the nodal displacement of the element

        for(int k = 0; k < 3; k++){
            for(int j = 0; j < sLen; j++){

                int idx = data.eNode[i][j]+k*nLen;
                ue(j+k*sLen) = u(idx);
            }
        }

        // Computes the elemental matrices and vectors

        elem[i].updateJ(shape3D);
        elem[i].updateF(shape3D,ue);
        matrix Kn = elem[i].selfKN(shape3D,data.LmR[i]);
        matrix Kl = elem[i].selfKL(shape3D,data.LmR[i]);
        darray Fx = elem[i].selfFX(shape3D,data.LmR[i]);
        elem[i].clean();

        // Inserts the elemental vector into the global B vector

        for(int k = 0; k < 3; k++){
            for(int j = 0; j < sLen; j++){

                int idx = data.eNode[i][j]+k*nLen;
                B(idx) -= Fx(j+k*sLen);
            }
        }
        
        // Inserts the elemental matrix into the global K matrix

        for(int m = 0; m < 3; m++){
            for(int n = 0; n < 3; n++){
                for(int j = 0; j < sLen; j++){
                    for(int k = 0; k < sLen; k++){

                        int row = data.eNode[i][j]+m*nLen;
                        int col = data.eNode[i][k]+n*nLen;

                        // Adds the element only in the upper triangle

                        if(row <= col){
                            
                            double val = Kn(j+m*sLen,k+n*sLen)+Kl(j+m*sLen,k+n*sLen);
                            K.coeffRef(row,col) += val;
                        }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------|
// Computes the total mass matrix M for global FEM    |
// ---------------------------------------------------|

void Mesh::totalM(sparse &M){

    int sLen = shape3D.N.rows();
    for(int i = 0; i < eLen; i++){

        // Computes the elemental M matrices

        elem[i].updateJ(shape3D);
        matrix Me = elem[i].selfM(shape3D,data.LmR[i][2]);
        elem[i].clean();

        // Inserts the elemental matrix into the global M matrix

        for(int m = 0; m < 3; m++){
            for(int n = 0; n < 3; n++){
                for(int j = 0; j < sLen; j++){
                    for(int k = 0; k < sLen; k++){

                        int row = data.eNode[i][j]+m*nLen;
                        int col = data.eNode[i][k]+n*nLen;

                        // Adds the element only in the upper triangle

                        if(row <= col){

                            double val = Me(j+m*sLen,k+n*sLen);
                            M.coeffRef(row,col) += val;
                        }
                    }
                }
            }
        }
    }
}

// -------------------------------------------------------------|
// Computes the vector of applied stresses B from Neumann BC    |
// -------------------------------------------------------------|

void Mesh::neumann(darray &B){

    int sLen = shape2D.N.cols();
    for(int i = 0; i < fLen; i++){

        // Coordinates of the nodes composing the faces

        vector<array3d> eXYZ(sLen);
        for(int j = 0; j < sLen; j++) eXYZ[j] = data.nXYZ[data.neuFace[i][j]];

        // Computes the elemental B vectors

        Face face(eXYZ,shape2D);
        darray Fe = face.selfFT(shape2D,data.neuVal[i]);

        // Inserts the elemental vectors into the global B vector

        for(int j = 0; j < sLen; j++){
            for(int k = 0; k < 3; k++){
                B(data.neuFace[i][j]+k*nLen) += Fe(j+k*sLen);
            }
        }
    }
}

// ----------------------------------------------------------------|
// Edits the matrix K and vector B for fixed nodal displacement    |
// ----------------------------------------------------------------|

void Mesh::dirichlet(sparse &K,darray &B){

    vector<ivector> row = math::sparsemap(K);

    // Gets the number of nodes in the dirichlet list

    for(int n = 0; n < 3; n++){
        int dLen = data.dirNode[n].size();

        // Gets the index of the node in the solution vector

        for(int i = 0; i < dLen; i++){
            int idx = data.dirNode[n][i]+n*nLen;

            // Substract the coefficients from K to B and cancels the row-colums

            for(int j : row[idx]){

                B(j) -= math::get(K,j,idx)*data.dirVal[n][i];
                math::symset(K,j,idx,0);
            }
            ivector().swap(row[idx]);
        }
    }

    // Places the dirichlet value the B vector and 1 in the diagonal of K

    for(int n = 0; n < 3; n++){
        for(int i = 0; i < data.dirNode[n].size(); i++){
            
            int idx = data.dirNode[n][i]+n*nLen;
            K.coeffRef(idx,idx) = 1;
            B(idx) = data.dirVal[n][i];
        }
    }
}

// ------------------------------------------------------------------|
// Edits the matrix K and vector B for coupled nodal displacement    |
// ------------------------------------------------------------------|

void Mesh::coupling(sparse &K,darray &B){

    vector<ivector> row = math::sparsemap(K);

    for(int n = 0; n < 3; n++){

        // Gets the number of nodes in the periodic list

        for(int i = 0; i < data.coupNode[n].size(); i++){
            int dLen = data.coupNode[n][i].size();

            // Gets the index of the curent and final coupled nodes

            for(int j = 0; j < dLen-1; j++){

                int idx1 = data.coupNode[n][i][j]+n*nLen;
                int idx2 = data.coupNode[n][i].back()+n*nLen;
                double fix = math::get(K,idx2,idx2)+math::get(K,idx1,idx1)+2*math::get(K,idx1,idx2);

                // Adds the current row to the final row of K if non-zero

                for(int k = 0; k < row[idx1].size(); k++){

                    int n = row[idx1][k];
                    double val = math::get(K,idx1,n);

                    // Updates the mapStruct only if an element is added to K

                    if(val != 0){
                        if(math::get(K,idx2,n) == 0){

                            row[n].push_back(idx2);
                            row[idx2].push_back(n);
                        }
                        math::symadd(K,idx2,n,val);
                    }
                }

                // Zeroes the column and row in K corresoponding to the coupled node

                for(int k : row[idx1]) math::symset(K,idx1,k,0);
                row[idx1] = ivector();
    
                // Corrects the diagonal elements in K and edits B

                K.coeffRef(idx2,idx2) = fix;
                K.coeffRef(idx1,idx1) = 1;
                B(idx2) += B(idx1);
                B(idx1) = 0;
            }
        }
    }
}

// --------------------------------------------------------------------|
// Edits the matrix K and vector B for a change of variable u => Δu    |
// --------------------------------------------------------------------|

void Mesh::delta(sparse &K,darray &B){

    vector<ivector> row = math::sparsemap(K);

    // Gets the number of nodes in the periodic list

    for(int n = 0; n < 3; n++){
        for(int i = 0; i < data.deltaNode[n].size(); i++){

            // Gets the index of the curent and final coupled nodes

            int idx1 = data.deltaNode[n][i].first+n*nLen;
            int idx2 = data.deltaNode[n][i].second+n*nLen;
            double fix = math::get(K,idx2,idx2)+math::get(K,idx1,idx1)+2*math::get(K,idx1,idx2);

            // Adds the current row to the final row of K if non-zero

            for(int k = 0; k < row[idx1].size(); k++){

                int n = row[idx1][k];
                double val = math::get(K,idx1,n);

                // Updates the mapStruct only if an element is added to K

                if(val != 0){
                    if(math::get(K,idx2,n) == 0){

                        row[n].push_back(idx2);
                        row[idx2].push_back(n);
                    }
                    math::symadd(K,idx2,n,val);
                }
            }

            // Corrects the diagonal elements in K and edits B

            K.coeffRef(idx2,idx2) = fix;
            B(idx2) += B(idx1);
        }
    }
}

// -----------------------------------------------------|
// Completes the solution for coupled or delta nodes    |
// -----------------------------------------------------|

void Mesh::complete(darray &u){

    for(int n = 0; n < 3; n++){

        // Copy the value of the last node into the coupled ones

        for(int i = 0; i < data.coupNode[n].size(); i++){
            for(int j = 0; j < data.coupNode[n][i].size(); j++){

                int idx1 = data.coupNode[n][i][j]+n*nLen;
                int idx2 = data.coupNode[n][i].back()+n*nLen;
                u[idx1] = u[idx2];
            }
        }

        // Comes back from the difference Δu to the actual displacement

        for(int i = 0; i < data.deltaNode[n].size(); i++){

            int idx1 = data.deltaNode[n][i].first+n*nLen;
            int idx2 = data.deltaNode[n][i].second+n*nLen;
            u[idx1] += u[idx2];
        }
    }
}

// ---------------------------------------------------------|
// Updates the coordinates of the node with the solution    |
// ---------------------------------------------------------|

void Mesh::update(darray &u){

    int sLen = shape3D.N.rows();

    // Adds the displacement to the corresponding node coordinate

    for(int i = 0; i < nLen; i++){
        for(int j = 0; j < 3; j++){
            data.nXYZ[i][j] += u[i+j*nLen];
        }
    }

    // Creates the list of elements with node coordinates

    for(int i = 0; i < eLen; i++){
        
        vector<array3d> eXYZ(sLen);
        for(int j = 0; j < sLen; j++) eXYZ[j] = data.nXYZ[data.eNode[i][j]];
        Elem hexa(eXYZ,data.eSurf[i]);
        elem[i] = hexa;
    }
}

// ---------------------------------------------|
// Averaged Von Mises stress in the elements    |
// ---------------------------------------------|

dvector Mesh::stress(darray &u,int svk){

    int sLen = shape3D.N.rows();
    darray ue(3*sLen);
    dvector VM(eLen);

    // Coordinates of the nodes of the element

    for(int i = 0; i < eLen; i++){

        // Stores the nodal displacement of the element

        for(int k = 0; k < 3; k++){
            for(int j = 0; j < sLen; j++){

                int idx = data.eNode[i][j]+k*nLen;
                ue(j+k*sLen) = u(idx);
            }
        }

        // Computes the averaged Von Mises stress

        if(svk){

            elem[i].updateJ(shape3D);
            elem[i].updateF(shape3D,ue);
            VM[i] = elem[i].stress(shape3D,data.LmR[i]);
            elem[i].clean();
        }
        else{

            elem[i].updateJ(shape3D);
            VM[i] = elem[i].stress(shape3D,data.LmR[i],ue);
            elem[i].clean();
        }
    }
    return VM;
}