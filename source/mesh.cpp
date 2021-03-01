#include "..\include\mesh.h"
using namespace std;

// -------------------------------------------------------------|
// Builds the elements list, shape functions and quadratures    |
// -------------------------------------------------------------|

Mesh::Mesh(meshStruct &mesh){
    
    quadStruct quad;
    this->mesh = mesh;
    nLen = mesh.nXYZ.size();
    eLen = mesh.eNode.size();
    fLen = mesh.neuFace.size();

    // Stores the quadrature rules and shape functions

    quad = math::legendre(3,mesh.order);
    shape3D = shape(quad);

    quad = math::legendre(2,mesh.order);
    vector<quadStruct> quadS(6,quad);
    shape2D = shape(quad);

    // Stores the quardature rules for surface FEM

    for(int i=0; i<quad.gLen; i++){
        for(int j=1; j<3; j++){

            quadS[j-1].gRST[i] = {pow(-1,j),quad.gRST[i][0],quad.gRST[i][1]};
            quadS[j+1].gRST[i] = {quad.gRST[i][0],pow(-1,j),quad.gRST[i][1]};
            quadS[j+3].gRST[i] = {quad.gRST[i][0],quad.gRST[i][1],pow(-1,j)};
        }
    }

    // Stores the shape functions for surface FEM

    for(int i=0; i<6; i++){
        shapeS[i] = shape(quadS[i]);
    }
}

// -------------------------------------------------------------|
// Returns the structure of 2D and 3D linear shape functions    |
// -------------------------------------------------------------|

shapeStruct Mesh::shape(quadStruct quad){

    shapeStruct shape;
    dvector node(mesh.order+1);
    int dim = quad.gRST[0].size();
    int nLen = 0.1+pow(mesh.order+1,dim);

    // Copies the quadrature in the shape structure

    shape.gLen = quad.gLen;
    shape.gRST = quad.gRST;
    shape.weight = quad.weight;

    // Memory allocation and 1D node list

    shape.dN.resize(dim);
    shape.N.setlength(nLen,quad.gLen);
    for(int i=0; i<=mesh.order; i++){node[i] = i*2.0/mesh.order-1;}
    for(int i=0; i<dim; i++){shape.dN[i].setlength(nLen,quad.gLen);}

    // Sets the coordinates of the Gauss points

    for(int i=0; i<shape.gLen; i++){

        vector<dvector> dN(dim);
        dvector N = math::lagrange(-1,node,quad.gRST[i]);
        for(int j=0; j<dim; j++){dN[j] = math::lagrange(j,node,quad.gRST[i]);}

        // Stores the shape functions evaluated at the Gauss points

        for(int j=0; j<nLen; j++){
            shape.N(j,i) = N[j];
            
            for(int k=0; k<dim; k++){
                shape.dN[k](j,i) = dN[k][j];
            }
        }
    }
    return shape;
}

// --------------------------------------------------------|
// Computes the total stiffness matrix K for global FEM    |
// --------------------------------------------------------|

sparse Mesh::totalK(){

    sparse K;
    int sLen = shape3D.N.cols();
    int size = 9*eLen*pow(mesh.order+1,6)/4;
    alglib::sparsecreate(3*nLen,3*nLen,size,K);

    // Coordinates of the nodes of the element

    for(int i=0; i<eLen; i++){
        
        vector<array3d> eXYZ(sLen);
        for(int j=0; j<sLen; j++){eXYZ[j] = mesh.nXYZ[mesh.eNode[i][j]];}

        // Computes the elemental K matrices

        Elem elem(eXYZ,mesh.eSurf[i]);
        matrix D = math::stiffness(mesh.EvR[i]);
        //matrix Ds = math::stiffness(mesh.EvS[i]);
        //matrix K2 = elem.selfKs(shapeS,Ds);
        matrix K1 = elem.selfK(shape3D,D);
        //math::add(1,1,K2,K1);

        // Inserts the elemental matrix into the global K matrix

        for(int j=0; j<sLen; j++){
            for(int k=0; k<sLen; k++){
                for(int m=0; m<3; m++){
                    for(int n=0; n<3; n++){

                        int row = mesh.eNode[i][j]+m*nLen;
                        int col = mesh.eNode[i][k]+n*nLen;

                        // Adds the element only in the upper triangle

                        if(row<=col){

                            double val = K1[j+m*sLen][k+n*sLen];
                            alglib::sparseadd(K,row,col,val);
                        }
                    }
                }
            }
        }
    }
    return K;
}






sparse Mesh::totalK2(){

    sparse K;
    int sLen = shape3D.N.cols();
    int size = 9*eLen*pow(mesh.order+1,6)/4;
    alglib::sparsecreate(9*nLen,9*nLen,size,K);

    // Coordinates of the nodes of the element

    for(int i=0; i<eLen; i++){
        
        vector<array3d> eXYZ(sLen);
        for(int j=0; j<sLen; j++){eXYZ[j] = mesh.nXYZ[mesh.eNode[i][j]];}

        // Computes the elemental K matrices

        Elem elem(eXYZ,mesh.eSurf[i]);
        matrix D = math::stiffness(mesh.EvR[i]);
        matrix Ka = elem.selfK(shape3D,D);



        D = math::couple(mesh.EvR[i][0],mesh.EvR[i][1],0);
        matrix Kb = elem.selfK(shape3D,D);
        matrix Ks = elem.selfM(shape3D,1);
        matrix Kg = elem.selfG(shape3D);

        // Inserts the elemental matrix into the global K matrix

        for(int j=0; j<sLen; j++){
            for(int k=0; k<sLen; k++){
                for(int m=0; m<3; m++){
                    for(int n=0; n<3; n++){

                        int row = mesh.eNode[i][j]+m*nLen;
                        int col = mesh.eNode[i][k]+n*nLen;

                        double val = Kg[j+m*sLen][k+n*sLen];
                        alglib::sparseadd(K,row,col+6*nLen,val); 

                        val = Ks[j+m*sLen][k+n*sLen];
                        alglib::sparseadd(K,row+3*nLen,col+6*nLen,val);

                        // Adds the element only in the upper triangle

                        if(row<=col){

                            val = Ka[j+m*sLen][k+n*sLen];
                            alglib::sparseadd(K,row,col,val);

                            val = Kb[j+m*sLen][k+n*sLen];
                            alglib::sparseadd(K,row+3*nLen,col+3*nLen,val);
                        }
                    }
                }
            }
        }
    }
    return K;
}





// ---------------------------------------------------|
// Computes the total mass matrix M for global FEM    |
// ---------------------------------------------------|

sparse Mesh::totalM(){

    sparse M;
    int sLen = shape3D.N.cols();
    int size = 9*eLen*pow(mesh.order+1,6)/4;
    alglib::sparsecreate(3*nLen,3*nLen,size,M);

    // Coordinates of the nodes of the element

    for(int i=0; i<eLen; i++){
        
        vector<array3d> eXYZ(sLen);
        for(int j=0; j<sLen; j++){eXYZ[j] = mesh.nXYZ[mesh.eNode[i][j]];}

        // Computes the elemental M matrices

        Elem elem(eXYZ,mesh.eSurf[i]);
        matrix M1 = elem.selfM(shape3D,mesh.EvR[i][2]);

        // Inserts the elemental matrix into the global M matrix

        for(int j=0; j<sLen; j++){
            for(int k=0; k<sLen; k++){
                for(int m=0; m<3; m++){
                    for(int n=0; n<3; n++){

                        int row = mesh.eNode[i][j]+m*nLen;
                        int col = mesh.eNode[i][k]+n*nLen;

                        // Adds the element only in the upper triangle

                        if(row<=col){

                            double val = M1[j+m*sLen][k+n*sLen];
                            alglib::sparseadd(M,row,col,val);
                        }
                    }
                }
            }
        }
    }
    return M;
}

// -------------------------------------------------------------|
// Computes the vector of applied stresses B from Neumann BC    |
// -------------------------------------------------------------|

darray Mesh::neumann(){

    darray B;
    int sLen = shape2D.N.cols();
    B.setlength(3*nLen);
    math::zero(B);

    // Coordinates of the nodes composing the faces

    for(int i=0; i<fLen; i++){

        vector<array3d> eXYZ(sLen);
        for(int j=0; j<sLen; j++){eXYZ[j] = mesh.nXYZ[mesh.neuFace[i][j]];}

        // Computes the elemental B vectors

        Face face(eXYZ,shape2D);
        darray B1 = face.selfB(shape2D,mesh.neuVal[i]);

        // Inserts the elemental vectors into the global B vector

        for(int j=0; j<sLen; j++){
            for(int k=0; k<3; k++){
                B(mesh.neuFace[i][j]+k*nLen) += B1(j+k*sLen);
            }
        }
    }
    return B;
}

// ----------------------------------------------------------------|
// Edits the matrix K and vector B for fixed nodal displacement    |
// ----------------------------------------------------------------|

void Mesh::dirichlet(sparse &K,darray &B){

    vector<ivector> row = math::sparsemap(K);

    // Gets the number of nodes in the dirichlet list

    for(int n=0; n<3; n++){
        int dLen = mesh.dirNode[n].size();

        // Gets the index of the node in the solution vector

        for(int i=0; i<dLen; i++){
            int idx = mesh.dirNode[n][i]+n*nLen;

            // Substract the coefficients from K to B and cancels the row-colums

            for(int j:row[idx]){

                B(j) -= math::get(K,j,idx)*mesh.dirVal[n][i];
                math::symset(K,j,idx,0);
            }
            row[idx].clear();
        }
    }

    // Places the dirichlet value the B vector and 1 in the diagonal of K

    for(int n=0; n<3; n++){
        for(int i=0; i<mesh.dirNode[n].size(); i++){
            
            int idx = mesh.dirNode[n][i]+n*nLen;
            alglib::sparseset(K,idx,idx,1);
            B(idx) = mesh.dirVal[n][i];
        }
    }
}

// ------------------------------------------------------------------|
// Edits the matrix K and vector B for coupled nodal displacement    |
// ------------------------------------------------------------------|

void Mesh::coupling(sparse &K,darray &B){

    vector<ivector> row = math::sparsemap(K);

    for(int n=0; n<3; n++){

        // Gets the number of nodes in the periodic list

        for(int i=0; i<mesh.coupNode[n].size(); i++){
            int dLen = mesh.coupNode[n][i].size();

            // Gets the index of the curent and final coupled nodes

            for(int j=0; j<dLen-1; j++){

                int idx1 = mesh.coupNode[n][i][j]+n*nLen;
                int idx2 = mesh.coupNode[n][i].back()+n*nLen;
                double fix = math::get(K,idx2,idx2)+math::get(K,idx1,idx1)+2*math::get(K,idx1,idx2);

                // Adds the current row to the final row of K if non-zero

                for(int k=0; k<row[idx1].size(); k++){

                    int n = row[idx1][k];
                    double val = math::get(K,idx1,n);

                    // Updates the mapStruct only if an element is added to K

                    if(val!=0){
                        if(math::get(K,idx2,n)==0){

                            row[n].push_back(idx2);
                            row[idx2].push_back(n);
                        }
                        math::symadd(K,idx2,n,val);
                    }
                }

                // Zeroes the column and row in K corresoponding to the coupled node

                for(int k:row[idx1]){math::symset(K,idx1,k,0);}
                row[idx1].clear();
    
                // Corrects the diagonal elements in K and edits B

                alglib::sparseset(K,idx2,idx2,fix);
                alglib::sparseset(K,idx1,idx1,1);
                B(idx2) += B(idx1);
                B(idx1) = 0;
            }
        }
    }
}

// -------------------------------------------------------------------|
// Edits the matrix K and vector B for a change of variable u = Δu    |
// -------------------------------------------------------------------|

void Mesh::delta(sparse &K,darray &B){

    vector<ivector> row = math::sparsemap(K);

    // Gets the number of nodes in the periodic list

    for(int n=0; n<3; n++){
        for(int i=0; i<mesh.deltaNode[n].size(); i++){

            // Gets the index of the curent and final coupled nodes

            int idx1 = mesh.deltaNode[n][i].first+n*nLen;
            int idx2 = mesh.deltaNode[n][i].second+n*nLen;
            double fix = math::get(K,idx2,idx2)+math::get(K,idx1,idx1)+2*math::get(K,idx1,idx2);

            // Adds the current row to the final row of K if non-zero

            for(int k=0; k<row[idx1].size(); k++){

                int n = row[idx1][k];
                double val = math::get(K,idx1,n);

                // Updates the mapStruct only if an element is added to K

                if(val!=0){
                    if(math::get(K,idx2,n)==0){

                        row[n].push_back(idx2);
                        row[idx2].push_back(n);
                    }
                    math::symadd(K,idx2,n,val);
                }
            }

            // Corrects the diagonal elements in K and edits B

            alglib::sparseset(K,idx2,idx2,fix);
            B(idx2) += B(idx1);
        }
    }
}

// -----------------------------------------------------|
// Completes the solution for coupled or delta nodes    |
// -----------------------------------------------------|

void Mesh::complete(darray &u){

    for(int n=0; n<3; n++){

        // Copy the value of the last node into the coupled ones

        for(int i=0; i<mesh.coupNode[n].size(); i++){
            for(int j=0; j<mesh.coupNode[n][i].size(); j++){

                int idx1 = mesh.coupNode[n][i][j]+n*nLen;
                int idx2 = mesh.coupNode[n][i].back()+n*nLen;
                u[idx1] = u[idx2];
            }
        }

        // Comes back from the difference Δu to the actual displacement

        for(int i=0; i<mesh.deltaNode[n].size(); i++){

            int idx1 = mesh.deltaNode[n][i].first+n*nLen;
            int idx2 = mesh.deltaNode[n][i].second+n*nLen;
            u[idx1] += u[idx2];
        }
    }
    
}

// ---------------------------------------------------------|
// Updates the coordinates of the node with the solution    |
// ---------------------------------------------------------|

void Mesh::update(darray &u){

    // Adds the displacement to the corresponding node coordinate

    for(int i=0; i<nLen; i++){
        for(int j=0; j<3; j++){
            mesh.nXYZ[i][j] += u[i+j*nLen];
        }
    }
}

// ----------------------------------------------------------|
// Computes the averaged Von Mises stress in each element    |
// ----------------------------------------------------------|

vector<darray> Mesh::stress(darray &u){

    vector<darray> sigma(eLen);
    int sLen = shape3D.N.cols();

    // Coordinates of the nodes of the element

    for(int i=0; i<eLen; i++){

        vector<array3d> eXYZ(sLen);
        for(int j=0; j<sLen; j++){eXYZ[j] = mesh.nXYZ[mesh.eNode[i][j]];}

        // Stores the displacement of the curent element nodes

        darray u1;
        u1.setlength(3*sLen);

        for(int j=0; j<sLen; j++){
            for(int k=0; k<3; k++){
                u1[j+k*sLen] = u(mesh.eNode[i][j]+k*nLen);
            }
        }

        // Computes the averaged Von Mises stress

        Elem elem(eXYZ,mesh.eSurf[i]);
        matrix D = math::stiffness(mesh.EvR[i]);
        sigma[i] = elem.stress(shape3D,D,u1);
    }
    return sigma;
}